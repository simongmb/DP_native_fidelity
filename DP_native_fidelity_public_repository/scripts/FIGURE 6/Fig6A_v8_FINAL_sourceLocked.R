# ============================================================
# Figure 6A v7 — DONOR-BLOCKED + CLEAN-SESSION
# Cross-validated coefficient prominence of non-programme genes
#
# FIX FROM v2:
# v2 used a percentile cutoff on tied coefficient ranks.
# In sparse elastic-net/lasso models, many coefficients are exactly
# zero and tied. That allowed large groups of zero-coefficient genes
# to be counted as "top 10%".
#
# v3 fixes this by:
#   - assigning ZERO prominence to exact-zero coefficients;
#   - ranking ONLY non-zero coefficients within each model;
#   - defining Top-10% membership as the strongest
#     min(ceil(10% of all predictors), number of non-zero genes)
#     non-zero coefficients;
#   - never filling the top-10% set with zero coefficients.
#
# This remains valid for:
#   ridge (all/most coefficients non-zero),
#   elastic net,
#   lasso.
#
# Interpretation:
# Repeated high prominence = stable model prominence.
# It does NOT establish causal regulation.
# Expression input is loaded reproducibly via scripts/load_subexpr_clean.R.
# ============================================================

required_pkgs <- c(
  "glmnet",
  "ggplot2",
  "readr",
  "dplyr",
  "tibble",
  "stringr",
  "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "\nMissing packages:\n",
      paste(missing_pkgs, collapse = ", "),
      "\n\nInstall with:\ninstall.packages(c(",
      paste0('"', missing_pkgs, '"', collapse = ", "),
      "))\n"
    )
  )
}

suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(scales)
})

EXPECTED_N_SAMPLES <- 110
EXPECTED_BACKGROUND_N <- 1553
EXPECTED_NATIVE_N <- 81
EXPECTED_ALT_N <- 99
EXPECTED_HELDOUT_N <- 1373

EXPECTED_REPEATS <- 20
EXPECTED_OUTER_K <- 10  # 10 held-out donors
EXPECTED_MODELS <- EXPECTED_REPEATS * EXPECTED_OUTER_K

TOP_N <- 20
TOP_FRACTION <- 0.10

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(
      dirname(
        normalizePath(
          sub("^--file=", "", file_arg[1])
        )
      )
    )
  }

  getwd()
}

start_dir <- get_script_dir()

candidate_roots <- unique(
  c(
    start_dir,
    getwd(),
    dirname(start_dir),
    dirname(dirname(start_dir)),
    dirname(dirname(dirname(start_dir))),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )
)

candidate_roots <- candidate_roots[
  dir.exists(candidate_roots)
]

score_root <- function(x) {
  s <- 0

  if (dir.exists(file.path(x, "figures"))) {
    s <- s + 2
  }

  if (dir.exists(file.path(x, "results"))) {
    s <- s + 2
  }

  if (
    length(
      list.files(
        x,
        pattern = "^Fig5C_outer_fold_metrics\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    ) > 0
  ) {
    s <- s + 10
  }

  if (
    length(
      list.files(
        x,
        pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    ) > 0
  ) {
    s <- s + 5
  }

  s
}

root_scores <- vapply(
  candidate_roots,
  score_root,
  numeric(1)
)

project_root <- candidate_roots[
  which.max(root_scores)
]

cat(
  "\nProject root:\n",
  project_root,
  "\n"
)

clean_gene <- function(x) {
  x <- toupper(
    stringr::str_trim(
      as.character(x)
    )
  )

  x[x == ""] <- NA_character_
  x
}

safe_read_csv <- function(path) {
  tryCatch(
    readr::read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    ),
    error = function(e) NULL
  )
}

pick_col <- function(df, candidates) {
  hit <- candidates[
    candidates %in% names(df)
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


# ------------------------------------------------------------
# Load clean-session expression matrix
# ------------------------------------------------------------

loader <- file.path(
  project_root,
  "scripts",
  "load_subexpr_clean.R"
)

if (!file.exists(loader)) {
  stop(
    paste0(
      "
QC STOP
Cannot find clean-session loader:
",
      loader,
      "
"
    )
  )
}

source(loader, local = FALSE)

if (!exists("subexpr", inherits = TRUE)) {
  stop(
    "
QC STOP
load_subexpr_clean.R did not create subexpr.
"
  )
}

expr <- get(
  "subexpr",
  inherits = TRUE
)

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

if (
  is.null(rownames(expr)) ||
  is.null(colnames(expr))
) {
  stop(
    "
QC STOP
subexpr must have gene row names and sample column names.
"
  )
}

rownames(expr) <- clean_gene(
  rownames(expr)
)

if (anyDuplicated(rownames(expr))) {
  stop(
    "
QC STOP
Duplicated gene symbols in subexpr.
"
  )
}

if (nrow(expr) != EXPECTED_BACKGROUND_N) {
  stop(
    paste0(
      "
QC STOP
Expected ",
      EXPECTED_BACKGROUND_N,
      " genes in clean-session subexpr, found ",
      nrow(expr),
      ".
"
    )
  )
}

if (ncol(expr) != EXPECTED_N_SAMPLES) {
  stop(
    paste0(
      "
QC STOP
Expected ",
      EXPECTED_N_SAMPLES,
      " samples, found ",
      ncol(expr),
      ".
"
    )
  )
}

cat(
  "
Clean-session expression matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples
"
)


# ------------------------------------------------------------
# Recover frozen background
# ------------------------------------------------------------
audit_files <- list.files(
  project_root,
  pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(audit_files) == 0) {
  stop(
    "\nQC STOP\nCannot find 02_gene_level_native_fidelity_audit.csv\n"
  )
}

audit <- safe_read_csv(audit_files[1])

if (is.null(audit)) {
  stop("\nQC STOP\nCannot read audit file.\n")
}

gene_col <- pick_col(
  audit,
  c(
    "GeneSymbol",
    "gene_symbol",
    "Gene",
    "gene",
    "SYMBOL",
    "symbol"
  )
)

if (is.na(gene_col)) {
  stop("\nQC STOP\nCannot identify gene column in audit file.\n")
}

background_genes <- unique(
  clean_gene(
    audit[[gene_col]]
  )
)

background_genes <- background_genes[
  !is.na(background_genes)
]

if (length(background_genes) != EXPECTED_BACKGROUND_N) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      EXPECTED_BACKGROUND_N,
      " background genes, found ",
      length(background_genes),
      ".\n"
    )
  )
}

if (
  length(
    setdiff(
      background_genes,
      rownames(expr)
    )
  ) > 0
) {
  stop(
    "\nQC STOP\nFrozen background is not fully present in subexpr.\n"
  )
}

# ------------------------------------------------------------
# Recover frozen 81/99 programmes
# ------------------------------------------------------------
all_csv <- list.files(
  project_root,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

read_exact_gene_list <- function(files, expected_n) {
  if (length(files) == 0) {
    return(character(0))
  }

  for (f in files) {
    tab <- safe_read_csv(f)

    if (
      is.null(tab) ||
      nrow(tab) == 0
    ) {
      next
    }

    gc <- pick_col(
      tab,
      c(
        "GeneSymbol",
        "gene_symbol",
        "Gene",
        "gene",
        "SYMBOL",
        "symbol"
      )
    )

    if (is.na(gc)) {
      next
    }

    genes <- unique(
      clean_gene(
        tab[[gc]]
      )
    )

    genes <- genes[
      !is.na(genes)
    ]

    if (length(genes) == expected_n) {
      cat(
        "\nRecovered n =",
        expected_n,
        "gene list from:\n",
        f,
        "\n"
      )

      return(genes)
    }
  }

  character(0)
}

native_files <- all_csv[
  grepl(
    "native.*81|81.*native|native_identity",
    basename(all_csv),
    ignore.case = TRUE
  )
]

alt_files <- all_csv[
  grepl(
    "alternative.*99|99.*alternative|alternative_state",
    basename(all_csv),
    ignore.case = TRUE
  )
]

native81 <- read_exact_gene_list(
  native_files,
  EXPECTED_NATIVE_N
)

alternative99 <- read_exact_gene_list(
  alt_files,
  EXPECTED_ALT_N
)

if (
  length(native81) != EXPECTED_NATIVE_N ||
  length(alternative99) != EXPECTED_ALT_N
) {
  exact_rows <- list()

  for (f in all_csv) {
    tab <- safe_read_csv(f)

    if (
      is.null(tab) ||
      nrow(tab) == 0
    ) {
      next
    }

    gc <- pick_col(
      tab,
      c(
        "GeneSymbol",
        "gene_symbol",
        "Gene",
        "gene",
        "SYMBOL",
        "symbol"
      )
    )

    if (is.na(gc)) {
      next
    }

    genes <- unique(
      clean_gene(
        tab[[gc]]
      )
    )

    genes <- genes[
      !is.na(genes)
    ]

    if (
      length(genes) %in%
        c(
          EXPECTED_NATIVE_N,
          EXPECTED_ALT_N
        )
    ) {
      exact_rows <- append(
        exact_rows,
        list(
          tibble::tibble(
            file = f,
            n = length(genes)
          )
        )
      )
    }
  }

  if (length(exact_rows) > 0) {
    exact_df <- dplyr::bind_rows(
      exact_rows
    ) %>%
      dplyr::distinct()

    if (length(native81) != EXPECTED_NATIVE_N) {
      c81 <- exact_df %>%
        dplyr::filter(
          n == EXPECTED_NATIVE_N
        ) %>%
        dplyr::mutate(
          score =
            10 *
              grepl(
                "FIGURE.?3|FIGURE.?4|Fig3|Fig4",
                file,
                ignore.case = TRUE
              ) +
            5 *
              grepl(
                "native|identity",
                file,
                ignore.case = TRUE
              )
        ) %>%
        dplyr::arrange(
          dplyr::desc(score)
        )

      native81 <- read_exact_gene_list(
        c81$file,
        EXPECTED_NATIVE_N
      )
    }

    if (length(alternative99) != EXPECTED_ALT_N) {
      c99 <- exact_df %>%
        dplyr::filter(
          n == EXPECTED_ALT_N
        ) %>%
        dplyr::mutate(
          score =
            10 *
              grepl(
                "FIGURE.?3|FIGURE.?4|Fig3|Fig4",
                file,
                ignore.case = TRUE
              ) +
            5 *
              grepl(
                "alternative|state",
                file,
                ignore.case = TRUE
              )
        ) %>%
        dplyr::arrange(
          dplyr::desc(score)
        )

      alternative99 <- read_exact_gene_list(
        c99$file,
        EXPECTED_ALT_N
      )
    }
  }
}

if (
  length(native81) != EXPECTED_NATIVE_N ||
  length(alternative99) != EXPECTED_ALT_N
) {
  stop(
    "\nQC STOP\nCould not safely recover frozen 81/99 programs.\n"
  )
}

program_genes <- unique(
  c(
    native81,
    alternative99
  )
)

if (
  length(program_genes) !=
    EXPECTED_NATIVE_N +
    EXPECTED_ALT_N
) {
  stop(
    "\nQC STOP\nUnexpected overlap between frozen programs.\n"
  )
}

# ------------------------------------------------------------
# Define non-programme predictors
# ------------------------------------------------------------
heldout_genes <- setdiff(
  background_genes,
  program_genes
)

if (length(heldout_genes) != EXPECTED_HELDOUT_N) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      EXPECTED_HELDOUT_N,
      " non-programme genes, found ",
      length(heldout_genes),
      ".\n"
    )
  )
}

x <- t(
  expr[
    heldout_genes,
    ,
    drop = FALSE
  ]
)

global_sd <- apply(
  x,
  2,
  sd,
  na.rm = TRUE
)

keep_global <- is.finite(global_sd) &
  global_sd > 0

x <- x[
  ,
  keep_global,
  drop = FALSE
]

if (anyNA(x)) {
  stop(
    "\nQC STOP\nMissing values in non-programme predictor matrix.\n"
  )
}

if (ncol(x) != EXPECTED_HELDOUT_N) {
  stop(
    paste0(
      "\nQC STOP\nExpected all ",
      EXPECTED_HELDOUT_N,
      " non-programme genes after variance QC, found ",
      ncol(x),
      ".\n"
    )
  )
}

# ------------------------------------------------------------
# Recover Figure 5C outer model settings
# ------------------------------------------------------------
metric_files <- list.files(
  project_root,
  pattern = "^Fig5C_outer_fold_metrics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(metric_files) == 0) {
  stop(
    "\nQC STOP\nCannot find Fig5C_outer_fold_metrics.csv\n"
  )
}

score_metric_file <- function(path) {
  s <- 0

  if (grepl("stateBalanceReconstruction_v12|Panel_C_stateBalanceReconstruction_v12", path, ignore.case = TRUE)) {
    s <- s + 200
  }

  if (grepl("donor_blocked_v8|Panel_C_donor", path, ignore.case = TRUE)) {
    s <- s + 100
  }

  if (grepl("Panel_C_v3", path, ignore.case = TRUE)) {
    s <- s + 30
  }

  if (grepl("Panel_C_v2", path, ignore.case = TRUE)) {
    s <- s + 20
  }

  if (grepl("FIGURE 5", path, ignore.case = TRUE)) {
    s <- s + 10
  }

  s
}

metric_scores <- vapply(
  metric_files,
  score_metric_file,
  numeric(1)
)

metric_file <- metric_files[
  which.max(metric_scores)
]

outer_metrics <- safe_read_csv(metric_file)

required_cols <- c(
  "Repeat",
  "Fold",
  "Test_donor",
  "alpha",
  "lambda"
)

if (
  is.null(outer_metrics) ||
  !all(
    required_cols %in%
      names(outer_metrics)
  )
) {
  stop(
    "\nQC STOP\nInvalid Figure 5C outer-fold metrics.\n"
  )
}

outer_metrics <- outer_metrics %>%
  dplyr::select(
    Repeat,
    Fold,
    Test_donor,
    alpha,
    lambda
  ) %>%
  dplyr::distinct()

if (nrow(outer_metrics) != EXPECTED_MODELS) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      EXPECTED_MODELS,
      " outer models, found ",
      nrow(outer_metrics),
      ".\n"
    )
  )
}

alpha_summary <- outer_metrics %>%
  dplyr::count(
    alpha,
    name = "N_models"
  ) %>%
  dplyr::arrange(alpha)

cat(
  "\nFigure 5C alpha distribution:\n"
)

print(
  alpha_summary,
  n = Inf
)

# ------------------------------------------------------------
# Recover Figure 5C target
# ------------------------------------------------------------
pred_files <- list.files(
  project_root,
  pattern = "^Fig5C_sample_predictions\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(pred_files) == 0) {
  stop(
    "\nQC STOP\nCannot find Fig5C_sample_predictions.csv\n"
  )
}

score_pred_file <- function(path) {
  s <- 0

  if (grepl("stateBalanceReconstruction_v12|Panel_C_stateBalanceReconstruction_v12", path, ignore.case = TRUE)) {
    s <- s + 200
  }

  if (grepl("donor_blocked_v8|Panel_C_donor", path, ignore.case = TRUE)) {
    s <- s + 100
  }

  if (grepl("Panel_C_v3", path, ignore.case = TRUE)) {
    s <- s + 30
  }

  if (grepl("Panel_C_v2", path, ignore.case = TRUE)) {
    s <- s + 20
  }

  if (grepl("FIGURE 5", path, ignore.case = TRUE)) {
    s <- s + 10
  }

  s
}

pred_scores <- vapply(
  pred_files,
  score_pred_file,
  numeric(1)
)

pred_file <- pred_files[
  which.max(pred_scores)
]

pred_tab <- safe_read_csv(pred_file)

if (
  is.null(pred_tab) ||
  !"Observed_balance" %in%
    names(pred_tab) ||
  nrow(pred_tab) != EXPECTED_N_SAMPLES
) {
  stop(
    "\nQC STOP\nInvalid Figure 5C prediction table.\n"
  )
}

y <- as.numeric(
  pred_tab$Observed_balance
)

expr_samples <- colnames(expr)

extract_donor <- function(sample_ids) {
  stringr::str_match(
    sample_ids,
    "^.*_([A-J])_[^_]+$"
  )[, 2]
}

donor_ids <- extract_donor(expr_samples)

if (anyNA(donor_ids) || length(unique(donor_ids)) != 10) {
  stop("\nQC STOP\nCould not recover exactly 10 donors from sample IDs.\n")
}

if (
  "Sample" %in%
    names(pred_tab) &&
  !is.null(expr_samples) &&
  setequal(
    as.character(
      pred_tab$Sample
    ),
    expr_samples
  )
) {
  y <- y[
    match(
      expr_samples,
      as.character(
        pred_tab$Sample
      )
    )
  ]

  target_alignment <-
    "Exact sample-name alignment"
} else {
  target_alignment <-
    "Figure 5C row order retained"
}

if (
  length(y) != EXPECTED_N_SAMPLES ||
  any(!is.finite(y))
) {
  stop(
    "\nQC STOP\nInvalid Figure 5C target vector.\n"
  )
}

# ------------------------------------------------------------
# Containers
# ------------------------------------------------------------
gene_names <- colnames(x)
P <- length(gene_names)

prominence_sum <- setNames(
  rep(0, P),
  gene_names
)

prominence_values <- setNames(
  vector("list", P),
  gene_names
)

top10_count <- setNames(
  rep(0L, P),
  gene_names
)

nonzero_count <- setNames(
  rep(0L, P),
  gene_names
)

positive_count <- setNames(
  rep(0L, P),
  gene_names
)

negative_count <- setNames(
  rep(0L, P),
  gene_names
)

abs_effect_sum <- setNames(
  rep(0, P),
  gene_names
)

effect_values <- setNames(
  vector("list", P),
  gene_names
)

model_qc_rows <- list()

TOP_K_MAX <- ceiling(
  TOP_FRACTION *
    P
)

# ------------------------------------------------------------
# Refit the same 200 outer models
# ------------------------------------------------------------
cat(
  "\nRefitting",
  EXPECTED_MODELS,
  "Figure 5C outer models...\n"
)

for (r in seq_len(EXPECTED_REPEATS)) {

  for (fold in seq_len(EXPECTED_OUTER_K)) {
    tune_row <- outer_metrics %>%
      dplyr::filter(
        Repeat == r,
        Fold == fold
      )

    if (nrow(tune_row) != 1) {
      stop(
        paste0(
          "\nQC STOP\nNo unique tuning row for repeat ",
          r,
          ", fold ",
          fold,
          ".\n"
        )
      )
    }

    test_donor <- as.character(
      tune_row$Test_donor
    )

    test_idx <- which(
      donor_ids == test_donor
    )

    train_idx <- which(
      donor_ids != test_donor
    )

    if (
      length(test_idx) == 0 ||
      any(donor_ids[test_idx] %in% donor_ids[train_idx])
    ) {
      stop(
        paste0(
          "\nQC STOP\nInvalid donor-blocked split for repeat ",
          r,
          ", fold ",
          fold,
          ".\n"
        )
      )
    }

    x_train <- x[
      train_idx,
      ,
      drop = FALSE
    ]

    y_train <- y[
      train_idx
    ]

    train_sd <- apply(
      x_train,
      2,
      sd
    )

    keep_train <- is.finite(train_sd) &
      train_sd > 0

    x_train2 <- x_train[
      ,
      keep_train,
      drop = FALSE
    ]

    train_sd2 <- train_sd[
      keep_train
    ]

    fit <- glmnet::glmnet(
      x = x_train2,
      y = y_train,
      family = "gaussian",
      alpha = as.numeric(
        tune_row$alpha
      ),
      lambda = as.numeric(
        tune_row$lambda
      ),
      standardize = TRUE,
      intercept = TRUE
    )

    beta <- as.matrix(
      coef(
        fit,
        s = as.numeric(
          tune_row$lambda
        )
      )
    )

    beta_vec <- beta[, 1]

    beta_vec <- beta_vec[
      names(beta_vec) !=
        "(Intercept)"
    ]

    std_effect <- beta_vec *
      train_sd2[
        names(beta_vec)
      ]

    full_effect <- setNames(
      rep(0, P),
      gene_names
    )

    full_effect[
      names(std_effect)
    ] <- std_effect

    abs_effect <- abs(
      full_effect
    )

    nz <- names(full_effect)[
      abs_effect >
        1e-12
    ]

    n_nz <- length(nz)

    # --------------------------------------------------------
    # Corrected prominence:
    # exact-zero genes always receive 0.
    # non-zero genes are ranked only among non-zero coefficients.
    #
    # Strongest non-zero coefficient = 1.
    # Weakest non-zero coefficient > 0.
    # --------------------------------------------------------
    prominence <- setNames(
      rep(0, P),
      gene_names
    )

    if (n_nz > 0) {
      # deterministic ordering:
      # first by descending |effect|, then gene symbol.
      ord_df <- tibble::tibble(
        GeneSymbol = nz,
        AbsEffect = abs_effect[nz]
      ) %>%
        dplyr::arrange(
          dplyr::desc(
            AbsEffect
          ),
          GeneSymbol
        )

      ordered_genes <- ord_df$GeneSymbol

      if (n_nz == 1) {
        prominence[
          ordered_genes
        ] <- 1
      } else {
        prominence[
          ordered_genes
        ] <- seq(
          from = 1,
          to = 1 / n_nz,
          length.out = n_nz
        )
      }

      # Exact top-10%:
      # use only non-zero genes;
      # NEVER fill with zero coefficients.
      k_this <- min(
        TOP_K_MAX,
        n_nz
      )

      top_genes <- ordered_genes[
        seq_len(k_this)
      ]

      top10_count[
        top_genes
      ] <- top10_count[
        top_genes
      ] +
        1L

      nonzero_count[
        nz
      ] <- nonzero_count[
        nz
      ] +
        1L

      pos <- nz[
        full_effect[nz] >
          0
      ]

      neg <- nz[
        full_effect[nz] <
          0
      ]

      if (length(pos) > 0) {
        positive_count[
          pos
        ] <- positive_count[
          pos
        ] +
          1L
      }

      if (length(neg) > 0) {
        negative_count[
          neg
        ] <- negative_count[
          neg
        ] +
          1L
      }

      abs_effect_sum[
        nz
      ] <- abs_effect_sum[
        nz
      ] +
        abs_effect[
          nz
        ]

      for (g in nz) {
        effect_values[[g]] <- c(
          effect_values[[g]],
          full_effect[g]
        )
      }
    }

    prominence_sum <- prominence_sum +
      prominence

    for (g in gene_names) {
      prominence_values[[g]] <- c(
        prominence_values[[g]],
        prominence[g]
      )
    }

    model_qc_rows <- append(
      model_qc_rows,
      list(
        tibble::tibble(
          Repeat = r,
          Fold = fold,
          alpha = as.numeric(
            tune_row$alpha
          ),
          lambda = as.numeric(
            tune_row$lambda
          ),
          N_train = length(train_idx),
          N_predictors = ncol(x_train2),
          N_nonzero = n_nz,
          N_top10 = ifelse(
            n_nz > 0,
            min(
              TOP_K_MAX,
              n_nz
            ),
            0
          ),
          Median_abs_standardized_effect =
            median(
              abs_effect
            ),
          Max_abs_standardized_effect =
            max(
              abs_effect
            )
        )
      )
    )
  }

  cat(
    "Completed repeat",
    r,
    "/",
    EXPECTED_REPEATS,
    "\n"
  )
}

model_qc <- dplyr::bind_rows(
  model_qc_rows
)

# ------------------------------------------------------------
# Aggregate metrics
# ------------------------------------------------------------
safe_median <- function(z) {
  if (length(z) == 0) {
    return(0)
  }

  median(z)
}

safe_q025 <- function(z) {
  if (length(z) == 0) {
    return(0)
  }

  unname(
    quantile(
      z,
      0.025
    )
  )
}

safe_q975 <- function(z) {
  if (length(z) == 0) {
    return(0)
  }

  unname(
    quantile(
      z,
      0.975
    )
  )
}

prominence_df <- tibble::tibble(
  GeneSymbol = gene_names,
  Mean_prominence =
    prominence_sum /
    EXPECTED_MODELS,
  Median_prominence =
    vapply(
      prominence_values,
      safe_median,
      numeric(1)
    ),
  Top10_count =
    as.integer(
      top10_count
    ),
  Top10_frequency =
    as.numeric(
      top10_count
    ) /
    EXPECTED_MODELS,
  Nonzero_count =
    as.integer(
      nonzero_count
    ),
  Nonzero_frequency =
    as.numeric(
      nonzero_count
    ) /
    EXPECTED_MODELS,
  Positive_count =
    as.integer(
      positive_count
    ),
  Negative_count =
    as.integer(
      negative_count
    ),
  Mean_abs_standardized_effect =
    ifelse(
      nonzero_count > 0,
      abs_effect_sum /
        nonzero_count,
      0
    ),
  Median_standardized_effect =
    vapply(
      effect_values,
      safe_median,
      numeric(1)
    ),
  Effect_q025 =
    vapply(
      effect_values,
      safe_q025,
      numeric(1)
    ),
  Effect_q975 =
    vapply(
      effect_values,
      safe_q975,
      numeric(1)
    )
) %>%
  dplyr::mutate(
    Sign_consistency =
      dplyr::case_when(
        Nonzero_count == 0 ~
          0,
        TRUE ~
          pmax(
            Positive_count,
            Negative_count
          ) /
          Nonzero_count
      ),
    Direction =
      dplyr::case_when(
        Median_standardized_effect > 0 ~
          "Native-aligned",
        Median_standardized_effect < 0 ~
          "Alternative-aligned",
        TRUE ~
          "Neutral"
      )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Mean_prominence
    ),
    dplyr::desc(
      Top10_frequency
    ),
    dplyr::desc(
      Sign_consistency
    )
  )

# ------------------------------------------------------------
# Select top 20
# ------------------------------------------------------------
plot_df <- prominence_df %>%
  dplyr::slice_head(
    n = TOP_N
  ) %>%
  dplyr::arrange(
    Mean_prominence
  )

plot_df$GeneSymbol <- factor(
  plot_df$GeneSymbol,
  levels = unique(
    plot_df$GeneSymbol
  )
)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 6",
  "Panel_A_prominence_donorBlocked_v7"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  prominence_df,
  file.path(
    out_dir,
    "Fig6A_all_nonprogramme_prominence_donorBlocked_v7.csv"
  )
)

readr::write_csv(
  plot_df,
  file.path(
    out_dir,
    "Fig6A_top20_prominence_donorBlocked_v7.csv"
  )
)

readr::write_csv(
  model_qc,
  file.path(
    out_dir,
    "Fig6A_outer_model_QC_donorBlocked_v7.csv"
  )
)

readr::write_csv(
  alpha_summary,
  file.path(
    out_dir,
    "Fig6A_alpha_distribution_donorBlocked_v7.csv"
  )
)

# ------------------------------------------------------------
# Plot -- Nature-style Figure 6A
# ------------------------------------------------------------
# Scientific calculations above are unchanged.
# This section only refines the visual layout for the final 2 x 2 Figure 6 assembly.

COL_NATIVE <- "#3D8C82"
COL_ALT <- "#C97857"
COL_NEUTRAL <- "#83949B"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"
COL_GRID <- "#E6ECEF"

p <- ggplot(
  plot_df,
  aes(
    x = Mean_prominence,
    y = GeneSymbol
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Mean_prominence,
      yend = GeneSymbol
    ),
    linewidth = 0.30,
    colour = "#D8E1E5"
  ) +
  geom_point(
    aes(
      size = Top10_frequency,
      colour = Direction,
      alpha = Sign_consistency
    ),
    stroke = 0
  ) +
  scale_colour_manual(
    values = c(
      "Native-aligned" = COL_NATIVE,
      "Alternative-aligned" = COL_ALT,
      "Neutral" = COL_NEUTRAL
    ),
    breaks = c(
      "Alternative-aligned",
      "Native-aligned"
    ),
    labels = c(
      "Alternative-aligned",
      "Native-aligned"
    ),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1.7, 3.6),
    limits = c(0, 1),
    breaks = c(0.60, 0.80, 1.00),
    labels = scales::label_percent(accuracy = 1),
    name = "Top-10% frequency"
  ) +
  scale_alpha_continuous(
    range = c(0.50, 1.00),
    limits = c(0, 1),
    breaks = c(0, 0.50, 1.00),
    labels = scales::label_percent(accuracy = 1),
    name = "Sign consistency"
  ) +
  scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.025))
  ) +
  labs(
    x = "Mean model prominence across 200 outer-CV models",
    y = NULL,
    title = "Recurrent non-programme genes show high model prominence",
    subtitle = "Top 20 non-programme genes across Figure 5C donor-blocked models"
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      override.aes = list(
        colour = "#20343D",
        alpha = 1
      )
    ),
    colour = guide_legend(
      order = 2,
      title.position = "top",
      override.aes = list(
        size = 2.2,
        alpha = 1
      )
    ),
    alpha = guide_legend(
      order = 3,
      title.position = "top",
      override.aes = list(
        size = 2.2,
        colour = "#20343D"
      )
    )
  ) +
  theme_classic(
    base_family = "Helvetica",
    base_size = 5.0
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      size = 6.0,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      lineheight = 1.0,
      margin = margin(b = 1.0, unit = "pt")
    ),
    plot.subtitle = element_text(
      size = 4.45,
      colour = COL_MUTED,
      hjust = 0.5,
      margin = margin(b = 2.2, unit = "pt")
    ),
    axis.title.x = element_text(
      size = 4.8,
      colour = COL_TEXT,
      margin = margin(t = 1.8, unit = "pt")
    ),
    axis.text.x = element_text(
      size = 4.35,
      colour = COL_MUTED
    ),
    axis.text.y = element_text(
      size = 4.55,
      colour = COL_TEXT,
      face = "bold",
      margin = margin(r = 1.8, unit = "pt")
    ),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.30,
      colour = COL_TEXT
    ),
    axis.ticks.x = element_line(
      linewidth = 0.25,
      colour = COL_TEXT
    ),
    axis.ticks.length = grid::unit(1.2, "pt"),
    panel.grid.major.x = element_line(
      linewidth = 0.22,
      colour = COL_GRID
    ),
    panel.grid.major.y = element_blank(),

    # Keep the three encodings, but make the right-side legend column much tighter.
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = grid::unit(1.6, "pt"),
    legend.margin = margin(0, 0, 0, 3, unit = "pt"),
    legend.title = element_text(
      size = 4.35,
      colour = COL_TEXT,
      margin = margin(b = 0.8, unit = "pt")
    ),
    legend.text = element_text(
      size = 4.15,
      colour = COL_MUTED
    ),
    legend.key.height = grid::unit(6.0, "pt"),
    legend.key.width = grid::unit(7.0, "pt"),

    plot.margin = margin(
      t = 2.2,
      r = 2.0,
      b = 2.0,
      l = 2.5,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# Nature-style export
# ------------------------------------------------------------
# One half of an ~183-mm double-column figure.
# A/B will be aligned to the same final panel dimensions.
FIG_W_MM <- 88
FIG_H_MM <- 62

out_dir_nature <- out_dir

dir.create(
  out_dir_nature,
  recursive = TRUE,
  showWarnings = FALSE
)

png_file <- file.path(
  out_dir_nature,
  "Fig6_PanelA_prominence_donorBlocked_v7.png"
)

pdf_file <- file.path(
  out_dir_nature,
  "Fig6_PanelA_prominence_donorBlocked_v7.pdf"
)

tiff_file <- file.path(
  out_dir_nature,
  "Fig6_PanelA_prominence_donorBlocked_v7.tiff"
)

ggsave(
  filename = png_file,
  plot = p,
  width = FIG_W_MM,
  height = FIG_H_MM,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# Use the base PDF device for macOS compatibility.
# This avoids cairo DLL failures seen in earlier figure export.
grDevices::pdf(
  file = pdf_file,
  width = FIG_W_MM / 25.4,
  height = FIG_H_MM / 25.4,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE,
  paper = "special"
)
print(p)
grDevices::dev.off()

if (!file.exists(pdf_file) ||
    is.na(file.info(pdf_file)$size) ||
    file.info(pdf_file)$size <= 0) {
  stop(
    paste0(
      "PDF export failed: ",
      pdf_file
    )
  )
}

ggsave(
  filename = tiff_file,
  plot = p,
  width = FIG_W_MM,
  height = FIG_H_MM,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

cat(
  "\nNature-style Figure 6A written to:\n",
  out_dir_nature,
  "\n"
)

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------
qc <- tibble::tibble(
  Metric = c(
    "Samples",
    "Frozen background genes",
    "Native programme genes",
    "Alternative programme genes",
    "Non-programme genes",
    "Figure 5C outer models",
    "Models with alpha = 0",
    "Models with alpha > 0",
    "Median non-zero coefficients per model",
    "Minimum non-zero coefficients per model",
    "Maximum non-zero coefficients per model",
    "Median top-set size per model",
    "Minimum top-set size per model",
    "Maximum top-set size per model",
    "Genes with top-10% frequency >= 0.50",
    "Genes with top-10% frequency >= 0.75",
    "Genes with mean prominence >= 0.50",
    "Genes with sign consistency >= 0.90"
  ),
  Value = c(
    EXPECTED_N_SAMPLES,
    length(background_genes),
    length(native81),
    length(alternative99),
    length(heldout_genes),
    nrow(model_qc),
    sum(
      model_qc$alpha == 0
    ),
    sum(
      model_qc$alpha > 0
    ),
    median(
      model_qc$N_nonzero
    ),
    min(
      model_qc$N_nonzero
    ),
    max(
      model_qc$N_nonzero
    ),
    median(
      model_qc$N_top10
    ),
    min(
      model_qc$N_top10
    ),
    max(
      model_qc$N_top10
    ),
    sum(
      prominence_df$Top10_frequency >=
        0.50
    ),
    sum(
      prominence_df$Top10_frequency >=
        0.75
    ),
    sum(
      prominence_df$Mean_prominence >=
        0.50
    ),
    sum(
      prominence_df$Sign_consistency >=
        0.90
    )
  )
)

readr::write_csv(
  qc,
  file.path(
    out_dir,
    "Fig6A_QC_summary_donorBlocked_v7.csv"
  )
)

readme <- c(
  "Figure 6A v7 — donor-blocked complete archive",
  "",
  "Stable model prominence of non-programme genes under donor-blocked cross-validation.",
  "",
  "Critical v3 correction:",
  "Zero coefficients are never ranked as model features.",
  "Only non-zero coefficients receive within-model prominence scores.",
  "Top-10% sets contain only non-zero genes and are capped at 10% of all predictors.",
  "",
  "This resolves the tied-zero inflation discovered in Figure 6A v2.",
  "",
  "Metrics:",
  "- Mean prominence: average within-model non-zero coefficient prominence; zero if not selected.",
  "- Top-10% frequency: fraction of the 200 outer models in which the gene is among the strongest non-zero coefficients.",
  "- Sign consistency: consistency of coefficient direction among models in which the gene has a non-zero coefficient.",
  "",
  "Interpretation:",
  "High prominence indicates stable model importance across leave-one-donor-out donor-blocked fits.",
  "It does not establish causal or regulatory activity.",
  "",
  paste0(
    "Figure 5C tuning source: ",
    normalizePath(
      metric_file
    )
  ),
  paste0(
    "Figure 5C target source: ",
    normalizePath(
      pred_file
    )
  ),
  paste0(
    "Target alignment: ",
    target_alignment
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig6_PanelA_donorBlocked_v7.txt"
  )
)

cat(
  "\n============================================\n"
)

cat(
  "FIGURE 6A v7 DONOR-BLOCKED COMPLETE\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Alpha distribution:\n"
)

print(
  alpha_summary,
  n = Inf
)

cat(
  "\nQC summary:\n"
)

print(
  qc,
  n = Inf
)

cat(
  "\nTop 20 predictive-prominence genes:\n"
)

print(
  prominence_df %>%
    dplyr::select(
      GeneSymbol,
      Mean_prominence,
      Top10_frequency,
      Nonzero_frequency,
      Sign_consistency,
      Median_standardized_effect,
      Direction
    ) %>%
    dplyr::slice_head(
      n = 20
    ),
  n = 20
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)

cat(
  "\nMain figure:\n",
  file.path(
    out_dir,
    "Fig6_PanelA_prominence_v3.pdf"
  ),
  "\n"
)
