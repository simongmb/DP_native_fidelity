# ============================================================
# Figure 5C v12 — DONOR-BLOCKED + CLEAN-SESSION
# Donor-blocked reconstruction of the prespecified anatomy-referenced state balance
#
# KEY CHANGE FROM v7:
#   Outer validation is leave-one-donor-out (10 donors: A-J).
#   All profiles from the same donor are always kept together.
#   Inner CV is also donor-grouped.
#   Bootstrap confidence intervals are donor-clustered.
#   Expression input is loaded reproducibly via scripts/load_subexpr_clean.R.
#
# Scientific question:
# Can the 1,373 non-programme genes reconstruct the prespecified Figure 4
# anatomy-referenced state balance without using any of the
# 81 + 99 defining programme genes as predictors?
#
# IMPORTANT:
# This is an INTERNAL conditional reconstruction analysis within GSE272333.
# The target axis is prespecified from the same atlas and is NOT externally validated.
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
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

MASTER_SEED <- 20261028
set.seed(MASTER_SEED)

# ------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------
EXPECTED_N_SAMPLES <- 110
EXPECTED_N_DONORS <- 10
EXPECTED_NATIVE_N <- 81
EXPECTED_ALT_N <- 99
EXPECTED_HELDOUT_N <- 1373

# Keep 20 repeats so downstream S5B remains compatible.
# The outer split is always donor-based; repeats mainly change
# the donor-grouped inner-CV tuning partitions.
N_REPEATS <- 20
INNER_K <- 5

ALPHA_GRID <- c(0, 0.25, 0.5, 0.75, 1)

# ------------------------------------------------------------
# 2. Locate project root
# ------------------------------------------------------------
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
    s <- s + 3
  }

  if (dir.exists(file.path(x, "results"))) {
    s <- s + 3
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

cat("\nProject root:\n", project_root, "\n")

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------
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

safe_cor <- function(x, y, method = "pearson") {
  suppressWarnings(
    cor(
      x,
      y,
      method = method,
      use = "complete.obs"
    )
  )
}

calc_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)

  obs <- obs[ok]
  pred <- pred[ok]

  pearson_r <- safe_cor(
    obs,
    pred,
    method = "pearson"
  )

  spearman_rho <- safe_cor(
    obs,
    pred,
    method = "spearman"
  )

  sse <- sum(
    (obs - pred)^2
  )

  sst <- sum(
    (obs - mean(obs))^2
  )

  r2 <- if (
    is.finite(sst) &&
    sst > 0
  ) {
    1 - sse / sst
  } else {
    NA_real_
  }

  mae <- mean(
    abs(obs - pred)
  )

  rmse <- sqrt(
    mean(
      (obs - pred)^2
    )
  )

  tibble::tibble(
    Pearson_r = pearson_r,
    Spearman_rho = spearman_rho,
    CV_R2 = r2,
    MAE = mae,
    RMSE = rmse
  )
}

# Extract donor from GSE272333 sample IDs such as:
# DP_B_6, HM1_A_3, ORS2_J_2
extract_donor <- function(sample_ids) {
  donor <- stringr::str_match(
    sample_ids,
    "^.*_([A-J])_[^_]+$"
  )[, 2]

  donor
}

# Build sample-level fold IDs while keeping each donor intact.
make_group_folds <- function(groups, k, seed) {
  groups <- as.character(groups)

  if (anyNA(groups) || any(groups == "")) {
    stop("\nQC STOP\nMissing donor labels in grouped CV.\n")
  }

  unique_groups <- unique(groups)

  if (length(unique_groups) < k) {
    stop(
      paste0(
        "\nQC STOP\n",
        "Need at least ",
        k,
        " donors for ",
        k,
        "-fold grouped inner CV, but only ",
        length(unique_groups),
        " are available.\n"
      )
    )
  }

  set.seed(seed)

  shuffled <- sample(
    unique_groups,
    length(unique_groups),
    replace = FALSE
  )

  donor_fold <- rep(
    seq_len(k),
    length.out = length(shuffled)
  )

  names(donor_fold) <- shuffled

  as.integer(
    donor_fold[groups]
  )
}

# ------------------------------------------------------------
# 4. Load clean-session expression matrix
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
Duplicated gene symbols detected in subexpr row names.
"
  )
}

if (nrow(expr) != 1553L) {
  stop(
    paste0(
      "
QC STOP
Expected exactly 1553 genes in clean-session subexpr, but found ",
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
Expected exactly ",
      EXPECTED_N_SAMPLES,
      " samples in subexpr, but found ",
      ncol(expr),
      ".
"
    )
  )
}

sample_ids <- colnames(expr)
donor_ids <- extract_donor(sample_ids)

if (anyNA(donor_ids)) {
  bad <- sample_ids[is.na(donor_ids)]

  stop(
    paste0(
      "
QC STOP
",
      "Could not infer donor from these sample IDs:
",
      paste(bad, collapse = ", "),
      "

Expected pattern like DP_B_6 or ORS2_J_2.
"
    )
  )
}

donors <- sort(unique(donor_ids))

if (length(donors) != EXPECTED_N_DONORS) {
  stop(
    paste0(
      "
QC STOP
",
      "Expected exactly ",
      EXPECTED_N_DONORS,
      " donors, but found ",
      length(donors),
      ": ",
      paste(donors, collapse = ", "),
      "
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

cat(
  "Donors:",
  paste(donors, collapse = ", "),
  "
"
)

cat(
  "Profiles per donor:
"
)

print(
  table(donor_ids)
)


# ------------------------------------------------------------
# 5. Recover frozen 81 / 99 program lists
# ------------------------------------------------------------
all_csv <- list.files(
  project_root,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

native_candidates <- all_csv[
  grepl(
    "native.*81|81.*native",
    basename(all_csv),
    ignore.case = TRUE
  )
]

alt_candidates <- all_csv[
  grepl(
    "alternative.*99|99.*alternative",
    basename(all_csv),
    ignore.case = TRUE
  )
]

read_gene_list <- function(files, expected_n) {
  if (length(files) == 0) {
    return(character(0))
  }

  for (f in files) {
    tab <- safe_read_csv(f)

    if (is.null(tab) || nrow(tab) == 0) {
      next
    }

    gcol <- pick_col(
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

    if (is.na(gcol)) {
      next
    }

    genes <- unique(
      clean_gene(
        tab[[gcol]]
      )
    )

    genes <- genes[
      !is.na(genes)
    ]

    if (length(genes) == expected_n) {
      cat(
        "\nRecovered",
        expected_n,
        "gene programmeme from:\n",
        f,
        "\n"
      )
      return(genes)
    }
  }

  character(0)
}

native81 <- read_gene_list(
  native_candidates,
  EXPECTED_NATIVE_N
)

alternative99 <- read_gene_list(
  alt_candidates,
  EXPECTED_ALT_N
)

# Fallback: search every CSV for exact-size gene lists
if (
  length(native81) != EXPECTED_NATIVE_N ||
  length(alternative99) != EXPECTED_ALT_N
) {
  exact_sets <- list()

  for (f in all_csv) {
    tab <- safe_read_csv(f)

    if (is.null(tab) || nrow(tab) == 0) {
      next
    }

    gcol <- pick_col(
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

    if (is.na(gcol)) {
      next
    }

    genes <- unique(
      clean_gene(
        tab[[gcol]]
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
      exact_sets[[length(exact_sets) + 1]] <-
        tibble::tibble(
          file = f,
          n = length(genes)
        )
    }
  }

  if (length(exact_sets) > 0) {
    exact_df <- dplyr::bind_rows(
      exact_sets
    ) %>%
      dplyr::distinct()

    if (length(native81) != EXPECTED_NATIVE_N) {
      cand <- exact_df %>%
        dplyr::filter(
          n == EXPECTED_NATIVE_N
        ) %>%
        dplyr::mutate(
          score =
            10 * grepl(
              "Fig4|FIGURE.?4|native|identity",
              file,
              ignore.case = TRUE
            ) +
            3 * grepl(
              "gene|program",
              file,
              ignore.case = TRUE
            )
        ) %>%
        dplyr::arrange(
          dplyr::desc(score)
        )

      if (nrow(cand) > 0) {
        native81 <- read_gene_list(
          cand$file,
          EXPECTED_NATIVE_N
        )
      }
    }

    if (length(alternative99) != EXPECTED_ALT_N) {
      cand <- exact_df %>%
        dplyr::filter(
          n == EXPECTED_ALT_N
        ) %>%
        dplyr::mutate(
          score =
            10 * grepl(
              "Fig4|FIGURE.?4|alternative|state",
              file,
              ignore.case = TRUE
            ) +
            3 * grepl(
              "gene|program",
              file,
              ignore.case = TRUE
            )
        ) %>%
        dplyr::arrange(
          dplyr::desc(score)
        )

      if (nrow(cand) > 0) {
        alternative99 <- read_gene_list(
          cand$file,
          EXPECTED_ALT_N
        )
      }
    }
  }
}

if (
  length(native81) != EXPECTED_NATIVE_N ||
  length(alternative99) != EXPECTED_ALT_N
) {
  stop(
    paste0(
      "\nQC STOP\n",
      "Could not safely recover the frozen 81 / 99 programme genes.\n",
      "Figure 5C will NOT redefine them.\n"
    )
  )
}

if (
  !all(native81 %in% rownames(expr)) ||
  !all(alternative99 %in% rownames(expr))
) {
  stop(
    "\nQC STOP\nProgram genes are not fully present in subexpr.\n"
  )
}

# ------------------------------------------------------------
# 6. Recover the frozen prespecified Figure 4 target
# ------------------------------------------------------------
# Reuse an already-exported per-sample Figure 4 table.
# Do not silently redefine the target from expression data.

find_target_table <- function(csv_files, n_expected) {
  candidates <- list()

  for (f in csv_files) {
    if (
      !grepl(
        "Fig4|FIGURE.?4|Panel.?A|Panel.?C|sample|trajectory|score|balance",
        f,
        ignore.case = TRUE
      )
    ) {
      next
    }

    tab <- safe_read_csv(f)

    if (
      is.null(tab) ||
      nrow(tab) != n_expected
    ) {
      next
    }

    nm <- names(tab)

    native_col <- nm[
      grepl(
        "native.*score|score.*native",
        nm,
        ignore.case = TRUE
      )
    ]

    alt_col <- nm[
      grepl(
        "alternative.*score|score.*alternative",
        nm,
        ignore.case = TRUE
      )
    ]

    balance_col <- nm[
      grepl(
        "balance|restoration.*axis|axis.*score",
        nm,
        ignore.case = TRUE
      )
    ]

    sample_col <- nm[
      grepl(
        "^sample$|sample.*id|profile.*id",
        nm,
        ignore.case = TRUE
      )
    ]

    score <- 0

    if (length(balance_col) > 0) {
      score <- score + 10
    }

    if (
      length(native_col) > 0 &&
      length(alt_col) > 0
    ) {
      score <- score + 8
    }

    if (length(sample_col) > 0) {
      score <- score + 4
    }

    if (score > 0) {
      candidates[[length(candidates) + 1]] <-
        list(
          file = f,
          tab = tab,
          native_col = if (length(native_col) > 0) native_col[1] else NA_character_,
          alt_col = if (length(alt_col) > 0) alt_col[1] else NA_character_,
          balance_col = if (length(balance_col) > 0) balance_col[1] else NA_character_,
          sample_col = if (length(sample_col) > 0) sample_col[1] else NA_character_,
          score = score
        )
    }
  }

  if (length(candidates) == 0) {
    return(NULL)
  }

  scores <- vapply(
    candidates,
    function(x) x$score,
    numeric(1)
  )

  candidates[[which.max(scores)]]
}

target_hit <- find_target_table(
  all_csv,
  EXPECTED_N_SAMPLES
)

target_source <- NA_character_

if (!is.null(target_hit)) {
  tt <- target_hit$tab

  # If the table contains sample IDs, explicitly align to expr columns.
  if (!is.na(target_hit$sample_col)) {
    tt_sample <- as.character(
      tt[[target_hit$sample_col]]
    )

    if (setequal(tt_sample, sample_ids)) {
      tt <- tt[
        match(sample_ids, tt_sample),
        ,
        drop = FALSE
      ]
    }
  }

  if (!is.na(target_hit$balance_col)) {
    y <- as.numeric(
      tt[[target_hit$balance_col]]
    )

    target_source <- paste0(
      "Frozen Figure 4 table: ",
      target_hit$file,
      " :: ",
      target_hit$balance_col
    )
  } else {
    y <- as.numeric(
      tt[[target_hit$native_col]]
    ) -
      as.numeric(
        tt[[target_hit$alt_col]]
      )

    target_source <- paste0(
      "Frozen Figure 4 table reconstructed from ",
      target_hit$native_col,
      " - ",
      target_hit$alt_col,
      ": ",
      target_hit$file
    )
  }

  cat(
    "\nTarget recovered from existing Figure 4 output:\n",
    target_source,
    "\n"
  )
} else {
  stop(
    paste0(
      "\nQC STOP\n",
      "No reusable frozen Figure 4 sample-score/state-balance table was found.\n",
      "Figure 5C will NOT reconstruct the target from expression data, because\n",
      "the prespecified Figure 4 target must be reused exactly for this analysis.\n"
    )
  )
}

if (
  length(y) != EXPECTED_N_SAMPLES ||
  any(!is.finite(y))
) {
  stop(
    "\nQC STOP\nInvalid prespecified anatomy-referenced state-balance target vector.\n"
  )
}

# ------------------------------------------------------------
# 7. Define non-programme predictor genes
# ------------------------------------------------------------
program_genes <- unique(
  c(
    native81,
    alternative99
  )
)

heldout_genes <- setdiff(
  rownames(expr),
  program_genes
)

audit_matches <- list.files(
  project_root,
  pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(audit_matches) > 0) {
  audit <- safe_read_csv(
    audit_matches[1]
  )

  if (!is.null(audit)) {
    gcol <- pick_col(
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

    if (!is.na(gcol)) {
      background_genes <- unique(
        clean_gene(
          audit[[gcol]]
        )
      )

      background_genes <- background_genes[
        !is.na(background_genes)
      ]

      background_in_expr <- intersect(
        background_genes,
        rownames(expr)
      )

      heldout_genes <- setdiff(
        background_in_expr,
        program_genes
      )
    }
  }
}

cat(
  "\nNon-programme predictor genes:",
  length(heldout_genes),
  "\n"
)

if (
  length(heldout_genes) != EXPECTED_HELDOUT_N
) {
  stop(
    paste0(
      "\nQC STOP\n",
      "Expected exactly ",
      EXPECTED_HELDOUT_N,
      " held-out non-programme genes, but recovered ",
      length(heldout_genes),
      ".\n",
      "Figure 5C has NOT been generated.\n"
    )
  )
}

leakage_genes <- intersect(
  heldout_genes,
  program_genes
)

if (length(leakage_genes) > 0) {
  stop(
    paste0(
      "\nQC STOP — predictor leakage detected.\n",
      "Program genes found among non-programme predictors:\n",
      paste(leakage_genes, collapse = ", "),
      "\n"
    )
  )
}

# Samples x genes
x <- t(
  expr[
    heldout_genes,
    ,
    drop = FALSE
  ]
)

if (
  nrow(x) != EXPECTED_N_SAMPLES ||
  ncol(x) != EXPECTED_HELDOUT_N
) {
  stop(
    "\nQC STOP\nUnexpected predictor matrix dimensions.\n"
  )
}

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

heldout_genes_model <- colnames(x)

cat(
  "Nonzero-variance non-programme predictors:",
  ncol(x),
  "\n"
)

if (ncol(x) < 1000) {
  stop(
    "\nQC STOP\nToo many non-programme predictors have zero/invalid variance.\n"
  )
}

if (anyNA(x)) {
  stop(
    paste0(
      "\nQC STOP\n",
      "Missing values detected in predictor matrix.\n",
      "Imputation is intentionally not performed silently.\n"
    )
  )
}

# ------------------------------------------------------------
# 8. DONOR-BLOCKED nested-CV elastic net
# ------------------------------------------------------------
pred_store <- matrix(
  NA_real_,
  nrow = EXPECTED_N_SAMPLES,
  ncol = N_REPEATS
)

rownames(pred_store) <- sample_ids
colnames(pred_store) <- paste0(
  "Repeat_",
  seq_len(N_REPEATS)
)

outer_fold_records <- list()
tuning_records <- list()
repeat_records <- list()

fit_inner_alpha <- function(
    x_train,
    y_train,
    donor_train,
    alpha_grid,
    inner_k,
    seed_base) {

  best_alpha <- NA_real_
  best_lambda <- NA_real_
  best_cvm <- Inf

  tuning_rows <- list()

  # CRITICAL: inner folds are built at DONOR level.
  inner_foldid <- make_group_folds(
    groups = donor_train,
    k = inner_k,
    seed = seed_base
  )

  # Hard QC: no donor can span multiple inner folds.
  donor_fold_check <- tapply(
    inner_foldid,
    donor_train,
    function(z) length(unique(z))
  )

  if (any(donor_fold_check != 1)) {
    stop(
      "\nQC STOP\nA donor was split across inner-CV folds.\n"
    )
  }

  for (a in alpha_grid) {
    set.seed(
      seed_base +
        round(a * 1000) +
        17
    )

    cvfit <- glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      family = "gaussian",
      alpha = a,
      foldid = inner_foldid,
      type.measure = "mse",
      standardize = TRUE,
      intercept = TRUE,
      nlambda = 100
    )

    i_min <- which.min(
      cvfit$cvm
    )

    this_cvm <- cvfit$cvm[
      i_min
    ]

    this_lambda <- cvfit$lambda[
      i_min
    ]

    tuning_rows[[length(tuning_rows) + 1]] <-
      tibble::tibble(
        alpha = a,
        lambda_min = this_lambda,
        inner_CV_MSE = this_cvm
      )

    if (
      is.finite(this_cvm) &&
      this_cvm < best_cvm
    ) {
      best_cvm <- this_cvm
      best_alpha <- a
      best_lambda <- this_lambda
    }
  }

  if (
    !is.finite(best_alpha) ||
    !is.finite(best_lambda)
  ) {
    stop(
      "\nInner CV failed to choose alpha/lambda.\n"
    )
  }

  list(
    alpha = best_alpha,
    lambda = best_lambda,
    mse = best_cvm,
    tuning = dplyr::bind_rows(
      tuning_rows
    )
  )
}

cat(
  "\nStarting DONOR-BLOCKED nested cross-validation...\n",
  "Outer CV: leave-one-donor-out across ",
  EXPECTED_N_DONORS,
  " donors x ",
  N_REPEATS,
  " repeats\n",
  "Inner CV: ",
  INNER_K,
  "-fold donor-grouped CV\n",
  "Alpha grid: ",
  paste(ALPHA_GRID, collapse = ", "),
  "\n\n",
  sep = ""
)

for (r in seq_len(N_REPEATS)) {

  # Every donor is held out exactly once per repeat.
  # Randomizing the order changes only the repeat bookkeeping;
  # the substantive outer split is leave-one-donor-out.
  set.seed(
    MASTER_SEED +
      r * 100
  )

  outer_donor_order <- sample(
    donors,
    length(donors),
    replace = FALSE
  )

  pred_r <- rep(
    NA_real_,
    EXPECTED_N_SAMPLES
  )

  for (fold in seq_along(outer_donor_order)) {
    test_donor <- outer_donor_order[fold]

    test_idx <- which(
      donor_ids == test_donor
    )

    train_idx <- which(
      donor_ids != test_donor
    )

    donor_train <- donor_ids[
      train_idx
    ]

    # Hard QC: held-out donor must be completely absent from training.
    if (
      any(donor_ids[test_idx] %in% donor_train)
    ) {
      stop(
        "\nQC STOP\nDonor leakage detected in outer CV.\n"
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

    x_test <- x[
      test_idx,
      ,
      drop = FALSE
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

    x_test2 <- x_test[
      ,
      keep_train,
      drop = FALSE
    ]

    if (ncol(x_train2) < 100) {
      stop(
        paste0(
          "\nQC STOP\n",
          "Too few nonzero-variance predictors in repeat ",
          r,
          ", held-out donor ",
          test_donor,
          ".\n"
        )
      )
    }

    tuned <- fit_inner_alpha(
      x_train = x_train2,
      y_train = y_train,
      donor_train = donor_train,
      alpha_grid = ALPHA_GRID,
      inner_k = INNER_K,
      seed_base =
        MASTER_SEED +
        r * 10000 +
        fold * 100
    )

    final_fit <- glmnet::glmnet(
      x = x_train2,
      y = y_train,
      family = "gaussian",
      alpha = tuned$alpha,
      lambda = tuned$lambda,
      standardize = TRUE,
      intercept = TRUE
    )

    pred_fold <- as.numeric(
      predict(
        final_fit,
        newx = x_test2,
        s = tuned$lambda
      )
    )

    pred_r[
      test_idx
    ] <- pred_fold

    fold_metric <- calc_metrics(
      y[test_idx],
      pred_fold
    ) %>%
      dplyr::mutate(
        Repeat = r,
        Fold = fold,
        Test_donor = test_donor,
        Train_donors = paste(
          sort(unique(donor_train)),
          collapse = ","
        ),
        N_train = length(train_idx),
        N_test = length(test_idx),
        N_predictors = ncol(x_train2),
        alpha = tuned$alpha,
        lambda = tuned$lambda,
        inner_CV_MSE = tuned$mse
      )

    outer_fold_records[[
      length(outer_fold_records) + 1
    ]] <- fold_metric

    tune_tab <- tuned$tuning %>%
      dplyr::mutate(
        Repeat = r,
        Fold = fold,
        Test_donor = test_donor,
        selected =
          alpha == tuned$alpha &
          abs(
            lambda_min -
              tuned$lambda
          ) < 1e-15
      )

    tuning_records[[
      length(tuning_records) + 1
    ]] <- tune_tab
  }

  if (any(!is.finite(pred_r))) {
    stop(
      paste0(
        "\nQC STOP\nMissing OOF predictions in repeat ",
        r,
        ".\n"
      )
    )
  }

  pred_store[
    ,
    r
  ] <- pred_r

  rep_metric <- calc_metrics(
    y,
    pred_r
  ) %>%
    dplyr::mutate(
      Repeat = r
    )

  repeat_records[[
    length(repeat_records) + 1
  ]] <- rep_metric

  cat(
    "Repeat",
    r,
    "/",
    N_REPEATS,
    " | Pearson r =",
    sprintf(
      "%.3f",
      rep_metric$Pearson_r
    ),
    " | Spearman rho =",
    sprintf(
      "%.3f",
      rep_metric$Spearman_rho
    ),
    " | CV R2 =",
    sprintf(
      "%.3f",
      rep_metric$CV_R2
    ),
    "\n"
  )
}

repeat_metrics <- dplyr::bind_rows(
  repeat_records
)

outer_fold_metrics <- dplyr::bind_rows(
  outer_fold_records
)

model_tuning <- dplyr::bind_rows(
  tuning_records
)

# ------------------------------------------------------------
# 9. Aggregate repeated OOF predictions
# ------------------------------------------------------------
pred_median <- apply(
  pred_store,
  1,
  median
)

pred_mean <- rowMeans(
  pred_store
)

pred_sd <- apply(
  pred_store,
  1,
  sd
)

aggregate_metrics <- calc_metrics(
  y,
  pred_median
)

prediction_df <- tibble::tibble(
  Sample = sample_ids,
  Donor = donor_ids,
  Observed_balance = y,
  Predicted_balance_median = pred_median,
  Predicted_balance_mean = pred_mean,
  Prediction_SD = pred_sd
)

# ------------------------------------------------------------
# 10. DONOR-CLUSTERED bootstrap CI
# ------------------------------------------------------------
# Resample donors, not individual profiles, so within-donor
# dependence is preserved.

set.seed(MASTER_SEED + 999)

B <- 5000

boot_metrics <- vector(
  "list",
  B
)

for (b in seq_len(B)) {

  sampled_donors <- sample(
    donors,
    size = length(donors),
    replace = TRUE
  )

  idx <- unlist(
    lapply(
      sampled_donors,
      function(d) which(donor_ids == d)
    ),
    use.names = FALSE
  )

  boot_metrics[[b]] <- calc_metrics(
    y[idx],
    pred_median[idx]
  )
}

boot_df <- dplyr::bind_rows(
  boot_metrics
)

ci_summary <- boot_df %>%
  dplyr::summarise(
    Pearson_low = quantile(
      Pearson_r,
      0.025,
      na.rm = TRUE
    ),
    Pearson_high = quantile(
      Pearson_r,
      0.975,
      na.rm = TRUE
    ),
    Spearman_low = quantile(
      Spearman_rho,
      0.025,
      na.rm = TRUE
    ),
    Spearman_high = quantile(
      Spearman_rho,
      0.975,
      na.rm = TRUE
    ),
    R2_low = quantile(
      CV_R2,
      0.025,
      na.rm = TRUE
    ),
    R2_high = quantile(
      CV_R2,
      0.975,
      na.rm = TRUE
    )
  )

# ------------------------------------------------------------
# 11. Plot — retain v7 visual style
# ------------------------------------------------------------
metric_text <- paste0(
  "Pearson r = ",
  sprintf("%.2f", aggregate_metrics$Pearson_r),
  " [95% donor-bootstrap CI ",
  sprintf("%.2f", ci_summary$Pearson_low),
  ", ",
  sprintf("%.2f", ci_summary$Pearson_high),
  "]\n",
  "Spearman rho = ",
  sprintf("%.2f", aggregate_metrics$Spearman_rho),
  " | CV R2 = ",
  sprintf("%.2f", aggregate_metrics$CV_R2),
  "\n",
  "Donor-blocked nested CV | MAE = ",
  sprintf("%.2f", aggregate_metrics$MAE)
)

lims <- range(
  c(
    prediction_df$Observed_balance,
    prediction_df$Predicted_balance_median
  ),
  finite = TRUE
)

pad <- diff(lims) * 0.055
if (!is.finite(pad) || pad == 0) pad <- 0.1
lims <- c(lims[1] - pad, lims[2] + pad)

p <- ggplot(
  prediction_df,
  aes(
    x = Observed_balance,
    y = Predicted_balance_median
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.32,
    colour = "#A7B2B7"
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.42,
    colour = "#456B78",
    fill = "#DDE7EA",
    alpha = 0.55
  ) +
  geom_point(
    size = 1.35,
    alpha = 0.82,
    colour = "#355F6B"
  ) +
  annotate(
    "text",
    x = lims[1] + 0.075 * diff(lims),
    y = lims[2] - 0.095 * diff(lims),
    label = metric_text,
    hjust = 0,
    vjust = 1,
    size = 1.46,
    family = "Helvetica",
    colour = "#20343D",
    lineheight = 1.24
  ) +
  coord_cartesian(
    xlim = lims,
    ylim = lims,
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    x = "Observed anatomy-referenced state balance",
    y = "Donor-blocked cross-validated reconstructed balance",
    title = "Donor-blocked reconstruction of the prespecified state balance",
    subtitle = "1,373 non-programme genes; donor-blocked model fitting"
  ) +
  theme_classic(
    base_family = "Helvetica",
    base_size = 5.0
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      size = 5.7,
      face = "bold",
      colour = "#20343D",
      hjust = 0.5,
      lineheight = 1.00,
      margin = margin(b = 1.0, unit = "pt")
    ),
    plot.subtitle = element_text(
      size = 4.25,
      colour = "#647B85",
      hjust = 0.5,
      margin = margin(b = 1.7, unit = "pt")
    ),
    axis.title.x = element_text(
      size = 5.0,
      colour = "#20343D",
      margin = margin(t = 1.5, unit = "pt")
    ),
    axis.title.y = element_text(
      size = 5.0,
      colour = "#20343D",
      margin = margin(r = 1.5, unit = "pt")
    ),
    axis.text = element_text(
      size = 4.45,
      colour = "#667D87"
    ),
    axis.line = element_line(
      linewidth = 0.32,
      colour = "#20343D"
    ),
    axis.ticks = element_line(
      linewidth = 0.28,
      colour = "#20343D"
    ),
    axis.ticks.length = grid::unit(1.15, "pt"),
    plot.margin = margin(
      t = 1.2,
      r = 3.2,
      b = 1.8,
      l = 3.0,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# 12. Output directory
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 5",
  "Panel_C_stateBalanceReconstruction_v12"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 13. Export figure
# ------------------------------------------------------------
FIG_W_MM <- 94
FIG_H_MM <- 56

ggsave(
  filename = file.path(
    out_dir,
    "Fig5_PanelC_stateBalanceReconstruction_v12.png"
  ),
  plot = p,
  width = FIG_W_MM,
  height = FIG_H_MM,
  units = "mm",
  dpi = 600,
  bg = "white"
)

pdf_file <- file.path(
  out_dir,
  "Fig5_PanelC_stateBalanceReconstruction_v12.pdf"
)

if (file.exists(pdf_file)) {
  file.remove(pdf_file)
}

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

if (
  !file.exists(pdf_file) ||
  is.na(file.info(pdf_file)$size) ||
  file.info(pdf_file)$size <= 0
) {
  stop(
    paste0(
      "PDF export failed: ",
      pdf_file
    )
  )
}

ggsave(
  filename = file.path(
    out_dir,
    "Fig5_PanelC_stateBalanceReconstruction_v12.tiff"
  ),
  plot = p,
  width = FIG_W_MM,
  height = FIG_H_MM,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ------------------------------------------------------------
# 14. Export data
# ------------------------------------------------------------
readr::write_csv(
  prediction_df,
  file.path(
    out_dir,
    "Fig5C_sample_predictions.csv"
  )
)

pred_matrix_df <- as.data.frame(
  pred_store
) %>%
  tibble::rownames_to_column(
    "Sample"
  )

readr::write_csv(
  pred_matrix_df,
  file.path(
    out_dir,
    "Fig5C_all_repeat_OOF_predictions.csv"
  )
)

readr::write_csv(
  repeat_metrics,
  file.path(
    out_dir,
    "Fig5C_repeat_metrics.csv"
  )
)

readr::write_csv(
  outer_fold_metrics,
  file.path(
    out_dir,
    "Fig5C_outer_fold_metrics.csv"
  )
)

readr::write_csv(
  model_tuning,
  file.path(
    out_dir,
    "Fig5C_model_tuning.csv"
  )
)

readr::write_csv(
  tibble::tibble(
    GeneSymbol = heldout_genes_model
  ),
  file.path(
    out_dir,
    "Fig5C_nonprogramme_predictor_genes.csv"
  )
)

readr::write_csv(
  boot_df,
  file.path(
    out_dir,
    "Fig5C_bootstrap_metric_CI_samples.csv"
  )
)

# Also save donor mapping for auditability.
readr::write_csv(
  tibble::tibble(
    Sample = sample_ids,
    Donor = donor_ids
  ),
  file.path(
    out_dir,
    "Fig5C_sample_donor_mapping.csv"
  )
)

# ------------------------------------------------------------
# 15. QC summary
# ------------------------------------------------------------
qc <- tibble::tibble(
  Metric = c(
    "Samples",
    "Donors",
    "Frozen native programme genes",
    "Frozen alternative programme genes",
    "Expected non-programme genes",
    "Recovered non-programme genes",
    "Nonzero-variance non-programme predictors",
    "Program genes leaking into predictors",
    "Outer CV design",
    "CV repeats",
    "Inner CV folds",
    "Inner CV grouping",
    "Bootstrap unit",
    "Aggregate Pearson r",
    "Aggregate Spearman rho",
    "Aggregate CV R2",
    "Aggregate MAE",
    "Aggregate RMSE"
  ),
  Value = c(
    EXPECTED_N_SAMPLES,
    length(donors),
    length(native81),
    length(alternative99),
    EXPECTED_HELDOUT_N,
    length(heldout_genes),
    ncol(x),
    length(leakage_genes),
    "leave-one-donor-out",
    N_REPEATS,
    INNER_K,
    "donor-grouped",
    "donor cluster",
    aggregate_metrics$Pearson_r,
    aggregate_metrics$Spearman_rho,
    aggregate_metrics$CV_R2,
    aggregate_metrics$MAE,
    aggregate_metrics$RMSE
  )
)

readr::write_csv(
  qc,
  file.path(
    out_dir,
    "Fig5C_QC_summary.csv"
  )
)

# ------------------------------------------------------------
# 16. README
# ------------------------------------------------------------
readme <- c(
  "Figure 5C v12 — donor-blocked conditional reconstruction",
  "",
  "Non-programme transcriptomic variation reconstructs the prespecified anatomy-referenced state balance.",
  "",
  "Critical change relative to the previous analysis:",
  "- All profiles from the same donor are kept together.",
  "- Outer model evaluation uses leave-one-donor-out cross-validation across 10 donors.",
  "- Inner hyperparameter tuning also uses donor-grouped folds.",
  "- Bootstrap confidence intervals resample donors rather than individual profiles.",
  "",
  "Question:",
  paste0(
    "Can the ",
    EXPECTED_HELDOUT_N,
    " non-programme genes predict the Figure 4 Native - Alternative balance?"
  ),
  "",
  "Predictor leakage control:",
  paste0(
    "- Native programme genes excluded from predictors: n = ",
    length(native81)
  ),
  paste0(
    "- Alternative programme genes excluded from predictors: n = ",
    length(alternative99)
  ),
  paste0(
    "- Program genes detected among predictors: n = ",
    length(leakage_genes)
  ),
  "",
  "Reconstruction design:",
  paste0(
    "- Profiles: n = ",
    EXPECTED_N_SAMPLES
  ),
  paste0(
    "- Donors: n = ",
    length(donors),
    " (",
    paste(donors, collapse = ", "),
    ")"
  ),
  paste0(
    "- Non-programme predictors: n = ",
    length(heldout_genes)
  ),
  paste0(
    "- Outer CV: leave-one-donor-out x ",
    N_REPEATS,
    " repeats"
  ),
  paste0(
    "- Inner CV: ",
    INNER_K,
    "-fold donor-grouped"
  ),
  paste0(
    "- Elastic-net alpha grid: ",
    paste(
      ALPHA_GRID,
      collapse = ", "
    )
  ),
  "- Lambda selected inside each outer training set.",
  "- Predictor standardization performed by glmnet within training fits.",
  "- No supervised genome-wide feature selection is performed before cross-validation.",
  "",
  "Target source:",
  target_source,
  "",
  "Primary aggregate metrics:",
  paste0(
    "- Pearson r = ",
    sprintf(
      "%.4f",
      aggregate_metrics$Pearson_r
    )
  ),
  paste0(
    "- Spearman rho = ",
    sprintf(
      "%.4f",
      aggregate_metrics$Spearman_rho
    )
  ),
  paste0(
    "- CV R2 = ",
    sprintf(
      "%.4f",
      aggregate_metrics$CV_R2
    )
  ),
  paste0(
    "- MAE = ",
    sprintf(
      "%.4f",
      aggregate_metrics$MAE
    )
  ),
  "",
  "Interpretation:",
  "This is an internal conditional reconstruction analysis within GSE272333.",
  "The donor-blocked design prevents profiles from the same donor",
  "from appearing in both outer training and test data.",
  "",
  "Main plotted reconstruction:",
  "Median out-of-fold reconstruction for each profile across repeated donor-blocked runs."
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig5_PanelC_v12_stateBalanceReconstruction.txt"
  )
)

# ------------------------------------------------------------
# 17. Console summary
# ------------------------------------------------------------
cat(
  "\n============================================\n"
)

cat(
  "Figure 5C donor-blocked conditional reconstruction completed successfully.\n"
)

cat(
  "============================================\n\n"
)

print(qc)

cat(
  "\n95% donor-cluster bootstrap CI:\n"
)

print(ci_summary)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)

cat(
  "\nMain figure:\n",
  pdf_file,
  "\n"
)
