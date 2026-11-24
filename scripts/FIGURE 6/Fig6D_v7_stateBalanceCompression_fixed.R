# ============================================================
# Figure 6D v7 — DONOR-BLOCKED, DESCRIPTIVE REPEAT COMPARISON
# Nested top-25 reconstruction performance under donor-blocked selection
#
# SCIENTIFIC FIX:
# Previous versions defined the top-25 marginal genes using all
# 110 samples before cross-validation. Because the outcome was
# used in that global feature selection, the resulting CV estimate
# was optimistic.
#
# v4 removes that leakage:
#
# For EACH outer training fold:
#   1) Use TRAINING SAMPLES ONLY.
#   2) Compute Spearman rho between each of the 1,373 non-programme genes
#      and the anatomy-referenced state balance.
#   3) Select the 25 genes with largest |rho| in that training fold.
#   4) Tune elastic net using inner CV inside that same outer training set.
#   5) Predict the untouched outer test fold.
#
# Comparator:
#   Full 1,373-gene model = frozen Figure 5C repeated nested-CV
#   performance, using the same 20 repeats and 10 outer folds.
#
# Main question:
# Can training-fold-selected top-25 genes retain most reconstruction
# performance when feature selection is confined to outer training folds?
#
# Outputs:
#   Fig6_PanelD_nestedTop25_v6.pdf/png/tiff
#   Fig6D_v6_repeat_metrics.csv
#   Fig6D_v6_aggregate_metrics.csv
#   Fig6D_v6_top25_selection_frequency.csv
#   Fig6D_v6_outer_fold_metrics.csv
#   Fig6D_QC_summary_v6.csv
# ============================================================

required_pkgs <- c(
  "glmnet",
  "ggplot2",
  "readr",
  "dplyr",
  "tibble",
  "stringr",
  "scales",
  "tidyr"
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
  library(tidyr)
})

EXPECTED_N_SAMPLES <- 110
EXPECTED_HELDOUT_N <- 1373

N_REPEATS <- 20
OUTER_K <- 10  # 10 held-out donors
INNER_K <- 5
TOP_N <- 25

ALPHA_GRID <- c(
  0,
  0.25,
  0.5,
  0.75,
  1
)

# ------------------------------------------------------------
# 1. Locate project root
# ------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)

  file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
  )

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

  if (
    length(
      list.files(
        x,
        pattern = "^Fig5C_repeat_metrics\\.csv$",
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
        pattern = "^Fig6B_all_nonprogramme_marginal_vs_multivariate\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    ) > 0
  ) {
    s <- s + 10
  }

  if (dir.exists(file.path(x, "figures"))) {
    s <- s + 2
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

# ------------------------------------------------------------
# 2. Helpers
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

make_group_folds <- function(groups, k, seed) {
  groups <- as.character(groups)
  if (anyNA(groups) || any(groups == "")) {
    stop("
QC STOP
Missing donor labels in grouped CV.
")
  }

  ug <- unique(groups)
  if (length(ug) < k) {
    stop("
QC STOP
Too few donors for grouped inner CV.
")
  }

  set.seed(seed)
  shuffled <- sample(ug, length(ug), replace = FALSE)
  donor_fold <- rep(seq_len(k), length.out = length(shuffled))
  names(donor_fold) <- shuffled
  as.integer(donor_fold[groups])
}

calc_metrics <- function(obs, pred) {
  ok <- is.finite(obs) &
    is.finite(pred)

  obs <- obs[ok]
  pred <- pred[ok]

  pearson_r <- suppressWarnings(
    cor(
      obs,
      pred,
      method = "pearson"
    )
  )

  spearman_rho <- suppressWarnings(
    cor(
      obs,
      pred,
      method = "spearman"
    )
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
    1 -
      sse /
      sst
  } else {
    NA_real_
  }

  mae <- mean(
    abs(
      obs -
        pred
    )
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

  inner_foldid <- make_group_folds(
    donor_train,
    inner_k,
    seed_base
  )

  donor_fold_check <- tapply(
    inner_foldid,
    donor_train,
    function(z) length(unique(z))
  )

  if (any(donor_fold_check != 1)) {
    stop("\nQC STOP\nA donor was split across inner folds.\n")
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

    if (
      is.finite(this_cvm) &&
      this_cvm <
        best_cvm
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
    mse = best_cvm
  )
}

# ------------------------------------------------------------
# 3. Load clean-session subexpr
# ------------------------------------------------------------
loader_file <- file.path(
  project_root,
  "scripts",
  "load_subexpr_clean.R"
)

if (!file.exists(loader_file)) {
  stop(
    paste0(
      "\nQC STOP\nCannot find clean-session loader:\n",
      loader_file,
      "\nRun/restore scripts/load_subexpr_clean.R first.\n"
    )
  )
}

source(loader_file, local = FALSE)

expr <- as.matrix(subexpr)
storage.mode(expr) <- "numeric"

if (nrow(expr) != 1553L || ncol(expr) != EXPECTED_N_SAMPLES) {
  stop(
    paste0(
      "\nQC STOP\nExpected subexpr = 1553 x ",
      EXPECTED_N_SAMPLES,
      "; found ",
      nrow(expr),
      " x ",
      ncol(expr),
      ".\n"
    )
  )
}

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("\nQC STOP\nsubexpr must contain gene row names and sample column names.\n")
}

rownames(expr) <- clean_gene(rownames(expr))

if (anyDuplicated(rownames(expr))) {
  stop("\nQC STOP\nDuplicated gene symbols in subexpr.\n")
}

cat(
  "\nLoaded clean-session subexpr: ",
  nrow(expr),
  " x ",
  ncol(expr),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Recover all 1,373 non-programme genes
# ------------------------------------------------------------
joint_files <- list.files(
  project_root,
  pattern = "^Fig6B_all_nonprogramme_marginal_vs_multivariate\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(joint_files) == 0) {
  stop(
    "\nQC STOP\nCannot find Fig6B_all_nonprogramme_marginal_vs_multivariate.csv\n"
  )
}

joint <- safe_read_csv(
  joint_files[1]
)

if (
  is.null(joint) ||
  !"GeneSymbol" %in%
    names(joint)
) {
  stop(
    "\nQC STOP\nInvalid Figure 6B non-programme table.\n"
  )
}

heldout_genes <- unique(
  clean_gene(
    joint$GeneSymbol
  )
)

heldout_genes <- heldout_genes[
  !is.na(heldout_genes)
]

if (
  length(heldout_genes) !=
    EXPECTED_HELDOUT_N
) {
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

if (
  length(
    setdiff(
      heldout_genes,
      rownames(expr)
    )
  ) > 0
) {
  stop(
    "\nQC STOP\nSome non-programme genes are missing from subexpr.\n"
  )
}

x_all <- t(
  expr[
    heldout_genes,
    ,
    drop = FALSE
  ]
)

if (anyNA(x_all)) {
  stop(
    "\nQC STOP\nMissing values in held-out expression matrix.\n"
  )
}

# ------------------------------------------------------------
# 5. Recover frozen target
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

  if (
    grepl(
      "stateBalanceReconstruction_v12|Panel_C_stateBalanceReconstruction_v12",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 100
  }

  if (
    grepl(
      "Panel_C_v3",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 30
  }

  if (
    grepl(
      "Panel_C_v2",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 20
  }

  if (
    grepl(
      "FIGURE 5",
      path,
      ignore.case = TRUE
    )
  ) {
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

pred_tab <- safe_read_csv(
  pred_file
)

if (
  is.null(pred_tab) ||
  !"Observed_balance" %in%
    names(pred_tab) ||
  nrow(pred_tab) !=
    EXPECTED_N_SAMPLES
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
donors <- sort(unique(donor_ids))

if (anyNA(donor_ids) || length(donors) != 10) {
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

  sample_ids <- expr_samples
  alignment_mode <-
    "Exact sample-name alignment"
} else {
  sample_ids <- if (
    !is.null(expr_samples)
  ) {
    expr_samples
  } else {
    paste0(
      "Sample_",
      seq_len(
        EXPECTED_N_SAMPLES
      )
    )
  }

  alignment_mode <-
    "Figure 5C row order retained"
}

if (
  length(y) !=
    EXPECTED_N_SAMPLES ||
  any(!is.finite(y))
) {
  stop(
    "\nQC STOP\nInvalid target vector.\n"
  )
}

# ------------------------------------------------------------
# 6. Recover frozen Figure 5C repeat metrics
# ------------------------------------------------------------
full_repeat_files <- list.files(
  project_root,
  pattern = "^Fig5C_repeat_metrics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(full_repeat_files) == 0) {
  stop(
    "\nQC STOP\nCannot find Fig5C_repeat_metrics.csv\n"
  )
}

score_full_file <- function(path) {
  s <- 0

  if (
    grepl(
      "stateBalanceReconstruction_v12|Panel_C_stateBalanceReconstruction_v12",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 100
  }

  if (
    grepl(
      "Panel_C_v3",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 30
  }

  if (
    grepl(
      "Panel_C_v2",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 20
  }

  if (
    grepl(
      "FIGURE 5",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 10
  }

  s
}

full_repeat_scores <- vapply(
  full_repeat_files,
  score_full_file,
  numeric(1)
)

full_repeat_file <- full_repeat_files[
  which.max(
    full_repeat_scores
  )
]

full_repeat <- safe_read_csv(
  full_repeat_file
)

if (
  is.null(full_repeat) ||
  !all(
    c(
      "Repeat",
      "Pearson_r",
      "Spearman_rho",
      "CV_R2",
      "MAE",
      "RMSE"
    ) %in%
      names(full_repeat)
  )
) {
  stop(
    "\nQC STOP\nInvalid Figure 5C repeat metrics.\n"
  )
}

if (
  nrow(full_repeat) !=
    N_REPEATS
) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      N_REPEATS,
      " Figure 5C repeats, found ",
      nrow(full_repeat),
      ".\n"
    )
  )
}

full_repeat <- full_repeat %>%
  dplyr::mutate(
    Model =
      "All 1,373 non-programme genes"
  )

# ------------------------------------------------------------
# 7. Nested top-25 outer-CV prediction
# ------------------------------------------------------------
pred_store <- matrix(
  NA_real_,
  nrow = EXPECTED_N_SAMPLES,
  ncol = N_REPEATS
)

rownames(pred_store) <- sample_ids

colnames(pred_store) <- paste0(
  "Repeat_",
  seq_len(
    N_REPEATS
  )
)

repeat_rows <- list()
fold_rows <- list()
selection_rows <- list()

cat(
  "\nStarting leakage-free nested top-25 prediction...\n"
)

for (r in seq_len(N_REPEATS)) {
  set.seed(
    20261028 +
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

  for (fold in seq_len(OUTER_K)) {
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

    if (any(donor_ids[test_idx] %in% donor_train)) {
      stop("\nQC STOP\nDonor leakage detected in outer CV.\n")
    }

    x_train_all <- x_all[
      train_idx,
      ,
      drop = FALSE
    ]

    x_test_all <- x_all[
      test_idx,
      ,
      drop = FALSE
    ]

    y_train <- y[
      train_idx
    ]

    # --------------------------------------------------------
    # LEAKAGE-FREE FEATURE SELECTION:
    # compute marginal rho using OUTER TRAINING data only
    # --------------------------------------------------------
    rho_train <- apply(
      x_train_all,
      2,
      function(z) {
        suppressWarnings(
          cor(
            z,
            y_train,
            method = "spearman",
            use = "complete.obs"
          )
        )
      }
    )

    rho_train[
      !is.finite(
        rho_train
      )
    ] <- 0

    rank_df <- tibble::tibble(
      GeneSymbol =
        names(
          rho_train
        ),
      Train_rho =
        as.numeric(
          rho_train
        ),
      Abs_train_rho =
        abs(
          as.numeric(
            rho_train
        )
      )
    ) %>%
      dplyr::arrange(
        dplyr::desc(
          Abs_train_rho
        ),
        GeneSymbol
      )

    selected_genes <- rank_df %>%
      dplyr::slice_head(
        n = TOP_N
      ) %>%
      dplyr::pull(
        GeneSymbol
      )

    if (
      length(
        selected_genes
      ) !=
        TOP_N
    ) {
      stop(
        paste0(
          "\nQC STOP\nCould not select exactly ",
          TOP_N,
          " genes in repeat ",
          r,
          ", fold ",
          fold,
          ".\n"
        )
      )
    }

    selection_rows <- append(
      selection_rows,
      list(
        rank_df %>%
          dplyr::slice_head(
            n = TOP_N
          ) %>%
          dplyr::mutate(
            Repeat = r,
            Fold = fold,
            Rank =
              dplyr::row_number()
          )
      )
    )

    x_train <- x_train_all[
      ,
      selected_genes,
      drop = FALSE
    ]

    x_test <- x_test_all[
      ,
      selected_genes,
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

    if (
      ncol(x_train2) <
        2
    ) {
      stop(
        "\nQC STOP\nToo few top-25 predictors after variance QC.\n"
      )
    }

    tuned <- fit_inner_alpha(
      x_train =
        x_train2,
      y_train =
        y_train,
      donor_train =
        donor_train,
      alpha_grid =
        ALPHA_GRID,
      inner_k =
        INNER_K,
      seed_base =
        20261028 +
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
      y[
        test_idx
      ],
      pred_fold
    ) %>%
      dplyr::mutate(
        Model =
          "Nested top-25 genes",
        Repeat =
          r,
        Fold =
          fold,
        N_features =
          ncol(
            x_train2
          ),
        alpha =
          tuned$alpha,
        lambda =
          tuned$lambda,
        inner_CV_MSE =
          tuned$mse
      )

    fold_rows <- append(
      fold_rows,
      list(
        fold_metric
      )
    )
  }

  if (
    any(
      !is.finite(
        pred_r
      )
    )
  ) {
    stop(
      paste0(
        "\nQC STOP\nMissing predictions in repeat ",
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
      Model =
        "Nested top-25 genes",
      Repeat =
        r,
      N_features =
        TOP_N
    )

  repeat_rows <- append(
    repeat_rows,
    list(
      rep_metric
    )
  )

  cat(
    "Repeat",
    r,
    "/",
    N_REPEATS,
    " | rho=",
    sprintf(
      "%.3f",
      rep_metric$Spearman_rho
    ),
    " | R2=",
    sprintf(
      "%.3f",
      rep_metric$CV_R2
    ),
    "\n"
  )
}

top25_repeat <- dplyr::bind_rows(
  repeat_rows
)

top25_fold <- dplyr::bind_rows(
  fold_rows
)

selection_all <- dplyr::bind_rows(
  selection_rows
)

# ------------------------------------------------------------
# 8. Aggregate nested top-25 predictions
# ------------------------------------------------------------
pred_median <- apply(
  pred_store,
  1,
  median
)

pred_mean <- rowMeans(
  pred_store
)

top25_aggregate <- calc_metrics(
  y,
  pred_median
) %>%
  dplyr::mutate(
    Model =
      "Nested top-25 genes",
    N_features =
      TOP_N
  )

# Frozen full-model aggregate from Figure 5C sample predictions
full_pred <- as.numeric(
  pred_tab$Predicted_balance_median
)

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
  full_pred <- full_pred[
    match(
      expr_samples,
      as.character(
        pred_tab$Sample
      )
    )
  ]
}

full_aggregate <- calc_metrics(
  y,
  full_pred
) %>%
  dplyr::mutate(
    Model =
      "All 1,373 non-programme genes",
    N_features =
      EXPECTED_HELDOUT_N
  )

aggregate_metrics <- dplyr::bind_rows(
  full_aggregate,
  top25_aggregate
)

# ------------------------------------------------------------
# 9. Repeat-level descriptive comparison
# ------------------------------------------------------------
repeat_compare <- dplyr::bind_rows(
  full_repeat %>%
    dplyr::select(
      Model,
      Repeat,
      Pearson_r,
      Spearman_rho,
      CV_R2,
      MAE,
      RMSE
    ),
  top25_repeat %>%
    dplyr::select(
      Model,
      Repeat,
      Pearson_r,
      Spearman_rho,
      CV_R2,
      MAE,
      RMSE
    )
)

# Repeat-level results are retained for descriptive stability assessment only.
# The 20 repeated CV runs reuse the same donors and profiles and therefore are
# not treated as independent experimental replicates for formal hypothesis testing.

# ------------------------------------------------------------
# 10. Selection stability
# ------------------------------------------------------------
selection_frequency <- selection_all %>%
  dplyr::group_by(
    GeneSymbol
  ) %>%
  dplyr::summarise(
    Selection_count =
      dplyr::n(),
    Selection_frequency =
      dplyr::n() /
      (
        N_REPEATS *
          OUTER_K
      ),
    Mean_abs_train_rho =
      mean(
        Abs_train_rho,
        na.rm = TRUE
      ),
    Median_rank =
      median(
        Rank,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Selection_frequency
    ),
    dplyr::desc(
      Mean_abs_train_rho
    )
  )

# ------------------------------------------------------------
# 11. Output directory
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 6",
  "Panel_D_nestedTop25_donorBlocked_v7_descriptive"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  repeat_compare,
  file.path(
    out_dir,
    "Fig6D_v7_repeat_metrics.csv"
  )
)

readr::write_csv(
  aggregate_metrics,
  file.path(
    out_dir,
    "Fig6D_v7_aggregate_metrics.csv"
  )
)

readr::write_csv(
  top25_fold,
  file.path(
    out_dir,
    "Fig6D_v7_outer_fold_metrics.csv"
  )
)

readr::write_csv(
  selection_all,
  file.path(
    out_dir,
    "Fig6D_v7_all_outer_fold_selected_genes.csv"
  )
)

readr::write_csv(
  selection_frequency,
  file.path(
    out_dir,
    "Fig6D_v7_top25_selection_frequency.csv"
  )
)

pred_df <- tibble::tibble(
  Sample =
    sample_ids,
  Observed_balance =
    y,
  Full_model_prediction =
    full_pred,
  Nested_top25_prediction =
    pred_median,
  Nested_top25_prediction_mean =
    pred_mean
)

readr::write_csv(
  pred_df,
  file.path(
    out_dir,
    "Fig6D_v7_sample_predictions.csv"
  )
)

# ------------------------------------------------------------
# Model labels used consistently in plotting and QC
# ------------------------------------------------------------
full_name <- "All 1,373 non-programme genes"
top25_name <- "Nested top-25 genes"

# ------------------------------------------------------------
# 12. Main plot -- Nature-style Figure 6D
# ------------------------------------------------------------
# Scientific analysis above is unchanged.
# Only visual layout is revised for the final Figure 6 assembly.

plot_order <- c(
  full_name,
  top25_name
)

repeat_compare$Model <- factor(
  repeat_compare$Model,
  levels = plot_order
)

aggregate_metrics$Model <- factor(
  aggregate_metrics$Model,
  levels = plot_order
)

label_df <- aggregate_metrics %>%
  dplyr::mutate(
    Label = paste0(
      "n=", N_features,
      "\n",
      "rho=", sprintf("%.2f", Spearman_rho),
      " | R2=", sprintf("%.2f", CV_R2)
    )
  )


COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"

# Use the repeat-level data range to keep the vertical plotting region compact.
rho_min <- min(
  repeat_compare$Spearman_rho,
  na.rm = TRUE
)

rho_max <- max(
  repeat_compare$Spearman_rho,
  na.rm = TRUE
)

p <- ggplot(
  repeat_compare,
  aes(
    x = Model,
    y = Spearman_rho
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    colour = "#D5DEE2"
  ) +
  geom_line(
    aes(group = Repeat),
    linewidth = 0.22,
    alpha = 0.22,
    colour = "#91A5AD"
  ) +
  geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    linewidth = 0.32,
    fill = "#E8EEF0",
    colour = "#526C76"
  ) +
  geom_jitter(
    width = 0.045,
    height = 0,
    size = 1.16,
    alpha = 0.68,
    colour = "#456B78"
  ) +
  geom_point(
    data = aggregate_metrics,
    aes(
      x = Model,
      y = Spearman_rho
    ),
    inherit.aes = FALSE,
    size = 2.35,
    shape = 21,
    stroke = 0.50,
    fill = "white",
    colour = COL_TEXT
  ) +
  geom_text(
    data = label_df,
    aes(
      x = Model,
      y = rho_max + 0.052,
      label = Label
    ),
    inherit.aes = FALSE,
    size = 1.68,
    family = "Helvetica",
    colour = COL_TEXT,
    lineheight = 1.08
  ) +
  coord_cartesian(
    ylim = c(
      rho_min - 0.020,
      rho_max + 0.090
    ),
    clip = "off"
  ) +
  scale_x_discrete(
    labels = c(
      "All 1,373\nnon-programme genes",
      "Nested top-25\ngenes"
    )
  ) +
  labs(
    x = NULL,
    y = "Out-of-fold Spearman correlation",
    title = "Nested top-25 selection retains most reconstruction performance",
    subtitle = paste0(
      "Top 25 genes re-selected within each outer training fold; ",
      N_REPEATS,
      " repeated donor-blocked ",
      OUTER_K,
      "-fold nested CV"
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
      margin = margin(b = 1.2, unit = "pt")
    ),
    axis.title.y = element_text(
      size = 4.75,
      colour = COL_TEXT,
      margin = margin(r = 1.8, unit = "pt")
    ),
    axis.text.x = element_text(
      size = 4.45,
      colour = COL_TEXT,
      lineheight = 1.0,
      margin = margin(t = 1.5, unit = "pt")
    ),
    axis.text.y = element_text(
      size = 4.25,
      colour = COL_MUTED
    ),
    axis.line = element_line(
      linewidth = 0.30,
      colour = COL_TEXT
    ),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(
      linewidth = 0.24,
      colour = COL_TEXT
    ),
    axis.ticks.length = grid::unit(1.1, "pt"),
    plot.margin = margin(
      t = 1.5,
      r = 2.5,
      b = 2.0,
      l = 2.5,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# Nature-style export
# ------------------------------------------------------------
# Match Figure 6C bottom-row geometry.
FIG_W_MM <- 88
FIG_H_MM <- 66

out_dir_nature <- file.path(
  project_root,
  "figures",
  "FIGURE 6",
  "Panel_D_Nature_v7_noWilcoxon"
)

dir.create(
  out_dir_nature,
  recursive = TRUE,
  showWarnings = FALSE
)

png_file <- file.path(
  out_dir_nature,
  "Fig6_PanelD_nestedTop25_Nature_v7.png"
)

pdf_file <- file.path(
  out_dir_nature,
  "Fig6_PanelD_nestedTop25_Nature_v7.pdf"
)

tiff_file <- file.path(
  out_dir_nature,
  "Fig6_PanelD_nestedTop25_Nature_v7.tiff"
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
  "\nNature-style Figure 6D written to:\n",
  out_dir_nature,
  "\n"
)

# ------------------------------------------------------------
# 13. QC summary
# ------------------------------------------------------------
full_row <- aggregate_metrics %>%
  dplyr::filter(
    Model ==
      full_name
  )

top25_row <- aggregate_metrics %>%
  dplyr::filter(
    Model ==
      top25_name
  )

qc <- tibble::tibble(
  Metric = c(
    "Samples",
    "Non-programme genes available for fold-wise selection",
    "Genes selected per outer fold",
    "Outer models",
    "Full-model aggregate Spearman rho",
    "Nested-top25 aggregate Spearman rho",
    "Spearman difference top25 - full",
    "Full-model aggregate CV R2",
    "Nested-top25 aggregate CV R2",
    "CV R2 difference top25 - full",
    "Unique genes selected at least once",
    "Genes selected in >= 50% outer folds",
    "Genes selected in >= 75% outer folds"
  ),
  Value = c(
    EXPECTED_N_SAMPLES,
    EXPECTED_HELDOUT_N,
    TOP_N,
    N_REPEATS *
      OUTER_K,
    full_row$Spearman_rho,
    top25_row$Spearman_rho,
    top25_row$Spearman_rho -
      full_row$Spearman_rho,
    full_row$CV_R2,
    top25_row$CV_R2,
    top25_row$CV_R2 -
      full_row$CV_R2,
    nrow(
      selection_frequency
    ),
    sum(
      selection_frequency$Selection_frequency >=
        0.50
    ),
    sum(
      selection_frequency$Selection_frequency >=
        0.75
    )
  )
)

readr::write_csv(
  qc,
  file.path(
    out_dir,
    "Fig6D_QC_summary_v7.csv"
  )
)

# ------------------------------------------------------------
# 14. README
# ------------------------------------------------------------
readme <- c(
  "Figure 6D v7",
  "",
  "Nested top-25 reconstruction using training-fold-only marginal selection.",
  "",
  "Critical correction:",
  "The top-25 genes are NOT selected once using all 110 samples.",
  "They are re-selected independently inside each outer training fold.",
  "Outer test samples therefore never contribute to gene selection or inner model tuning.",
  "",
  "Procedure per outer fold:",
  "1. Compute Spearman rho for all 1,373 non-programme genes using training samples only.",
  "2. Select the 25 genes with largest absolute training-fold rho.",
  "3. Tune elastic-net alpha/lambda using inner 5-fold CV within the training fold.",
  "4. Reconstruct the untouched outer test fold.",
  "",
  "Comparator:",
  "Frozen Figure 5C all-1,373-gene repeated nested-CV model.",
  "",
  "Repeat-level inference note:",
  "The 20 repeated CV runs reuse the same 10 donors and 110 profiles.",
  "They are therefore used to describe stability/consistency rather than as independent replicates for formal hypothesis testing.",
  "No paired Wilcoxon P value is reported in v6.",
  "",
  "Interpretation:",
  "If nested top-25 performance remains strong, this supports substantial compression of state-balance information.",
  "The comparison is descriptive; no claim of complete leakage-free discovery is made because the target axis was prespecified from the same atlas.",
  "",
  paste0(
    "Full-model repeat source: ",
    normalizePath(
      full_repeat_file
    )
  ),
  paste0(
    "Target source: ",
    normalizePath(
      pred_file
    )
  ),
  paste0(
    "Sample alignment: ",
    alignment_mode
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig6_PanelD_v7.txt"
  )
)

cat(
  "\n============================================\n"
)

cat(
  "FIGURE 6D v7 COMPLETE\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Aggregate metrics:\n"
)

print(
  aggregate_metrics,
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
  "\nMost frequently selected genes:\n"
)

print(
  selection_frequency %>%
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
    "Fig6_PanelD_nestedTop25_v6.pdf"
  ),
  "\n"
)
