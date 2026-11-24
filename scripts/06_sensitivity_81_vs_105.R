# ============================================================
# 06_sensitivity_81_vs_105.R
#
# Sensitivity analysis for the six proliferation-associated genes
# removed from the manuscript-final 99-gene Alternative-state programme.
#
# Primary comparison:
#   Native-identity = 81 genes
#   Alternative-primary = 99 genes
#   Alternative-sensitivity = original threshold-defined 105 genes
#
# Outputs:
#   1) Native81 vs Alternative105 programme-score Spearman rho
#   2) Primary 81/99 balance vs sensitivity 81/105 balance Spearman rho
#   3) Number of significant non-programme genes for the 81/105 axis
#   4) Donor-blocked nested-CV Spearman rho
#   5) Donor-blocked nested-CV R2
#
# This is a sensitivity analysis only. It does NOT redefine the
# manuscript-final frozen programme.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(glmnet)
})

MASTER_SEED <- 20261028L
N_REPEATS <- 20L
ALPHA_GRID <- c(0, 0.25, 0.5, 0.75, 1)

EXPECTED_NATIVE <- 81L
EXPECTED_ALT99 <- 99L
EXPECTED_ALT105 <- 105L
EXPECTED_BACKGROUND <- 1553L
EXPECTED_SAMPLES <- 110L
EXPECTED_PREDICTORS_105 <- 1367L   # 1553 - 81 - 105

# ------------------------------------------------------------
# 1. Project root
# ------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

find_project_root <- function(start = get_script_dir()) {
  x <- normalizePath(start, winslash = "/", mustWork = TRUE)
  cand <- unique(c(
    x, getwd(),
    dirname(x), dirname(dirname(x)), dirname(dirname(dirname(x))),
    dirname(getwd()), dirname(dirname(getwd()))
  ))
  for (z in cand) {
    if (dir.exists(file.path(z, "scripts")) &&
        dir.exists(file.path(z, "results"))) {
      return(z)
    }
  }
  stop("Cannot identify project root containing both scripts/ and results/.")
}

project_root <- find_project_root()
cat("Project root: ", project_root, "\n", sep = "")

# ------------------------------------------------------------
# 2. Load clean-session subexpr (1553 x 110)
# ------------------------------------------------------------

loader_candidates <- c(
  file.path(project_root, "scripts", "load_subexpr_clean.R"),
  file.path(project_root, "load_subexpr_clean.R")
)
loader <- loader_candidates[file.exists(loader_candidates)][1]
if (is.na(loader)) stop("Cannot find scripts/load_subexpr_clean.R")

source(loader, local = FALSE)

expr_mat <- as.matrix(subexpr)
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- toupper(trimws(rownames(expr_mat)))

if (nrow(expr_mat) != EXPECTED_BACKGROUND ||
    ncol(expr_mat) != EXPECTED_SAMPLES) {
  stop("QC FAIL: subexpr must be 1553 x 110.")
}

# ------------------------------------------------------------
# 3. Read frozen programme files and 105-gene audit branch
# ------------------------------------------------------------

frozen_dir <- file.path(project_root, "results", "FROZEN_PROGRAMS")

native_file <- file.path(
  frozen_dir, "Native_identity_restoration_81_FROZEN.csv"
)
alt99_file <- file.path(
  frozen_dir, "Alternative_state_rescue_99_FROZEN.csv"
)
alt105_file <- file.path(
  frozen_dir, "Alternative_state_pre_exclusion_105_AUDIT.csv"
)

for (f in c(native_file, alt99_file, alt105_file)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

read_genes <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (!"GeneSymbol" %in% names(x)) stop("GeneSymbol missing: ", path)
  unique(toupper(trimws(as.character(x$GeneSymbol))))
}

native81 <- read_genes(native_file)
alt99 <- read_genes(alt99_file)
alt105 <- read_genes(alt105_file)

stopifnot(
  length(native81) == EXPECTED_NATIVE,
  length(alt99) == EXPECTED_ALT99,
  length(alt105) == EXPECTED_ALT105,
  all(alt99 %in% alt105),
  length(setdiff(alt105, alt99)) == 6L,
  length(intersect(native81, alt105)) == 0L
)

cat("Native programme:      ", length(native81), "\n", sep = "")
cat("Alternative primary:   ", length(alt99), "\n", sep = "")
cat("Alternative sensitivity:", length(alt105), "\n", sep = "")
cat("Six retained genes:    ",
    paste(sort(setdiff(alt105, alt99)), collapse = ", "), "\n\n", sep = "")

# ------------------------------------------------------------
# 4. Direction signs from final native-fidelity audit
# ------------------------------------------------------------

audit_hits <- list.files(
  project_root,
  pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
audit_hits <- audit_hits[!grepl("__MACOSX", audit_hits, fixed = TRUE)]
if (length(audit_hits) == 0L) stop("Cannot find 02_gene_level_native_fidelity_audit.csv")
if (length(audit_hits) > 1L) {
  res <- audit_hits[grepl("/results/", audit_hits)]
  if (length(res) == 1L) audit_hits <- res
}
audit_file <- audit_hits[1]

audit <- readr::read_csv(audit_file, show_col_types = FALSE)
required <- c("GeneSymbol", "Rescue_type")
if (!all(required %in% names(audit))) {
  stop("Audit table must contain GeneSymbol and Rescue_type.")
}

direction_tbl <- audit %>%
  transmute(
    GeneSymbol = toupper(trimws(as.character(GeneSymbol))),
    Rescue_type = as.character(Rescue_type),
    Direction_sign = case_when(
      grepl("up", Rescue_type, ignore.case = TRUE) ~ 1,
      grepl("down", Rescue_type, ignore.case = TRUE) ~ -1,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(is.finite(Direction_sign)) %>%
  distinct(GeneSymbol, .keep_all = TRUE)

direction_map <- setNames(direction_tbl$Direction_sign, direction_tbl$GeneSymbol)

all_needed <- unique(c(native81, alt99, alt105))
if (!all(all_needed %in% names(direction_map))) {
  stop("QC FAIL: missing restoration-direction signs for one or more programme genes.")
}
if (!all(all_needed %in% rownames(expr_mat))) {
  stop("QC FAIL: one or more programme genes are absent from subexpr.")
}

# ------------------------------------------------------------
# 5. Direction-aligned z-scored expression and programme scores
# ------------------------------------------------------------

row_z <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  s <- apply(m, 1, sd, na.rm = TRUE)
  if (any(!is.finite(s) | s == 0)) {
    stop("Zero/non-finite gene SD encountered.")
  }
  sweep(sweep(m, 1, mu, "-"), 1, s, "/")
}

z <- row_z(expr_mat)
aligned <- z

signs_all <- direction_map[rownames(aligned)]
if (any(!is.finite(signs_all))) {
  stop("QC FAIL: direction sign absent for one or more 1,553 genes.")
}
aligned <- aligned * as.numeric(signs_all)

score_program <- function(genes) {
  colMeans(aligned[genes, , drop = FALSE])
}

native_score <- score_program(native81)
alt99_score <- score_program(alt99)
alt105_score <- score_program(alt105)

balance99 <- native_score - alt99_score
balance105 <- native_score - alt105_score

rho_prog105 <- suppressWarnings(
  cor(native_score, alt105_score, method = "spearman", use = "complete.obs")
)
rho_axis_99_105 <- suppressWarnings(
  cor(balance99, balance105, method = "spearman", use = "complete.obs")
)
pearson_axis_99_105 <- suppressWarnings(
  cor(balance99, balance105, method = "pearson", use = "complete.obs")
)

# ------------------------------------------------------------
# 6. Held-out non-programme association for 81/105
# ------------------------------------------------------------

program105 <- unique(c(native81, alt105))
predictor_genes105 <- setdiff(rownames(expr_mat), program105)

if (length(predictor_genes105) != EXPECTED_PREDICTORS_105) {
  stop(
    "QC FAIL: expected 1,367 non-programme predictors for 81/105; found ",
    length(predictor_genes105)
  )
}

assoc_rho <- vapply(
  predictor_genes105,
  function(g) suppressWarnings(
    cor(expr_mat[g, ], balance105, method = "spearman", use = "complete.obs")
  ),
  numeric(1)
)

# Approximate two-sided Spearman p values via cor.test, matching manuscript logic.
assoc_p <- vapply(
  predictor_genes105,
  function(g) {
    suppressWarnings(
      tryCatch(
        cor.test(
          as.numeric(expr_mat[g, ]),
          balance105,
          method = "spearman",
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      )
    )
  },
  numeric(1)
)

assoc_fdr <- p.adjust(assoc_p, method = "BH")
n_sig105 <- sum(is.finite(assoc_fdr) & assoc_fdr < 0.05)

assoc_tbl <- tibble(
  GeneSymbol = predictor_genes105,
  Spearman_rho = assoc_rho,
  P_value = assoc_p,
  BH_FDR = assoc_fdr,
  Significant_FDR_0_05 = assoc_fdr < 0.05
) %>%
  arrange(BH_FDR, desc(abs(Spearman_rho)))

# ------------------------------------------------------------
# 7. Donor metadata from sample IDs
# ------------------------------------------------------------

sample_names <- colnames(expr_mat)
parts <- strsplit(sample_names, "_", fixed = TRUE)
if (any(lengths(parts) < 2L)) {
  stop("Cannot parse donor from sample names. Expected COMPARTMENT_DONOR_...")
}
donor_ids <- vapply(parts, `[`, character(1), 2)
donors <- sort(unique(donor_ids))

if (length(donors) != 10L) {
  stop("QC FAIL: expected 10 donors; found ", length(donors))
}

# Predictor matrix: samples x genes, unscaled here.
# glmnet standardizes inside each training fit.
X <- t(expr_mat[predictor_genes105, , drop = FALSE])
Y <- as.numeric(balance105)

# ------------------------------------------------------------
# 8. Helpers for donor-grouped nested CV
# ------------------------------------------------------------

make_inner_foldid <- function(train_donors, donor_train, seed) {
  set.seed(seed)
  shuffled <- sample(train_donors, length(train_donors), replace = FALSE)
  donor_fold <- setNames(rep(1:5, length.out = length(shuffled)), shuffled)
  foldid <- unname(donor_fold[donor_train])

  # Hard QC: every donor belongs to exactly one inner fold.
  chk <- tapply(foldid, donor_train, function(x) length(unique(x)))
  if (any(chk != 1L)) stop("Inner donor leakage detected.")
  foldid
}

r2_score <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

metrics <- function(obs, pred) {
  c(
    Pearson = suppressWarnings(cor(obs, pred, method = "pearson")),
    Spearman = suppressWarnings(cor(obs, pred, method = "spearman")),
    R2 = r2_score(obs, pred),
    MAE = mean(abs(obs - pred)),
    RMSE = sqrt(mean((obs - pred)^2))
  )
}

# ------------------------------------------------------------
# 9. 20 repeated donor-blocked nested-CV runs
# ------------------------------------------------------------

all_repeat_predictions <- matrix(
  NA_real_,
  nrow = length(Y),
  ncol = N_REPEATS,
  dimnames = list(sample_names, paste0("repeat_", seq_len(N_REPEATS)))
)

repeat_metrics <- vector("list", N_REPEATS)

for (rep_i in seq_len(N_REPEATS)) {
  cat("Repeat ", rep_i, "/", N_REPEATS, "\n", sep = "")

  pred <- rep(NA_real_, length(Y))

  for (outer_i in seq_along(donors)) {
    test_donor <- donors[outer_i]
    test_idx <- donor_ids == test_donor
    train_idx <- !test_idx

    donor_train <- donor_ids[train_idx]
    train_donors <- sort(unique(donor_train))

    if (test_donor %in% train_donors) {
      stop("Outer donor leakage detected.")
    }

    x_train <- X[train_idx, , drop = FALSE]
    y_train <- Y[train_idx]
    x_test <- X[test_idx, , drop = FALSE]

    foldid <- make_inner_foldid(
      train_donors = train_donors,
      donor_train = donor_train,
      seed = MASTER_SEED + rep_i * 1000L + outer_i
    )

    best_mse <- Inf
    best_fit <- NULL
    best_alpha <- NA_real_
    best_lambda <- NA_real_

    for (a in ALPHA_GRID) {
      set.seed(MASTER_SEED + rep_i * 10000L + outer_i * 100L + round(a * 10))
      cvfit <- glmnet::cv.glmnet(
        x = x_train,
        y = y_train,
        family = "gaussian",
        alpha = a,
        nlambda = 100,
        foldid = foldid,
        standardize = TRUE,
        intercept = TRUE,
        type.measure = "mse"
      )

      j <- which.min(cvfit$cvm)
      mse <- cvfit$cvm[j]

      if (is.finite(mse) && mse < best_mse) {
        best_mse <- mse
        best_fit <- cvfit
        best_alpha <- a
        best_lambda <- cvfit$lambda[j]
      }
    }

    if (is.null(best_fit)) stop("No valid inner-CV model found.")

    pred[test_idx] <- as.numeric(
      predict(
        best_fit,
        newx = x_test,
        s = best_lambda
      )
    )
  }

  if (any(!is.finite(pred))) stop("Missing/non-finite OOF predictions.")

  all_repeat_predictions[, rep_i] <- pred

  mm <- metrics(Y, pred)
  repeat_metrics[[rep_i]] <- tibble(
    Repeat = rep_i,
    Pearson_r = unname(mm["Pearson"]),
    Spearman_rho = unname(mm["Spearman"]),
    CV_R2 = unname(mm["R2"]),
    MAE = unname(mm["MAE"]),
    RMSE = unname(mm["RMSE"])
  )
}

repeat_metrics <- bind_rows(repeat_metrics)

aggregate_pred <- apply(all_repeat_predictions, 1, median, na.rm = TRUE)
agg <- metrics(Y, aggregate_pred)

# ------------------------------------------------------------
# 10. Output
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "results",
  "SENSITIVITY_81_vs_105"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_tbl <- tibble(
  Metric = c(
    "Native81_vs_Alternative105_programme_Spearman",
    "Balance99_vs_Balance105_Spearman",
    "Balance99_vs_Balance105_Pearson",
    "Nonprogramme_predictors_81_105",
    "Significant_nonprogramme_genes_BH_FDR_lt_0.05",
    "Aggregate_CV_Pearson",
    "Aggregate_CV_Spearman",
    "Aggregate_CV_R2",
    "Aggregate_CV_MAE",
    "Aggregate_CV_RMSE",
    "Median_repeat_CV_Spearman",
    "Median_repeat_CV_R2"
  ),
  Value = c(
    rho_prog105,
    rho_axis_99_105,
    pearson_axis_99_105,
    length(predictor_genes105),
    n_sig105,
    unname(agg["Pearson"]),
    unname(agg["Spearman"]),
    unname(agg["R2"]),
    unname(agg["MAE"]),
    unname(agg["RMSE"]),
    median(repeat_metrics$Spearman_rho),
    median(repeat_metrics$CV_R2)
  )
)

scores_tbl <- tibble(
  Sample = sample_names,
  Donor = donor_ids,
  Native81_score = native_score,
  Alternative99_score = alt99_score,
  Alternative105_score = alt105_score,
  Balance_81_99 = balance99,
  Balance_81_105 = balance105,
  Aggregate_OOF_prediction_81_105 = aggregate_pred
)

readr::write_csv(
  summary_tbl,
  file.path(out_dir, "81_vs_105_sensitivity_summary.csv")
)
readr::write_csv(
  scores_tbl,
  file.path(out_dir, "81_vs_105_sample_scores.csv")
)
readr::write_csv(
  assoc_tbl,
  file.path(out_dir, "81_vs_105_nonprogramme_associations.csv")
)
readr::write_csv(
  repeat_metrics,
  file.path(out_dir, "81_vs_105_nestedCV_repeat_metrics.csv")
)

pred_tbl <- as.data.frame(all_repeat_predictions) %>%
  rownames_to_column("Sample")
readr::write_csv(
  pred_tbl,
  file.path(out_dir, "81_vs_105_all_OOF_predictions.csv")
)

cat("\n============================================================\n")
cat("81 vs 105 SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Native81 vs Alt105 programme rho:      %.6f\n", rho_prog105))
cat(sprintf("Balance99 vs Balance105 Spearman:      %.6f\n", rho_axis_99_105))
cat(sprintf("Balance99 vs Balance105 Pearson:       %.6f\n", pearson_axis_99_105))
cat(sprintf("Non-programme predictors (81/105):     %d\n", length(predictor_genes105)))
cat(sprintf("BH-FDR<0.05 non-programme genes:       %d\n", n_sig105))
cat(sprintf("Aggregate CV Spearman:                 %.6f\n", agg["Spearman"]))
cat(sprintf("Aggregate CV R2:                       %.6f\n", agg["R2"]))
cat(sprintf("Median repeat CV Spearman:             %.6f\n", median(repeat_metrics$Spearman_rho)))
cat(sprintf("Median repeat CV R2:                   %.6f\n", median(repeat_metrics$CV_R2)))
cat("\nOutput directory:\n", out_dir, "\n", sep = "")
