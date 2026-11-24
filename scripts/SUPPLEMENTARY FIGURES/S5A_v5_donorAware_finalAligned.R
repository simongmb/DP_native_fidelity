# ============================================================
# Supplementary Figure S5A v5 — DONOR-AWARE
# Fixed-reconstruction permutation test for donor-blocked Figure 5C
#
# Scientific question:
# Is the donor-blocked state-balance reconstruction observed in
# Figure 5C stronger than expected after disrupting within-donor
# label alignment while preserving donor-specific outcome structure?
#
# Input:
#   Fig5C_sample_predictions.csv from the finalized Figure 5C state-balance reconstruction
#
# Primary statistic:
#   Spearman correlation between the observed anatomy-referenced
#   state balance and donor-blocked reconstructed balance
#
# Null model:
#   10,000 WITHIN-DONOR permutations of observed balance.
#   Donor-blocked reconstructions are held fixed.
#
# Why:
#   This preserves each donor's outcome distribution and avoids
#   treating all 110 profiles as freely exchangeable.
#
# IMPORTANT:
# - No model is refit here.
# - This is a fixed-reconstruction conditional permutation test.
# - It is NOT a full-pipeline outcome permutation.
# ============================================================

required_pkgs <- c(
  "ggplot2",
  "readr",
  "dplyr",
  "tibble",
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
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(scales)
})

MASTER_SEED <- 20261028
set.seed(MASTER_SEED)

N_PERM <- 10000
EXPECTED_N <- 110
EXPECTED_DONORS <- 10

# ------------------------------------------------------------
# Locate project root / latest donor-blocked Fig5C output
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
    dirname(dirname(dirname(start_dir))),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )
)

candidate_roots <- candidate_roots[
  dir.exists(candidate_roots)
]

prediction_files <- character(0)

for (root in candidate_roots) {
  prediction_files <- c(
    prediction_files,
    list.files(
      root,
      pattern = "^Fig5C_sample_predictions\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
  )
}

prediction_files <- unique(prediction_files)

if (length(prediction_files) == 0) {
  stop(
    "\nQC STOP\nCannot find Fig5C_sample_predictions.csv.\n"
  )
}

score_file <- function(path) {
  s <- 0
  if (grepl("stateBalanceReconstruction_v12|Panel_C_stateBalanceReconstruction_v12", path, ignore.case = TRUE)) s <- s + 200
  if (grepl("donor_blocked_v8", path, ignore.case = TRUE)) s <- s + 100
  if (grepl("Panel_C_donor", path, ignore.case = TRUE)) s <- s + 50
  if (grepl("FIGURE 5", path, ignore.case = TRUE)) s <- s + 10
  s
}

scores <- vapply(
  prediction_files,
  score_file,
  numeric(1)
)

prediction_file <- prediction_files[
  which.max(scores)
]

cat(
  "\nUsing Figure 5C prediction file:\n",
  prediction_file,
  "\n"
)

pred0 <- readr::read_csv(
  prediction_file,
  show_col_types = FALSE
)

required_exact <- c(
  "Sample",
  "Donor",
  "Observed_balance",
  "Predicted_balance_median"
)

if (
  !all(
    required_exact %in%
      names(pred0)
  )
) {
  stop(
    paste0(
      "\nQC STOP\nRequired columns missing.\n",
      "Need: ",
      paste(required_exact, collapse = ", "),
      "\nAvailable: ",
      paste(names(pred0), collapse = ", "),
      "\n"
    )
  )
}

pred <- pred0 %>%
  dplyr::transmute(
    Sample = as.character(Sample),
    Donor = as.character(Donor),
    Observed = as.numeric(Observed_balance),
    Predicted = as.numeric(Predicted_balance_median)
  )

if (nrow(pred) != EXPECTED_N) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      EXPECTED_N,
      " profiles but found ",
      nrow(pred),
      ".\n"
    )
  )
}

if (dplyr::n_distinct(pred$Donor) != EXPECTED_DONORS) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      EXPECTED_DONORS,
      " donors but found ",
      dplyr::n_distinct(pred$Donor),
      ".\n"
    )
  )
}

if (
  any(!is.finite(pred$Observed)) ||
  any(!is.finite(pred$Predicted))
) {
  stop(
    "\nQC STOP\nNon-finite observed or predicted values detected.\n"
  )
}

# ------------------------------------------------------------
# Observed statistics
# ------------------------------------------------------------
observed_spearman <- suppressWarnings(
  cor(
    pred$Observed,
    pred$Predicted,
    method = "spearman"
  )
)

observed_pearson <- suppressWarnings(
  cor(
    pred$Observed,
    pred$Predicted,
    method = "pearson"
  )
)

sse_obs <- sum(
  (pred$Observed - pred$Predicted)^2
)

sst_obs <- sum(
  (pred$Observed - mean(pred$Observed))^2
)

observed_r2 <- 1 - sse_obs / sst_obs

observed_mae <- mean(
  abs(
    pred$Observed - pred$Predicted
  )
)

cat(
  "\nObserved donor-blocked Figure 5C reconstruction metrics:\n",
  "Spearman rho = ",
  sprintf("%.4f", observed_spearman),
  "\nPearson r = ",
  sprintf("%.4f", observed_pearson),
  "\nCV R2 = ",
  sprintf("%.4f", observed_r2),
  "\nMAE = ",
  sprintf("%.4f", observed_mae),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Donor-aware permutation null
# ------------------------------------------------------------
cat(
  "\nRunning ",
  format(N_PERM, big.mark = ","),
  " within-donor permutations...\n",
  sep = ""
)

null_rho <- numeric(N_PERM)
null_pearson <- numeric(N_PERM)
null_r2 <- numeric(N_PERM)

donor_levels <- sort(
  unique(
    pred$Donor
  )
)

set.seed(MASTER_SEED + 101)

for (b in seq_len(N_PERM)) {

  perm_obs <- pred$Observed

  for (d in donor_levels) {
    idx <- which(
      pred$Donor == d
    )

    perm_obs[idx] <- sample(
      pred$Observed[idx],
      length(idx),
      replace = FALSE
    )
  }

  null_rho[b] <- suppressWarnings(
    cor(
      perm_obs,
      pred$Predicted,
      method = "spearman"
    )
  )

  null_pearson[b] <- suppressWarnings(
    cor(
      perm_obs,
      pred$Predicted,
      method = "pearson"
    )
  )

  sse_perm <- sum(
    (perm_obs - pred$Predicted)^2
  )

  sst_perm <- sum(
    (perm_obs - mean(perm_obs))^2
  )

  null_r2[b] <- 1 -
    sse_perm /
    sst_perm
}

null_df <- tibble::tibble(
  Permutation = seq_len(N_PERM),
  Spearman_rho = null_rho,
  Pearson_r = null_pearson,
  R2 = null_r2
)

# Positive-direction empirical P values
p_spearman <- (
  sum(
    null_rho >= observed_spearman,
    na.rm = TRUE
  ) + 1
) / (
  N_PERM + 1
)

p_pearson <- (
  sum(
    null_pearson >= observed_pearson,
    na.rm = TRUE
  ) + 1
) / (
  N_PERM + 1
)

p_r2 <- (
  sum(
    null_r2 >= observed_r2,
    na.rm = TRUE
  ) + 1
) / (
  N_PERM + 1
)

null_q <- quantile(
  null_rho,
  c(0.025, 0.5, 0.975),
  na.rm = TRUE
)

summary_df <- tibble::tibble(
  Metric = c(
    "Observed_Spearman_rho",
    "Observed_Pearson_r",
    "Observed_CV_R2",
    "Observed_MAE",
    "Permutation_P_Spearman",
    "Permutation_P_Pearson",
    "Permutation_P_R2",
    "Null_Spearman_2.5pct",
    "Null_Spearman_median",
    "Null_Spearman_97.5pct",
    "N_profiles",
    "N_donors",
    "N_permutations"
  ),
  Value = c(
    observed_spearman,
    observed_pearson,
    observed_r2,
    observed_mae,
    p_spearman,
    p_pearson,
    p_r2,
    null_q[[1]],
    null_q[[2]],
    null_q[[3]],
    nrow(pred),
    length(donor_levels),
    N_PERM
  )
)

print(summary_df)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plot_df <- null_df %>%
  dplyr::filter(
    is.finite(Spearman_rho)
  )

p <- ggplot(
  plot_df,
  aes(
    x = Spearman_rho
  )
) +
  geom_histogram(
    bins = 55,
    fill = "#DDE7EA",
    colour = "white",
    linewidth = 0.18
  ) +
  geom_vline(
    xintercept = observed_spearman,
    colour = "#B55D43",
    linewidth = 0.65
  ) +
  annotate(
    "text",
    x = observed_spearman,
    y = Inf,
    label = paste0(
      "Observed rho = ",
      sprintf("%.3f", observed_spearman),
      "\nP = ",
      format.pval(
        p_spearman,
        digits = 2,
        eps = 1 / (N_PERM + 1)
      )
    ),
    hjust = 1.02,
    vjust = 1.35,
    size = 2.4,
    colour = "#7D3B2A"
  ) +
  labs(
    title = "Donor-aware permutation test",
    subtitle = "Observed labels permuted within donor; reconstructions held fixed",
    x = "Permuted Spearman rho",
    y = "Count"
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 8.5
    ),
    plot.subtitle = element_text(
      size = 7
    )
  )

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------
project_root <- dirname(
  dirname(
    dirname(
      prediction_file
    )
  )
)

if (
  !dir.exists(
    file.path(
      project_root,
      "figures"
    )
  )
) {
  project_root <- start_dir
}

out_dir <- file.path(
  project_root,
  "figures",
  "SUPPLEMENTARY FIGURES",
  "S5_Panel_A_predictionPermutation_donorAware_v5"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(
    out_dir,
    "S5_Panel_A_predictionPermutation_donorAware_v5.pdf"
  ),
  p,
  width = 100,
  height = 68,
  units = "mm",
  device = "pdf"
)

ggsave(
  file.path(
    out_dir,
    "S5_Panel_A_predictionPermutation_donorAware_v5.png"
  ),
  p,
  width = 100,
  height = 68,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "S5_Panel_A_predictionPermutation_donorAware_v5.tiff"
  ),
  p,
  width = 100,
  height = 68,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

readr::write_csv(
  null_df,
  file.path(
    out_dir,
    "S5A_permutation_null.csv"
  )
)

readr::write_csv(
  summary_df,
  file.path(
    out_dir,
    "S5A_permutation_summary.csv"
  )
)

readr::write_csv(
  pred,
  file.path(
    out_dir,
    "S5A_input_predictions.csv"
  )
)

readme <- c(
  "Supplementary Figure S5A v5 — donor-aware fixed-reconstruction permutation",
  "",
  "Input:",
  prediction_file,
  "",
  "Null model:",
  "Observed anatomy-referenced state balance is permuted within donor.",
  "Donor-blocked out-of-fold reconstructions are held fixed.",
  "",
  "Purpose:",
  "Test whether reconstruction/label concordance exceeds",
  "that expected after preserving donor-specific outcome structure.",
  "",
  "This is NOT a full-pipeline refitting permutation."
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_S5A_donorAware_v5.txt"
  )
)

cat(
  "\nS5A donor-aware permutation completed.\n",
  "Output directory:\n",
  out_dir,
  "\n",
  sep = ""
)