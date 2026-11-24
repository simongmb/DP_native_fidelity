# ============================================================
# 07_donor_aware_association_sensitivity.R
#
# Donor-aware sensitivity analysis for sample-level associations
# in Figures 4 and 5.
#
# Rationale:
# GSE272333 contains 110 follicular-compartment profiles from 10 donors.
# Sample-level Spearman rho is retained as the effect-size/ranking statistic,
# but nominal independence of all 110 profiles is not assumed for inference.
#
# This script adds:
#   A) within-donor permutation P for Native vs Alternative programme scores
#   B) within-donor permutation P and BH FDR for all 1,373 held-out
#      gene-vs-restoration-balance associations
#
# No programme membership is changed.
# No figure is overwritten.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

MASTER_SEED <- 20261102L
N_PERM <- 5000L

EXPECTED_BACKGROUND <- 1553L
EXPECTED_SAMPLES <- 110L
EXPECTED_DONORS <- 10L
EXPECTED_NATIVE <- 81L
EXPECTED_ALT <- 99L
EXPECTED_HELDOUT <- 1373L

# ------------------------------------------------------------
# 1. Locate project root
# ------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

find_project_root <- function(start = get_script_dir()) {
  x <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(
    x, getwd(),
    dirname(x), dirname(dirname(x)), dirname(dirname(dirname(x))),
    dirname(getwd()), dirname(dirname(getwd()))
  ))
  for (z in candidates) {
    if (dir.exists(file.path(z, "scripts")) &&
        dir.exists(file.path(z, "results"))) {
      return(z)
    }
  }
  stop("Cannot identify project root containing scripts/ and results/.")
}

project_root <- find_project_root()

cat("Project root: ", project_root, "\n", sep = "")

# ------------------------------------------------------------
# 2. Clean-session subexpr
# ------------------------------------------------------------

loader <- file.path(project_root, "scripts", "load_subexpr_clean.R")
if (!file.exists(loader)) {
  stop("Cannot find clean-session loader:\n", loader)
}

source(loader, local = FALSE)

expr <- as.matrix(subexpr)
storage.mode(expr) <- "numeric"
rownames(expr) <- toupper(trimws(rownames(expr)))

if (nrow(expr) != EXPECTED_BACKGROUND ||
    ncol(expr) != EXPECTED_SAMPLES) {
  stop(
    "QC FAIL: expected subexpr 1553 x 110; found ",
    nrow(expr), " x ", ncol(expr), "."
  )
}

if (is.null(colnames(expr)) || is.null(rownames(expr))) {
  stop("QC FAIL: subexpr must have gene and sample names.")
}

# ------------------------------------------------------------
# 3. Frozen programmes
# ------------------------------------------------------------

frozen_dir <- file.path(project_root, "results", "FROZEN_PROGRAMS")

native_file <- file.path(
  frozen_dir, "Native_identity_restoration_81_FROZEN.csv"
)
alt_file <- file.path(
  frozen_dir, "Alternative_state_rescue_99_FROZEN.csv"
)

for (f in c(native_file, alt_file)) {
  if (!file.exists(f)) stop("Missing frozen programme file:\n", f)
}

read_genes <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (!"GeneSymbol" %in% names(x)) stop("GeneSymbol missing in ", path)
  unique(toupper(trimws(as.character(x$GeneSymbol))))
}

native81 <- read_genes(native_file)
alt99 <- read_genes(alt_file)

stopifnot(
  length(native81) == EXPECTED_NATIVE,
  length(alt99) == EXPECTED_ALT,
  length(intersect(native81, alt99)) == 0L,
  all(native81 %in% rownames(expr)),
  all(alt99 %in% rownames(expr))
)

heldout <- setdiff(rownames(expr), c(native81, alt99))

if (length(heldout) != EXPECTED_HELDOUT) {
  stop(
    "QC FAIL: expected 1,373 held-out genes; found ",
    length(heldout), "."
  )
}

# ------------------------------------------------------------
# 4. Restoration-direction signs from canonical audit
# ------------------------------------------------------------

audit_file <- file.path(
  project_root,
  "results",
  "DP_native_fidelity_reconstruction",
  "02_gene_level_native_fidelity_audit.csv"
)

if (!file.exists(audit_file)) {
  stop("Cannot find canonical native-fidelity audit:\n", audit_file)
}

audit <- readr::read_csv(audit_file, show_col_types = FALSE)

if (!all(c("GeneSymbol", "Rescue_type") %in% names(audit))) {
  stop("Audit must contain GeneSymbol and Rescue_type.")
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

direction_map <- setNames(
  direction_tbl$Direction_sign,
  direction_tbl$GeneSymbol
)

program_genes <- c(native81, alt99)

if (!all(program_genes %in% names(direction_map))) {
  stop("QC FAIL: missing direction signs for programme genes.")
}

# ------------------------------------------------------------
# 5. Reconstruct primary programme scores and balance
# ------------------------------------------------------------

row_z <- function(m) {
  mu <- rowMeans(m)
  s <- apply(m, 1, sd)
  if (any(!is.finite(s) | s == 0)) {
    stop("QC FAIL: zero/non-finite SD in programme expression.")
  }
  sweep(sweep(m, 1, mu, "-"), 1, s, "/")
}

program_z <- row_z(expr[program_genes, , drop = FALSE])
program_sign <- direction_map[rownames(program_z)]
aligned <- program_z * as.numeric(program_sign)

native_score <- colMeans(aligned[native81, , drop = FALSE])
alt_score <- colMeans(aligned[alt99, , drop = FALSE])
balance <- native_score - alt_score

rho_programme <- suppressWarnings(
  cor(native_score, alt_score, method = "spearman")
)

# ------------------------------------------------------------
# 6. Recover donor IDs
# ------------------------------------------------------------

sample_ids <- colnames(expr)

# Same donor parsing convention used in donor-blocked Figure 5/6 scripts.
extract_donor <- function(ids) {
  out <- stringr::str_match(ids, "^.*_([A-J])_[^_]+$")[, 2]

  # Fallback for IDs formatted as COMPARTMENT_DONOR_...
  if (anyNA(out)) {
    parts <- strsplit(ids, "_", fixed = TRUE)
    if (all(lengths(parts) >= 2L)) {
      candidate <- vapply(parts, `[`, character(1), 2)
      if (length(unique(candidate)) == EXPECTED_DONORS) {
        out <- candidate
      }
    }
  }
  out
}

donor <- extract_donor(sample_ids)

if (anyNA(donor) || length(unique(donor)) != EXPECTED_DONORS) {
  stop(
    "QC FAIL: could not recover exactly 10 donor IDs from sample names.\n",
    paste(sample_ids, collapse = "\n")
  )
}

cat("Donors: ", paste(sort(unique(donor)), collapse = ", "), "\n", sep = "")

# ------------------------------------------------------------
# 7. Helpers: ranks and within-donor permutations
# ------------------------------------------------------------

rank_standardize <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  r <- r - mean(r, na.rm = TRUE)
  denom <- sqrt(sum(r^2, na.rm = TRUE))
  if (!is.finite(denom) || denom == 0) {
    stop("Cannot rank-standardize constant/non-finite vector.")
  }
  r / denom
}

permute_within_donor <- function(x, donor) {
  out <- x
  for (d in unique(donor)) {
    idx <- which(donor == d)
    out[idx] <- sample(x[idx], length(idx), replace = FALSE)
  }
  out
}

# ------------------------------------------------------------
# 8. Figure 4 donor-aware programme-score permutation
# ------------------------------------------------------------

set.seed(MASTER_SEED)

obs_prog <- rho_programme
b_prog <- 0L

for (b in seq_len(N_PERM)) {
  alt_perm <- permute_within_donor(alt_score, donor)
  rho_perm <- suppressWarnings(
    cor(native_score, alt_perm, method = "spearman")
  )
  if (is.finite(rho_perm) && abs(rho_perm) >= abs(obs_prog)) {
    b_prog <- b_prog + 1L
  }
}

p_prog_blocked <- (b_prog + 1) / (N_PERM + 1)

cat(
  sprintf(
    "Figure 4 programme correlation: rho = %.6f; within-donor permutation P = %.6g\n",
    obs_prog, p_prog_blocked
  )
)

# ------------------------------------------------------------
# 9. Figure 5 donor-aware held-out gene associations
# ------------------------------------------------------------

# Observed rho: same sample-level Spearman statistic used in the manuscript.
held_expr <- expr[heldout, , drop = FALSE]

obs_rho <- vapply(
  seq_len(nrow(held_expr)),
  function(i) suppressWarnings(
    cor(held_expr[i, ], balance, method = "spearman")
  ),
  numeric(1)
)

# Pre-rank and unit-standardize each gene once.
gene_rank_std <- t(apply(
  held_expr,
  1,
  rank_standardize
))

# Count null correlations at least as extreme as observed.
extreme_count <- integer(length(heldout))

set.seed(MASTER_SEED + 1L)

for (b in seq_len(N_PERM)) {
  bal_perm <- permute_within_donor(balance, donor)
  bal_rank_std <- rank_standardize(bal_perm)

  # Because both vectors are centered and unit norm after rank transform,
  # dot product equals Spearman correlation.
  rho_null <- as.numeric(gene_rank_std %*% bal_rank_std)

  extreme_count <- extreme_count +
    as.integer(abs(rho_null) >= abs(obs_rho))

  if (b %% 500L == 0L) {
    cat("Held-out permutation ", b, "/", N_PERM, "\n", sep = "")
  }
}

p_blocked <- (extreme_count + 1) / (N_PERM + 1)
fdr_blocked <- p.adjust(p_blocked, method = "BH")

n_sig_blocked <- sum(fdr_blocked < 0.05, na.rm = TRUE)

assoc_tbl <- tibble(
  GeneSymbol = heldout,
  Sample_level_Spearman_rho = obs_rho,
  Within_donor_empirical_P = p_blocked,
  Within_donor_BH_FDR = fdr_blocked,
  Significant_within_donor_FDR_0_05 = fdr_blocked < 0.05
) %>%
  arrange(Within_donor_BH_FDR, desc(abs(Sample_level_Spearman_rho)))

# ------------------------------------------------------------
# 10. Output
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "results",
  "DONOR_AWARE_ASSOCIATION_SENSITIVITY"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_tbl <- tibble(
  Analysis = c(
    "Figure 4 Native-vs-Alternative programme scores",
    "Figure 5 non-programme gene associations"
  ),
  Observed_effect = c(
    obs_prog,
    NA_real_
  ),
  Donor_aware_result = c(
    p_prog_blocked,
    n_sig_blocked
  ),
  Result_type = c(
    paste0("Within-donor empirical P (", N_PERM, " permutations)"),
    "Genes with within-donor BH FDR < 0.05"
  ),
  Primary_manuscript_reference = c(
    "Spearman rho = -0.81",
    "921 of 1,373 genes at sample-level BH FDR < 0.05"
  )
)

readr::write_csv(
  summary_tbl,
  file.path(out_dir, "donor_aware_sensitivity_summary.csv")
)

readr::write_csv(
  assoc_tbl,
  file.path(out_dir, "Fig5_nonprogramme_donor_aware_associations.csv")
)

programme_tbl <- tibble(
  Sample = sample_ids,
  Donor = donor,
  Native_score = native_score,
  Alternative_score = alt_score,
  Restoration_balance = balance
)

readr::write_csv(
  programme_tbl,
  file.path(out_dir, "Fig4_programme_scores_with_donor.csv")
)

cat("\n============================================================\n")
cat("DONOR-AWARE ASSOCIATION SENSITIVITY COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Programme score Spearman rho:          %.6f\n", obs_prog))
cat(sprintf("Within-donor permutation P:            %.6g\n", p_prog_blocked))
cat(sprintf("Non-programme genes tested:           %d\n", length(heldout)))
cat(sprintf("Primary sample-level BH-FDR count:     %d\n", 921L))
cat(sprintf("Donor-aware BH-FDR<0.05 count:         %d\n", n_sig_blocked))
cat(sprintf("Permutations:                           %d\n", N_PERM))
cat("\nOutput directory:\n", out_dir, "\n", sep = "")
