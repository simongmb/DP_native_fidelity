# ============================================================
# 04_freeze_restoration_programmes.R
#
# Freeze the manuscript-final 81-gene Native-identity restoration
# programme and 99-gene Alternative-state rescue programme from the
# finalized 1,553-gene native-fidelity audit table.
#
# This script is intentionally independent of the R global environment.
# It does NOT use x, tt, screenshots, or historical programme files.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

RESCUE_Q3   <- 1.091014036
FIDELITY_Q1 <- -0.149005173
FIDELITY_Q3 <- 0.912572450

PROLIFERATION_EXCLUSIONS <- c(
  "ATAD2", "BIRC5", "CCNA2", "CENPA", "CIP2A", "MCM10"
)

EXPECTED_BACKGROUND <- 1553L
EXPECTED_HIGH_RESCUE <- 388L
EXPECTED_NATIVE <- 81L
EXPECTED_LOW_PRE <- 105L
EXPECTED_ALTERNATIVE <- 99L

find_project_root <- function(start = getwd()) {
  x <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(x, dirname(x), dirname(dirname(x)), dirname(dirname(dirname(x)))))
  for (z in candidates) {
    if (dir.exists(file.path(z, "results")) || dir.exists(file.path(z, "scripts"))) {
      return(z)
    }
  }
  x
}

find_unique_file <- function(root, pattern, label) {
  hits <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("__MACOSX", hits, fixed = TRUE)]
  if (length(hits) == 0L) stop("Cannot find ", label, " under project root: ", root)
  if (length(hits) > 1L) {
    stop(
      "Found multiple candidate files for ", label,
      ". Keep only the finalized file or pass a clean project tree.\n",
      paste(hits, collapse = "\n")
    )
  }
  hits[[1]]
}

project_root <- find_project_root()
audit_file <- find_unique_file(
  project_root,
  "^02_gene_level_native_fidelity_audit\\.csv$",
  "final native-fidelity audit"
)

cat("Project root: ", project_root, "\n", sep = "")
cat("Audit file:   ", audit_file, "\n\n", sep = "")

audit <- readr::read_csv(audit_file, show_col_types = FALSE)

required <- c(
  "GeneSymbol",
  "Rescue_fraction",
  "Rescue_type",
  "native_fidelity_direction_aligned",
  "anatomical_consistency_n_of_7"
)
missing <- setdiff(required, names(audit))
if (length(missing) > 0L) {
  stop("Audit table is missing required columns: ", paste(missing, collapse = ", "))
}

audit <- audit %>%
  mutate(
    GeneSymbol = toupper(trimws(as.character(GeneSymbol))),
    Rescue_fraction = as.numeric(Rescue_fraction),
    native_fidelity_direction_aligned = as.numeric(native_fidelity_direction_aligned),
    anatomical_consistency_n_of_7 = as.integer(anatomical_consistency_n_of_7)
  )

if (nrow(audit) != EXPECTED_BACKGROUND || dplyr::n_distinct(audit$GeneSymbol) != EXPECTED_BACKGROUND) {
  stop("QC FAIL: expected 1,553 unique analyzable genes; found ", nrow(audit), " rows and ", dplyr::n_distinct(audit$GeneSymbol), " unique symbols.")
}

# Exact manuscript-final inequalities. Note the strict >/< thresholds.
high_rescue <- audit %>%
  filter(Rescue_fraction > RESCUE_Q3)

native81 <- audit %>%
  filter(
    Rescue_fraction > RESCUE_Q3,
    native_fidelity_direction_aligned > FIDELITY_Q3,
    anatomical_consistency_n_of_7 >= 5L
  )

low_pre105 <- audit %>%
  filter(
    Rescue_fraction > RESCUE_Q3,
    native_fidelity_direction_aligned < FIDELITY_Q1,
    anatomical_consistency_n_of_7 <= 2L
  )

alternative99 <- low_pre105 %>%
  filter(!GeneSymbol %in% PROLIFERATION_EXCLUSIONS)

# Hard QC against manuscript-final counts.
stopifnot(
  nrow(high_rescue) == EXPECTED_HIGH_RESCUE,
  nrow(native81) == EXPECTED_NATIVE,
  nrow(low_pre105) == EXPECTED_LOW_PRE,
  nrow(alternative99) == EXPECTED_ALTERNATIVE,
  length(intersect(native81$GeneSymbol, alternative99$GeneSymbol)) == 0L
)

# Direction composition QC.
native_direction <- table(native81$Rescue_type)
alt_direction <- table(alternative99$Rescue_type)

if (unname(native_direction["Restored_upward"]) != 55L ||
    unname(native_direction["Restored_downward"]) != 26L) {
  stop("QC FAIL: Native programme direction composition is not 55 up / 26 down.")
}
if (unname(alt_direction["Restored_upward"]) != 30L ||
    unname(alt_direction["Restored_downward"]) != 69L) {
  stop("QC FAIL: Alternative programme direction composition is not 30 up / 69 down.")
}

# Verify all six exclusions are present in the low-fidelity pre-exclusion branch
# and absent from the final Alternative programme.
if (!all(PROLIFERATION_EXCLUSIONS %in% low_pre105$GeneSymbol)) {
  stop("QC FAIL: one or more predefined proliferation exclusions are absent from the 105-gene pre-exclusion branch.")
}
if (any(PROLIFERATION_EXCLUSIONS %in% alternative99$GeneSymbol)) {
  stop("QC FAIL: one or more predefined proliferation exclusions remain in the 99-gene Alternative programme.")
}

out_dir <- file.path(project_root, "results", "FROZEN_PROGRAMS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Full annotated programme tables are useful for audit/provenance.
readr::write_csv(
  native81 %>% arrange(GeneSymbol),
  file.path(out_dir, "Native_identity_restoration_81_FROZEN.csv")
)
readr::write_csv(
  alternative99 %>% arrange(GeneSymbol),
  file.path(out_dir, "Alternative_state_rescue_99_FROZEN.csv")
)
readr::write_csv(
  audit %>% arrange(GeneSymbol),
  file.path(out_dir, "Analyzable_background_1553_FROZEN.csv")
)
readr::write_csv(
  low_pre105 %>% arrange(GeneSymbol),
  file.path(out_dir, "Alternative_state_pre_exclusion_105_AUDIT.csv")
)

qc <- tibble::tibble(
  Item = c(
    "Analyzable background",
    "High rescue (> Q3)",
    "Native-identity restoration",
    "Low-fidelity pre-exclusion",
    "Alternative-state rescue",
    "Programme overlap",
    "Native restored-upward",
    "Native restored-downward",
    "Alternative restored-upward",
    "Alternative restored-downward",
    "Excluded proliferation genes removed"
  ),
  Value = c(
    nrow(audit),
    nrow(high_rescue),
    nrow(native81),
    nrow(low_pre105),
    nrow(alternative99),
    length(intersect(native81$GeneSymbol, alternative99$GeneSymbol)),
    unname(native_direction["Restored_upward"]),
    unname(native_direction["Restored_downward"]),
    unname(alt_direction["Restored_upward"]),
    unname(alt_direction["Restored_downward"]),
    sum(PROLIFERATION_EXCLUSIONS %in% low_pre105$GeneSymbol & !PROLIFERATION_EXCLUSIONS %in% alternative99$GeneSymbol)
  )
)
readr::write_csv(qc, file.path(out_dir, "Frozen_program_QC.csv"))

cat("============================================\n")
cat("FROZEN PROGRAMMES CREATED SUCCESSFULLY\n")
cat("============================================\n")
cat("Analyzable background:       ", nrow(audit), "\n", sep = "")
cat("High rescue:                 ", nrow(high_rescue), "\n", sep = "")
cat("Native-identity programme:   ", nrow(native81), " (55 up / 26 down)\n", sep = "")
cat("Low-fidelity pre-exclusion:  ", nrow(low_pre105), "\n", sep = "")
cat("Alternative-state programme: ", nrow(alternative99), " (30 up / 69 down)\n", sep = "")
cat("Output directory:            ", out_dir, "\n", sep = "")
