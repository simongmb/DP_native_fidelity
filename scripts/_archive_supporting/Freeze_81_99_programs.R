# ============================================================
# Freeze_81_99_programs.R
#
# Purpose:
# Permanently freeze the two Figure 3 programs into
# simple standardized CSV files.
#
# Native-identity restoration = 81 genes
# Alternative-state rescue    = 99 genes
#
# This script intentionally avoids automatic project-wide
# searching except where absolutely necessary.
# ============================================================


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})


# ------------------------------------------------------------
# 1. Helper
# ------------------------------------------------------------

clean_gene <- function(x) {
  
  x <- as.character(x)
  
  x <- stringr::str_trim(x)
  
  x[x == ""] <- NA_character_
  
  x <- toupper(x)
  
  unique(
    x[
      !is.na(x)
    ]
  )
}


# ------------------------------------------------------------
# 2. Locate project root
# ------------------------------------------------------------

project_root <- getwd()


# If currently inside scripts/ or another subfolder,
# walk upward until results/ and figures/ are found.

candidate_roots <- unique(
  c(
    project_root,
    dirname(project_root),
    dirname(dirname(project_root)),
    dirname(dirname(dirname(project_root)))
  )
)


valid_root <- candidate_roots[
  vapply(
    candidate_roots,
    function(z) {
      dir.exists(file.path(z, "results")) &&
        dir.exists(file.path(z, "figures"))
    },
    logical(1)
  )
]


if (length(valid_root) == 0) {
  
  stop(
    "\nCannot identify project root containing results/ and figures/.\n"
  )
}


project_root <- valid_root[1]


cat(
  "\nProject root:\n",
  project_root,
  "\n"
)


# ------------------------------------------------------------
# 3. Locate background audit
# ------------------------------------------------------------

audit_files <- list.files(
  project_root,
  pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)


if (length(audit_files) == 0) {
  
  stop(
    "\nCannot find 02_gene_level_native_fidelity_audit.csv\n"
  )
}


audit_file <- audit_files[1]


audit <- readr::read_csv(
  audit_file,
  show_col_types = FALSE
)


gene_candidates <- c(
  "GeneSymbol",
  "gene_symbol",
  "Gene",
  "gene",
  "symbol",
  "SYMBOL"
)


gene_col <- gene_candidates[
  gene_candidates %in% names(audit)
][1]


if (is.na(gene_col)) {
  
  stop(
    "\nCannot identify gene-symbol column in audit table.\n"
  )
}


background_genes <- clean_gene(
  audit[[gene_col]]
)


if (length(background_genes) != 1553) {
  
  stop(
    paste0(
      "\nQC FAIL: expected 1553 background genes, found ",
      length(background_genes),
      ".\n"
    )
  )
}


cat(
  "\nBackground QC: 1553 genes OK\n"
)


# ------------------------------------------------------------
# 4. Recover native 81
#
# First use an already exported native list if available.
# ------------------------------------------------------------

native_files <- list.files(
  project_root,
  pattern = "native.*81.*\\.csv$|81.*native.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


native81 <- character(0)


for (f in native_files) {
  
  tmp <- tryCatch(
    readr::read_csv(
      f,
      show_col_types = FALSE
    ),
    error = function(e) NULL
  )
  
  
  if (is.null(tmp)) next
  
  
  gc <- gene_candidates[
    gene_candidates %in% names(tmp)
  ][1]
  
  
  if (is.na(gc)) next
  
  
  g <- clean_gene(
    tmp[[gc]]
  )
  
  
  if (length(g) == 81) {
    
    native81 <- g
    
    cat(
      "\nNative 81 recovered from:\n",
      f,
      "\n"
    )
    
    break
  }
}


# ------------------------------------------------------------
# 5. Fallback for native 81:
# use tt in current Environment if it is the 81-gene table
# ------------------------------------------------------------

if (
  length(native81) != 81 &&
  exists(
    "tt",
    envir = .GlobalEnv
  )
) {
  
  tt_obj <- get(
    "tt",
    envir = .GlobalEnv
  )
  
  
  if (
    is.data.frame(tt_obj) &&
    nrow(tt_obj) == 81
  ) {
    
    gc <- gene_candidates[
      gene_candidates %in%
        names(tt_obj)
    ][1]
    
    
    if (!is.na(gc)) {
      
      native81 <- clean_gene(
        tt_obj[[gc]]
      )
      
      cat(
        "\nNative 81 recovered from current object: tt\n"
      )
    }
  }
}


if (length(native81) != 81) {
  
  stop(
    paste0(
      "\nCould not recover the frozen native program exactly.\n",
      "Expected 81 genes, found ",
      length(native81),
      ".\n"
    )
  )
}


# ------------------------------------------------------------
# 6. Recover alternative 99 DIRECTLY from current object x
#
# Your current Environment shows:
# x = chr [1:99]
# ------------------------------------------------------------

if (
  !exists(
    "x",
    envir = .GlobalEnv
  )
) {
  
  stop(
    paste0(
      "\nObject 'x' is not present in the current R Environment.\n",
      "The screenshot shows x as a character vector of length 99,\n",
      "so run this script in the same current RStudio session.\n"
    )
  )
}


x_obj <- get(
  "x",
  envir = .GlobalEnv
)


alternative99 <- clean_gene(
  x_obj
)


cat(
  "\nObject x contains:",
  length(alternative99),
  "unique genes\n"
)


if (length(alternative99) != 99) {
  
  stop(
    paste0(
      "\nQC FAIL: object x must contain exactly 99 unique genes.\n",
      "Current unique count = ",
      length(alternative99),
      "\n"
    )
  )
}


# ------------------------------------------------------------
# 7. Strong biological QC for the 99 set
# ------------------------------------------------------------

expected_alt_markers <- c(
  "MYL9",
  "MYLK",
  "LMOD1",
  "FLNA",
  "VCL",
  "FERMT2"
)


excluded_proliferation <- c(
  "ATAD2",
  "BIRC5",
  "CCNA2",
  "CENPA",
  "CIP2A",
  "MCM10"
)


marker_present <- expected_alt_markers %in%
  alternative99


excluded_present <- excluded_proliferation %in%
  alternative99


cat(
  "\nAlternative-state marker QC:\n"
)


print(
  tibble::tibble(
    Gene = expected_alt_markers,
    Present = marker_present
  )
)


cat(
  "\nExcluded proliferation-gene QC:\n"
)


print(
  tibble::tibble(
    Gene = excluded_proliferation,
    Present = excluded_present
  )
)


if (!all(marker_present)) {
  
  stop(
    paste0(
      "\nQC FAIL:\n",
      "Object x has 99 genes but does NOT contain all expected ",
      "alternative-state representative genes.\n",
      "Therefore I will not assume x is the frozen 99 program.\n"
    )
  )
}


if (any(excluded_present)) {
  
  stop(
    paste0(
      "\nQC FAIL:\n",
      "Object x still contains one or more of the six excluded ",
      "proliferation-associated genes.\n"
    )
  )
}


# ------------------------------------------------------------
# 8. Background membership QC
# ------------------------------------------------------------

native_missing <- setdiff(
  native81,
  background_genes
)


alt_missing <- setdiff(
  alternative99,
  background_genes
)


if (length(native_missing) > 0) {
  
  stop(
    paste0(
      "\nNative genes missing from 1553 background:\n",
      paste(
        native_missing,
        collapse = ", "
      ),
      "\n"
    )
  )
}


if (length(alt_missing) > 0) {
  
  stop(
    paste0(
      "\nAlternative genes missing from 1553 background:\n",
      paste(
        alt_missing,
        collapse = ", "
      ),
      "\n"
    )
  )
}


# ------------------------------------------------------------
# 9. Check overlap
# ------------------------------------------------------------

program_overlap <- intersect(
  native81,
  alternative99
)


cat(
  "\n81/99 overlap:",
  length(program_overlap),
  "genes\n"
)


if (length(program_overlap) > 0) {
  
  cat(
    paste(
      program_overlap,
      collapse = ", "
    ),
    "\n"
  )
}


# ------------------------------------------------------------
# 10. Save STANDARD frozen files
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "results",
  "FROZEN_PROGRAMS"
)


dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


native_file <- file.path(
  out_dir,
  "Native_identity_restoration_81_FROZEN.csv"
)


alternative_file <- file.path(
  out_dir,
  "Alternative_state_rescue_99_FROZEN.csv"
)


background_file <- file.path(
  out_dir,
  "Analyzable_background_1553_FROZEN.csv"
)


readr::write_csv(
  
  tibble::tibble(
    GeneSymbol = sort(native81)
  ),
  
  native_file
)


readr::write_csv(
  
  tibble::tibble(
    GeneSymbol = sort(alternative99)
  ),
  
  alternative_file
)


readr::write_csv(
  
  tibble::tibble(
    GeneSymbol = sort(background_genes)
  ),
  
  background_file
)


# ------------------------------------------------------------
# 11. Save QC table
# ------------------------------------------------------------

qc <- tibble::tibble(
  
  Item = c(
    "Analyzable background",
    "Native-identity restoration",
    "Alternative-state rescue",
    "Program overlap",
    "Expected alternative markers present",
    "Excluded proliferation genes remaining"
  ),
  
  Value = c(
    length(background_genes),
    length(native81),
    length(alternative99),
    length(program_overlap),
    sum(marker_present),
    sum(excluded_present)
  ),
  
  Expected = c(
    "1553",
    "81",
    "99",
    "report only",
    "6",
    "0"
  )
)


readr::write_csv(
  qc,
  file.path(
    out_dir,
    "Frozen_program_QC.csv"
  )
)


# ------------------------------------------------------------
# 12. Final console output
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)


cat(
  "FROZEN PROGRAM FILES CREATED SUCCESSFULLY\n"
)


cat(
  "============================================\n"
)


cat(
  "\nBackground =",
  length(background_genes),
  "\n"
)


cat(
  "Native =",
  length(native81),
  "\n"
)


cat(
  "Alternative =",
  length(alternative99),
  "\n"
)


cat(
  "\nSaved to:\n",
  out_dir,
  "\n\n"
)


cat(
  native_file,
  "\n"
)


cat(
  alternative_file,
  "\n"
)


cat(
  background_file,
  "\n"
)