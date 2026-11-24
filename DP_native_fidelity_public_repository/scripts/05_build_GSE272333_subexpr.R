# ============================================================
# 05_build_GSE272333_subexpr.R
# Clean-session reconstruction of the expression matrix used by Figures 4–6.
# Correct Figure 4–6 background: 1,553 genes x 110 profiles.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

EXPECTED_STRICT <- 1571L
EXPECTED_ANALYZABLE <- 1553L
EXPECTED_SAMPLES <- 110L

find_project_root <- function(start = getwd()) {
  x <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(x, dirname(x), dirname(dirname(x)), dirname(dirname(dirname(x)))))
  for (z in candidates) {
    if (dir.exists(file.path(z, "results")) || dir.exists(file.path(z, "scripts"))) return(z)
  }
  x
}

find_unique <- function(root, patterns, label) {
  hits <- character(0)
  for (pat in patterns) {
    hits <- c(hits, list.files(root, pattern = pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
  }
  hits <- unique(hits[file.exists(hits) & !grepl("__MACOSX", hits, fixed = TRUE)])
  if (length(hits) == 0L) stop("Cannot find ", label, ".")
  if (length(hits) > 1L) stop("Multiple candidate files found for ", label, ":\n", paste(hits, collapse = "\n"))
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

read_any <- function(path) {
  if (grepl("\\.(tsv|tab|txt)(\\.gz)?$", path, ignore.case = TRUE)) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
}

strip_ensembl_version <- function(x) sub("\\.\\d+$", "", as.character(x))

project_root <- find_project_root()
strict_file <- find_unique(project_root, c("^GSE44765_STRICT_rescue_genes\\.csv$"), "GSE44765 strict-rescue table")
mapping_file <- find_unique(project_root, c("^GSE272333_gene_ID_mapping\\.csv$", "gene_ID_mapping_audit\\.csv$"), "GSE272333 Ensembl-to-HGNC mapping table")
counts_file <- find_unique(project_root, c(
  "GSE272333.*normalized.*count.*\\.tab\\.gz$",
  "GSE272333.*normalized.*count.*\\.tsv\\.gz$",
  "GSE272333.*normalized.*count.*\\.txt\\.gz$",
  "GSE272333.*normalized.*count.*\\.csv\\.gz$"
), "GSE272333 normalized-count matrix")

strict <- read_any(strict_file)
if (!"GeneSymbol" %in% names(strict)) stop("Strict-rescue table lacks GeneSymbol.")
strict_symbols <- unique(toupper(trimws(as.character(strict$GeneSymbol))))
strict_symbols <- strict_symbols[!is.na(strict_symbols) & strict_symbols != ""]
if (length(strict_symbols) != EXPECTED_STRICT) stop("QC FAIL: expected 1,571 strict-rescue genes; found ", length(strict_symbols))

mapping <- read_any(mapping_file)
for (nm in names(mapping)) {
  low <- tolower(nm)
  if (low %in% c("ensembl_id","ensembl","gene_id","ensembl_gene_id")) names(mapping)[names(mapping) == nm] <- "Ensembl_ID"
  if (low %in% c("approved_symbol","genesymbol","gene_symbol","symbol","hgnc_symbol")) names(mapping)[names(mapping) == nm] <- "Approved_Symbol"
}
if (!all(c("Ensembl_ID","Approved_Symbol") %in% names(mapping))) stop("Mapping table must contain Ensembl_ID and Approved_Symbol.")
mapping <- mapping %>%
  transmute(Ensembl_ID = strip_ensembl_version(Ensembl_ID), Approved_Symbol = toupper(trimws(as.character(Approved_Symbol)))) %>%
  filter(!is.na(Ensembl_ID), !is.na(Approved_Symbol), Approved_Symbol != "") %>%
  distinct(Ensembl_ID, .keep_all = TRUE)

counts <- read_any(counts_file)
if (ncol(counts) < 2L) stop("Normalized-count matrix has <2 columns.")
names(counts)[1] <- "Ensembl_ID"
counts$Ensembl_ID <- strip_ensembl_version(counts$Ensembl_ID)
sample_cols <- setdiff(names(counts), "Ensembl_ID")
if (length(sample_cols) != EXPECTED_SAMPLES) stop("QC FAIL: expected 110 atlas profiles; found ", length(sample_cols))

parts <- strsplit(sample_cols, "_", fixed = TRUE)
if (any(lengths(parts) < 2L)) stop("Cannot parse compartment/donor from one or more sample IDs.")
sample_meta <- tibble(sample = sample_cols, compartment = vapply(parts, `[`, character(1), 1), donor = vapply(parts, `[`, character(1), 2))

x <- counts %>% left_join(mapping, by = "Ensembl_ID") %>% filter(Approved_Symbol %in% strict_symbols)
for (s in sample_cols) x[[s]] <- suppressWarnings(as.numeric(x[[s]]))

x_gene <- x %>%
  group_by(Approved_Symbol) %>%
  summarise(across(all_of(sample_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(Approved_Symbol)

if (nrow(x_gene) != EXPECTED_ANALYZABLE) stop("QC FAIL: expected 1,553 analyzable genes; found ", nrow(x_gene))

mat <- as.matrix(x_gene[, sample_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- x_gene$Approved_Symbol
subexpr <- log2(mat + 1)

if (anyDuplicated(rownames(subexpr))) stop("QC FAIL: duplicated gene symbols.")
if (nrow(subexpr) != 1553L || ncol(subexpr) != 110L) stop("QC FAIL: subexpr dimensions must be 1553 x 110.")
if (any(!is.finite(subexpr))) stop("QC FAIL: non-finite values in subexpr.")

out_dir <- file.path(project_root, "results", "GSE272333_expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(subexpr, file.path(out_dir, "subexpr_1553_log2.rds"), compress = "xz")
readr::write_csv(data.frame(GeneSymbol = rownames(subexpr), subexpr, check.names = FALSE), file.path(out_dir, "subexpr_1553_log2.csv.gz"))
readr::write_csv(sample_meta, file.path(out_dir, "subexpr_sample_metadata.csv"))
readr::write_csv(tibble(
  item = c("strict_rescue_genes","subexpr_genes","subexpr_samples","unique_gene_symbols","unique_donors","compartments"),
  observed = c(length(strict_symbols), nrow(subexpr), ncol(subexpr), length(unique(rownames(subexpr))), length(unique(sample_meta$donor)), paste(sort(unique(sample_meta$compartment)), collapse = ";")),
  expected = c("1571","1553","110","1553","10","CTS;DC;DP;HM1;HM2;IRS1;ORS1;ORS2")
), file.path(out_dir, "subexpr_QC.csv"))

cat("QC PASS\nsubexpr dimensions: ", nrow(subexpr), " x ", ncol(subexpr), "\n", sep = "")
