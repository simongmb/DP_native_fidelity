# Shared clean-session loader for Figures 4–6.
find_project_root_for_subexpr <- function(start = getwd()) {
  x <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(x, dirname(x), dirname(dirname(x)), dirname(dirname(dirname(x)))))
  for (z in candidates) {
    if (dir.exists(file.path(z, "results")) || dir.exists(file.path(z, "scripts"))) return(z)
  }
  x
}

project_root_subexpr <- find_project_root_for_subexpr()
subexpr_rds <- file.path(project_root_subexpr, "results", "GSE272333_expression", "subexpr_1553_log2.rds")
if (!file.exists(subexpr_rds)) {
  builder_candidates <- c(
    file.path(project_root_subexpr, "scripts", "05_build_GSE272333_subexpr.R"),
    file.path(project_root_subexpr, "05_build_GSE272333_subexpr.R")
  )
  builder_file <- builder_candidates[file.exists(builder_candidates)][1]
  if (is.na(builder_file)) stop("subexpr RDS is absent and 05_build_GSE272333_subexpr.R cannot be found. Expected: ", subexpr_rds)
  source(builder_file, local = FALSE)
}
subexpr <- readRDS(subexpr_rds)
subexpr <- as.matrix(subexpr)
storage.mode(subexpr) <- "numeric"
rownames(subexpr) <- toupper(trimws(rownames(subexpr)))
if (nrow(subexpr) != 1553L) stop("QC FAIL: subexpr must contain exactly 1,553 genes; found ", nrow(subexpr))
if (ncol(subexpr) != 110L) stop("QC FAIL: subexpr must contain exactly 110 profiles; found ", ncol(subexpr))
if (anyDuplicated(rownames(subexpr))) stop("QC FAIL: duplicated gene symbols in subexpr.")
if (is.null(colnames(subexpr))) stop("QC FAIL: subexpr has no sample names.")
cat("Loaded clean-session subexpr: ", nrow(subexpr), " genes x ", ncol(subexpr), " profiles\n", sep = "")
