############################################################
# GSE44765 preprocessing + QC + PCA + correlation heatmap
# HairRegeneration_AI project
#
# INPUT:
#   data_raw/GSE44765/*.CEL.gz
#
# OUTPUT:
#   results/
#   figures/
#
# Workflow:
#   CEL -> RMA -> QC -> PCA -> correlation heatmap
############################################################


############################################################
# 0. Start clean
############################################################

options(stringsAsFactors = FALSE)

# Close any graphics devices left open by previous failed runs
graphics.off()

cat("\n")
cat("====================================================\n")
cat("   GSE44765 preprocessing started\n")
cat("====================================================\n\n")

cat("Working directory:\n")
cat(getwd(), "\n\n")


############################################################
# 1. Create folders
############################################################

dir.create(
  "data_raw/GSE44765",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "scripts",
  showWarnings = FALSE
)

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)


############################################################
# 2. Install required CRAN packages
############################################################

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "R.utils"
)

for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    cat("Installing CRAN package:", pkg, "\n")
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}


############################################################
# 3. Install Bioconductor manager
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


############################################################
# 4. Install required Bioconductor packages
############################################################

bioc_packages <- c(
  "affy",
  "limma",
  "AnnotationDbi",
  "hgu133plus2.db",
  "hgu133plus2cdf"
)

for (pkg in bioc_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    cat("Installing Bioconductor package:", pkg, "\n")
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


############################################################
# 5. Load packages
############################################################

suppressPackageStartupMessages({
  
  library(affy)
  
  library(limma)
  
  library(AnnotationDbi)
  
  library(hgu133plus2.db)
  
  library(hgu133plus2cdf)
  
  library(ggplot2)
  
  library(ggrepel)
  
  library(pheatmap)
  
  library(R.utils)
  
})

cat("\nAll required packages loaded successfully.\n\n")


############################################################
# 6. Locate CEL files
############################################################

cel_dir <- "data_raw/GSE44765"

cel_files <- list.files(
  path = cel_dir,
  pattern = "\\.CEL(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat(
  "Number of CEL files detected:",
  length(cel_files),
  "\n\n"
)


############################################################
# 7. Check that exactly 18 CEL files exist
############################################################

if (length(cel_files) != 18) {
  
  cat("\nDetected files:\n")
  
  print(
    basename(cel_files)
  )
  
  stop(
    paste0(
      "\nERROR:\n",
      "Expected 18 CEL/CEL.gz files, but found ",
      length(cel_files),
      ".\n\n",
      "Please place all 18 GSE44765 CEL.gz files in:\n",
      "data_raw/GSE44765/\n"
    )
  )
}


############################################################
# 8. Extract DP number from filename
############################################################

extract_dp_number <- function(x) {
  
  x <- basename(x)
  
  x <- sub(
    ".*_DP([0-9]+)\\.CEL(\\.gz)?$",
    "\\1",
    x,
    ignore.case = TRUE
  )
  
  as.integer(x)
}


############################################################
# 9. Sort DP1 -> DP18 correctly
############################################################

dp_number <- extract_dp_number(cel_files)

if (any(is.na(dp_number))) {
  
  stop(
    "\nERROR: Could not extract DP numbers from filenames.\n"
  )
}

cel_files <- cel_files[
  order(dp_number)
]

dp_number <- extract_dp_number(cel_files)

cat("CEL files after sorting:\n\n")

print(
  basename(cel_files)
)

cat("\n")

if (!identical(dp_number, 1:18)) {
  
  stop(
    "\nERROR: Files are not exactly DP1 to DP18.\n"
  )
}


############################################################
# 10. Build sample metadata
############################################################

meta <- data.frame(
  
  sample = paste0(
    "DP",
    1:18
  ),
  
  donor = rep(
    c(
      "D5",
      "D6",
      "D7"
    ),
    times = 6
  ),
  
  group = rep(
    c(
      "Intact",
      "P0",
      "P1",
      "P3_2D",
      "P5_2D",
      "P3_3D"
    ),
    each = 3
  ),
  
  stringsAsFactors = FALSE
)

rownames(meta) <- meta$sample


############################################################
# 11. Set factor order
############################################################

meta$group <- factor(
  
  meta$group,
  
  levels = c(
    "Intact",
    "P0",
    "P1",
    "P3_2D",
    "P5_2D",
    "P3_3D"
  )
)

meta$donor <- factor(
  
  meta$donor,
  
  levels = c(
    "D5",
    "D6",
    "D7"
  )
)


############################################################
# 12. Print metadata
############################################################

cat("\nSample metadata:\n\n")

print(meta)

cat("\n")


############################################################
# 13. Save metadata
############################################################

write.csv(
  meta,
  file = "results/GSE44765_sample_metadata.csv",
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 14. Try reading CEL.gz directly
############################################################

cat("\nReading CEL files...\n")
cat("This may take several minutes.\n\n")

raw_affy <- tryCatch(
  
  {
    
    ReadAffy(
      filenames = cel_files
    )
    
  },
  
  error = function(e) {
    
    cat("\nDirect CEL.gz reading failed.\n")
    cat("Error message:\n")
    cat(conditionMessage(e), "\n\n")
    
    NULL
  }
)


############################################################
# 15. If needed, automatically decompress CEL.gz
############################################################

if (is.null(raw_affy)) {
  
  cat("Trying automatic decompression...\n\n")
  
  temp_cel_dir <- file.path(
    tempdir(),
    "GSE44765_CEL_uncompressed"
  )
  
  dir.create(
    temp_cel_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  uncompressed_files <- character(
    length(cel_files)
  )
  
  for (i in seq_along(cel_files)) {
    
    input_file <- cel_files[i]
    
    output_name <- sub(
      "\\.gz$",
      "",
      basename(input_file),
      ignore.case = TRUE
    )
    
    output_file <- file.path(
      temp_cel_dir,
      output_name
    )
    
    R.utils::gunzip(
      filename = input_file,
      destname = output_file,
      overwrite = TRUE,
      remove = FALSE
    )
    
    uncompressed_files[i] <- output_file
  }
  
  raw_affy <- ReadAffy(
    filenames = uncompressed_files
  )
}


############################################################
# 16. Rename raw samples DP1 -> DP18
############################################################

sampleNames(raw_affy) <- paste0(
  "DP",
  1:18
)

cat("\nSamples successfully loaded:\n\n")

print(
  sampleNames(raw_affy)
)

cat("\n")


############################################################
# 17. Save raw AffyBatch object
############################################################

saveRDS(
  raw_affy,
  file = "results/GSE44765_raw_AffyBatch.rds"
)


############################################################
# 18. RMA normalization
############################################################

cat("\n====================================================\n")
cat("Running RMA normalization...\n")
cat("====================================================\n\n")

eset <- rma(
  raw_affy
)

expr <- exprs(
  eset
)

colnames(expr) <- paste0(
  "DP",
  1:18
)


############################################################
# 19. Check expression matrix
############################################################

cat("\nRMA completed.\n\n")

cat("Expression matrix dimensions:\n")

print(
  dim(expr)
)

cat("\nExpression value summary:\n")

print(
  summary(
    as.vector(expr)
  )
)

cat("\n")

if (ncol(expr) != 18) {
  
  stop(
    "\nERROR: RMA output does not contain 18 samples.\n"
  )
}


############################################################
# 20. Match metadata to expression matrix
############################################################

stopifnot(
  identical(
    colnames(expr),
    rownames(meta)
  )
)

cat(
  "Expression matrix and metadata match correctly.\n\n"
)


############################################################
# 21. Check for NA / infinite values
############################################################

cat(
  "Number of NA values in expression matrix:",
  sum(is.na(expr)),
  "\n"
)

cat(
  "Number of infinite values:",
  sum(!is.finite(expr)),
  "\n\n"
)

if (
  any(is.na(expr)) ||
  any(!is.finite(expr))
) {
  
  stop(
    "\nERROR: Expression matrix contains NA or infinite values.\n"
  )
}


############################################################
# 22. Save RMA expression matrix
############################################################

saveRDS(
  expr,
  file = "results/GSE44765_RMA_probe_expression.rds"
)

saveRDS(
  eset,
  file = "results/GSE44765_RMA_eset.rds"
)


############################################################
# 23. Save compressed CSV expression matrix
############################################################

expr_connection <- gzfile(
  "results/GSE44765_RMA_probe_expression.csv.gz",
  open = "wt"
)

write.csv(
  expr,
  file = expr_connection,
  quote = FALSE
)

close(
  expr_connection
)


############################################################
# 24. RMA boxplot
############################################################

pdf(
  "figures/GSE44765_RMA_boxplot.pdf",
  width = 11,
  height = 6
)

boxplot(
  expr,
  las = 2,
  outline = FALSE,
  main = "GSE44765 after RMA normalization",
  ylab = "log2 expression"
)

dev.off()


############################################################
# 25. Expression density plot
############################################################

pdf(
  "figures/GSE44765_RMA_density.pdf",
  width = 9,
  height = 6
)

limma::plotDensities(
  expr,
  legend = FALSE,
  main = "GSE44765 expression density after RMA"
)

dev.off()


############################################################
# 26. Calculate variance of all probes
############################################################

probe_variance <- apply(
  expr,
  1,
  var,
  na.rm = TRUE
)

probe_variance <- sort(
  probe_variance,
  decreasing = TRUE
)


############################################################
# 27. Select top 5000 variable probes
############################################################

n_top <- min(
  5000,
  length(probe_variance)
)

top_probes <- names(
  probe_variance[
    seq_len(n_top)
  ]
)

cat(
  "\nNumber of probes used for PCA:",
  length(top_probes),
  "\n\n"
)


############################################################
# 28. PCA
############################################################

pca <- prcomp(
  
  t(
    expr[
      top_probes,
      ,
      drop = FALSE
    ]
  ),
  
  center = TRUE,
  
  scale. = FALSE
)


############################################################
# 29. PCA variance explained
############################################################

variance_explained <- 100 *
  pca$sdev^2 /
  sum(
    pca$sdev^2
  )

cat(
  "Variance explained by first five PCs:\n"
)

print(
  round(
    variance_explained[1:5],
    2
  )
)

cat("\n")


############################################################
# 30. Build PCA dataframe
############################################################

pca_df <- data.frame(
  
  sample = rownames(
    pca$x
  ),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  stringsAsFactors = FALSE
)

pca_df$group <- meta[
  pca_df$sample,
  "group"
]

pca_df$donor <- meta[
  pca_df$sample,
  "donor"
]


############################################################
# 31. PCA figure
############################################################

p_pca <- ggplot(
  
  pca_df,
  
  aes(
    x = PC1,
    y = PC2,
    color = group,
    shape = donor,
    label = sample
  )
  
) +
  
  geom_point(
    size = 4
  ) +
  
  ggrepel::geom_text_repel(
    size = 3.2,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  
  labs(
    
    title =
      "GSE44765: PCA after RMA normalization",
    
    subtitle =
      "Top 5,000 most variable probes",
    
    x = paste0(
      "PC1 (",
      round(
        variance_explained[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        variance_explained[2],
        1
      ),
      "%)"
    ),
    
    color =
      "Condition",
    
    shape =
      "Donor"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    legend.position = "right"
  )


############################################################
# 32. Display PCA
############################################################

print(
  p_pca
)


############################################################
# 33. Save PCA PDF
############################################################

ggsave(
  
  filename =
    "figures/GSE44765_PCA_RMA.pdf",
  
  plot =
    p_pca,
  
  width =
    8.5,
  
  height =
    6.5
)


############################################################
# 34. Save PCA PNG
############################################################

ggsave(
  
  filename =
    "figures/GSE44765_PCA_RMA.png",
  
  plot =
    p_pca,
  
  width =
    8.5,
  
  height =
    6.5,
  
  dpi =
    300
)


############################################################
# 35. Save PCA coordinates
############################################################

write.csv(
  
  pca_df,
  
  file =
    "results/GSE44765_PCA_coordinates.csv",
  
  row.names =
    FALSE,
  
  quote =
    FALSE
)


############################################################
# 36. Save PCA variance
############################################################

pca_variance_df <- data.frame(
  
  PC = paste0(
    "PC",
    seq_along(
      variance_explained
    )
  ),
  
  VariancePercent =
    variance_explained
)

write.csv(
  
  pca_variance_df,
  
  file =
    "results/GSE44765_PCA_variance.csv",
  
  row.names =
    FALSE,
  
  quote =
    FALSE
)


############################################################
# 37. Calculate sample correlation matrix
############################################################

sample_cor <- cor(
  
  expr[
    top_probes,
    ,
    drop = FALSE
  ],
  
  method =
    "pearson",
  
  use =
    "pairwise.complete.obs"
)


############################################################
# 38. Check correlation matrix
############################################################

cat("\nCorrelation matrix dimensions:\n")

print(
  dim(sample_cor)
)

cat("\nCorrelation range:\n")

print(
  range(
    sample_cor,
    na.rm = TRUE
  )
)

cat(
  "\nNumber of NA values in correlation matrix:",
  sum(
    is.na(sample_cor)
  ),
  "\n\n"
)

if (
  any(
    is.na(sample_cor)
  )
) {
  
  stop(
    "\nERROR: Correlation matrix contains NA values.\n"
  )
}


############################################################
# 39. Create annotation for heatmap
############################################################

annotation_df <- data.frame(
  
  Group =
    as.character(
      meta$group
    ),
  
  Donor =
    as.character(
      meta$donor
    ),
  
  stringsAsFactors =
    FALSE
)

rownames(
  annotation_df
) <- rownames(
  meta
)


############################################################
# 40. Check annotation
############################################################

annotation_df <- annotation_df[
  colnames(sample_cor),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(annotation_df),
    colnames(sample_cor)
  )
)

cat(
  "Heatmap annotation successfully matched.\n\n"
)


############################################################
# 41. Remove broken old heatmap files if present
############################################################

heatmap_pdf <-
  "figures/GSE44765_sample_correlation_heatmap.pdf"

heatmap_png <-
  "figures/GSE44765_sample_correlation_heatmap.png"

if (
  file.exists(
    heatmap_pdf
  )
) {
  
  file.remove(
    heatmap_pdf
  )
}

if (
  file.exists(
    heatmap_png
  )
) {
  
  file.remove(
    heatmap_png
  )
}


############################################################
# 42. Generate correlation heatmap PDF
#
# IMPORTANT:
# pheatmap now writes directly to filename.
# This avoids the previous 0-byte PDF problem.
############################################################

pheatmap(
  
  sample_cor,
  
  annotation_col =
    annotation_df,
  
  annotation_row =
    annotation_df,
  
  border_color =
    NA,
  
  fontsize =
    9,
  
  fontsize_row =
    8,
  
  fontsize_col =
    8,
  
  main =
    "GSE44765 sample correlation",
  
  filename =
    heatmap_pdf,
  
  width =
    10,
  
  height =
    9
)


############################################################
# 43. Generate correlation heatmap PNG
############################################################

pheatmap(
  
  sample_cor,
  
  annotation_col =
    annotation_df,
  
  annotation_row =
    annotation_df,
  
  border_color =
    NA,
  
  fontsize =
    9,
  
  fontsize_row =
    8,
  
  fontsize_col =
    8,
  
  main =
    "GSE44765 sample correlation",
  
  filename =
    heatmap_png,
  
  width =
    10,
  
  height =
    9
)


############################################################
# 44. Verify heatmap output files
############################################################

pdf_size <- file.info(
  heatmap_pdf
)$size

png_size <- file.info(
  heatmap_png
)$size

cat(
  "\nHeatmap PDF size:",
  pdf_size,
  "bytes\n"
)

cat(
  "Heatmap PNG size:",
  png_size,
  "bytes\n\n"
)

if (
  is.na(pdf_size) ||
  pdf_size == 0
) {
  
  stop(
    "\nERROR: Heatmap PDF is empty.\n"
  )
}

if (
  is.na(png_size) ||
  png_size == 0
) {
  
  stop(
    "\nERROR: Heatmap PNG is empty.\n"
  )
}


############################################################
# 45. Basic QC summary
############################################################

qc_summary <- data.frame(
  
  Sample =
    colnames(expr),
  
  Group =
    as.character(
      meta[
        colnames(expr),
        "group"
      ]
    ),
  
  Donor =
    as.character(
      meta[
        colnames(expr),
        "donor"
      ]
    ),
  
  Median =
    apply(
      expr,
      2,
      median
    ),
  
  Mean =
    colMeans(
      expr
    ),
  
  SD =
    apply(
      expr,
      2,
      sd
    ),
  
  Min =
    apply(
      expr,
      2,
      min
    ),
  
  Max =
    apply(
      expr,
      2,
      max
    )
)


############################################################
# 46. Save QC summary
############################################################

write.csv(
  
  qc_summary,
  
  file =
    "results/GSE44765_basic_QC_summary.csv",
  
  row.names =
    FALSE,
  
  quote =
    FALSE
)


############################################################
# 47. Save correlation matrix
############################################################

write.csv(
  
  sample_cor,
  
  file =
    "results/GSE44765_sample_correlation_matrix.csv",
  
  quote =
    FALSE
)


############################################################
# 48. Save session information
############################################################

sink(
  "results/GSE44765_sessionInfo.txt"
)

sessionInfo()

sink()


############################################################
# 49. Save all important objects together
############################################################

save(
  
  expr,
  meta,
  pca,
  pca_df,
  variance_explained,
  sample_cor,
  
  file =
    "results/GSE44765_analysis_workspace.RData"
)


############################################################
# 50. Final report
############################################################

cat("\n")
cat("====================================================\n")
cat("   GSE44765 preprocessing completed successfully\n")
cat("====================================================\n\n")

cat("Expression matrix dimensions:\n")

print(
  dim(expr)
)

cat("\n")


cat("Expression summary:\n")

print(
  summary(
    as.vector(expr)
  )
)

cat("\n")


cat("First five PC variance percentages:\n")

print(
  round(
    variance_explained[1:5],
    2
  )
)

cat("\n")


cat("Correlation range:\n")

print(
  range(
    sample_cor
  )
)

cat("\n")


cat("Number of NA values in correlation matrix:\n")

print(
  sum(
    is.na(sample_cor)
  )
)

cat("\n")


cat("Heatmap PDF size:\n")

print(
  file.info(
    heatmap_pdf
  )$size
)

cat("\n")


cat("Heatmap PNG size:\n")

print(
  file.info(
    heatmap_png
  )$size
)

cat("\n")


cat("Main output files:\n\n")

cat(
  "results/GSE44765_RMA_probe_expression.rds\n"
)

cat(
  "results/GSE44765_RMA_probe_expression.csv.gz\n"
)

cat(
  "results/GSE44765_sample_metadata.csv\n"
)

cat(
  "results/GSE44765_PCA_coordinates.csv\n"
)

cat(
  "results/GSE44765_PCA_variance.csv\n"
)

cat(
  "results/GSE44765_sample_correlation_matrix.csv\n"
)

cat(
  "results/GSE44765_basic_QC_summary.csv\n"
)

cat(
  "results/GSE44765_analysis_workspace.RData\n"
)

cat("\nFigures:\n\n")

cat(
  "figures/GSE44765_RMA_boxplot.pdf\n"
)

cat(
  "figures/GSE44765_RMA_density.pdf\n"
)

cat(
  "figures/GSE44765_PCA_RMA.pdf\n"
)

cat(
  "figures/GSE44765_PCA_RMA.png\n"
)

cat(
  "figures/GSE44765_sample_correlation_heatmap.pdf\n"
)

cat(
  "figures/GSE44765_sample_correlation_heatmap.png\n"
)

cat("\n")


cat("IMPORTANT:\n")

cat(
  "RMA values are already log2 transformed.\n"
)

cat(
  "DO NOT apply log2(expr + 1) again.\n"
)

cat("\n")

cat(
  "Analysis finished successfully.\n"
)

cat("====================================================\n")