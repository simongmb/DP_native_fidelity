############################################################
# 06_GSE66663_AGA_baseline_validation.R
#
# COMPLETE ONE-CLICK VERSION
#
# Purpose:
#   Validate the predefined GSE44765 strict rescue signature
#   in the GSE66663 BAB/BAN AGA dermal-papilla cell-line model.
#
# Biological question:
#   Does non-balding DP (BAN) show a more rescue-like
#   transcriptional state than balding DP (BAB)?
#
# IMPORTANT LIMITATION:
#   BAB and BAN are immortalized DP cell lines.
#   The three replicates per group are experimental/technical
#   replicates, NOT six independent patients.
#
# Main external contrast:
#   BAN - BAB
#
# Expected direction if the GSE44765 rescue signature is
# concordant with the non-balding state:
#   1) predefined rescue score: BAN > BAB
#   2) Restored_upward genes: enriched UP in BAN
#   3) Restored_downward genes: enriched DOWN in BAN
#   4) Rescue_logFC vs BAN-BAB logFC: positive correlation
#
# Input folder:
#   data_raw/GSE66663/
#
# Required:
#   GSE66663_non-normalized_data.txt.gz
#   one or both GPL10558 annotation files
#
# Optional:
#   six IDAT.gz files may remain in the folder.
#   They are NOT required for this first validation.
#
# Required discovery signature:
#   results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv
#
# Outputs:
#   results/GSE66663_AGA/
#   figures/GSE66663_AGA/
############################################################


############################################################
# 0. Basic settings
############################################################

graphics.off()
options(stringsAsFactors = FALSE)
set.seed(66663)

cat("\n")
cat("============================================================\n")
cat(" GSE66663 AGA BASELINE VALIDATION STARTED\n")
cat("============================================================\n\n")

cat("Working directory:\n")
cat(getwd(), "\n\n")


############################################################
# 1. Paths
############################################################

raw_dir <- "data_raw/GSE66663"

results_dir <- "results/GSE66663_AGA"

figures_dir <- "figures/GSE66663_AGA"

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 2. Locate discovery signature robustly
############################################################

preferred_signature <- file.path(
  "results",
  "GSE44765_rescue",
  "GSE44765_STRICT_rescue_genes.csv"
)

if (file.exists(preferred_signature)) {
  
  discovery_signature_file <- preferred_signature
  
} else {
  
  signature_candidates <- list.files(
    path = "results",
    pattern = "^GSE44765_STRICT_rescue_genes\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(signature_candidates) == 0) {
    
    stop(
      paste0(
        "\nERROR:\n",
        "Cannot find GSE44765_STRICT_rescue_genes.csv under results/.\n",
        "Please confirm that the GSE44765 rescue analysis has already been run.\n"
      )
    )
  }
  
  discovery_signature_file <- signature_candidates[1]
}

cat(
  "Discovery signature file:\n",
  discovery_signature_file,
  "\n\n"
)


############################################################
# 3. Install/load packages
############################################################

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap"
)

for (pkg in cran_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace(
  "limma",
  quietly = TRUE
)) {
  
  BiocManager::install(
    "limma",
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  
  library(limma)
  
  library(ggplot2)
  
  library(ggrepel)
  
  library(pheatmap)
  
})

cat("Packages loaded successfully.\n\n")


############################################################
# 4. Locate GSE66663 non-normalized matrix
############################################################

non_normalized_candidates <- list.files(
  path = raw_dir,
  pattern = "^GSE66663_non-normalized_data\\.txt\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(non_normalized_candidates) != 1) {
  
  stop(
    paste0(
      "\nERROR:\n",
      "Expected exactly one file named\n",
      "GSE66663_non-normalized_data.txt.gz\n",
      "inside:\n",
      raw_dir,
      "\n"
    )
  )
}

non_normalized_file <- non_normalized_candidates[1]

cat(
  "Non-normalized matrix:\n",
  non_normalized_file,
  "\n\n"
)


############################################################
# 5. Read GSE66663 non-normalized data
############################################################

raw_table <- read.delim(
  non_normalized_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

if (!"ILMN_ID" %in% colnames(raw_table)) {
  
  stop(
    "ERROR: ILMN_ID column is missing from the GSE66663 matrix."
  )
}

cat(
  "Rows in supplied non-normalized matrix:",
  nrow(raw_table),
  "\n"
)

cat(
  "Columns in supplied non-normalized matrix:",
  ncol(raw_table),
  "\n\n"
)


############################################################
# 6. Expected six signal columns
############################################################

expected_signal_columns <- c(
  "AVG_Signal.BAB_repA",
  "AVG_Signal.BAB_repB",
  "AVG_Signal.BAB_repC",
  "AVG_Signal.BAN_repA",
  "AVG_Signal.BAN_repB",
  "AVG_Signal.BAN_repC"
)

missing_signal_columns <- setdiff(
  expected_signal_columns,
  colnames(raw_table)
)

if (length(missing_signal_columns) > 0) {
  
  stop(
    paste0(
      "ERROR: missing expected signal columns:\n",
      paste(
        missing_signal_columns,
        collapse = "\n"
      )
    )
  )
}


############################################################
# 7. Build sample metadata
############################################################

meta <- data.frame(
  
  GSM = c(
    "GSM1627296",
    "GSM1627297",
    "GSM1627298",
    "GSM1627299",
    "GSM1627300",
    "GSM1627301"
  ),
  
  Sample = c(
    "BAB_repA",
    "BAB_repB",
    "BAB_repC",
    "BAN_repA",
    "BAN_repB",
    "BAN_repC"
  ),
  
  Group = c(
    "BAB",
    "BAB",
    "BAB",
    "BAN",
    "BAN",
    "BAN"
  ),
  
  SignalColumn = expected_signal_columns,
  
  stringsAsFactors = FALSE
)

meta$Group <- factor(
  meta$Group,
  levels = c(
    "BAB",
    "BAN"
  )
)

rownames(meta) <- meta$Sample

write.csv(
  meta,
  file.path(
    results_dir,
    "GSE66663_sample_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat("Sample design:\n")
print(meta)
cat("\n")


############################################################
# 8. Build raw intensity matrix
############################################################

raw_signal <- as.matrix(
  raw_table[
    ,
    expected_signal_columns,
    drop = FALSE
  ]
)

storage.mode(raw_signal) <- "numeric"

rownames(raw_signal) <- raw_table$ILMN_ID

colnames(raw_signal) <- meta$Sample

if (
  any(!is.finite(raw_signal)) ||
  any(raw_signal <= 0)
) {
  
  stop(
    "ERROR: non-positive or non-finite signal values detected."
  )
}


############################################################
# 9. Save raw intensity QC
############################################################

raw_qc <- data.frame(
  
  Sample = colnames(raw_signal),
  
  Group = as.character(meta$Group),
  
  MinSignal = apply(
    raw_signal,
    2,
    min,
    na.rm = TRUE
  ),
  
  MedianSignal = apply(
    raw_signal,
    2,
    median,
    na.rm = TRUE
  ),
  
  MeanSignal = apply(
    raw_signal,
    2,
    mean,
    na.rm = TRUE
  ),
  
  MaxSignal = apply(
    raw_signal,
    2,
    max,
    na.rm = TRUE
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  raw_qc,
  file.path(
    results_dir,
    "GSE66663_raw_signal_QC.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 10. Log2 transformation + quantile normalization
#
# GEO states that processed expression was generated using
# log2 transformation and quantile normalization.
############################################################

log2_signal <- log2(raw_signal)

expr_probe <- limma::normalizeBetweenArrays(
  log2_signal,
  method = "quantile"
)

saveRDS(
  expr_probe,
  file.path(
    results_dir,
    "GSE66663_probe_level_log2_quantile_expression.rds"
  )
)


############################################################
# 11. Helper: read only the [Probes] section of GPL10558
############################################################

read_gpl10558_probes <- function(filename) {
  
  con <- gzfile(
    filename,
    open = "rt"
  )
  
  on.exit(
    close(con),
    add = TRUE
  )
  
  all_lines <- readLines(
    con,
    warn = FALSE
  )
  
  probe_section <- which(
    trimws(all_lines) == "[Probes]"
  )
  
  control_section <- which(
    trimws(all_lines) == "[Controls]"
  )
  
  if (
    length(probe_section) == 0 ||
    length(control_section) == 0
  ) {
    
    stop(
      paste0(
        "ERROR: cannot locate [Probes]/[Controls] in ",
        basename(filename)
      )
    )
  }
  
  probe_start <- probe_section[1] + 1
  
  control_start <- control_section[
    control_section > probe_section[1]
  ][1]
  
  probe_end <- control_start - 1
  
  probe_lines <- all_lines[
    probe_start:probe_end
  ]
  
  temp_file <- tempfile(
    fileext = ".txt"
  )
  
  writeLines(
    probe_lines,
    temp_file
  )
  
  probe_table <- read.delim(
    temp_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  
  unlink(temp_file)
  
  required_columns <- c(
    "Probe_Id",
    "Symbol"
  )
  
  if (!all(
    required_columns %in%
    colnames(probe_table)
  )) {
    
    stop(
      paste0(
        "ERROR: annotation file ",
        basename(filename),
        " lacks Probe_Id/Symbol columns."
      )
    )
  }
  
  probe_table
}


############################################################
# 12. Locate GPL10558 annotation candidates
############################################################

annotation_candidates <- list.files(
  path = raw_dir,
  pattern = "^GPL10558.*\\.txt\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(annotation_candidates) == 0) {
  
  stop(
    paste0(
      "\nERROR:\n",
      "No GPL10558 annotation text file found in:\n",
      raw_dir,
      "\n"
    )
  )
}

cat(
  "GPL10558 annotation candidates:",
  length(annotation_candidates),
  "\n"
)


############################################################
# 13. Automatically choose annotation with best probe match
############################################################

annotation_tables <- list()

annotation_match_rows <- list()

annotation_counter <- 1

for (annotation_file in annotation_candidates) {
  
  current_annotation <- read_gpl10558_probes(
    annotation_file
  )
  
  n_overlap <- length(
    intersect(
      rownames(expr_probe),
      current_annotation$Probe_Id
    )
  )
  
  annotation_tables[
    basename(annotation_file)
  ] <- list(
    current_annotation
  )
  
  annotation_match_rows[
    annotation_counter
  ] <- list(
    data.frame(
      
      AnnotationFile = basename(annotation_file),
      
      AnnotationProbeRows = nrow(current_annotation),
      
      ExpressionProbeRows = nrow(expr_probe),
      
      MatchingProbeIDs = n_overlap,
      
      MatchFraction = n_overlap / nrow(expr_probe),
      
      stringsAsFactors = FALSE
    )
  )
  
  annotation_counter <- annotation_counter + 1
}

annotation_match_table <- do.call(
  rbind,
  annotation_match_rows
)

annotation_match_table <- annotation_match_table[
  order(
    -annotation_match_table$MatchingProbeIDs
  ),
  ,
  drop = FALSE
]

write.csv(
  annotation_match_table,
  file.path(
    results_dir,
    "GSE66663_annotation_probe_match.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

best_annotation_name <- annotation_match_table$AnnotationFile[1]

annotation <- getElement(
  annotation_tables,
  best_annotation_name
)

cat("\nAnnotation match:\n")
print(annotation_match_table)

cat(
  "\nSelected annotation:",
  best_annotation_name,
  "\n\n"
)


############################################################
# 14. Map probes to gene symbols
############################################################

annotation_small <- annotation[
  ,
  c(
    "Probe_Id",
    "Symbol"
  ),
  drop = FALSE
]

annotation_small$Probe_Id <- as.character(
  annotation_small$Probe_Id
)

annotation_small$Symbol <- trimws(
  as.character(
    annotation_small$Symbol
  )
)

annotation_small <- annotation_small[
  !is.na(annotation_small$Symbol) &
    annotation_small$Symbol != "",
  ,
  drop = FALSE
]

annotation_small <- annotation_small[
  !duplicated(annotation_small$Probe_Id),
  ,
  drop = FALSE
]

probe_to_symbol <- annotation_small$Symbol

names(probe_to_symbol) <- annotation_small$Probe_Id

common_probe_ids <- intersect(
  rownames(expr_probe),
  names(probe_to_symbol)
)

expr_mapped <- expr_probe[
  common_probe_ids,
  ,
  drop = FALSE
]

gene_symbol <- unname(
  probe_to_symbol[
    common_probe_ids
  ]
)

cat(
  "Expression probes:",
  nrow(expr_probe),
  "\n"
)

cat(
  "Mapped probes with non-empty symbols:",
  nrow(expr_mapped),
  "\n"
)


############################################################
# 15. Collapse multiple probes to gene level
############################################################

expr_gene <- limma::avereps(
  expr_mapped,
  ID = gene_symbol
)

expr_gene <- expr_gene[
  !is.na(rownames(expr_gene)) &
    rownames(expr_gene) != "",
  ,
  drop = FALSE
]

cat(
  "Gene-level features after collapsing probes:",
  nrow(expr_gene),
  "\n\n"
)

saveRDS(
  expr_gene,
  file.path(
    results_dir,
    "GSE66663_gene_level_expression.rds"
  )
)


############################################################
# 16. Gene-level expression CSV
############################################################

expression_output <- data.frame(
  GeneSymbol = rownames(expr_gene),
  expr_gene,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expression_connection <- gzfile(
  file.path(
    results_dir,
    "GSE66663_gene_level_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  expression_output,
  expression_connection,
  row.names = FALSE,
  quote = FALSE
)

close(expression_connection)


############################################################
# 17. PCA using top variable genes
############################################################

gene_variance <- apply(
  expr_gene,
  1,
  var,
  na.rm = TRUE
)

gene_variance <- sort(
  gene_variance,
  decreasing = TRUE
)

n_pca_genes <- min(
  5000,
  length(gene_variance)
)

pca_genes <- names(
  gene_variance[
    seq_len(n_pca_genes)
  ]
)

pca <- prcomp(
  t(
    expr_gene[
      pca_genes,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- 100 *
  pca$sdev^2 /
  sum(pca$sdev^2)

pca_df <- data.frame(
  
  Sample = rownames(pca$x),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  stringsAsFactors = FALSE
)

pca_df$Group <- meta[
  pca_df$Sample,
  "Group"
]

p_pca <- ggplot(
  
  pca_df,
  
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    label = Sample
  )
  
) +
  
  geom_point(
    size = 4
  ) +
  
  ggrepel::geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  labs(
    
    title =
      "GSE66663: BAB versus BAN dermal papilla cell lines",
    
    subtitle =
      "Top 5,000 variable genes; log2 + quantile normalization",
    
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
    
    color = "Cell line"
  ) +
  
  theme_classic(
    base_size = 13
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66663_PCA.png"
  ),
  p_pca,
  width = 8,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66663_PCA.pdf"
  ),
  p_pca,
  width = 8,
  height = 6.5
)


############################################################
# 18. Sample correlation heatmap
############################################################

sample_correlation <- cor(
  expr_gene[
    pca_genes,
    ,
    drop = FALSE
  ],
  method = "pearson"
)

annotation_col <- data.frame(
  Group = as.character(meta$Group),
  stringsAsFactors = FALSE
)

rownames(annotation_col) <- meta$Sample

pheatmap(
  sample_correlation,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  border_color = NA,
  main = "GSE66663 sample correlation",
  filename = file.path(
    figures_dir,
    "GSE66663_sample_correlation_heatmap.png"
  ),
  width = 7.5,
  height = 7
)


############################################################
# 19. Differential expression: BAN - BAB
############################################################

design <- model.matrix(
  ~ 0 + Group,
  data = meta
)

colnames(design) <- sub(
  "^Group",
  "",
  colnames(design)
)

contrast_matrix <- limma::makeContrasts(
  BAN_vs_BAB = BAN - BAB,
  levels = design
)

fit <- limma::lmFit(
  expr_gene,
  design
)

fit2 <- limma::contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- limma::eBayes(
  fit2,
  trend = TRUE,
  robust = TRUE
)

de_all <- limma::topTable(
  fit2,
  coef = "BAN_vs_BAB",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

de_all$GeneSymbol <- rownames(de_all)

de_all <- de_all[
  ,
  c(
    "GeneSymbol",
    setdiff(
      colnames(de_all),
      "GeneSymbol"
    )
  ),
  drop = FALSE
]

write.csv(
  de_all,
  file.path(
    results_dir,
    "GSE66663_BAN_vs_BAB_all_genes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 20. Load frozen GSE44765 strict rescue signature
############################################################

signature <- read.csv(
  discovery_signature_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_signature_columns <- c(
  "GeneSymbol",
  "Rescue_type",
  "Rescue_logFC"
)

if (!all(
  required_signature_columns %in%
  colnames(signature)
)) {
  
  stop(
    paste0(
      "ERROR: discovery signature must contain:\n",
      paste(
        required_signature_columns,
        collapse = ", "
      )
    )
  )
}

signature$GeneSymbol <- trimws(
  as.character(
    signature$GeneSymbol
  )
)

signature <- signature[
  !is.na(signature$GeneSymbol) &
    signature$GeneSymbol != "",
  ,
  drop = FALSE
]

signature <- signature[
  !duplicated(signature$GeneSymbol),
  ,
  drop = FALSE
]


############################################################
# 21. Signature overlap
############################################################

signature_overlap <- signature[
  signature$GeneSymbol %in%
    rownames(expr_gene),
  ,
  drop = FALSE
]

up_genes <- signature_overlap$GeneSymbol[
  signature_overlap$Rescue_type ==
    "Restored_upward"
]

down_genes <- signature_overlap$GeneSymbol[
  signature_overlap$Rescue_type ==
    "Restored_downward"
]

signature_coverage <- nrow(signature_overlap) /
  nrow(signature)

write.csv(
  signature_overlap,
  file.path(
    results_dir,
    "GSE44765_signature_overlap_in_GSE66663.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Discovery strict signature genes:",
  nrow(signature),
  "\n"
)

cat(
  "Signature genes represented in GSE66663:",
  nrow(signature_overlap),
  "\n"
)

cat(
  "Signature coverage:",
  round(
    100 * signature_coverage,
    1
  ),
  "%\n"
)

cat(
  "Restored-upward overlap:",
  length(up_genes),
  "\n"
)

cat(
  "Restored-downward overlap:",
  length(down_genes),
  "\n\n"
)


############################################################
# 22. Predefined rescue score
#
# Score = mean Z(restored-upward) - mean Z(restored-downward)
#
# The score is descriptive at the replicate level.
# Because BAB/BAN replicates are not independent patients,
# do not interpret replicate-level P values as patient-level
# evidence.
############################################################

signature_expr <- expr_gene[
  signature_overlap$GeneSymbol,
  ,
  drop = FALSE
]

signature_z <- t(
  scale(
    t(signature_expr)
  )
)

signature_z[
  !is.finite(signature_z)
] <- 0

up_score <- colMeans(
  signature_z[
    up_genes,
    ,
    drop = FALSE
  ]
)

down_score <- colMeans(
  signature_z[
    down_genes,
    ,
    drop = FALSE
  ]
)

rescue_score <- up_score - down_score

score_df <- data.frame(
  
  Sample = colnames(expr_gene),
  
  Group = as.character(
    meta[
      colnames(expr_gene),
      "Group"
    ]
  ),
  
  Up_score = as.numeric(
    up_score[
      colnames(expr_gene)
    ]
  ),
  
  Down_score = as.numeric(
    down_score[
      colnames(expr_gene)
    ]
  ),
  
  Rescue_signature_score = as.numeric(
    rescue_score[
      colnames(expr_gene)
    ]
  ),
  
  stringsAsFactors = FALSE
)

score_df$Group <- factor(
  score_df$Group,
  levels = c(
    "BAB",
    "BAN"
  )
)

write.csv(
  score_df,
  file.path(
    results_dir,
    "GSE66663_predefined_rescue_scores.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 23. Group-level score summary
############################################################

score_summary <- aggregate(
  Rescue_signature_score ~ Group,
  data = score_df,
  FUN = function(x) {
    c(
      Mean = mean(x),
      SD = sd(x),
      Median = median(x)
    )
  }
)

score_summary_flat <- data.frame(
  
  Group = score_summary$Group,
  
  Mean = score_summary$Rescue_signature_score[
    ,
    "Mean"
  ],
  
  SD = score_summary$Rescue_signature_score[
    ,
    "SD"
  ],
  
  Median = score_summary$Rescue_signature_score[
    ,
    "Median"
  ],
  
  stringsAsFactors = FALSE
)

write.csv(
  score_summary_flat,
  file.path(
    results_dir,
    "GSE66663_predefined_rescue_score_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

bab_mean <- mean(
  score_df$Rescue_signature_score[
    score_df$Group == "BAB"
  ]
)

ban_mean <- mean(
  score_df$Rescue_signature_score[
    score_df$Group == "BAN"
  ]
)

ban_minus_bab_score <- ban_mean - bab_mean


############################################################
# 24. Rescue score plot
############################################################

p_score <- ggplot(
  
  score_df,
  
  aes(
    x = Group,
    y = Rescue_signature_score,
    color = Group,
    label = Sample
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_point(
    size = 4
  ) +
  
  ggrepel::geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  
  labs(
    
    title =
      "Predefined GSE44765 rescue signature in GSE66663",
    
    subtitle =
      "Higher score = more GSE44765 rescue-like direction",
    
    x = NULL,
    
    y =
      "Predefined rescue signature score"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66663_predefined_rescue_score.png"
  ),
  p_score,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66663_predefined_rescue_score.pdf"
  ),
  p_score,
  width = 7,
  height = 6
)


############################################################
# 25. CAMERA validation of the two frozen signature arms
############################################################

signature_index <- list(
  
  Restored_upward = which(
    rownames(expr_gene) %in%
      up_genes
  ),
  
  Restored_downward = which(
    rownames(expr_gene) %in%
      down_genes
  )
)

camera_results <- limma::camera(
  
  y = expr_gene,
  
  index = signature_index,
  
  design = design,
  
  contrast = contrast_matrix[
    ,
    "BAN_vs_BAB"
  ],
  
  inter.gene.cor = 0.01
)

camera_results$GeneSet <- rownames(
  camera_results
)

camera_results$Contrast <- "BAN_vs_BAB"

camera_results <- camera_results[
  ,
  c(
    "Contrast",
    "GeneSet",
    setdiff(
      colnames(camera_results),
      c(
        "Contrast",
        "GeneSet"
      )
    )
  ),
  drop = FALSE
]

write.csv(
  camera_results,
  file.path(
    results_dir,
    "GSE66663_CAMERA_predefined_signature.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 26. Cross-study effect concordance
#
# Discovery:
#   GSE44765 Rescue_logFC = P3_3D - P3_2D
#
# External:
#   GSE66663 logFC = BAN - BAB
############################################################

external_effect <- de_all[
  ,
  c(
    "GeneSymbol",
    "logFC",
    "t",
    "P.Value",
    "adj.P.Val"
  ),
  drop = FALSE
]

colnames(external_effect) <- c(
  "GeneSymbol",
  "External_logFC_BAN_minus_BAB",
  "External_t",
  "External_P",
  "External_FDR"
)

discovery_effect <- signature_overlap[
  ,
  c(
    "GeneSymbol",
    "Rescue_type",
    "Rescue_logFC"
  ),
  drop = FALSE
]

merged_effects <- merge(
  discovery_effect,
  external_effect,
  by = "GeneSymbol"
)

valid_effect <- is.finite(
  merged_effects$Rescue_logFC
) &
  is.finite(
    merged_effects$External_logFC_BAN_minus_BAB
  )

merged_effects$Directional_concordance <- sign(
  merged_effects$Rescue_logFC
) ==
  sign(
    merged_effects$External_logFC_BAN_minus_BAB
  )

merged_effects$Expected_signed_external_logFC <- sign(
  merged_effects$Rescue_logFC
) *
  merged_effects$External_logFC_BAN_minus_BAB

spearman_rho <- suppressWarnings(
  cor(
    merged_effects$Rescue_logFC[
      valid_effect
    ],
    merged_effects$External_logFC_BAN_minus_BAB[
      valid_effect
    ],
    method = "spearman"
  )
)

directional_fraction <- mean(
  merged_effects$Directional_concordance[
    valid_effect
  ],
  na.rm = TRUE
)

median_expected_signed <- median(
  merged_effects$Expected_signed_external_logFC[
    valid_effect
  ],
  na.rm = TRUE
)

mean_expected_signed <- mean(
  merged_effects$Expected_signed_external_logFC[
    valid_effect
  ],
  na.rm = TRUE
)

concordance_summary <- data.frame(
  
  Contrast = "BAN_vs_BAB",
  
  N_signature_genes = sum(valid_effect),
  
  Spearman_rho = spearman_rho,
  
  Directional_concordance = directional_fraction,
  
  Median_expected_signed_external_logFC =
    median_expected_signed,
  
  Mean_expected_signed_external_logFC =
    mean_expected_signed,
  
  stringsAsFactors = FALSE
)

write.csv(
  merged_effects,
  file.path(
    results_dir,
    "GSE66663_cross_study_gene_effects.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  concordance_summary,
  file.path(
    results_dir,
    "GSE66663_cross_study_concordance_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 27. Cross-study concordance plot
############################################################

plot_effects <- merged_effects[
  valid_effect,
  ,
  drop = FALSE
]

p_concordance <- ggplot(
  
  plot_effects,
  
  aes(
    x = Rescue_logFC,
    y = External_logFC_BAN_minus_BAB
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.35
  ) +
  
  geom_point(
    alpha = 0.35,
    size = 1.2
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  
  labs(
    
    title =
      "Cross-study effect concordance: GSE44765 vs GSE66663",
    
    subtitle = paste0(
      "Spearman rho = ",
      round(
        spearman_rho,
        3
      ),
      "; directional concordance = ",
      round(
        100 * directional_fraction,
        1
      ),
      "%"
    ),
    
    x =
      "GSE44765 rescue log2FC (P3 3D - P3 2D)",
    
    y =
      "GSE66663 log2FC (BAN - BAB)"
  ) +
  
  theme_classic(
    base_size = 13
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66663_cross_study_effect_concordance.png"
  ),
  p_concordance,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66663_cross_study_effect_concordance.pdf"
  ),
  p_concordance,
  width = 7.5,
  height = 6.5
)


############################################################
# 28. Predefined signature heatmap
#
# IMPORTANT:
# Genes are selected using discovery information only.
# GSE66663 performance is NOT used to choose genes.
############################################################

if (
  "Combined_FDR" %in%
  colnames(signature_overlap)
) {
  
  signature_ranked <- signature_overlap[
    order(
      signature_overlap$Combined_FDR,
      -abs(
        signature_overlap$Rescue_logFC
      ),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
} else {
  
  signature_ranked <- signature_overlap[
    order(
      -abs(
        signature_overlap$Rescue_logFC
      )
    ),
    ,
    drop = FALSE
  ]
}

n_heatmap_genes <- min(
  50,
  nrow(signature_ranked)
)

heatmap_genes <- unique(
  signature_ranked$GeneSymbol[
    seq_len(n_heatmap_genes)
  ]
)

heatmap_expr <- expr_gene[
  heatmap_genes,
  ,
  drop = FALSE
]

heatmap_z <- t(
  scale(
    t(heatmap_expr)
  )
)

heatmap_z[
  !is.finite(heatmap_z)
] <- 0

pheatmap(
  heatmap_z,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  main =
    "Predefined GSE44765 strict rescue signature in GSE66663",
  filename = file.path(
    figures_dir,
    "GSE66663_predefined_signature_heatmap.png"
  ),
  width = 8.5,
  height = 11
)

pheatmap(
  heatmap_z,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  main =
    "Predefined GSE44765 strict rescue signature in GSE66663",
  filename = file.path(
    figures_dir,
    "GSE66663_predefined_signature_heatmap.pdf"
  ),
  width = 8.5,
  height = 11
)


############################################################
# 29. Signature-arm mean external effects
############################################################

arm_summary_list <- list()

arm_names <- c(
  "Restored_upward",
  "Restored_downward"
)

for (i in seq_along(arm_names)) {
  
  current_arm <- arm_names[i]
  
  values <- merged_effects$External_logFC_BAN_minus_BAB[
    merged_effects$Rescue_type ==
      current_arm
  ]
  
  arm_summary_list[
    i
  ] <- list(
    data.frame(
      
      Rescue_type = current_arm,
      
      N = length(values),
      
      Mean_external_logFC = mean(
        values,
        na.rm = TRUE
      ),
      
      Median_external_logFC = median(
        values,
        na.rm = TRUE
      ),
      
      Fraction_external_up = mean(
        values > 0,
        na.rm = TRUE
      ),
      
      Fraction_external_down = mean(
        values < 0,
        na.rm = TRUE
      ),
      
      stringsAsFactors = FALSE
    )
  )
}

arm_summary <- do.call(
  rbind,
  arm_summary_list
)

write.csv(
  arm_summary,
  file.path(
    results_dir,
    "GSE66663_signature_arm_external_effect_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 30. Human-readable analysis summary
############################################################

summary_file <- file.path(
  results_dir,
  "ANALYSIS_SUMMARY_GSE66663.txt"
)

sink(summary_file)

cat(
  "GSE66663 AGA baseline validation\n"
)

cat(
  "====================================\n\n"
)

cat(
  "Biological model:\n"
)

cat(
  "BAB = immortalized balding dermal papilla cell line\n"
)

cat(
  "BAN = immortalized non-balding dermal papilla cell line\n"
)

cat(
  "Replicates are not independent patients.\n\n"
)

cat(
  "Samples analyzed: 6\n"
)

cat(
  "BAB replicates: 3\n"
)

cat(
  "BAN replicates: 3\n\n"
)

cat(
  "Supplied expression probes:",
  nrow(raw_signal),
  "\n"
)

cat(
  "Selected GPL10558 annotation:",
  best_annotation_name,
  "\n"
)

cat(
  "Annotation probe matches:",
  annotation_match_table$MatchingProbeIDs[1],
  "/",
  nrow(raw_signal),
  "\n"
)

cat(
  "Mapped gene-level features:",
  nrow(expr_gene),
  "\n\n"
)

cat(
  "Discovery strict rescue genes:",
  nrow(signature),
  "\n"
)

cat(
  "Signature genes represented in GSE66663:",
  nrow(signature_overlap),
  "\n"
)

cat(
  "Signature coverage:",
  round(
    100 * signature_coverage,
    1
  ),
  "%\n"
)

cat(
  "Restored-upward overlap:",
  length(up_genes),
  "\n"
)

cat(
  "Restored-downward overlap:",
  length(down_genes),
  "\n\n"
)

cat(
  "Predefined rescue score by sample:\n"
)

print(score_df)

cat("\n")

cat(
  "Predefined rescue score by group:\n"
)

print(score_summary_flat)

cat("\n")

cat(
  "BAN - BAB mean score difference:",
  round(
    ban_minus_bab_score,
    4
  ),
  "\n\n"
)

cat(
  "CAMERA predefined signature validation:\n"
)

print(camera_results)

cat("\n")

cat(
  "Cross-study effect concordance:\n"
)

print(concordance_summary)

cat("\n")

cat(
  "Signature-arm external effect summary:\n"
)

print(arm_summary)

cat("\n")

cat(
  "Interpretation rule:\n"
)

cat(
  "A positive BAN-BAB score difference, positive cross-study correlation,\n"
)

cat(
  "Restored_upward enrichment UP in BAN, and Restored_downward enrichment DOWN\n"
)

cat(
  "in BAN would support concordance between the GSE44765 rescue direction\n"
)

cat(
  "and the non-balding DP cell-line state.\n\n"
)

cat(
  "IMPORTANT LIMITATION:\n"
)

cat(
  "Do not interpret the three BAB/BAN replicates as six independent patients.\n"
)

cat(
  "This dataset provides cell-line-model evidence, not patient-cohort validation.\n"
)

sink()


############################################################
# 31. Save R workspace and session information
############################################################

save(
  
  raw_table,
  raw_signal,
  expr_probe,
  annotation_match_table,
  annotation,
  expr_gene,
  meta,
  pca,
  de_all,
  signature,
  signature_overlap,
  score_df,
  score_summary_flat,
  camera_results,
  concordance_summary,
  merged_effects,
  arm_summary,
  
  file = file.path(
    results_dir,
    "GSE66663_AGA_validation_workspace.RData"
  )
)

sink(
  file.path(
    results_dir,
    "sessionInfo_GSE66663.txt"
  )
)

sessionInfo()

sink()


############################################################
# 32. Final console report
############################################################

cat("\n")
cat("============================================================\n")
cat(" GSE66663 AGA BASELINE VALIDATION COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "Selected annotation:",
  best_annotation_name,
  "\n"
)

cat(
  "Gene-level features:",
  nrow(expr_gene),
  "\n"
)

cat(
  "Signature coverage:",
  nrow(signature_overlap),
  "/",
  nrow(signature),
  " (",
  round(
    100 * signature_coverage,
    1
  ),
  "%)\n",
  sep = ""
)

cat("\nRescue score summary:\n")
print(score_summary_flat)

cat(
  "\nBAN - BAB mean rescue-score difference:",
  round(
    ban_minus_bab_score,
    4
  ),
  "\n"
)

cat("\nCAMERA:\n")
print(camera_results)

cat("\nCross-study concordance:\n")
print(concordance_summary)

cat("\nMain figures:\n")

cat(
  file.path(
    figures_dir,
    "GSE66663_PCA.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66663_predefined_rescue_score.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66663_cross_study_effect_concordance.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66663_predefined_signature_heatmap.png"
  ),
  "\n"
)

cat("\nMain summary:\n")
cat(summary_file, "\n")

cat("\nIMPORTANT:\n")
cat(
  "GSE66663 is a BAB/BAN cell-line model, not a patient cohort.\n"
)

cat("============================================================\n")
