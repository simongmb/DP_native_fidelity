############################################################
# GSE293743 independent external validation
# COMPLETE ONE-CLICK VERSION
#
# Project:
#   HairRegeneration_AI
#
# Purpose:
#   Validate the GSE44765-derived strict rescue signature
#   in an entirely independent Clariom S dataset.
#
# Primary samples:
#   15 human dermal papilla cell arrays
#   5 conditions x 3 replicate blocks:
#     Control
#     DPAC
#     CHIR
#     FGF9
#     CHIR_FGF9
#
# GSM8889398 (Fibro) is intentionally excluded from the
# PRIMARY validation because it is a different cell type.
# It can be used later as a supplementary outgroup.
#
# Required raw data:
#   data_raw/GSE293743/*.CEL.gz
#
# Required discovery result:
#   results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv
#
# Main workflow:
#   Clariom S CEL -> oligo RMA -> gene annotation -> QC
#   -> predefined signature score
#   -> paired/block-aware limma contrasts
#   -> CAMERA gene-set validation
#   -> cross-study effect concordance
#
############################################################


############################################################
# 0. Start clean
############################################################

graphics.off()
options(stringsAsFactors = FALSE)
set.seed(293743)

cat("\n")
cat("============================================================\n")
cat(" GSE293743 independent external validation started\n")
cat("============================================================\n\n")
cat("Working directory:\n")
cat(getwd(), "\n\n")


############################################################
# 1. Paths
############################################################

raw_dir <-
  "data_raw/GSE293743"

discovery_signature_file <-
  "results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv"

results_dir <-
  "results/GSE293743_external"

figures_dir <-
  "figures/GSE293743_external"

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
# 2. Install required CRAN packages
############################################################

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "R.utils"
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


############################################################
# 3. Install Bioconductor manager
############################################################

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


############################################################
# 4. Install Bioconductor packages
############################################################

bioc_packages <- c(
  "oligo",
  "pd.clariom.s.human",
  "clariomshumantranscriptcluster.db",
  "AnnotationDbi",
  "limma"
)

for (pkg in bioc_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
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
  
  library(oligo)
  
  library(pd.clariom.s.human)
  
  library(clariomshumantranscriptcluster.db)
  
  library(AnnotationDbi)
  
  library(limma)
  
  library(ggplot2)
  
  library(ggrepel)
  
  library(pheatmap)
  
  library(R.utils)
  
})

cat("All required packages loaded successfully.\n\n")


############################################################
# 6. Check discovery signature exists
############################################################

if (!file.exists(
  discovery_signature_file
)) {
  
  stop(
    paste0(
      "\nERROR:\nCannot find discovery signature:\n",
      discovery_signature_file,
      "\nPlease run GSE44765 Script 02 first.\n"
    )
  )
}


############################################################
# 7. Expected GSE293743 sample structure
############################################################

expected_gsm <- c(
  
  # Replicate block 1
  "GSM8889383",
  "GSM8889384",
  "GSM8889385",
  "GSM8889386",
  "GSM8889387",
  
  # Replicate block 2
  "GSM8889388",
  "GSM8889389",
  "GSM8889390",
  "GSM8889391",
  "GSM8889392",
  
  # Replicate block 3
  "GSM8889393",
  "GSM8889394",
  "GSM8889395",
  "GSM8889396",
  "GSM8889397",
  
  # Fibroblast outgroup
  "GSM8889398"
)


############################################################
# 8. Locate CEL/CEL.gz files
############################################################

all_cel_files <- list.files(
  path = raw_dir,
  pattern = "\\.CEL(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat(
  "CEL/CEL.gz files detected:",
  length(all_cel_files),
  "\n"
)

if (length(all_cel_files) == 0) {
  
  stop(
    paste0(
      "\nERROR:\nNo CEL files were found in:\n",
      raw_dir,
      "\n\nCopy all 16 GSE293743 CEL.gz files into that folder.\n"
    )
  )
}


############################################################
# 9. Map each GSM accession to exactly one CEL file
############################################################

find_file_for_gsm <- function(
    gsm,
    file_vector
) {
  
  hits <- file_vector[
    grepl(
      gsm,
      basename(file_vector),
      fixed = TRUE
    )
  ]
  
  if (length(hits) != 1) {
    
    stop(
      paste0(
        "\nERROR: expected exactly one file for ",
        gsm,
        " but found ",
        length(hits),
        ".\n"
      )
    )
  }
  
  hits
}


ordered_files <- vapply(
  expected_gsm,
  find_file_for_gsm,
  file_vector = all_cel_files,
  FUN.VALUE = character(1)
)

names(ordered_files) <-
  expected_gsm


############################################################
# 10. Save raw-file manifest
############################################################

manifest <- data.frame(
  
  GSM =
    expected_gsm,
  
  File =
    basename(
      ordered_files
    ),
  
  SizeBytes =
    as.numeric(
      file.info(
        ordered_files
      )$size
    ),
  
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(
    results_dir,
    "GSE293743_raw_file_manifest.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat("\nRaw-file manifest:\n")
print(manifest)
cat("\n")


############################################################
# 11. Build metadata
############################################################

condition_order <- c(
  "Control",
  "DPAC",
  "CHIR",
  "FGF9",
  "CHIR_FGF9"
)

meta_dpc <- data.frame(
  
  GSM =
    expected_gsm[1:15],
  
  sample =
    c(
      paste0(
        condition_order,
        "_R1"
      ),
      paste0(
        condition_order,
        "_R2"
      ),
      paste0(
        condition_order,
        "_R3"
      )
    ),
  
  condition =
    rep(
      condition_order,
      times = 3
    ),
  
  replicate =
    rep(
      c(
        "R1",
        "R2",
        "R3"
      ),
      each = 5
    ),
  
  stringsAsFactors = FALSE
)

rownames(meta_dpc) <-
  meta_dpc$GSM

meta_dpc$condition <- factor(
  meta_dpc$condition,
  levels = condition_order
)

meta_dpc$replicate <- factor(
  meta_dpc$replicate,
  levels = c(
    "R1",
    "R2",
    "R3"
  )
)

write.csv(
  meta_dpc,
  file.path(
    results_dir,
    "GSE293743_DPC_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat("Primary DPC metadata:\n")
print(meta_dpc)
cat("\n")

cat(
  "NOTE: GSM8889398 Fibro is excluded from the PRIMARY validation.\n\n"
)


############################################################
# 12. Select the 15 DPC CEL files
############################################################

dpc_files <- ordered_files[
  meta_dpc$GSM
]


############################################################
# 13. Read Clariom S CEL files
#
# First try compressed CEL.gz directly.
# If that fails, automatically gunzip to a temp folder.
############################################################

cat("Reading 15 DPC Clariom S CEL files...\n")
cat("This can take several minutes.\n\n")

raw_data <- tryCatch(
  
  {
    
    oligo::read.celfiles(
      filenames =
        dpc_files
    )
    
  },
  
  error = function(e) {
    
    cat(
      "Direct reading of CEL.gz failed.\n"
    )
    
    cat(
      "Error message:\n",
      conditionMessage(e),
      "\n\n"
    )
    
    NULL
  }
)


############################################################
# 14. Automatic decompression fallback
############################################################

if (is.null(raw_data)) {
  
  cat(
    "Trying automatic decompression...\n\n"
  )
  
  temp_cel_dir <- file.path(
    tempdir(),
    "GSE293743_CEL_uncompressed"
  )
  
  dir.create(
    temp_cel_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  uncompressed_files <- character(
    length(dpc_files)
  )
  
  for (i in seq_along(dpc_files)) {
    
    input_file <-
      dpc_files[i]
    
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
    
    if (grepl(
      "\\.gz$",
      input_file,
      ignore.case = TRUE
    )) {
      
      R.utils::gunzip(
        filename =
          input_file,
        destname =
          output_file,
        overwrite =
          TRUE,
        remove =
          FALSE
      )
      
    } else {
      
      file.copy(
        from =
          input_file,
        to =
          output_file,
        overwrite =
          TRUE
      )
    }
    
    uncompressed_files[i] <-
      output_file
  }
  
  raw_data <- oligo::read.celfiles(
    filenames =
      uncompressed_files
  )
}


############################################################
# 15. Rename raw samples using GSM accessions
############################################################

Biobase::sampleNames(
  raw_data
) <- meta_dpc$GSM

cat(
  "Raw Clariom S object class:\n"
)

print(
  class(raw_data)
)

cat("\n")


############################################################
# 16. RMA normalization
#
# Clariom S is handled by oligo.
# Try target='core'; fall back to standard RMA if the
# installed object class does not expose a target argument.
############################################################

cat(
  "Running RMA normalization...\n"
)

eset <- tryCatch(
  
  {
    
    oligo::rma(
      raw_data,
      target = "core"
    )
    
  },
  
  error = function(e) {
    
    cat(
      "RMA with target='core' was not accepted by this object class.\n"
    )
    
    cat(
      "Retrying standard oligo::rma(raw_data)...\n\n"
    )
    
    oligo::rma(
      raw_data
    )
  }
)


############################################################
# 17. Extract transcript-cluster expression
############################################################

expr_tc <- Biobase::exprs(
  eset
)

colnames(expr_tc) <-
  meta_dpc$GSM

cat(
  "RMA transcript-cluster matrix dimensions:\n"
)

print(
  dim(expr_tc)
)

cat(
  "\nRMA value summary:\n"
)

print(
  summary(
    as.vector(expr_tc)
  )
)

cat("\n")

if (
  any(is.na(expr_tc)) ||
  any(!is.finite(expr_tc))
) {
  
  stop(
    "ERROR: RMA expression matrix contains NA/Inf."
  )
}


############################################################
# 18. Save transcript-cluster matrix
############################################################

saveRDS(
  expr_tc,
  file.path(
    results_dir,
    "GSE293743_RMA_transcriptcluster_expression.rds"
  )
)


############################################################
# 19. Annotate Clariom S transcript clusters
############################################################

cat(
  "Annotating Clariom S transcript clusters to gene symbols...\n"
)

available_keytypes <- AnnotationDbi::keytypes(
  clariomshumantranscriptcluster.db
)

if (!"PROBEID" %in% available_keytypes) {
  
  cat(
    "Available key types:\n"
  )
  
  print(
    available_keytypes
  )
  
  stop(
    "ERROR: PROBEID key type not available in annotation package."
  )
}


gene_symbol <- AnnotationDbi::mapIds(
  
  clariomshumantranscriptcluster.db,
  
  keys =
    rownames(expr_tc),
  
  column =
    "SYMBOL",
  
  keytype =
    "PROBEID",
  
  multiVals =
    "first"
)


############################################################
# 20. Annotation quality check
############################################################

n_mapped <- sum(
  !is.na(gene_symbol) &
    gene_symbol != ""
)

mapping_fraction <-
  n_mapped /
  length(gene_symbol)

cat(
  "Mapped transcript clusters:",
  n_mapped,
  "/",
  length(gene_symbol),
  " (",
  round(
    100 * mapping_fraction,
    1
  ),
  "%)\n\n",
  sep = ""
)

if (n_mapped < 5000) {
  
  stop(
    paste0(
      "ERROR: only ",
      n_mapped,
      " transcript clusters mapped to gene symbols. ",
      "Please stop and send the console output before continuing."
    )
  )
}


############################################################
# 21. Save annotation table
############################################################

annotation_table <- data.frame(
  
  TranscriptCluster =
    rownames(expr_tc),
  
  GeneSymbol =
    as.character(
      gene_symbol
    ),
  
  stringsAsFactors = FALSE
)

write.csv(
  annotation_table,
  file.path(
    results_dir,
    "GSE293743_transcriptcluster_annotation.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 22. Collapse transcript clusters to one gene
############################################################

valid_annotation <-
  !is.na(gene_symbol) &
  gene_symbol != ""

expr_gene <- limma::avereps(
  
  expr_tc[
    valid_annotation,
    ,
    drop = FALSE
  ],
  
  ID =
    as.character(
      gene_symbol[
        valid_annotation
      ]
    )
)

cat(
  "Gene-level expression matrix dimensions:\n"
)

print(
  dim(expr_gene)
)

cat("\n")


############################################################
# 23. Save gene-level matrix
############################################################

saveRDS(
  expr_gene,
  file.path(
    results_dir,
    "GSE293743_RMA_gene_expression.rds"
  )
)

gene_connection <- gzfile(
  file.path(
    results_dir,
    "GSE293743_RMA_gene_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  expr_gene,
  gene_connection,
  quote = FALSE
)

close(
  gene_connection
)


############################################################
# 24. QC: RMA boxplot
############################################################

pdf(
  file.path(
    figures_dir,
    "GSE293743_RMA_boxplot.pdf"
  ),
  width = 10,
  height = 6
)

boxplot(
  expr_gene,
  las = 2,
  outline = FALSE,
  main =
    "GSE293743 after RMA normalization",
  ylab =
    "log2 expression"
)

dev.off()


############################################################
# 25. QC: density plot
############################################################

pdf(
  file.path(
    figures_dir,
    "GSE293743_RMA_density.pdf"
  ),
  width = 9,
  height = 6
)

limma::plotDensities(
  expr_gene,
  legend = FALSE,
  main =
    "GSE293743 expression density after RMA"
)

dev.off()


############################################################
# 26. PCA on top 5000 variable genes
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
    seq_len(
      n_pca_genes
    )
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

variance_explained <-
  100 *
  pca$sdev^2 /
  sum(
    pca$sdev^2
  )

pca_df <- data.frame(
  
  GSM =
    rownames(
      pca$x
    ),
  
  PC1 =
    pca$x[, 1],
  
  PC2 =
    pca$x[, 2],
  
  stringsAsFactors = FALSE
)

pca_df$condition <-
  meta_dpc[
    pca_df$GSM,
    "condition"
  ]

pca_df$replicate <-
  meta_dpc[
    pca_df$GSM,
    "replicate"
  ]


############################################################
# 27. PCA figure
############################################################

p_pca <- ggplot(
  
  pca_df,
  
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    shape = replicate,
    label = GSM
  )
  
) +
  
  geom_point(
    size = 4
  ) +
  
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  labs(
    
    title =
      "GSE293743 independent validation: PCA",
    
    subtitle =
      "Top 5,000 most variable genes; 15 DPC samples",
    
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
      "Replicate"
  ) +
  
  theme_classic(
    base_size = 13
  )


ggsave(
  file.path(
    figures_dir,
    "GSE293743_PCA.png"
  ),
  p_pca,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE293743_PCA.pdf"
  ),
  p_pca,
  width = 9,
  height = 7
)


############################################################
# 28. Sample correlation heatmap
############################################################

sample_cor <- cor(
  expr_gene[
    pca_genes,
    ,
    drop = FALSE
  ],
  method = "pearson"
)

annotation_col <- data.frame(
  
  Condition =
    as.character(
      meta_dpc$condition
    ),
  
  Replicate =
    as.character(
      meta_dpc$replicate
    ),
  
  stringsAsFactors = FALSE
)

rownames(annotation_col) <-
  meta_dpc$GSM

pheatmap(
  sample_cor,
  annotation_col =
    annotation_col,
  annotation_row =
    annotation_col,
  border_color = NA,
  main =
    "GSE293743 sample correlation",
  filename =
    file.path(
      figures_dir,
      "GSE293743_sample_correlation_heatmap.png"
    ),
  width = 10,
  height = 9
)

pheatmap(
  sample_cor,
  annotation_col =
    annotation_col,
  annotation_row =
    annotation_col,
  border_color = NA,
  main =
    "GSE293743 sample correlation",
  filename =
    file.path(
      figures_dir,
      "GSE293743_sample_correlation_heatmap.pdf"
    ),
  width = 10,
  height = 9
)


############################################################
# 29. Load PREDEFINED GSE44765 strict rescue signature
############################################################

signature <- read.csv(
  discovery_signature_file,
  stringsAsFactors = FALSE
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
      "ERROR: discovery signature lacks required columns: ",
      paste(
        required_signature_columns,
        collapse = ", "
      )
    )
  )
}


############################################################
# 30. Define upward/downward signature genes
############################################################

signature$GeneSymbol <-
  as.character(
    signature$GeneSymbol
  )

signature <- signature[
  !is.na(signature$GeneSymbol) &
    signature$GeneSymbol != "",
  ,
  drop = FALSE
]

signature <- signature[
  !duplicated(
    signature$GeneSymbol
  ),
  ,
  drop = FALSE
]

up_genes_all <- signature$GeneSymbol[
  signature$Rescue_type ==
    "Restored_upward"
]

down_genes_all <- signature$GeneSymbol[
  signature$Rescue_type ==
    "Restored_downward"
]

up_genes_overlap <- intersect(
  up_genes_all,
  rownames(expr_gene)
)

down_genes_overlap <- intersect(
  down_genes_all,
  rownames(expr_gene)
)

signature_overlap <- signature[
  signature$GeneSymbol %in%
    rownames(expr_gene),
  ,
  drop = FALSE
]

write.csv(
  signature_overlap,
  file.path(
    results_dir,
    "GSE44765_signature_overlap_in_GSE293743.csv"
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
  "Overlapping strict signature genes:",
  nrow(signature_overlap),
  "\n"
)

cat(
  "Upward genes available:",
  length(up_genes_overlap),
  "/",
  length(up_genes_all),
  "\n"
)

cat(
  "Downward genes available:",
  length(down_genes_overlap),
  "/",
  length(down_genes_all),
  "\n\n"
)


############################################################
# 31. Require adequate cross-platform signature coverage
############################################################

signature_coverage <-
  nrow(signature_overlap) /
  nrow(signature)

if (signature_coverage < 0.50) {
  
  stop(
    paste0(
      "ERROR: only ",
      round(
        100 * signature_coverage,
        1
      ),
      "% of the discovery signature is represented. ",
      "Please stop and inspect annotation/platform mapping."
    )
  )
}


############################################################
# 32. Calculate gene-wise Z scores within external cohort
#
# IMPORTANT:
# The signature was defined in GSE44765.
# No GSE293743 labels are used for feature selection here.
############################################################

signature_genes_overlap <- unique(
  c(
    up_genes_overlap,
    down_genes_overlap
  )
)

signature_expr <- expr_gene[
  signature_genes_overlap,
  ,
  drop = FALSE
]

signature_z <- t(
  scale(
    t(
      signature_expr
    )
  )
)

signature_z[
  !is.finite(
    signature_z
  )
] <- 0


############################################################
# 33. Compute predefined signature scores
#
# High total score =
# more "GSE44765 3D-rescue-like / Intact-like" direction.
############################################################

up_score <- colMeans(
  signature_z[
    up_genes_overlap,
    ,
    drop = FALSE
  ]
)

down_score <- colMeans(
  signature_z[
    down_genes_overlap,
    ,
    drop = FALSE
  ]
)

total_score <-
  up_score -
  down_score


############################################################
# 34. Build score table
############################################################

score_df <- data.frame(
  
  GSM =
    colnames(expr_gene),
  
  condition =
    meta_dpc[
      colnames(expr_gene),
      "condition"
    ],
  
  replicate =
    meta_dpc[
      colnames(expr_gene),
      "replicate"
    ],
  
  Up_score =
    as.numeric(
      up_score[
        colnames(expr_gene)
      ]
    ),
  
  Down_score =
    as.numeric(
      down_score[
        colnames(expr_gene)
      ]
    ),
  
  Rescue_signature_score =
    as.numeric(
      total_score[
        colnames(expr_gene)
      ]
    ),
  
  stringsAsFactors = FALSE
)

score_df$condition <- factor(
  score_df$condition,
  levels = condition_order
)

score_df$replicate <- factor(
  score_df$replicate,
  levels = c(
    "R1",
    "R2",
    "R3"
  )
)

write.csv(
  score_df,
  file.path(
    results_dir,
    "GSE293743_predefined_signature_scores.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 35. Signature score trajectory
############################################################

p_score <- ggplot(
  
  score_df,
  
  aes(
    x = condition,
    y = Rescue_signature_score,
    group = replicate,
    shape = replicate
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_line(
    alpha = 0.65
  ) +
  
  geom_point(
    size = 3.3
  ) +
  
  labs(
    
    title =
      "Independent validation of the GSE44765 rescue signature",
    
    subtitle =
      "Predefined signature; each line represents one replicate block",
    
    x =
      NULL,
    
    y =
      "Predefined rescue signature score",
    
    shape =
      "Replicate"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      )
  )


ggsave(
  file.path(
    figures_dir,
    "GSE293743_rescue_signature_trajectory.png"
  ),
  p_score,
  width = 9,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE293743_rescue_signature_trajectory.pdf"
  ),
  p_score,
  width = 9,
  height = 6.5
)


############################################################
# 36. Mean score summary by condition
############################################################

score_summary_list <- lapply(
  
  condition_order,
  
  function(cond) {
    
    x <- score_df$Rescue_signature_score[
      score_df$condition ==
        cond
    ]
    
    data.frame(
      Condition =
        cond,
      N =
        length(x),
      Mean =
        mean(x),
      SD =
        sd(x),
      Median =
        median(x),
      stringsAsFactors =
        FALSE
    )
  }
)

score_summary <- do.call(
  rbind,
  score_summary_list
)

write.csv(
  score_summary,
  file.path(
    results_dir,
    "GSE293743_signature_score_summary_by_condition.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 37. Paired/block-aware differential expression model
#
# replicate is included to remove block-specific variation.
############################################################

design <- model.matrix(
  ~ 0 + condition + replicate,
  data = meta_dpc
)

cat(
  "External-validation design matrix rank:",
  qr(design)$rank,
  "/",
  ncol(design),
  "\n"
)

if (qr(design)$rank != ncol(design)) {
  
  stop(
    "ERROR: external-validation design matrix is not full rank."
  )
}

fit <- limma::lmFit(
  expr_gene,
  design
)


############################################################
# 38. Treatment vs Control contrasts
############################################################

contrast_matrix <- limma::makeContrasts(
  
  DPAC_vs_Control =
    conditionDPAC -
    conditionControl,
  
  CHIR_vs_Control =
    conditionCHIR -
    conditionControl,
  
  FGF9_vs_Control =
    conditionFGF9 -
    conditionControl,
  
  CHIR_FGF9_vs_Control =
    conditionCHIR_FGF9 -
    conditionControl,
  
  levels =
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


############################################################
# 39. Save external DE tables
############################################################

contrast_names <- colnames(
  contrast_matrix
)

external_de_list <- list()

for (contrast_name in contrast_names) {
  
  tt <- limma::topTable(
    fit2,
    coef =
      contrast_name,
    number =
      Inf,
    adjust.method =
      "BH",
    sort.by =
      "none"
  )
  
  tt$GeneSymbol <-
    rownames(tt)
  
  tt <- tt[
    ,
    c(
      "GeneSymbol",
      setdiff(
        colnames(tt),
        "GeneSymbol"
      )
    ),
    drop = FALSE
  ]
  
  external_de_list[contrast_name] <- list(tt)
  
  write.csv(
    tt,
    file.path(
      results_dir,
      paste0(
        "DE_",
        contrast_name,
        "_all_genes.csv"
      )
    ),
    row.names = FALSE,
    quote = FALSE
  )
}


############################################################
# 40. CAMERA validation of predefined gene sets
#
# Expected:
#   Restored_upward -> Up in a rescue-like intervention
#   Restored_downward -> Down in a rescue-like intervention
############################################################

signature_index <- list(
  
  Restored_upward =
    which(
      rownames(expr_gene) %in%
        up_genes_overlap
    ),
  
  Restored_downward =
    which(
      rownames(expr_gene) %in%
        down_genes_overlap
    )
)

camera_results_list <- list()

for (j in seq_len(
  ncol(
    contrast_matrix
  )
)) {
  
  contrast_name <-
    colnames(
      contrast_matrix
    )[j]
  
  camera_result <- limma::camera(
    
    y =
      expr_gene,
    
    index =
      signature_index,
    
    design =
      design,
    
    contrast =
      contrast_matrix[, j],
    
    inter.gene.cor =
      0.01
  )
  
  camera_result$GeneSet <-
    rownames(
      camera_result
    )
  
  camera_result$Contrast <-
    contrast_name
  
  camera_result <- camera_result[
    ,
    c(
      "Contrast",
      "GeneSet",
      setdiff(
        colnames(camera_result),
        c(
          "Contrast",
          "GeneSet"
        )
      )
    ),
    drop = FALSE
  ]
  
  camera_results_list[contrast_name] <- list(camera_result)
}

camera_results <- do.call(
  rbind,
  camera_results_list
)

rownames(camera_results) <-
  NULL

write.csv(
  camera_results,
  file.path(
    results_dir,
    "GSE293743_CAMERA_predefined_signature_validation.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 41. Cross-study effect concordance
#
# Compare the discovery rescue effect from GSE44765:
#   P3_3D - P3_2D
#
# against each independent intervention effect:
#   Treatment - Control
############################################################

discovery_effect <- signature[
  ,
  c(
    "GeneSymbol",
    "Rescue_logFC",
    "Rescue_type"
  ),
  drop = FALSE
]

concordance_rows <- list()
merged_effect_tables <- list()

counter <- 1

for (contrast_name in contrast_names) {
  
  external_tt <- getElement(external_de_list, contrast_name)
  
  external_small <- external_tt[
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
  
  colnames(
    external_small
  ) <- c(
    "GeneSymbol",
    "External_logFC",
    "External_t",
    "External_P",
    "External_FDR"
  )
  
  merged <- merge(
    discovery_effect,
    external_small,
    by = "GeneSymbol"
  )
  
  merged$Expected_signed_external_logFC <-
    sign(
      merged$Rescue_logFC
    ) *
    merged$External_logFC
  
  merged$Concordant_direction <-
    sign(
      merged$Rescue_logFC
    ) ==
    sign(
      merged$External_logFC
    )
  
  valid_for_cor <-
    is.finite(
      merged$Rescue_logFC
    ) &
    is.finite(
      merged$External_logFC
    )
  
  rho <- suppressWarnings(
    cor(
      merged$Rescue_logFC[
        valid_for_cor
      ],
      merged$External_logFC[
        valid_for_cor
      ],
      method = "spearman"
    )
  )
  
  concordance_fraction <- mean(
    merged$Concordant_direction[
      valid_for_cor
    ],
    na.rm = TRUE
  )
  
  median_signed_effect <- median(
    merged$Expected_signed_external_logFC[
      valid_for_cor
    ],
    na.rm = TRUE
  )
  
  mean_signed_effect <- mean(
    merged$Expected_signed_external_logFC[
      valid_for_cor
    ],
    na.rm = TRUE
  )
  
  concordance_rows[counter] <- list(data.frame(
    
    Contrast =
      contrast_name,
    
    N_signature_genes =
      sum(
        valid_for_cor
      ),
    
    Spearman_rho =
      rho,
    
    Directional_concordance =
      concordance_fraction,
    
    Median_expected_signed_external_logFC =
      median_signed_effect,
    
    Mean_expected_signed_external_logFC =
      mean_signed_effect,
    
    stringsAsFactors =
      FALSE
  ))
  
  merged$Contrast <-
    contrast_name
  
  merged_effect_tables[contrast_name] <- list(merged)
  
  counter <-
    counter + 1
}

concordance_summary <- do.call(
  rbind,
  concordance_rows
)

write.csv(
  concordance_summary,
  file.path(
    results_dir,
    "GSE293743_cross_study_concordance_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

all_merged_effects <- do.call(
  rbind,
  merged_effect_tables
)

rownames(all_merged_effects) <-
  NULL

write.csv(
  all_merged_effects,
  file.path(
    results_dir,
    "GSE293743_cross_study_gene_effects.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 42. Cross-study effect concordance plot
############################################################

plot_effects <- all_merged_effects[
  is.finite(
    all_merged_effects$Rescue_logFC
  ) &
    is.finite(
      all_merged_effects$External_logFC
    ),
  ,
  drop = FALSE
]

plot_effects$Signature_direction <-
  ifelse(
    plot_effects$Rescue_type ==
      "Restored_upward",
    "Restored upward",
    "Restored downward"
  )

p_concordance <- ggplot(
  
  plot_effects,
  
  aes(
    x = Rescue_logFC,
    y = External_logFC
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
    size = 1.0
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  
  facet_wrap(
    ~ Contrast,
    ncol = 2
  ) +
  
  labs(
    
    title =
      "Cross-study effect concordance",
    
    subtitle =
      "Discovery effect: GSE44765 P3 3D - P3 2D; external effect: treatment - control",
    
    x =
      "GSE44765 discovery rescue log2FC",
    
    y =
      "GSE293743 external log2FC"
  ) +
  
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    figures_dir,
    "GSE293743_cross_study_effect_concordance.png"
  ),
  p_concordance,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE293743_cross_study_effect_concordance.pdf"
  ),
  p_concordance,
  width = 10,
  height = 8
)


############################################################
# 43. Discovery-ranked signature heatmap
#
# Genes are ranked ONLY by GSE44765 discovery statistics,
# not by GSE293743 performance.
############################################################

if ("Combined_FDR" %in%
    colnames(signature)) {
  
  signature_for_heatmap <- signature[
    order(
      signature$Combined_FDR,
      -abs(
        signature$Rescue_logFC
      ),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
} else {
  
  signature_for_heatmap <- signature[
    order(
      -abs(
        signature$Rescue_logFC
      )
    ),
    ,
    drop = FALSE
  ]
}

signature_for_heatmap <- signature_for_heatmap[
  signature_for_heatmap$GeneSymbol %in%
    rownames(expr_gene),
  ,
  drop = FALSE
]

n_heatmap_genes <- min(
  50,
  nrow(
    signature_for_heatmap
  )
)

heatmap_genes <- unique(
  signature_for_heatmap$GeneSymbol[
    seq_len(
      n_heatmap_genes
    )
  ]
)

heatmap_expr <- expr_gene[
  heatmap_genes,
  ,
  drop = FALSE
]

heatmap_z <- t(
  scale(
    t(
      heatmap_expr
    )
  )
)

heatmap_z[
  !is.finite(
    heatmap_z
  )
] <- 0

pheatmap(
  heatmap_z,
  annotation_col =
    annotation_col,
  cluster_rows =
    TRUE,
  cluster_cols =
    TRUE,
  show_rownames =
    TRUE,
  fontsize_row =
    7,
  border_color =
    NA,
  main =
    "Predefined GSE44765 strict rescue signature in GSE293743",
  filename =
    file.path(
      figures_dir,
      "GSE293743_predefined_signature_heatmap.png"
    ),
  width = 11,
  height = 11
)

pheatmap(
  heatmap_z,
  annotation_col =
    annotation_col,
  cluster_rows =
    TRUE,
  cluster_cols =
    TRUE,
  show_rownames =
    TRUE,
  fontsize_row =
    7,
  border_color =
    NA,
  main =
    "Predefined GSE44765 strict rescue signature in GSE293743",
  filename =
    file.path(
      figures_dir,
      "GSE293743_predefined_signature_heatmap.pdf"
    ),
  width = 11,
  height = 11
)


############################################################
# 44. Simple block-matched score differences
############################################################

control_scores <- score_df[
  score_df$condition ==
    "Control",
  c(
    "replicate",
    "Rescue_signature_score"
  )
]

colnames(control_scores)[2] <-
  "Control_score"

score_difference_rows <- list()

counter <- 1

for (cond in condition_order[
  condition_order !=
  "Control"
]) {
  
  treatment_scores <- score_df[
    score_df$condition ==
      cond,
    c(
      "replicate",
      "Rescue_signature_score"
    )
  ]
  
  colnames(
    treatment_scores
  )[2] <-
    "Treatment_score"
  
  paired <- merge(
    control_scores,
    treatment_scores,
    by = "replicate"
  )
  
  paired$Difference <-
    paired$Treatment_score -
    paired$Control_score
  
  score_difference_rows[counter] <- list(data.frame(
    
    Condition =
      cond,
    
    Mean_difference =
      mean(
        paired$Difference
      ),
    
    Median_difference =
      median(
        paired$Difference
      ),
    
    All_3_blocks_positive =
      all(
        paired$Difference > 0
      ),
    
    R1_difference =
      paired$Difference[
        paired$replicate ==
          "R1"
      ],
    
    R2_difference =
      paired$Difference[
        paired$replicate ==
          "R2"
      ],
    
    R3_difference =
      paired$Difference[
        paired$replicate ==
          "R3"
      ],
    
    stringsAsFactors =
      FALSE
  ))
  
  counter <-
    counter + 1
}

score_differences <- do.call(
  rbind,
  score_difference_rows
)

write.csv(
  score_differences,
  file.path(
    results_dir,
    "GSE293743_block_matched_score_differences.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 45. Save design and contrasts
############################################################

write.csv(
  design,
  file.path(
    results_dir,
    "GSE293743_design_matrix.csv"
  ),
  quote = FALSE
)

write.csv(
  contrast_matrix,
  file.path(
    results_dir,
    "GSE293743_contrast_matrix.csv"
  ),
  quote = FALSE
)


############################################################
# 46. Human-readable analysis summary
############################################################

summary_file <- file.path(
  results_dir,
  "ANALYSIS_SUMMARY_GSE293743.txt"
)

sink(
  summary_file
)

cat(
  "GSE293743 independent external validation\n"
)

cat(
  "==========================================\n\n"
)

cat(
  "Platform: Affymetrix Clariom S Human\n"
)

cat(
  "Primary analysis samples: 15 DPC arrays\n"
)

cat(
  "Conditions: Control, DPAC, CHIR, FGF9, CHIR_FGF9\n"
)

cat(
  "Replicate blocks: R1, R2, R3\n"
)

cat(
  "Fibroblast GSM8889398 excluded from primary validation.\n\n"
)

cat(
  "RMA transcript clusters:",
  nrow(expr_tc),
  "\n"
)

cat(
  "Mapped transcript clusters:",
  n_mapped,
  "\n"
)

cat(
  "Gene-level features:",
  nrow(expr_gene),
  "\n\n"
)

cat(
  "Discovery strict signature genes:",
  nrow(signature),
  "\n"
)

cat(
  "Signature genes represented in external platform:",
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
  "Upward signature overlap:",
  length(up_genes_overlap),
  "\n"
)

cat(
  "Downward signature overlap:",
  length(down_genes_overlap),
  "\n\n"
)

cat(
  "Predefined signature score by condition:\n"
)

for (i in seq_len(
  nrow(score_summary)
)) {
  
  cat(
    score_summary$Condition[i],
    ": mean=",
    round(
      score_summary$Mean[i],
      4
    ),
    ", SD=",
    round(
      score_summary$SD[i],
      4
    ),
    "\n",
    sep = ""
  )
}

cat("\n")

cat(
  "Block-matched treatment - control score differences:\n"
)

for (i in seq_len(
  nrow(score_differences)
)) {
  
  cat(
    score_differences$Condition[i],
    ": mean difference=",
    round(
      score_differences$Mean_difference[i],
      4
    ),
    ", all 3 blocks positive=",
    score_differences$All_3_blocks_positive[i],
    "\n",
    sep = ""
  )
}

cat("\n")

cat(
  "CAMERA predefined gene-set validation:\n"
)

print(
  camera_results
)

cat("\n")

cat(
  "Cross-study effect concordance:\n"
)

print(
  concordance_summary
)

sink()


############################################################
# 47. Save workspace
############################################################

save(
  
  expr_tc,
  expr_gene,
  meta_dpc,
  signature,
  signature_overlap,
  score_df,
  score_summary,
  score_differences,
  design,
  contrast_matrix,
  fit2,
  camera_results,
  concordance_summary,
  all_merged_effects,
  
  file =
    file.path(
      results_dir,
      "GSE293743_external_validation_workspace.RData"
    )
)


############################################################
# 48. Save session information
############################################################

sink(
  file.path(
    results_dir,
    "sessionInfo_GSE293743.txt"
  )
)

sessionInfo()

sink()


############################################################
# 49. Final report
############################################################

cat("\n")
cat("============================================================\n")
cat(" GSE293743 EXTERNAL VALIDATION COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "Gene-level expression features:",
  nrow(expr_gene),
  "\n"
)

cat(
  "Strict discovery signature coverage:",
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

cat("\nPredefined signature score means:\n")
print(
  score_summary
)

cat("\nBlock-matched score differences:\n")
print(
  score_differences
)

cat("\nCAMERA validation:\n")
print(
  camera_results
)

cat("\nCross-study concordance:\n")
print(
  concordance_summary
)

cat("\nMain figures:\n")
cat(
  file.path(
    figures_dir,
    "GSE293743_PCA.png"
  ),
  "\n"
)
cat(
  file.path(
    figures_dir,
    "GSE293743_rescue_signature_trajectory.png"
  ),
  "\n"
)
cat(
  file.path(
    figures_dir,
    "GSE293743_cross_study_effect_concordance.png"
  ),
  "\n"
)
cat(
  file.path(
    figures_dir,
    "GSE293743_predefined_signature_heatmap.png"
  ),
  "\n"
)

cat("\nMain summary:\n")
cat(
  summary_file,
  "\n"
)

cat("\nIMPORTANT INTERPRETATION RULE:\n")
cat(
  "This is a true external validation: GSE293743 labels were not used to select discovery genes.\n"
)
cat(
  "Do not alter the GSE44765 signature based on GSE293743 performance.\n"
)

cat("============================================================\n")
