############################################################
# 09_GSE184501_GSE184866_2D3D_AGA_validation.R
#
# COMPLETE ONE-CLICK VERSION
#
# Scientific purpose
# ------------------
# Compare the BAB/BAN transcriptional contrast under:
#
#   GSE184501: 2D culture + 10 nM DHT, 24 h
#   GSE184866: 3D culture + 10 nM DHT, 24 h
#
# IMPORTANT:
#   The two studies use different RNA-seq experiments/platforms
#   and different processed-data formats.
#
# Therefore this script NEVER merges their expression matrices
# and NEVER treats "2D versus 3D" as a direct sample-level
# contrast.
#
# Instead, each dataset is analyzed independently, and only
# within-study BAB/BAN EFFECTS are compared across datasets.
#
# Primary frozen biological reference:
#   GSE44765 strict rescue signature
#
# Primary question:
#   Does the BAN-versus-BAB effect align with the GSE44765
#   rescue direction in 2D, and is that alignment weaker,
#   preserved, or reversed in the independent 3D dataset?
#
# Secondary mechanistic questions:
#   1. What happens to the predefined GSE44765 mitotic/
#      cell-cycle rescue module?
#   2. What happens to the exploratory DHT anti-rescue
#      metabolic sentinel genes identified in GSE66664?
#
# Statistical interpretation:
#   All BAB/BAN replicates here are experimental/technical
#   replicates of cell-line models, NOT independent patients.
#
# Required input files
# --------------------
# data_raw/GSE184501/
#   GSE184501_2D_cuffnorm_fpkm_quantilenormalised.txt.gz
#
# data_raw/GSE184866/
#   GSE184866_3D_star_count.txt.gz
#
# Required discovery signature
# ----------------------------
# results/**/GSE44765_STRICT_rescue_genes.csv
#
# Optional mechanistic module from script 08
# ------------------------------------------
# results/GSE66664_mechanism/
#   GSE44765_predefined_mitotic_cellcycle_module_genes.csv
#
# Outputs
# -------
# results/GSE184501_GSE184866/
# figures/GSE184501_GSE184866/
############################################################


############################################################
# 0. Basic settings
############################################################

graphics.off()

options(
  stringsAsFactors = FALSE
)

set.seed(
  184501866
)

cat("\n")
cat("============================================================\n")
cat(" GSE184501 + GSE184866 2D/3D AGA VALIDATION STARTED\n")
cat("============================================================\n\n")

cat(
  "Working directory:\n",
  getwd(),
  "\n\n"
)


############################################################
# 1. Paths
############################################################

raw_2d_dir <-
  "data_raw/GSE184501"

raw_3d_dir <-
  "data_raw/GSE184866"

results_dir <-
  "results/GSE184501_GSE184866"

figures_dir <-
  "figures/GSE184501_GSE184866"

dir.create(
  raw_2d_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  raw_3d_dir,
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
# 2. Install/load packages
############################################################

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap"
)

for (
  pkg in cran_packages
) {
  
  if (
    !requireNamespace(
      pkg,
      quietly = TRUE
    )
  ) {
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

if (
  !requireNamespace(
    "BiocManager",
    quietly = TRUE
  )
) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_packages <- c(
  "edgeR",
  "limma",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

for (
  pkg in bioc_packages
) {
  
  if (
    !requireNamespace(
      pkg,
      quietly = TRUE
    )
  ) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

suppressPackageStartupMessages({
  
  library(
    edgeR
  )
  
  library(
    limma
  )
  
  library(
    org.Hs.eg.db
  )
  
  library(
    AnnotationDbi
  )
  
  library(
    ggplot2
  )
  
  library(
    ggrepel
  )
  
  library(
    pheatmap
  )
  
})

cat(
  "Packages loaded successfully.\n\n"
)


############################################################
# 3. Locate input files robustly
############################################################

file_2d_candidates <- list.files(
  path = raw_2d_dir,
  pattern = "^GSE184501_2D_cuffnorm_fpkm_quantilenormalised.*\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (
  length(
    file_2d_candidates
  ) == 0
) {
  
  stop(
    paste0(
      "\nERROR:\nNo GSE184501 2D file found under:\n",
      raw_2d_dir,
      "\n"
    )
  )
}

file_2d <- file_2d_candidates[
  which.max(
    file.info(
      file_2d_candidates
    )$size
  )
]

file_3d_candidates <- list.files(
  path = raw_3d_dir,
  pattern = "^GSE184866_3D_star_count.*\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (
  length(
    file_3d_candidates
  ) == 0
) {
  
  stop(
    paste0(
      "\nERROR:\nNo GSE184866 3D count file found under:\n",
      raw_3d_dir,
      "\n"
    )
  )
}

file_3d <- file_3d_candidates[
  which.max(
    file.info(
      file_3d_candidates
    )$size
  )
]

cat(
  "2D input:\n",
  file_2d,
  "\n\n"
)

cat(
  "3D input:\n",
  file_3d,
  "\n\n"
)


############################################################
# 4. Locate frozen GSE44765 strict rescue signature
############################################################

preferred_signature <- file.path(
  "results",
  "GSE44765_rescue",
  "GSE44765_STRICT_rescue_genes.csv"
)

if (
  file.exists(
    preferred_signature
  )
) {
  
  signature_file <-
    preferred_signature
  
} else {
  
  signature_candidates <- list.files(
    path = "results",
    pattern = "^GSE44765_STRICT_rescue_genes\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (
    length(
      signature_candidates
    ) == 0
  ) {
    
    stop(
      paste0(
        "\nERROR:\n",
        "Cannot find GSE44765_STRICT_rescue_genes.csv under results/.\n"
      )
    )
  }
  
  signature_file <-
    signature_candidates[1]
}

cat(
  "Frozen discovery signature:\n",
  signature_file,
  "\n\n"
)


############################################################
# 5. Load frozen signature
############################################################

signature <- read.csv(
  signature_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_signature_columns <- c(
  "GeneSymbol",
  "Rescue_type",
  "Rescue_logFC"
)

if (
  !all(
    required_signature_columns %in%
    colnames(
      signature
    )
  )
) {
  
  stop(
    paste0(
      "ERROR: signature lacks required columns:\n",
      paste(
        setdiff(
          required_signature_columns,
          colnames(
            signature
          )
        ),
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
  !is.na(
    signature$GeneSymbol
  ) &
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

cat(
  "Frozen GSE44765 strict rescue genes:",
  nrow(
    signature
  ),
  "\n\n"
)


############################################################
# 6. Read GSE184501 2D quantile-normalized FPKM file
############################################################

raw_2d <- read.delim(
  file_2d,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

required_2d_columns <- c(
  "gene_short_name",
  "nearest_ref_id",
  "locus",
  "BAN_repA",
  "BAN_repB",
  "BAN_repC",
  "BAB_repA",
  "BAB_repB",
  "BAB_repC"
)

if (
  !all(
    required_2d_columns %in%
    colnames(
      raw_2d
    )
  )
) {
  
  stop(
    paste0(
      "ERROR: GSE184501 file lacks expected columns:\n",
      paste(
        setdiff(
          required_2d_columns,
          colnames(
            raw_2d
          )
        ),
        collapse = ", "
      )
    )
  )
}

cat(
  "GSE184501 transcript-level rows:",
  nrow(
    raw_2d
  ),
  "\n"
)

cat(
  "GSE184501 unique gene_short_name values:",
  length(
    unique(
      raw_2d$gene_short_name
    )
  ),
  "\n\n"
)


############################################################
# 7. Build GSE184501 2D gene-level expression
#
# The supplied file is transcript-level quantile-normalized
# FPKM on a linear scale.
#
# Primary aggregation:
#   sum transcript-level FPKM within gene_short_name,
#   then log2(gene FPKM + 1).
#
# No second quantile normalization is applied.
############################################################

sample_cols_2d <- c(
  "BAN_repA",
  "BAN_repB",
  "BAN_repC",
  "BAB_repA",
  "BAB_repB",
  "BAB_repC"
)

gene_name_2d <- trimws(
  as.character(
    raw_2d$gene_short_name
  )
)

valid_gene_2d <- !is.na(
  gene_name_2d
) &
  gene_name_2d != "" &
  gene_name_2d != "-"

fpkm_2d_transcript <- as.matrix(
  raw_2d[
    valid_gene_2d,
    sample_cols_2d,
    drop = FALSE
  ]
)

storage.mode(
  fpkm_2d_transcript
) <- "numeric"

if (
  any(
    !is.finite(
      fpkm_2d_transcript
    )
  ) ||
  any(
    fpkm_2d_transcript < 0
  )
) {
  
  stop(
    "ERROR: invalid GSE184501 FPKM values detected."
  )
}

gene_name_2d_valid <- gene_name_2d[
  valid_gene_2d
]

fpkm_2d_gene <- rowsum(
  fpkm_2d_transcript,
  group =
    gene_name_2d_valid,
  reorder = FALSE
)

expr_2d <- log2(
  fpkm_2d_gene + 1
)

expr_2d <- expr_2d[
  apply(
    expr_2d,
    1,
    function(x) {
      any(
        is.finite(
          x
        )
      ) &&
        max(
          x,
          na.rm = TRUE
        ) > 0
    }
  ),
  ,
  drop = FALSE
]

cat(
  "GSE184501 gene-level features after aggregation:",
  nrow(
    expr_2d
  ),
  "\n\n"
)

saveRDS(
  expr_2d,
  file.path(
    results_dir,
    "GSE184501_2D_gene_level_log2_FPKM.rds"
  )
)


############################################################
# 8. GSE184501 metadata
############################################################

meta_2d <- data.frame(
  
  Sample =
    colnames(
      expr_2d
    ),
  
  Group =
    c(
      "BAN",
      "BAN",
      "BAN",
      "BAB",
      "BAB",
      "BAB"
    ),
  
  Culture =
    "2D",
  
  DHT =
    "10nM_24h",
  
  stringsAsFactors =
    FALSE
)

meta_2d$Group <- factor(
  meta_2d$Group,
  levels = c(
    "BAB",
    "BAN"
  )
)

rownames(
  meta_2d
) <- meta_2d$Sample

write.csv(
  meta_2d,
  file.path(
    results_dir,
    "GSE184501_2D_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 9. GSE184501 PCA
############################################################

variance_2d <- apply(
  expr_2d,
  1,
  var,
  na.rm = TRUE
)

variance_2d <- sort(
  variance_2d,
  decreasing = TRUE
)

n_pca_2d <- min(
  5000,
  length(
    variance_2d
  )
)

pca_genes_2d <- names(
  variance_2d[
    seq_len(
      n_pca_2d
    )
  ]
)

pca_2d <- prcomp(
  t(
    expr_2d[
      pca_genes_2d,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

pca_2d_var <- 100 *
  pca_2d$sdev^2 /
  sum(
    pca_2d$sdev^2
  )

pca_2d_df <- data.frame(
  
  Sample =
    rownames(
      pca_2d$x
    ),
  
  PC1 =
    pca_2d$x[, 1],
  
  PC2 =
    pca_2d$x[, 2],
  
  stringsAsFactors =
    FALSE
)

pca_2d_df$Group <- meta_2d[
  pca_2d_df$Sample,
  "Group"
]

p_2d_pca <- ggplot(
  
  pca_2d_df,
  
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
      "GSE184501: 2D BAB/BAN after 10 nM DHT for 24 h",
    
    subtitle =
      "Author-supplied quantile-normalized FPKM; gene-level log2(FPKM+1)",
    
    x = paste0(
      "PC1 (",
      round(
        pca_2d_var[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        pca_2d_var[2],
        1
      ),
      "%)"
    ),
    
    color =
      "Cell line"
  ) +
  
  theme_classic(
    base_size = 13
  )

ggsave(
  file.path(
    figures_dir,
    "GSE184501_2D_PCA.png"
  ),
  p_2d_pca,
  width = 8,
  height = 6.5,
  dpi = 300
)


############################################################
# 10. GSE184501 differential expression: BAN - BAB
############################################################

design_2d <- model.matrix(
  ~ 0 + Group,
  data = meta_2d
)

colnames(
  design_2d
) <- sub(
  "^Group",
  "",
  colnames(
    design_2d
  )
)

contrast_2d <- limma::makeContrasts(
  BAN_vs_BAB = BAN - BAB,
  levels = design_2d
)

fit_2d <- limma::lmFit(
  expr_2d,
  design_2d
)

fit_2d <- limma::contrasts.fit(
  fit_2d,
  contrast_2d
)

fit_2d <- limma::eBayes(
  fit_2d,
  trend = TRUE,
  robust = TRUE
)

de_2d <- limma::topTable(
  fit_2d,
  coef = "BAN_vs_BAB",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

de_2d$GeneSymbol <- rownames(
  de_2d
)

de_2d <- de_2d[
  ,
  c(
    "GeneSymbol",
    setdiff(
      colnames(
        de_2d
      ),
      "GeneSymbol"
    )
  ),
  drop = FALSE
]

write.csv(
  de_2d,
  file.path(
    results_dir,
    "GSE184501_2D_BAN_vs_BAB_all_genes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 11. Read GSE184866 3D raw count matrix
############################################################

raw_3d <- read.delim(
  file_3d,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

expected_3d_columns <- c(
  "gene_name",
  "3D-Bald_repA",
  "3D-Bald_repB",
  "3D-Bald_repC",
  "3D-NonBald_repA",
  "3D-NonBald_repB",
  "3D-NonBald_repC"
)

if (
  !all(
    expected_3d_columns %in%
    colnames(
      raw_3d
    )
  )
) {
  
  stop(
    paste0(
      "ERROR: GSE184866 file lacks expected columns:\n",
      paste(
        setdiff(
          expected_3d_columns,
          colnames(
            raw_3d
          )
        ),
        collapse = ", "
      )
    )
  )
}

cat(
  "GSE184866 Ensembl rows:",
  nrow(
    raw_3d
  ),
  "\n\n"
)


############################################################
# 12. Map GSE184866 Ensembl IDs to gene symbols
############################################################

sample_cols_3d <- c(
  "3D-Bald_repA",
  "3D-Bald_repB",
  "3D-Bald_repC",
  "3D-NonBald_repA",
  "3D-NonBald_repB",
  "3D-NonBald_repC"
)

ensembl_with_version <- as.character(
  raw_3d$gene_name
)

ensembl_id <- sub(
  "\\..*$",
  "",
  ensembl_with_version
)

counts_3d_raw <- as.matrix(
  raw_3d[
    ,
    sample_cols_3d,
    drop = FALSE
  ]
)

storage.mode(
  counts_3d_raw
) <- "numeric"

if (
  any(
    !is.finite(
      counts_3d_raw
    )
  ) ||
  any(
    counts_3d_raw < 0
  )
) {
  
  stop(
    "ERROR: invalid GSE184866 count values detected."
  )
}

symbol_map_3d <- AnnotationDbi::mapIds(
  
  org.Hs.eg.db,
  
  keys =
    unique(
      ensembl_id
    ),
  
  column =
    "SYMBOL",
  
  keytype =
    "ENSEMBL",
  
  multiVals =
    "first"
)

mapped_symbol_3d <- unname(
  symbol_map_3d[
    ensembl_id
  ]
)

valid_map_3d <- !is.na(
  mapped_symbol_3d
) &
  mapped_symbol_3d != ""

counts_3d_mapped <- counts_3d_raw[
  valid_map_3d,
  ,
  drop = FALSE
]

symbol_3d_valid <- mapped_symbol_3d[
  valid_map_3d
]

counts_3d_gene <- rowsum(
  counts_3d_mapped,
  group =
    symbol_3d_valid,
  reorder = FALSE
)

cat(
  "GSE184866 mapped gene symbols:",
  nrow(
    counts_3d_gene
  ),
  "\n\n"
)

saveRDS(
  counts_3d_gene,
  file.path(
    results_dir,
    "GSE184866_3D_gene_level_counts.rds"
  )
)


############################################################
# 13. GSE184866 metadata
############################################################

meta_3d <- data.frame(
  
  Sample =
    colnames(
      counts_3d_gene
    ),
  
  Group =
    c(
      "BAB",
      "BAB",
      "BAB",
      "BAN",
      "BAN",
      "BAN"
    ),
  
  Culture =
    "3D",
  
  DHT =
    "10nM_24h",
  
  stringsAsFactors =
    FALSE
)

meta_3d$Group <- factor(
  meta_3d$Group,
  levels = c(
    "BAB",
    "BAN"
  )
)

rownames(
  meta_3d
) <- meta_3d$Sample

write.csv(
  meta_3d,
  file.path(
    results_dir,
    "GSE184866_3D_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 14. edgeR/TMM + voom for GSE184866
############################################################

dge_3d <- edgeR::DGEList(
  counts =
    counts_3d_gene,
  group =
    meta_3d$Group
)

design_3d <- model.matrix(
  ~ 0 + Group,
  data = meta_3d
)

colnames(
  design_3d
) <- sub(
  "^Group",
  "",
  colnames(
    design_3d
  )
)

keep_3d <- edgeR::filterByExpr(
  dge_3d,
  design =
    design_3d
)

dge_3d <- dge_3d[
  keep_3d,
  ,
  keep.lib.sizes = FALSE
]

dge_3d <- edgeR::calcNormFactors(
  dge_3d,
  method = "TMM"
)

voom_3d <- limma::voom(
  dge_3d,
  design =
    design_3d,
  plot = FALSE
)

expr_3d <- voom_3d$E

saveRDS(
  expr_3d,
  file.path(
    results_dir,
    "GSE184866_3D_TMM_voom_expression.rds"
  )
)

cat(
  "GSE184866 genes retained after filterByExpr:",
  nrow(
    expr_3d
  ),
  "\n\n"
)


############################################################
# 15. GSE184866 PCA
############################################################

variance_3d <- apply(
  expr_3d,
  1,
  var,
  na.rm = TRUE
)

variance_3d <- sort(
  variance_3d,
  decreasing = TRUE
)

n_pca_3d <- min(
  5000,
  length(
    variance_3d
  )
)

pca_genes_3d <- names(
  variance_3d[
    seq_len(
      n_pca_3d
    )
  ]
)

pca_3d <- prcomp(
  t(
    expr_3d[
      pca_genes_3d,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

pca_3d_var <- 100 *
  pca_3d$sdev^2 /
  sum(
    pca_3d$sdev^2
  )

pca_3d_df <- data.frame(
  
  Sample =
    rownames(
      pca_3d$x
    ),
  
  PC1 =
    pca_3d$x[, 1],
  
  PC2 =
    pca_3d$x[, 2],
  
  stringsAsFactors =
    FALSE
)

pca_3d_df$Group <- meta_3d[
  pca_3d_df$Sample,
  "Group"
]

p_3d_pca <- ggplot(
  
  pca_3d_df,
  
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
      "GSE184866: 3D BAB/BAN after 10 nM DHT for 24 h",
    
    subtitle =
      "Raw gene counts; TMM normalization + voom",
    
    x = paste0(
      "PC1 (",
      round(
        pca_3d_var[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        pca_3d_var[2],
        1
      ),
      "%)"
    ),
    
    color =
      "Cell line"
  ) +
  
  theme_classic(
    base_size = 13
  )

ggsave(
  file.path(
    figures_dir,
    "GSE184866_3D_PCA.png"
  ),
  p_3d_pca,
  width = 8,
  height = 6.5,
  dpi = 300
)


############################################################
# 16. GSE184866 differential expression: BAN - BAB
############################################################

contrast_3d <- limma::makeContrasts(
  BAN_vs_BAB = BAN - BAB,
  levels = design_3d
)

fit_3d <- limma::lmFit(
  voom_3d,
  design_3d
)

fit_3d <- limma::contrasts.fit(
  fit_3d,
  contrast_3d
)

fit_3d <- limma::eBayes(
  fit_3d,
  trend = TRUE,
  robust = TRUE
)

de_3d <- limma::topTable(
  fit_3d,
  coef = "BAN_vs_BAB",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

de_3d$GeneSymbol <- rownames(
  de_3d
)

de_3d <- de_3d[
  ,
  c(
    "GeneSymbol",
    setdiff(
      colnames(
        de_3d
      ),
      "GeneSymbol"
    )
  ),
  drop = FALSE
]

write.csv(
  de_3d,
  file.path(
    results_dir,
    "GSE184866_3D_BAN_vs_BAB_all_genes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 17. Helper: predefined rescue score
############################################################

make_rescue_score <- function(
    expression_matrix,
    signature_table,
    dataset_label,
    metadata
) {
  
  overlap <- signature_table[
    signature_table$GeneSymbol %in%
      rownames(
        expression_matrix
      ),
    ,
    drop = FALSE
  ]
  
  up_genes <- overlap$GeneSymbol[
    overlap$Rescue_type ==
      "Restored_upward"
  ]
  
  down_genes <- overlap$GeneSymbol[
    overlap$Rescue_type ==
      "Restored_downward"
  ]
  
  if (
    length(
      up_genes
    ) == 0 ||
    length(
      down_genes
    ) == 0
  ) {
    
    stop(
      paste0(
        "ERROR: insufficient signature overlap in ",
        dataset_label
      )
    )
  }
  
  z <- t(
    scale(
      t(
        expression_matrix[
          overlap$GeneSymbol,
          ,
          drop = FALSE
        ]
      )
    )
  )
  
  z[
    !is.finite(
      z
    )
  ] <- 0
  
  up_score <- colMeans(
    z[
      up_genes,
      ,
      drop = FALSE
    ]
  )
  
  down_score <- colMeans(
    z[
      down_genes,
      ,
      drop = FALSE
    ]
  )
  
  total_score <-
    up_score -
    down_score
  
  score_df <- data.frame(
    
    Dataset =
      dataset_label,
    
    Sample =
      colnames(
        expression_matrix
      ),
    
    Group =
      as.character(
        metadata[
          colnames(
            expression_matrix
          ),
          "Group"
        ]
      ),
    
    Up_score =
      as.numeric(
        up_score[
          colnames(
            expression_matrix
          )
        ]
      ),
    
    Down_score =
      as.numeric(
        down_score[
          colnames(
            expression_matrix
          )
        ]
      ),
    
    Rescue_signature_score =
      as.numeric(
        total_score[
          colnames(
            expression_matrix
          )
        ]
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  list(
    overlap =
      overlap,
    up_genes =
      up_genes,
    down_genes =
      down_genes,
    score_df =
      score_df
  )
}


############################################################
# 18. Rescue scores in 2D and 3D separately
############################################################

score_2d_object <- make_rescue_score(
  expression_matrix =
    expr_2d,
  signature_table =
    signature,
  dataset_label =
    "GSE184501_2D",
  metadata =
    meta_2d
)

score_3d_object <- make_rescue_score(
  expression_matrix =
    expr_3d,
  signature_table =
    signature,
  dataset_label =
    "GSE184866_3D",
  metadata =
    meta_3d
)

score_2d <- score_2d_object$score_df

score_3d <- score_3d_object$score_df

score_all <- rbind(
  score_2d,
  score_3d
)

score_all$Group <- factor(
  score_all$Group,
  levels = c(
    "BAB",
    "BAN"
  )
)

write.csv(
  score_all,
  file.path(
    results_dir,
    "GSE184501_GSE184866_predefined_rescue_scores.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 19. Descriptive score summaries
############################################################

score_summary_rows <- list()

score_counter <- 1

for (
  current_dataset in unique(
    score_all$Dataset
  )
) {
  
  x <- score_all[
    score_all$Dataset ==
      current_dataset,
    ,
    drop = FALSE
  ]
  
  bab <- x$Rescue_signature_score[
    x$Group ==
      "BAB"
  ]
  
  ban <- x$Rescue_signature_score[
    x$Group ==
      "BAN"
  ]
  
  pooled_sd <- sqrt(
    (
      (
        length(
          bab
        ) - 1
      ) *
        var(
          bab
        ) +
        (
          length(
            ban
          ) - 1
        ) *
        var(
          ban
        )
    ) /
      (
        length(
          bab
        ) +
          length(
            ban
          ) -
          2
      )
  )
  
  cohen_d <- (
    mean(
      ban
    ) -
      mean(
        bab
      )
  ) /
    pooled_sd
  
  score_summary_rows[
    score_counter
  ] <- list(
    data.frame(
      
      Dataset =
        current_dataset,
      
      BAB_mean =
        mean(
          bab
        ),
      
      BAN_mean =
        mean(
          ban
        ),
      
      BAN_minus_BAB =
        mean(
          ban
        ) -
        mean(
          bab
        ),
      
      Descriptive_Cohen_d =
        cohen_d,
      
      stringsAsFactors =
        FALSE
    )
  )
  
  score_counter <-
    score_counter + 1
}

score_summary <- do.call(
  rbind,
  score_summary_rows
)

write.csv(
  score_summary,
  file.path(
    results_dir,
    "GSE184501_GSE184866_rescue_score_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 20. Rescue score figure
#
# IMPORTANT:
# Scores are gene-wise standardized WITHIN each dataset.
# Therefore compare BAB vs BAN within each panel.
# Do not interpret the absolute vertical level across panels.
############################################################

p_scores <- ggplot(
  
  score_all,
  
  aes(
    x = Group,
    y = Rescue_signature_score,
    color = Group,
    label = Sample
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  
  geom_point(
    size = 3.5
  ) +
  
  ggrepel::geom_text_repel(
    size = 2.8,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.55,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ Dataset,
    nrow = 1,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "Predefined GSE44765 rescue signature in 2D and 3D AGA cell-line models",
    
    subtitle =
      "Scores standardized separately within each dataset; compare BAB vs BAN within panel",
    
    x =
      NULL,
    
    y =
      "Predefined rescue signature score"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position =
      "none"
  )

ggsave(
  file.path(
    figures_dir,
    "GSE184501_GSE184866_predefined_rescue_scores.png"
  ),
  p_scores,
  width = 11,
  height = 6,
  dpi = 300
)


############################################################
# 21. CAMERA helper
############################################################

run_camera_for_signature <- function(
    y_object,
    expression_rownames,
    design,
    contrast_vector,
    overlap_table,
    dataset_label
) {
  
  up_genes <- overlap_table$GeneSymbol[
    overlap_table$Rescue_type ==
      "Restored_upward"
  ]
  
  down_genes <- overlap_table$GeneSymbol[
    overlap_table$Rescue_type ==
      "Restored_downward"
  ]
  
  index_list <- list(
    
    Restored_upward =
      which(
        expression_rownames %in%
          up_genes
      ),
    
    Restored_downward =
      which(
        expression_rownames %in%
          down_genes
      )
  )
  
  result <- limma::camera(
    
    y =
      y_object,
    
    index =
      index_list,
    
    design =
      design,
    
    contrast =
      contrast_vector,
    
    inter.gene.cor =
      0.01
  )
  
  result$GeneSet <- rownames(
    result
  )
  
  result$Dataset <-
    dataset_label
  
  result$Contrast <-
    "BAN_vs_BAB"
  
  result <- result[
    ,
    c(
      "Dataset",
      "Contrast",
      "GeneSet",
      setdiff(
        colnames(
          result
        ),
        c(
          "Dataset",
          "Contrast",
          "GeneSet"
        )
      )
    ),
    drop = FALSE
  ]
  
  result
}


############################################################
# 22. CAMERA in 2D and 3D
############################################################

camera_2d <- run_camera_for_signature(
  y_object =
    expr_2d,
  expression_rownames =
    rownames(
      expr_2d
    ),
  design =
    design_2d,
  contrast_vector =
    contrast_2d[
      ,
      "BAN_vs_BAB"
    ],
  overlap_table =
    score_2d_object$overlap,
  dataset_label =
    "GSE184501_2D"
)

camera_3d <- run_camera_for_signature(
  y_object =
    voom_3d,
  expression_rownames =
    rownames(
      expr_3d
    ),
  design =
    design_3d,
  contrast_vector =
    contrast_3d[
      ,
      "BAN_vs_BAB"
    ],
  overlap_table =
    score_3d_object$overlap,
  dataset_label =
    "GSE184866_3D"
)

camera_all <- rbind(
  camera_2d,
  camera_3d
)

write.csv(
  camera_all,
  file.path(
    results_dir,
    "GSE184501_GSE184866_CAMERA_predefined_signature.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 23. Cross-study alignment helper
#
# Positive aligned effect:
#   BAN-BAB moves in the same direction as GSE44765 rescue.
############################################################

make_alignment_table <- function(
    de_table,
    signature_table,
    dataset_label
) {
  
  small <- de_table[
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
    small
  ) <- c(
    "GeneSymbol",
    "External_logFC",
    "External_t",
    "External_P",
    "External_FDR"
  )
  
  merged <- merge(
    signature_table[
      ,
      c(
        "GeneSymbol",
        "Rescue_type",
        "Rescue_logFC"
      ),
      drop = FALSE
    ],
    small,
    by =
      "GeneSymbol"
  )
  
  merged$Dataset <-
    dataset_label
  
  merged$Directional_concordance <- sign(
    merged$Rescue_logFC
  ) ==
    sign(
      merged$External_logFC
    )
  
  merged$RescueAligned_logFC <- sign(
    merged$Rescue_logFC
  ) *
    merged$External_logFC
  
  merged$RescueAligned_t <- sign(
    merged$Rescue_logFC
  ) *
    merged$External_t
  
  merged
}


############################################################
# 24. Build alignment tables
############################################################

alignment_2d <- make_alignment_table(
  de_table =
    de_2d,
  signature_table =
    signature,
  dataset_label =
    "GSE184501_2D"
)

alignment_3d <- make_alignment_table(
  de_table =
    de_3d,
  signature_table =
    signature,
  dataset_label =
    "GSE184866_3D"
)

alignment_all <- rbind(
  alignment_2d,
  alignment_3d
)

write.csv(
  alignment_all,
  file.path(
    results_dir,
    "GSE184501_GSE184866_signature_gene_alignment.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 25. Cross-study concordance summary
############################################################

summarize_alignment <- function(
    x
) {
  
  valid <- is.finite(
    x$Rescue_logFC
  ) &
    is.finite(
      x$External_logFC
    ) &
    is.finite(
      x$External_t
    )
  
  data.frame(
    
    Dataset =
      x$Dataset[1],
    
    N_signature_genes =
      sum(
        valid
      ),
    
    Spearman_rho =
      suppressWarnings(
        cor(
          x$Rescue_logFC[
            valid
          ],
          x$External_logFC[
            valid
          ],
          method =
            "spearman"
        )
      ),
    
    Directional_concordance =
      mean(
        x$Directional_concordance[
          valid
        ],
        na.rm = TRUE
      ),
    
    Mean_rescue_aligned_logFC =
      mean(
        x$RescueAligned_logFC[
          valid
        ],
        na.rm = TRUE
      ),
    
    Median_rescue_aligned_logFC =
      median(
        x$RescueAligned_logFC[
          valid
        ],
        na.rm = TRUE
      ),
    
    Mean_rescue_aligned_t =
      mean(
        x$RescueAligned_t[
          valid
        ],
        na.rm = TRUE
      ),
    
    Median_rescue_aligned_t =
      median(
        x$RescueAligned_t[
          valid
        ],
        na.rm = TRUE
      ),
    
    stringsAsFactors =
      FALSE
  )
}

alignment_summary <- rbind(
  summarize_alignment(
    alignment_2d
  ),
  summarize_alignment(
    alignment_3d
  )
)

write.csv(
  alignment_summary,
  file.path(
    results_dir,
    "GSE184501_GSE184866_cross_study_concordance_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 26. Cross-study concordance figure
############################################################

p_concordance <- ggplot(
  
  alignment_all,
  
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
    size = 1.1
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ Dataset,
    nrow = 1,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "Alignment of BAB/BAN effects with the GSE44765 rescue direction",
    
    subtitle =
      "External effect = BAN - BAB within each dataset",
    
    x =
      "GSE44765 rescue log2FC (P3 3D - P3 2D)",
    
    y =
      "External BAN - BAB log2FC"
  ) +
  
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    figures_dir,
    "GSE184501_GSE184866_cross_study_concordance.png"
  ),
  p_concordance,
  width = 11,
  height = 6.5,
  dpi = 300
)


############################################################
# 27. Compare the SAME signature genes across 2D and 3D
#
# This is an effect-level comparison only.
# It is NOT a direct 2D-versus-3D sample-level test.
############################################################

compare_2d <- alignment_2d[
  ,
  c(
    "GeneSymbol",
    "Rescue_type",
    "External_logFC",
    "External_t",
    "RescueAligned_logFC",
    "RescueAligned_t"
  ),
  drop = FALSE
]

colnames(
  compare_2d
)[
  colnames(
    compare_2d
  ) != "GeneSymbol" &
    colnames(
      compare_2d
    ) != "Rescue_type"
] <- paste0(
  colnames(
    compare_2d
  )[
    colnames(
      compare_2d
    ) != "GeneSymbol" &
      colnames(
        compare_2d
      ) != "Rescue_type"
  ],
  "_2D"
)

compare_3d <- alignment_3d[
  ,
  c(
    "GeneSymbol",
    "Rescue_type",
    "External_logFC",
    "External_t",
    "RescueAligned_logFC",
    "RescueAligned_t"
  ),
  drop = FALSE
]

colnames(
  compare_3d
)[
  colnames(
    compare_3d
  ) != "GeneSymbol" &
    colnames(
      compare_3d
    ) != "Rescue_type"
] <- paste0(
  colnames(
    compare_3d
  )[
    colnames(
      compare_3d
    ) != "GeneSymbol" &
      colnames(
        compare_3d
      ) != "Rescue_type"
  ],
  "_3D"
)

effect_compare <- merge(
  compare_2d,
  compare_3d,
  by = c(
    "GeneSymbol",
    "Rescue_type"
  )
)

effect_compare$AlignmentChange_3D_minus_2D <-
  effect_compare$RescueAligned_t_3D -
  effect_compare$RescueAligned_t_2D

effect_compare$AbsBANBAB_t_change_3D_minus_2D <-
  abs(
    effect_compare$External_t_3D
  ) -
  abs(
    effect_compare$External_t_2D
  )

write.csv(
  effect_compare,
  file.path(
    results_dir,
    "GSE184501_GSE184866_effect_level_2D3D_comparison.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 28. 2D-versus-3D effect architecture summary
############################################################

effect_architecture_summary <- data.frame(
  
  N_common_signature_genes =
    nrow(
      effect_compare
    ),
  
  Spearman_2D_vs_3D_BAN_BAB_logFC =
    suppressWarnings(
      cor(
        effect_compare$External_logFC_2D,
        effect_compare$External_logFC_3D,
        method = "spearman",
        use = "complete.obs"
      )
    ),
  
  Spearman_2D_vs_3D_BAN_BAB_t =
    suppressWarnings(
      cor(
        effect_compare$External_t_2D,
        effect_compare$External_t_3D,
        method = "spearman",
        use = "complete.obs"
      )
    ),
  
  Median_alignment_change_3D_minus_2D =
    median(
      effect_compare$AlignmentChange_3D_minus_2D,
      na.rm = TRUE
    ),
  
  Mean_alignment_change_3D_minus_2D =
    mean(
      effect_compare$AlignmentChange_3D_minus_2D,
      na.rm = TRUE
    ),
  
  Fraction_genes_with_weaker_abs_t_in_3D =
    mean(
      abs(
        effect_compare$External_t_3D
      ) <
        abs(
          effect_compare$External_t_2D
        ),
      na.rm = TRUE
    ),
  
  stringsAsFactors =
    FALSE
)

write.csv(
  effect_architecture_summary,
  file.path(
    results_dir,
    "GSE184501_GSE184866_effect_architecture_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 29. 2D vs 3D rescue-aligned t scatter
############################################################

p_aligned_t <- ggplot(
  
  effect_compare,
  
  aes(
    x = RescueAligned_t_2D,
    y = RescueAligned_t_3D,
    color = Rescue_type
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
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.35,
    linetype = 2
  ) +
  
  geom_point(
    alpha = 0.4,
    size = 1.2
  ) +
  
  labs(
    
    title =
      "Effect-level comparison of the rescue signature in 2D and 3D",
    
    subtitle =
      "Positive values mean BAN-BAB aligns with the GSE44765 rescue direction",
    
    x =
      "GSE184501 2D rescue-aligned moderated t",
    
    y =
      "GSE184866 3D rescue-aligned moderated t",
    
    color =
      "Rescue arm"
  ) +
  
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    figures_dir,
    "GSE184501_GSE184866_rescue_aligned_t_2D_vs_3D.png"
  ),
  p_aligned_t,
  width = 8,
  height = 7,
  dpi = 300
)


############################################################
# 30. Alignment summary by rescue arm
############################################################

arm_summary_rows <- list()

arm_counter <- 1

for (
  current_dataset in c(
    "2D",
    "3D"
  )
) {
  
  current_alignment <- if (
    current_dataset == "2D"
  ) {
    alignment_2d
  } else {
    alignment_3d
  }
  
  for (
    current_arm in c(
      "Restored_upward",
      "Restored_downward"
    )
  ) {
    
    x <- current_alignment[
      current_alignment$Rescue_type ==
        current_arm,
      ,
      drop = FALSE
    ]
    
    arm_summary_rows[
      arm_counter
    ] <- list(
      data.frame(
        
        Culture =
          current_dataset,
        
        Rescue_type =
          current_arm,
        
        N =
          nrow(
            x
          ),
        
        Mean_aligned_t =
          mean(
            x$RescueAligned_t,
            na.rm = TRUE
          ),
        
        Median_aligned_t =
          median(
            x$RescueAligned_t,
            na.rm = TRUE
          ),
        
        Fraction_aligned =
          mean(
            x$RescueAligned_t > 0,
            na.rm = TRUE
          ),
        
        Mean_abs_external_t =
          mean(
            abs(
              x$External_t
            ),
            na.rm = TRUE
          ),
        
        stringsAsFactors =
          FALSE
      )
    )
    
    arm_counter <-
      arm_counter + 1
  }
}

arm_summary <- do.call(
  rbind,
  arm_summary_rows
)

write.csv(
  arm_summary,
  file.path(
    results_dir,
    "GSE184501_GSE184866_rescue_arm_alignment_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 31. Optional predefined 27-gene mitotic/cell-cycle module
############################################################

module_candidates <- list.files(
  path = "results",
  pattern = "^GSE44765_predefined_mitotic_cellcycle_module_genes\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

mitotic_module_genes <- character(0)

if (
  length(
    module_candidates
  ) > 0
) {
  
  module_table <- read.csv(
    module_candidates[1],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  module_gene_column <- NULL
  
  if (
    "SYMBOL" %in%
    colnames(
      module_table
    )
  ) {
    
    module_gene_column <-
      "SYMBOL"
    
  } else if (
    "GeneSymbol" %in%
    colnames(
      module_table
    )
  ) {
    
    module_gene_column <-
      "GeneSymbol"
  }
  
  if (
    !is.null(
      module_gene_column
    )
  ) {
    
    mitotic_module_genes <- unique(
      trimws(
        as.character(
          module_table[
            ,
            module_gene_column
          ]
        )
      )
    )
    
    mitotic_module_genes <- mitotic_module_genes[
      !is.na(
        mitotic_module_genes
      ) &
        mitotic_module_genes != ""
    ]
  }
}


############################################################
# 32. Module summary helper
############################################################

summarize_gene_panel <- function(
    alignment_table,
    genes,
    panel_name,
    culture_label
) {
  
  x <- alignment_table[
    alignment_table$GeneSymbol %in%
      genes,
    ,
    drop = FALSE
  ]
  
  data.frame(
    
    Panel =
      panel_name,
    
    Culture =
      culture_label,
    
    N_present =
      nrow(
        x
      ),
    
    Mean_rescue_aligned_t =
      if (
        nrow(
          x
        ) > 0
      ) {
        mean(
          x$RescueAligned_t,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    
    Median_rescue_aligned_t =
      if (
        nrow(
          x
        ) > 0
      ) {
        median(
          x$RescueAligned_t,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    
    Fraction_rescue_aligned =
      if (
        nrow(
          x
        ) > 0
      ) {
        mean(
          x$RescueAligned_t > 0,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    
    Mean_abs_BAN_BAB_t =
      if (
        nrow(
          x
        ) > 0
      ) {
        mean(
          abs(
            x$External_t
          ),
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    
    stringsAsFactors =
      FALSE
  )
}


############################################################
# 33. Predefined mitotic module and exploratory metabolic panel
############################################################

metabolic_sentinel_genes <- c(
  "ASNS",
  "PHGDH",
  "PSAT1",
  "SLC7A5",
  "GOT1",
  "SHMT2",
  "MTHFD2"
)

panel_summary_list <- list()

panel_counter <- 1

if (
  length(
    mitotic_module_genes
  ) > 0
) {
  
  panel_summary_list[
    panel_counter
  ] <- list(
    summarize_gene_panel(
      alignment_2d,
      mitotic_module_genes,
      "Predefined_GSE44765_mitotic_cellcycle",
      "2D"
    )
  )
  
  panel_counter <-
    panel_counter + 1
  
  panel_summary_list[
    panel_counter
  ] <- list(
    summarize_gene_panel(
      alignment_3d,
      mitotic_module_genes,
      "Predefined_GSE44765_mitotic_cellcycle",
      "3D"
    )
  )
  
  panel_counter <-
    panel_counter + 1
}

panel_summary_list[
  panel_counter
] <- list(
  summarize_gene_panel(
    alignment_2d,
    metabolic_sentinel_genes,
    "Exploratory_GSE66664_metabolic_sentinel",
    "2D"
  )
)

panel_counter <-
  panel_counter + 1

panel_summary_list[
  panel_counter
] <- list(
  summarize_gene_panel(
    alignment_3d,
    metabolic_sentinel_genes,
    "Exploratory_GSE66664_metabolic_sentinel",
    "3D"
  )
)

panel_summary <- do.call(
  rbind,
  panel_summary_list
)

write.csv(
  panel_summary,
  file.path(
    results_dir,
    "GSE184501_GSE184866_mechanistic_panel_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 34. Mechanistic panel figure
############################################################

if (
  nrow(
    panel_summary
  ) > 0
) {
  
  p_panels <- ggplot(
    
    panel_summary,
    
    aes(
      x = Culture,
      y = Mean_rescue_aligned_t,
      group = Panel,
      shape = Panel
    )
    
  ) +
    
    geom_hline(
      yintercept = 0,
      linewidth = 0.35
    ) +
    
    geom_line(
      linewidth = 0.7
    ) +
    
    geom_point(
      size = 3.2
    ) +
    
    labs(
      
      title =
        "Mechanistic panel alignment in 2D and 3D BAB/BAN models",
      
      subtitle =
        "Positive = BAN-BAB aligns with the GSE44765 rescue direction",
      
      x =
        "Culture dataset",
      
      y =
        "Mean rescue-aligned moderated t",
      
      shape =
        "Gene panel"
    ) +
    
    theme_classic(
      base_size = 12
    )
  
  ggsave(
    file.path(
      figures_dir,
      "GSE184501_GSE184866_mechanistic_panel_alignment.png"
    ),
    p_panels,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )
}


############################################################
# 35. Discovery-ranked heatmap gene list
#
# Gene selection uses GSE44765 discovery information only.
############################################################

common_signature_genes <- intersect(
  score_2d_object$overlap$GeneSymbol,
  score_3d_object$overlap$GeneSymbol
)

signature_common <- signature[
  signature$GeneSymbol %in%
    common_signature_genes,
  ,
  drop = FALSE
]

if (
  "Combined_FDR" %in%
  colnames(
    signature_common
  )
) {
  
  signature_common <- signature_common[
    order(
      signature_common$Combined_FDR,
      -abs(
        signature_common$Rescue_logFC
      ),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
} else {
  
  signature_common <- signature_common[
    order(
      -abs(
        signature_common$Rescue_logFC
      )
    ),
    ,
    drop = FALSE
  ]
}

n_heatmap_genes <- min(
  50,
  nrow(
    signature_common
  )
)

heatmap_genes <- signature_common$GeneSymbol[
  seq_len(
    n_heatmap_genes
  )
]


############################################################
# 36. Separate heatmaps for 2D and 3D
############################################################

heatmap_2d <- expr_2d[
  heatmap_genes,
  ,
  drop = FALSE
]

heatmap_2d_z <- t(
  scale(
    t(
      heatmap_2d
    )
  )
)

heatmap_2d_z[
  !is.finite(
    heatmap_2d_z
  )
] <- 0

ann_2d <- data.frame(
  Group =
    as.character(
      meta_2d$Group
    ),
  stringsAsFactors = FALSE
)

rownames(
  ann_2d
) <- meta_2d$Sample

pheatmap(
  heatmap_2d_z,
  annotation_col =
    ann_2d,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 7,
  main =
    "GSE184501 2D: predefined GSE44765 rescue genes",
  filename =
    file.path(
      figures_dir,
      "GSE184501_2D_predefined_signature_heatmap.png"
    ),
  width = 8.5,
  height = 11
)

heatmap_3d <- expr_3d[
  heatmap_genes,
  ,
  drop = FALSE
]

heatmap_3d_z <- t(
  scale(
    t(
      heatmap_3d
    )
  )
)

heatmap_3d_z[
  !is.finite(
    heatmap_3d_z
  )
] <- 0

ann_3d <- data.frame(
  Group =
    as.character(
      meta_3d$Group
    ),
  stringsAsFactors = FALSE
)

rownames(
  ann_3d
) <- meta_3d$Sample

pheatmap(
  heatmap_3d_z,
  annotation_col =
    ann_3d,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 7,
  main =
    "GSE184866 3D: predefined GSE44765 rescue genes",
  filename =
    file.path(
      figures_dir,
      "GSE184866_3D_predefined_signature_heatmap.png"
    ),
  width = 8.5,
  height = 11
)


############################################################
# 37. Input/QC summary
############################################################

input_qc <- data.frame(
  
  Dataset = c(
    "GSE184501_2D",
    "GSE184866_3D"
  ),
  
  InputRows = c(
    nrow(
      raw_2d
    ),
    nrow(
      raw_3d
    )
  ),
  
  GeneLevelRows = c(
    nrow(
      expr_2d
    ),
    nrow(
      expr_3d
    )
  ),
  
  SignatureGenesPresent = c(
    nrow(
      score_2d_object$overlap
    ),
    nrow(
      score_3d_object$overlap
    )
  ),
  
  SignatureCoverage = c(
    nrow(
      score_2d_object$overlap
    ) /
      nrow(
        signature
      ),
    nrow(
      score_3d_object$overlap
    ) /
      nrow(
        signature
      )
  ),
  
  stringsAsFactors =
    FALSE
)

write.csv(
  input_qc,
  file.path(
    results_dir,
    "GSE184501_GSE184866_input_QC_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 38. Human-readable summary
############################################################

summary_file <- file.path(
  results_dir,
  "ANALYSIS_SUMMARY_GSE184501_GSE184866.txt"
)

sink(
  summary_file
)

cat(
  "GSE184501 + GSE184866 2D/3D AGA validation\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "IMPORTANT DESIGN NOTE:\n"
)

cat(
  "GSE184501 and GSE184866 are separate RNA-seq experiments.\n"
)

cat(
  "Their expression matrices were NOT merged.\n"
)

cat(
  "All 2D-versus-3D statements are based on comparison of\n"
)

cat(
  "within-dataset BAN-versus-BAB effects, not direct sample-level\n"
)

cat(
  "2D-versus-3D testing.\n\n"
)

cat(
  "Both datasets contain BAB and BAN cell-line model replicates,\n"
)

cat(
  "not independent patient cohorts.\n\n"
)

cat(
  "Input/QC summary:\n"
)

print(
  input_qc
)

cat("\n")

cat(
  "Predefined rescue-score summary:\n"
)

print(
  score_summary
)

cat("\n")

cat(
  "CAMERA predefined signature results:\n"
)

print(
  camera_all
)

cat("\n")

cat(
  "Cross-study alignment with GSE44765 rescue effect:\n"
)

print(
  alignment_summary
)

cat("\n")

cat(
  "2D-versus-3D effect architecture summary:\n"
)

print(
  effect_architecture_summary
)

cat("\n")

cat(
  "Rescue-arm alignment summary:\n"
)

print(
  arm_summary
)

cat("\n")

cat(
  "Mechanistic panel summary:\n"
)

print(
  panel_summary
)

cat("\n")

cat(
  "Interpretation guide:\n"
)

cat(
  "1. Positive rescue-aligned effects mean BAN-BAB points in the\n"
)

cat(
  "   same direction as the GSE44765 P3_3D-P3_2D rescue effect.\n"
)

cat(
  "2. A smaller positive alignment in the 3D dataset, compared with\n"
)

cat(
  "   the 2D dataset, is consistent with attenuation of BAB/BAN\n"
)

cat(
  "   rescue-axis divergence in 3D, but does NOT prove a causal\n"
)

cat(
  "   3D effect because dataset/platform is confounded with culture.\n"
)

cat(
  "3. A sign reversal in 3D indicates a qualitatively different\n"
)

cat(
  "   BAB/BAN effect architecture rather than simple attenuation.\n"
)

cat(
  "4. The mitotic/cell-cycle module is predefined from GSE44765.\n"
)

cat(
  "5. The metabolic sentinel panel is exploratory and was motivated\n"
)

cat(
  "   by the independent GSE66664 DHT anti-rescue analysis.\n"
)

cat(
  "6. Do not interpret replicate-level separation as patient-level\n"
)

cat(
  "   clinical validation.\n"
)

sink()


############################################################
# 39. Save workspace
############################################################

save(
  
  signature,
  raw_2d,
  expr_2d,
  meta_2d,
  pca_2d,
  de_2d,
  raw_3d,
  counts_3d_gene,
  dge_3d,
  voom_3d,
  expr_3d,
  meta_3d,
  pca_3d,
  de_3d,
  score_2d_object,
  score_3d_object,
  score_all,
  score_summary,
  camera_all,
  alignment_2d,
  alignment_3d,
  alignment_all,
  alignment_summary,
  effect_compare,
  effect_architecture_summary,
  arm_summary,
  mitotic_module_genes,
  metabolic_sentinel_genes,
  panel_summary,
  input_qc,
  
  file = file.path(
    results_dir,
    "GSE184501_GSE184866_workspace.RData"
  )
)


############################################################
# 40. Session information
############################################################

sink(
  file.path(
    results_dir,
    "sessionInfo_GSE184501_GSE184866.txt"
  )
)

sessionInfo()

sink()


############################################################
# 41. Final console report
############################################################

cat("\n")
cat("============================================================\n")
cat(" GSE184501 + GSE184866 VALIDATION COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "2D gene-level features:",
  nrow(
    expr_2d
  ),
  "\n"
)

cat(
  "3D genes after edgeR filtering:",
  nrow(
    expr_3d
  ),
  "\n"
)

cat(
  "2D signature coverage:",
  nrow(
    score_2d_object$overlap
  ),
  "/",
  nrow(
    signature
  ),
  " (",
  round(
    100 *
      nrow(
        score_2d_object$overlap
      ) /
      nrow(
        signature
      ),
    1
  ),
  "%)\n",
  sep = ""
)

cat(
  "3D signature coverage:",
  nrow(
    score_3d_object$overlap
  ),
  "/",
  nrow(
    signature
  ),
  " (",
  round(
    100 *
      nrow(
        score_3d_object$overlap
      ) /
      nrow(
        signature
      ),
    1
  ),
  "%)\n",
  sep = ""
)

cat("\nRescue-score summary:\n")
print(
  score_summary
)

cat("\nCAMERA:\n")
print(
  camera_all
)

cat("\nGSE44765 alignment:\n")
print(
  alignment_summary
)

cat("\nEffect architecture:\n")
print(
  effect_architecture_summary
)

cat("\nMechanistic panels:\n")
print(
  panel_summary
)

cat("\nMain figures:\n")

cat(
  file.path(
    figures_dir,
    "GSE184501_2D_PCA.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE184866_3D_PCA.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE184501_GSE184866_predefined_rescue_scores.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE184501_GSE184866_cross_study_concordance.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE184501_GSE184866_rescue_aligned_t_2D_vs_3D.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE184501_GSE184866_mechanistic_panel_alignment.png"
  ),
  "\n"
)

cat("\nSummary:\n")
cat(
  summary_file,
  "\n"
)

cat("\nIMPORTANT:\n")
cat(
  "Do not directly interpret the two datasets as a randomized 2D-vs-3D experiment.\n"
)

cat("============================================================\n")
