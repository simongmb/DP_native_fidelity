############################################################
# GSE149189 independent developmental / cell-type validation
# COMPLETE ONE-CLICK VERSION
#
# Project:
#   HairRegeneration_AI
#
# Required raw files in:
#   data_raw/GSE149189/
#
# Expected files:
#   GSM4491915_Fetus_DP_ELW24.txt.gz
#   GSM4491916_Fetus_DP_ELW26.txt.gz
#   GSM4491917_Fetus_DSC_ELW25.txt.gz
#   GSM4491918_Fetus_DSC_ELW27.txt.gz
#   GSM4491919_Fetus_IFD_ELW28.txt.gz
#   GSM4491920_Fetus_IFD_ELW29.txt.gz
#   GSM4491921_Adult_DP_ELW30.txt.gz
#   GSM4491922_Adult_DP_ELW33.txt.gz
#
# Optional processed file:
#   GSE149189_CPM_FetusAdult_16623Gene.txt.gz
#
# Required discovery signature:
#   results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv
#
# Main questions:
#   1. Does the predefined GSE44765 rescue signature differ
#      between fetal DP and adult DP?
#   2. Is the signature preferentially associated with DP
#      compared with fetal DSC / IFD?
#   3. Are external gene-level effects concordant with the
#      discovery rescue effect from GSE44765?
#
# IMPORTANT:
#   GSE149189 is used only for external validation.
#   It is NOT used to redefine or tune the discovery signature.
############################################################


############################################################
# 0. Basic settings
############################################################

graphics.off()
options(stringsAsFactors = FALSE)
set.seed(149189)

cat("\n")
cat("============================================================\n")
cat(" GSE149189 external validation started\n")
cat("============================================================\n\n")
cat("Working directory:\n")
cat(getwd(), "\n\n")


############################################################
# 1. Paths
############################################################

raw_dir <-
  "data_raw/GSE149189"

discovery_signature_file <-
  "results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv"

results_dir <-
  "results/GSE149189_external"

figures_dir <-
  "figures/GSE149189_external"

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
  "edgeR",
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
  
  library(edgeR)
  
  library(limma)
  
  library(ggplot2)
  
  library(ggrepel)
  
  library(pheatmap)
  
})

cat("Packages loaded successfully.\n\n")


############################################################
# 6. Check discovery signature
############################################################

if (!file.exists(
  discovery_signature_file
)) {
  
  stop(
    paste0(
      "\nERROR:\nCannot find:\n",
      discovery_signature_file,
      "\n\nRun GSE44765 rescue analysis first.\n"
    )
  )
}


############################################################
# 7. Expected sample metadata
############################################################

meta <- data.frame(
  
  GSM = c(
    "GSM4491915",
    "GSM4491916",
    "GSM4491917",
    "GSM4491918",
    "GSM4491919",
    "GSM4491920",
    "GSM4491921",
    "GSM4491922"
  ),
  
  sample = c(
    "Fetus_DP_ELW24",
    "Fetus_DP_ELW26",
    "Fetus_DSC_ELW25",
    "Fetus_DSC_ELW27",
    "Fetus_IFD_ELW28",
    "Fetus_IFD_ELW29",
    "Adult_DP_ELW30",
    "Adult_DP_ELW33"
  ),
  
  group = c(
    "Fetus_DP",
    "Fetus_DP",
    "Fetus_DSC",
    "Fetus_DSC",
    "Fetus_IFD",
    "Fetus_IFD",
    "Adult_DP",
    "Adult_DP"
  ),
  
  stringsAsFactors = FALSE
)

rownames(meta) <-
  meta$GSM

meta$group <- factor(
  meta$group,
  levels = c(
    "Adult_DP",
    "Fetus_DP",
    "Fetus_DSC",
    "Fetus_IFD"
  )
)

write.csv(
  meta,
  file.path(
    results_dir,
    "GSE149189_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 8. Locate raw count files
############################################################

all_count_files <- list.files(
  path = raw_dir,
  pattern = "^GSM44919[0-9]+.*\\.txt\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat(
  "Candidate raw count files detected:",
  length(all_count_files),
  "\n"
)


############################################################
# 9. Find exactly one file per GSM
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
        "\nERROR: expected exactly one raw count file for ",
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
  
  meta$GSM,
  
  find_file_for_gsm,
  
  file_vector =
    all_count_files,
  
  FUN.VALUE =
    character(1)
)

names(ordered_files) <-
  meta$GSM


############################################################
# 10. Read one two-column count file
############################################################

read_count_file <- function(
    filename
) {
  
  x <- read.delim(
    filename,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )
  
  if (ncol(x) != 2) {
    
    stop(
      paste0(
        "ERROR: ",
        basename(filename),
        " does not contain exactly two columns."
      )
    )
  }
  
  colnames(x) <- c(
    "GeneSymbol",
    "Count"
  )
  
  x$GeneSymbol <-
    as.character(
      x$GeneSymbol
    )
  
  x$Count <-
    suppressWarnings(
      as.numeric(
        x$Count
      )
    )
  
  if (
    any(is.na(x$Count)) ||
    any(x$Count < 0)
  ) {
    
    stop(
      paste0(
        "ERROR: invalid counts detected in ",
        basename(filename)
      )
    )
  }
  
  if (anyDuplicated(
    x$GeneSymbol
  ) > 0) {
    
    stop(
      paste0(
        "ERROR: duplicated gene names in ",
        basename(filename)
      )
    )
  }
  
  x
}


############################################################
# 11. Read all 8 count files
############################################################

count_tables <- lapply(
  ordered_files,
  read_count_file
)

reference_genes <-
  getElement(
    count_tables,
    1
  )$GeneSymbol


############################################################
# 12. Confirm identical gene order in all files
############################################################

for (i in seq_along(
  count_tables
)) {
  
  current_genes <-
    getElement(
      count_tables,
      i
    )$GeneSymbol
  
  if (!identical(
    reference_genes,
    current_genes
  )) {
    
    stop(
      paste0(
        "ERROR: gene order differs in ",
        basename(
          ordered_files[i]
        ),
        "."
      )
    )
  }
}


############################################################
# 13. Build raw count matrix
############################################################

count_matrix <- do.call(
  
  cbind,
  
  lapply(
    count_tables,
    function(x) {
      x$Count
    }
  )
)

rownames(count_matrix) <-
  reference_genes

colnames(count_matrix) <-
  meta$GSM

storage.mode(
  count_matrix
) <- "integer"


############################################################
# 14. Basic raw-count QC
############################################################

library_sizes <- colSums(
  count_matrix
)

raw_qc <- data.frame(
  
  GSM =
    colnames(
      count_matrix
    ),
  
  Sample =
    meta[
      colnames(count_matrix),
      "sample"
    ],
  
  Group =
    as.character(
      meta[
        colnames(count_matrix),
        "group"
      ]
    ),
  
  LibrarySize =
    as.numeric(
      library_sizes
    ),
  
  DetectedGenes =
    colSums(
      count_matrix > 0
    ),
  
  stringsAsFactors = FALSE
)

write.csv(
  raw_qc,
  file.path(
    results_dir,
    "GSE149189_raw_count_QC.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Raw count matrix dimensions:\n"
)

print(
  dim(count_matrix)
)

cat(
  "\nLibrary sizes:\n"
)

print(
  library_sizes
)

cat("\n")


############################################################
# 15. Save raw count matrix
############################################################

saveRDS(
  count_matrix,
  file.path(
    results_dir,
    "GSE149189_raw_count_matrix.rds"
  )
)


############################################################
# 16. Build edgeR DGEList
############################################################

dge <- edgeR::DGEList(
  counts =
    count_matrix,
  group =
    meta$group
)


############################################################
# 17. Filter low-expression genes
############################################################

design_filter <- model.matrix(
  ~ 0 + group,
  data = meta
)

keep <- edgeR::filterByExpr(
  dge,
  design =
    design_filter
)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

cat(
  "Genes retained after filterByExpr:",
  nrow(dge),
  "\n\n"
)


############################################################
# 18. TMM normalization
############################################################

dge <- edgeR::calcNormFactors(
  dge,
  method = "TMM"
)


############################################################
# 19. Log2 CPM for QC and visualization
############################################################

log_cpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 2
)

saveRDS(
  log_cpm,
  file.path(
    results_dir,
    "GSE149189_TMM_logCPM.rds"
  )
)

logcpm_connection <- gzfile(
  file.path(
    results_dir,
    "GSE149189_TMM_logCPM.csv.gz"
  ),
  open = "wt"
)

write.csv(
  log_cpm,
  logcpm_connection,
  quote = FALSE
)

close(
  logcpm_connection
)


############################################################
# 20. Optional cross-check against supplied CPM file
############################################################

cpm_candidates <- list.files(
  path = raw_dir,
  pattern = "GSE149189_CPM.*\\.txt\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(cpm_candidates) >= 1) {
  
  supplied_cpm <- read.delim(
    cpm_candidates[1],
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  if (
    "refGene" %in%
    colnames(
      supplied_cpm
    )
  ) {
    
    rownames(
      supplied_cpm
    ) <-
      supplied_cpm$refGene
    
    supplied_cpm$refGene <-
      NULL
    
    sample_name_map <- c(
      GSM4491915 = "Fetus_DP_ELW24",
      GSM4491916 = "Fetus_DP_ELW26",
      GSM4491917 = "Fetus_DSC_ELW25",
      GSM4491918 = "Fetus_DSC_ELW27",
      GSM4491919 = "Fetus_IFD_ELW28",
      GSM4491920 = "Fetus_IFD_ELW29",
      GSM4491921 = "Adult_DP_ELW30",
      GSM4491922 = "Adult_DP_ELW33"
    )
    
    cpm_check_rows <- list()
    cpm_counter <- 1
    
    for (gsm in meta$GSM) {
      
      supplied_name <-
        unname(
          sample_name_map[gsm]
        )
      
      common_genes <- intersect(
        rownames(
          supplied_cpm
        ),
        rownames(
          count_matrix
        )
      )
      
      raw_standard_cpm <-
        count_matrix[
          common_genes,
          gsm
        ] /
        sum(
          count_matrix[
            ,
            gsm
          ]
        ) *
        1e6
      
      if (
        supplied_name %in%
        colnames(
          supplied_cpm
        )
      ) {
        
        supplied_values <-
          supplied_cpm[
            common_genes,
            supplied_name
          ]
        
        cpm_check_rows[
          cpm_counter
        ] <- list(
          data.frame(
            
            GSM =
              gsm,
            
            SuppliedColumn =
              supplied_name,
            
            NCommonGenes =
              length(
                common_genes
              ),
            
            PearsonCorrelation =
              suppressWarnings(
                cor(
                  raw_standard_cpm,
                  supplied_values,
                  method = "pearson",
                  use = "complete.obs"
                )
              ),
            
            SpearmanCorrelation =
              suppressWarnings(
                cor(
                  raw_standard_cpm,
                  supplied_values,
                  method = "spearman",
                  use = "complete.obs"
                )
              ),
            
            stringsAsFactors =
              FALSE
          )
        )
        
        cpm_counter <-
          cpm_counter + 1
      }
    }
    
    if (length(
      cpm_check_rows
    ) > 0) {
      
      cpm_check <- do.call(
        rbind,
        cpm_check_rows
      )
      
      write.csv(
        cpm_check,
        file.path(
          results_dir,
          "GSE149189_supplied_CPM_crosscheck.csv"
        ),
        row.names = FALSE,
        quote = FALSE
      )
    }
  }
}


############################################################
# 21. PCA using top 5000 variable genes
############################################################

gene_variance <- apply(
  log_cpm,
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
  length(
    gene_variance
  )
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
    log_cpm[
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

pca_df$Sample <-
  meta[
    pca_df$GSM,
    "sample"
  ]

pca_df$Group <-
  meta[
    pca_df$GSM,
    "group"
  ]


############################################################
# 22. PCA plot
############################################################

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
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  labs(
    
    title =
      "GSE149189: developmental and cell-type structure",
    
    subtitle =
      "TMM-normalized RNA-seq; top 5,000 variable genes",
    
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
      "Group"
  ) +
  
  theme_classic(
    base_size = 13
  )

ggsave(
  file.path(
    figures_dir,
    "GSE149189_PCA.png"
  ),
  p_pca,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE149189_PCA.pdf"
  ),
  p_pca,
  width = 9,
  height = 7
)


############################################################
# 23. Sample-correlation heatmap
############################################################

sample_cor <- cor(
  log_cpm[
    pca_genes,
    ,
    drop = FALSE
  ],
  method = "pearson"
)

heatmap_annotation <- data.frame(
  
  Group =
    as.character(
      meta$group
    ),
  
  stringsAsFactors = FALSE
)

rownames(
  heatmap_annotation
) <- meta$GSM

pheatmap(
  sample_cor,
  annotation_col =
    heatmap_annotation,
  annotation_row =
    heatmap_annotation,
  border_color = NA,
  main =
    "GSE149189 sample correlation",
  filename =
    file.path(
      figures_dir,
      "GSE149189_sample_correlation_heatmap.png"
    ),
  width = 9,
  height = 8
)


############################################################
# 24. Load predefined discovery signature
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
  colnames(
    signature
  )
)) {
  
  stop(
    "ERROR: discovery signature lacks required columns."
  )
}

signature$GeneSymbol <-
  as.character(
    signature$GeneSymbol
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


############################################################
# 25. Signature overlap
############################################################

signature_overlap <- signature[
  signature$GeneSymbol %in%
    rownames(
      log_cpm
    ),
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

signature_coverage <-
  nrow(
    signature_overlap
  ) /
  nrow(
    signature
  )

write.csv(
  signature_overlap,
  file.path(
    results_dir,
    "GSE44765_signature_overlap_in_GSE149189.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Discovery signature genes:",
  nrow(signature),
  "\n"
)

cat(
  "Signature overlap after RNA-seq filtering:",
  nrow(signature_overlap),
  "\n"
)

cat(
  "Signature coverage:",
  round(
    100 *
      signature_coverage,
    1
  ),
  "%\n\n"
)


############################################################
# 26. Predefined score across all 8 samples
############################################################

signature_expression <- log_cpm[
  signature_overlap$GeneSymbol,
  ,
  drop = FALSE
]

signature_z <- t(
  scale(
    t(
      signature_expression
    )
  )
)

signature_z[
  !is.finite(
    signature_z
  )
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

total_score <-
  up_score -
  down_score

score_df <- data.frame(
  
  GSM =
    colnames(
      log_cpm
    ),
  
  Sample =
    meta[
      colnames(log_cpm),
      "sample"
    ],
  
  Group =
    meta[
      colnames(log_cpm),
      "group"
    ],
  
  Up_score =
    as.numeric(
      up_score[
        colnames(
          log_cpm
        )
      ]
    ),
  
  Down_score =
    as.numeric(
      down_score[
        colnames(
          log_cpm
        )
      ]
    ),
  
  Rescue_signature_score =
    as.numeric(
      total_score[
        colnames(
          log_cpm
        )
      ]
    ),
  
  stringsAsFactors = FALSE
)

score_df$Group <- factor(
  score_df$Group,
  levels = c(
    "Adult_DP",
    "Fetus_DP",
    "Fetus_DSC",
    "Fetus_IFD"
  )
)

write.csv(
  score_df,
  file.path(
    results_dir,
    "GSE149189_predefined_signature_scores.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 27. Signature score plot
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
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  labs(
    
    title =
      "Predefined GSE44765 rescue signature in GSE149189",
    
    subtitle =
      "No GSE149189 labels were used to define the signature",
    
    x =
      NULL,
    
    y =
      "Predefined rescue signature score"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    legend.position =
      "none",
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figures_dir,
    "GSE149189_predefined_signature_scores.png"
  ),
  p_score,
  width = 9,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE149189_predefined_signature_scores.pdf"
  ),
  p_score,
  width = 9,
  height = 6.5
)


############################################################
# 28. Voom-limma model
############################################################

design <- model.matrix(
  ~ 0 + group,
  data = meta
)

if (qr(
  design
)$rank != ncol(
  design
)) {
  
  stop(
    "ERROR: design matrix is not full rank."
  )
}

voom_object <- limma::voom(
  dge,
  design =
    design,
  plot = FALSE
)

fit <- limma::lmFit(
  voom_object,
  design
)


############################################################
# 29. External-validation contrasts
#
# PRIMARY:
#   Fetus_DP vs Adult_DP
#
# SECONDARY specificity:
#   Fetus_DP vs Fetus_DSC
#   Fetus_DP vs Fetus_IFD
#   Fetus_DP vs average fetal non-DP compartment
############################################################

contrast_matrix <- limma::makeContrasts(
  
  Fetus_DP_vs_Adult_DP =
    groupFetus_DP -
    groupAdult_DP,
  
  Fetus_DP_vs_Fetus_DSC =
    groupFetus_DP -
    groupFetus_DSC,
  
  Fetus_DP_vs_Fetus_IFD =
    groupFetus_DP -
    groupFetus_IFD,
  
  Fetus_DP_vs_Fetal_NonDP =
    groupFetus_DP -
    (
      groupFetus_DSC +
        groupFetus_IFD
    ) / 2,
  
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
# 30. Save DE tables
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
    rownames(
      tt
    )
  
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
  
  external_de_list[
    contrast_name
  ] <- list(
    tt
  )
  
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
# 31. CAMERA predefined gene-set validation
############################################################

signature_index <- list(
  
  Restored_upward =
    which(
      rownames(
        voom_object$E
      ) %in%
        up_genes
    ),
  
  Restored_downward =
    which(
      rownames(
        voom_object$E
      ) %in%
        down_genes
    )
)

camera_list <- list()
camera_counter <- 1

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
      voom_object,
    
    index =
      signature_index,
    
    design =
      design,
    
    contrast =
      contrast_matrix[
        ,
        j
      ],
    
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
        colnames(
          camera_result
        ),
        c(
          "Contrast",
          "GeneSet"
        )
      )
    ),
    drop = FALSE
  ]
  
  camera_list[
    camera_counter
  ] <- list(
    camera_result
  )
  
  camera_counter <-
    camera_counter + 1
}

camera_results <- do.call(
  rbind,
  camera_list
)

rownames(
  camera_results
) <- NULL

write.csv(
  camera_results,
  file.path(
    results_dir,
    "GSE149189_CAMERA_predefined_signature_validation.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 32. Cross-study gene-level concordance
############################################################

discovery_effect <- signature_overlap[
  ,
  c(
    "GeneSymbol",
    "Rescue_logFC",
    "Rescue_type"
  ),
  drop = FALSE
]

concordance_rows <- list()
all_merged_rows <- list()
concordance_counter <- 1

for (contrast_name in contrast_names) {
  
  external_tt <-
    getElement(
      external_de_list,
      contrast_name
    )
  
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
    by =
      "GeneSymbol"
  )
  
  valid <-
    is.finite(
      merged$Rescue_logFC
    ) &
    is.finite(
      merged$External_logFC
    )
  
  merged$Concordant_direction <-
    sign(
      merged$Rescue_logFC
    ) ==
    sign(
      merged$External_logFC
    )
  
  merged$Expected_signed_external_logFC <-
    sign(
      merged$Rescue_logFC
    ) *
    merged$External_logFC
  
  rho <- suppressWarnings(
    cor(
      merged$Rescue_logFC[
        valid
      ],
      merged$External_logFC[
        valid
      ],
      method =
        "spearman"
    )
  )
  
  directional_concordance <- mean(
    merged$Concordant_direction[
      valid
    ],
    na.rm = TRUE
  )
  
  concordance_rows[
    concordance_counter
  ] <- list(
    data.frame(
      
      Contrast =
        contrast_name,
      
      N_signature_genes =
        sum(
          valid
        ),
      
      Spearman_rho =
        rho,
      
      Directional_concordance =
        directional_concordance,
      
      Median_expected_signed_external_logFC =
        median(
          merged$Expected_signed_external_logFC[
            valid
          ],
          na.rm = TRUE
        ),
      
      Mean_expected_signed_external_logFC =
        mean(
          merged$Expected_signed_external_logFC[
            valid
          ],
          na.rm = TRUE
        ),
      
      stringsAsFactors =
        FALSE
    )
  )
  
  merged$Contrast <-
    contrast_name
  
  all_merged_rows[
    concordance_counter
  ] <- list(
    merged
  )
  
  concordance_counter <-
    concordance_counter + 1
}

concordance_summary <- do.call(
  rbind,
  concordance_rows
)

all_merged_effects <- do.call(
  rbind,
  all_merged_rows
)

rownames(
  all_merged_effects
) <- NULL

write.csv(
  concordance_summary,
  file.path(
    results_dir,
    "GSE149189_cross_study_concordance_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  all_merged_effects,
  file.path(
    results_dir,
    "GSE149189_cross_study_gene_effects.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 33. Cross-study concordance figure
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
    ncol = 2,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "GSE149189 cross-study effect concordance",
    
    subtitle =
      "Discovery effect: GSE44765 P3 3D - P3 2D",
    
    x =
      "GSE44765 discovery rescue log2FC",
    
    y =
      "GSE149189 external log2FC"
  ) +
  
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    figures_dir,
    "GSE149189_cross_study_effect_concordance.png"
  ),
  p_concordance,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE149189_cross_study_effect_concordance.pdf"
  ),
  p_concordance,
  width = 10,
  height = 8
)


############################################################
# 34. Predefined signature heatmap
#
# Genes are ranked ONLY by GSE44765 discovery statistics.
############################################################

if (
  "Combined_FDR" %in%
  colnames(
    signature_overlap
  )
) {
  
  signature_for_heatmap <- signature_overlap[
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
  
  signature_for_heatmap <- signature_overlap[
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

heatmap_expr <- log_cpm[
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
    heatmap_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  main =
    "Predefined GSE44765 strict rescue signature in GSE149189",
  filename =
    file.path(
      figures_dir,
      "GSE149189_predefined_signature_heatmap.png"
    ),
  width = 10,
  height = 11
)

pheatmap(
  heatmap_z,
  annotation_col =
    heatmap_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  main =
    "Predefined GSE44765 strict rescue signature in GSE149189",
  filename =
    file.path(
      figures_dir,
      "GSE149189_predefined_signature_heatmap.pdf"
    ),
  width = 10,
  height = 11
)


############################################################
# 35. Group-level signature score summary
############################################################

group_levels <- levels(
  score_df$Group
)

score_summary_list <- lapply(
  
  group_levels,
  
  function(current_group) {
    
    values <- score_df$Rescue_signature_score[
      score_df$Group ==
        current_group
    ]
    
    data.frame(
      
      Group =
        current_group,
      
      N =
        length(
          values
        ),
      
      Mean =
        mean(
          values
        ),
      
      SD =
        sd(
          values
        ),
      
      Median =
        median(
          values
        ),
      
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
    "GSE149189_signature_score_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 36. Save design and contrast matrices
############################################################

write.csv(
  design,
  file.path(
    results_dir,
    "GSE149189_design_matrix.csv"
  ),
  quote = FALSE
)

write.csv(
  contrast_matrix,
  file.path(
    results_dir,
    "GSE149189_contrast_matrix.csv"
  ),
  quote = FALSE
)


############################################################
# 37. Human-readable summary
############################################################

summary_file <- file.path(
  results_dir,
  "ANALYSIS_SUMMARY_GSE149189.txt"
)

sink(
  summary_file
)

cat(
  "GSE149189 external validation\n"
)

cat(
  "====================================\n\n"
)

cat(
  "Raw samples analyzed: 8\n"
)

cat(
  "Groups: Adult_DP, Fetus_DP, Fetus_DSC, Fetus_IFD\n\n"
)

cat(
  "Raw genes:",
  nrow(
    count_matrix
  ),
  "\n"
)

cat(
  "Genes retained after edgeR filterByExpr:",
  nrow(
    dge
  ),
  "\n\n"
)

cat(
  "Discovery strict rescue genes:",
  nrow(
    signature
  ),
  "\n"
)

cat(
  "Signature genes represented after filtering:",
  nrow(
    signature_overlap
  ),
  "\n"
)

cat(
  "Signature coverage:",
  round(
    100 *
      signature_coverage,
    1
  ),
  "%\n"
)

cat(
  "Restored-upward genes represented:",
  length(
    up_genes
  ),
  "\n"
)

cat(
  "Restored-downward genes represented:",
  length(
    down_genes
  ),
  "\n\n"
)

cat(
  "Predefined signature score by group:\n"
)

print(
  score_summary
)

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

cat("\n")

cat(
  "Interpretation hierarchy:\n"
)

cat(
  "1. Fetus_DP_vs_Adult_DP is the PRIMARY developmental comparison.\n"
)

cat(
  "2. Fetus_DP_vs_Fetus_DSC and Fetus_DP_vs_Fetus_IFD are SECONDARY cell-type specificity comparisons.\n"
)

cat(
  "3. Do not redefine the GSE44765 signature based on this dataset.\n"
)

sink()


############################################################
# 38. Save workspace
############################################################

save(
  
  count_matrix,
  dge,
  log_cpm,
  meta,
  pca,
  signature,
  signature_overlap,
  score_df,
  score_summary,
  design,
  contrast_matrix,
  fit2,
  camera_results,
  concordance_summary,
  all_merged_effects,
  
  file =
    file.path(
      results_dir,
      "GSE149189_external_validation_workspace.RData"
    )
)


############################################################
# 39. Save session information
############################################################

sink(
  file.path(
    results_dir,
    "sessionInfo_GSE149189.txt"
  )
)

sessionInfo()

sink()


############################################################
# 40. Final console report
############################################################

cat("\n")
cat("============================================================\n")
cat(" GSE149189 EXTERNAL VALIDATION COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "Genes retained:",
  nrow(
    dge
  ),
  "\n"
)

cat(
  "Signature coverage:",
  nrow(
    signature_overlap
  ),
  "/",
  nrow(
    signature
  ),
  " (",
  round(
    100 *
      signature_coverage,
    1
  ),
  "%)\n",
  sep = ""
)

cat("\nSignature score summary:\n")
print(
  score_summary
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
    "GSE149189_PCA.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE149189_predefined_signature_scores.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE149189_cross_study_effect_concordance.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE149189_predefined_signature_heatmap.png"
  ),
  "\n"
)

cat("\nMain summary:\n")
cat(
  summary_file,
  "\n"
)

cat("\nIMPORTANT:\n")
cat(
  "GSE149189 is an external validation dataset only.\n"
)
cat(
  "Do not tune the discovery signature using these results.\n"
)
cat("============================================================\n")
