############################################################
# GSE44765
# Dermal papilla transcriptional rescue analysis
#
# COMPLETE REVISED VERSION
#
# Required input:
# results/GSE44765_RMA_probe_expression.rds
#
# Output:
# results/GSE44765_rescue/
# figures/GSE44765_rescue/
############################################################


############################################################
# 0. Start clean
############################################################

graphics.off()

options(stringsAsFactors = FALSE)

cat("\n")
cat("========================================================\n")
cat(" GSE44765 rescue-signature analysis started\n")
cat("========================================================\n\n")

cat("Working directory:\n")
cat(getwd(), "\n\n")


############################################################
# 1. Fixed analysis thresholds
############################################################

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.50

TOP_HEATMAP_GENES <- 40
TOP_LABEL_GENES <- 12
TOP_TRAJECTORY_GENES <- 12


############################################################
# 2. Output folders
############################################################

dir.create(
  "results/GSE44765_rescue",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/GSE44765_rescue",
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Install CRAN packages
############################################################

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap"
)

for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}


############################################################
# 4. Install Bioconductor
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


############################################################
# 5. Install Bioconductor packages
############################################################

bioc_packages <- c(
  "limma",
  "AnnotationDbi",
  "hgu133plus2.db"
)

for (pkg in bioc_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


############################################################
# 6. Load packages
############################################################

suppressPackageStartupMessages({
  
  library(limma)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  
})

cat("Packages loaded successfully.\n\n")


############################################################
# 7. Check input file
############################################################

expression_file <-
  "results/GSE44765_RMA_probe_expression.rds"

if (!file.exists(expression_file)) {
  
  stop(
    paste0(
      "\nERROR:\n",
      "Cannot find:\n",
      expression_file,
      "\n\n",
      "Please run Script 01 first.\n"
    )
  )
}


############################################################
# 8. Load RMA matrix
############################################################

expr_probe <- readRDS(
  expression_file
)

cat("Probe-level matrix dimensions:\n")
print(dim(expr_probe))
cat("\n")


############################################################
# 9. Check samples
############################################################

expected_samples <- paste0(
  "DP",
  1:18
)

if (!all(
  expected_samples %in% colnames(expr_probe)
)) {
  
  stop(
    "\nERROR: DP1-DP18 were not all detected.\n"
  )
}

expr_probe <- expr_probe[
  ,
  expected_samples,
  drop = FALSE
]


############################################################
# 10. Build metadata
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
# 11. Factor ordering
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
# 12. Match metadata
############################################################

meta <- meta[
  colnames(expr_probe),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(meta),
    colnames(expr_probe)
  )
)

cat("Metadata matched successfully.\n\n")


############################################################
# 13. Probe-to-gene annotation
############################################################

cat("Mapping probes to gene symbols...\n\n")

probe_ids <- rownames(expr_probe)

gene_symbol <- AnnotationDbi::mapIds(
  
  hgu133plus2.db,
  
  keys = probe_ids,
  
  column = "SYMBOL",
  
  keytype = "PROBEID",
  
  multiVals = "first"
)


############################################################
# 14. Save probe annotation
############################################################

annotation_table <- data.frame(
  
  ProbeID = probe_ids,
  
  GeneSymbol = as.character(
    gene_symbol
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  
  annotation_table,
  
  file =
    "results/GSE44765_rescue/GSE44765_probe_annotation.csv",
  
  row.names = FALSE,
  
  quote = FALSE
)


############################################################
# 15. Remove unannotated probes
############################################################

valid_probe <-
  !is.na(gene_symbol) &
  gene_symbol != ""

expr_annotated <- expr_probe[
  valid_probe,
  ,
  drop = FALSE
]

gene_symbol_valid <-
  as.character(
    gene_symbol[
      valid_probe
    ]
  )

cat(
  "Annotated probes:",
  nrow(expr_annotated),
  "\n"
)


############################################################
# 16. Collapse probes to genes
############################################################

expr_gene <- limma::avereps(
  
  expr_annotated,
  
  ID = gene_symbol_valid
)

cat(
  "Unique genes:",
  nrow(expr_gene),
  "\n\n"
)


############################################################
# 17. Check matrix
############################################################

if (
  any(is.na(expr_gene)) ||
  any(!is.finite(expr_gene))
) {
  
  stop(
    "\nERROR: Gene matrix contains NA/Inf.\n"
  )
}


############################################################
# 18. Save gene-level expression matrix
############################################################

saveRDS(
  
  expr_gene,
  
  file =
    "results/GSE44765_rescue/GSE44765_gene_expression_RMA.rds"
)

con <- gzfile(
  
  "results/GSE44765_rescue/GSE44765_gene_expression_RMA.csv.gz",
  
  open = "wt"
)

write.csv(
  expr_gene,
  con,
  quote = FALSE
)

close(con)


############################################################
# 19. Design matrix
#
# Donor is explicitly included in the model.
############################################################

design <- model.matrix(
  
  ~ 0 + group + donor,
  
  data = meta
)

cat("Design matrix:\n")
print(design)
cat("\n")


############################################################
# 20. Check design rank
############################################################

if (
  qr(design)$rank != ncol(design)
) {
  
  stop(
    "\nERROR: Design matrix is not full rank.\n"
  )
}


############################################################
# 21. Fit limma model
############################################################

fit <- lmFit(
  expr_gene,
  design
)


############################################################
# 22. Contrasts
############################################################

contrast_matrix <- makeContrasts(
  
  Intact_vs_P3_2D =
    groupIntact - groupP3_2D,
  
  P3_3D_vs_P3_2D =
    groupP3_3D - groupP3_2D,
  
  P3_3D_vs_Intact =
    groupP3_3D - groupIntact,
  
  levels = design
)


############################################################
# 23. Contrast fitting
############################################################

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- eBayes(
  fit2,
  trend = TRUE,
  robust = TRUE
)


############################################################
# 24. Extract differential-expression tables
############################################################

tt_loss <- topTable(
  
  fit2,
  
  coef = "Intact_vs_P3_2D",
  
  number = Inf,
  
  adjust.method = "BH",
  
  sort.by = "none"
)

tt_rescue <- topTable(
  
  fit2,
  
  coef = "P3_3D_vs_P3_2D",
  
  number = Inf,
  
  adjust.method = "BH",
  
  sort.by = "none"
)

tt_final <- topTable(
  
  fit2,
  
  coef = "P3_3D_vs_Intact",
  
  number = Inf,
  
  adjust.method = "BH",
  
  sort.by = "none"
)


############################################################
# 25. Add gene symbols
############################################################

tt_loss$GeneSymbol <- rownames(tt_loss)
tt_rescue$GeneSymbol <- rownames(tt_rescue)
tt_final$GeneSymbol <- rownames(tt_final)


############################################################
# 26. Reorder columns
############################################################

move_gene_first <- function(x) {
  
  x[
    ,
    c(
      "GeneSymbol",
      setdiff(
        colnames(x),
        "GeneSymbol"
      )
    ),
    drop = FALSE
  ]
}

tt_loss <- move_gene_first(tt_loss)
tt_rescue <- move_gene_first(tt_rescue)
tt_final <- move_gene_first(tt_final)


############################################################
# 27. DEG classification function
############################################################

classify_de <- function(
    table,
    fdr_cutoff = FDR_CUTOFF,
    logfc_cutoff = LOGFC_CUTOFF
) {
  
  status <- rep(
    "Not_significant",
    nrow(table)
  )
  
  status[
    table$adj.P.Val < fdr_cutoff &
      table$logFC >= logfc_cutoff
  ] <- "Up"
  
  status[
    table$adj.P.Val < fdr_cutoff &
      table$logFC <= -logfc_cutoff
  ] <- "Down"
  
  status
}


############################################################
# 28. Add DEG status
############################################################

tt_loss$Status <-
  classify_de(tt_loss)

tt_rescue$Status <-
  classify_de(tt_rescue)

tt_final$Status <-
  classify_de(tt_final)


############################################################
# 29. Save DE tables
############################################################

write.csv(
  
  tt_loss,
  
  "results/GSE44765_rescue/DE_Intact_vs_P3_2D_all_genes.csv",
  
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  
  tt_rescue,
  
  "results/GSE44765_rescue/DE_P3_3D_vs_P3_2D_all_genes.csv",
  
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  
  tt_final,
  
  "results/GSE44765_rescue/DE_P3_3D_vs_Intact_all_genes.csv",
  
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 30. Print DEG counts
############################################################

cat("\n")
cat("========================================================\n")
cat(" Differential expression summary\n")
cat("========================================================\n\n")

cat("Intact vs P3_2D:\n")
print(table(tt_loss$Status))

cat("\nP3_3D vs P3_2D:\n")
print(table(tt_rescue$Status))

cat("\nP3_3D vs Intact:\n")
print(table(tt_final$Status))

cat("\n")


############################################################
# 31. Prepare merge tables
############################################################

loss_small <- tt_loss[
  ,
  c(
    "GeneSymbol",
    "logFC",
    "AveExpr",
    "P.Value",
    "adj.P.Val"
  )
]

colnames(loss_small) <- c(
  "GeneSymbol",
  "Loss_logFC",
  "Loss_AveExpr",
  "Loss_P",
  "Loss_FDR"
)


rescue_small <- tt_rescue[
  ,
  c(
    "GeneSymbol",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )
]

colnames(rescue_small) <- c(
  "GeneSymbol",
  "Rescue_logFC",
  "Rescue_P",
  "Rescue_FDR"
)


final_small <- tt_final[
  ,
  c(
    "GeneSymbol",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )
]

colnames(final_small) <- c(
  "GeneSymbol",
  "Final_vs_Intact_logFC",
  "Final_vs_Intact_P",
  "Final_vs_Intact_FDR"
)


############################################################
# 32. Merge contrast results
############################################################

rescue_table <- merge(
  loss_small,
  rescue_small,
  by = "GeneSymbol",
  all = TRUE
)

rescue_table <- merge(
  rescue_table,
  final_small,
  by = "GeneSymbol",
  all = TRUE
)


############################################################
# 33. Define significant change during 2D culture
############################################################

rescue_table$Altered_in_2D <-
  
  rescue_table$Loss_FDR < FDR_CUTOFF &
  
  abs(
    rescue_table$Loss_logFC
  ) >= LOGFC_CUTOFF


############################################################
# 34. Define significant response to 3D
############################################################

rescue_table$Changed_in_3D <-
  
  rescue_table$Rescue_FDR < FDR_CUTOFF &
  
  abs(
    rescue_table$Rescue_logFC
  ) >= LOGFC_CUTOFF


############################################################
# 35. Direction consistency
############################################################

rescue_table$Same_direction <-
  
  sign(
    rescue_table$Loss_logFC
  ) ==
  
  sign(
    rescue_table$Rescue_logFC
  )


############################################################
# 36. Is 3D closer to Intact?
############################################################

rescue_table$Closer_to_Intact <-
  
  abs(
    rescue_table$Final_vs_Intact_logFC
  ) <
  
  abs(
    rescue_table$Loss_logFC
  )


############################################################
# 37. Rescue fraction
############################################################

rescue_table$Rescue_fraction <-
  
  rescue_table$Rescue_logFC /
  rescue_table$Loss_logFC


############################################################
# 38. Directional rescue
############################################################

rescue_table$Directional_rescue <-
  
  rescue_table$Altered_in_2D &
  
  rescue_table$Same_direction &
  
  rescue_table$Closer_to_Intact


############################################################
# 39. Strict rescue
############################################################

rescue_table$Strict_rescue <-
  
  rescue_table$Altered_in_2D &
  
  rescue_table$Changed_in_3D &
  
  rescue_table$Same_direction &
  
  rescue_table$Closer_to_Intact


############################################################
# 40. Rescue type
############################################################

rescue_table$Rescue_type <-
  "No_rescue"

rescue_table$Rescue_type[
  rescue_table$Directional_rescue &
    rescue_table$Loss_logFC > 0
] <- "Restored_upward"

rescue_table$Rescue_type[
  rescue_table$Directional_rescue &
    rescue_table$Loss_logFC < 0
] <- "Restored_downward"


############################################################
# 41. Combined FDR
############################################################

rescue_table$Combined_FDR <- pmax(
  
  rescue_table$Loss_FDR,
  
  rescue_table$Rescue_FDR,
  
  na.rm = FALSE
)


############################################################
# 42. Sort table
############################################################

rescue_table <- rescue_table[
  order(
    rescue_table$Combined_FDR,
    -abs(
      rescue_table$Rescue_logFC
    ),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


############################################################
# 43. Strict and directional tables
############################################################

strict_rescue <- rescue_table[
  !is.na(rescue_table$Strict_rescue) &
    rescue_table$Strict_rescue,
  ,
  drop = FALSE
]

directional_rescue <- rescue_table[
  !is.na(rescue_table$Directional_rescue) &
    rescue_table$Directional_rescue,
  ,
  drop = FALSE
]


############################################################
# 44. Save rescue tables
############################################################

write.csv(
  
  rescue_table,
  
  "results/GSE44765_rescue/GSE44765_rescue_analysis_all_genes.csv",
  
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  
  strict_rescue,
  
  "results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv",
  
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  
  directional_rescue,
  
  "results/GSE44765_rescue/GSE44765_directional_rescue_genes.csv",
  
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 45. Rescue summary
############################################################

cat("\n")
cat("========================================================\n")
cat(" Rescue summary\n")
cat("========================================================\n\n")

cat(
  "Genes altered in P3 2D vs Intact:",
  sum(
    rescue_table$Altered_in_2D,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Genes significantly changed by 3D:",
  sum(
    rescue_table$Changed_in_3D,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Directional rescue genes:",
  nrow(directional_rescue),
  "\n"
)

cat(
  "Strict rescue genes:",
  nrow(strict_rescue),
  "\n\n"
)


############################################################
# 46. Volcano plotting function
############################################################

make_volcano <- function(
    tt,
    title_text,
    output_prefix
) {
  
  temp <- tt
  
  temp$minus_log10_FDR <-
    
    -log10(
      pmax(
        temp$adj.P.Val,
        1e-300
      )
    )
  
  temp$Plot_status <-
    classify_de(temp)
  
  sig <- temp[
    temp$Plot_status != "Not_significant",
    ,
    drop = FALSE
  ]
  
  sig <- sig[
    order(
      sig$adj.P.Val
    ),
    ,
    drop = FALSE
  ]
  
  label_n <- min(
    TOP_LABEL_GENES,
    nrow(sig)
  )
  
  if (label_n > 0) {
    
    label_data <- sig[
      seq_len(label_n),
      ,
      drop = FALSE
    ]
    
  } else {
    
    label_data <- sig[
      0,
      ,
      drop = FALSE
    ]
  }
  
  
  p <- ggplot(
    
    temp,
    
    aes(
      x = logFC,
      y = minus_log10_FDR,
      color = Plot_status
    )
    
  ) +
    
    geom_point(
      alpha = 0.65,
      size = 1.4
    ) +
    
    geom_vline(
      xintercept = c(
        -LOGFC_CUTOFF,
        LOGFC_CUTOFF
      ),
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept =
        -log10(FDR_CUTOFF),
      linetype = "dashed"
    ) +
    
    ggrepel::geom_text_repel(
      
      data = label_data,
      
      aes(
        x = logFC,
        y = minus_log10_FDR,
        label = GeneSymbol
      ),
      
      inherit.aes = FALSE,
      
      size = 3,
      
      max.overlaps = Inf
    ) +
    
    labs(
      
      title = title_text,
      
      subtitle = paste0(
        "FDR < ",
        FDR_CUTOFF,
        ", |log2FC| >= ",
        LOGFC_CUTOFF
      ),
      
      x = "log2 fold change",
      
      y = "-log10(FDR)",
      
      color = "Status"
    ) +
    
    theme_classic(
      base_size = 13
    )
  
  
  ggsave(
    
    paste0(
      "figures/GSE44765_rescue/",
      output_prefix,
      ".pdf"
    ),
    
    p,
    
    width = 8,
    height = 6
  )
  
  ggsave(
    
    paste0(
      "figures/GSE44765_rescue/",
      output_prefix,
      ".png"
    ),
    
    p,
    
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(p)
}


############################################################
# 47. Volcano plots
############################################################

volcano_loss <- make_volcano(
  
  tt_loss,
  
  "Intact DP vs P3 2D",
  
  "Volcano_Intact_vs_P3_2D"
)

volcano_rescue <- make_volcano(
  
  tt_rescue,
  
  "P3 3D spheroid vs P3 2D",
  
  "Volcano_P3_3D_vs_P3_2D"
)


############################################################
# 48. Rescue landscape dataframe
############################################################

rescue_plot_data <- rescue_table[
  
  is.finite(
    rescue_table$Loss_logFC
  ) &
    
    is.finite(
      rescue_table$Rescue_logFC
    ),
  
  ,
  drop = FALSE
]

rescue_plot_data$Plot_class <-
  "Other"

rescue_plot_data$Plot_class[
  !is.na(
    rescue_plot_data$Directional_rescue
  ) &
    rescue_plot_data$Directional_rescue
] <- "Directional rescue"

rescue_plot_data$Plot_class[
  !is.na(
    rescue_plot_data$Strict_rescue
  ) &
    rescue_plot_data$Strict_rescue
] <- "Strict rescue"


############################################################
# 49. Select rescue genes to label
############################################################

if (nrow(strict_rescue) > 0) {
  
  label_rescue <- strict_rescue
  
} else {
  
  label_rescue <- directional_rescue
}

label_rescue <- label_rescue[
  is.finite(label_rescue$Loss_logFC) &
    is.finite(label_rescue$Rescue_logFC),
  ,
  drop = FALSE
]

label_rescue <- label_rescue[
  order(
    label_rescue$Combined_FDR,
    -abs(
      label_rescue$Rescue_logFC
    ),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

label_n <- min(
  TOP_LABEL_GENES,
  nrow(label_rescue)
)

if (label_n > 0) {
  
  label_rescue <- label_rescue[
    seq_len(label_n),
    ,
    drop = FALSE
  ]
  
} else {
  
  label_rescue <- rescue_plot_data[
    0,
    ,
    drop = FALSE
  ]
}


############################################################
# 50. Rescue landscape plot
#
# IMPORTANT FIX:
# The text layer uses inherit.aes = FALSE.
# Therefore it does NOT require Plot_class.
############################################################

p_rescue_scatter <- ggplot(
  
  rescue_plot_data,
  
  aes(
    x = Loss_logFC,
    y = Rescue_logFC,
    color = Plot_class
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  
  geom_point(
    alpha = 0.55,
    size = 1.6
  ) +
  
  ggrepel::geom_text_repel(
    
    data = label_rescue,
    
    aes(
      x = Loss_logFC,
      y = Rescue_logFC,
      label = GeneSymbol
    ),
    
    inherit.aes = FALSE,
    
    size = 3,
    
    max.overlaps = Inf
  ) +
  
  labs(
    
    title =
      "Transcriptional rescue landscape",
    
    subtitle =
      "Dashed diagonal indicates approximate complete restoration toward Intact DP",
    
    x =
      "Intact − P3 2D log2FC",
    
    y =
      "P3 3D − P3 2D log2FC",
    
    color =
      "Classification"
  ) +
  
  theme_classic(
    base_size = 13
  )


############################################################
# 51. Save rescue landscape
############################################################

ggsave(
  
  "figures/GSE44765_rescue/Rescue_landscape.pdf",
  
  p_rescue_scatter,
  
  width = 8,
  height = 7
)

ggsave(
  
  "figures/GSE44765_rescue/Rescue_landscape.png",
  
  p_rescue_scatter,
  
  width = 8,
  height = 7,
  dpi = 300
)


############################################################
# 52. Select heatmap candidates
############################################################

if (nrow(strict_rescue) > 0) {
  
  heatmap_candidates <- strict_rescue
  
  cat(
    "\nUsing STRICT rescue genes for heatmap.\n"
  )
  
} else if (nrow(directional_rescue) > 0) {
  
  heatmap_candidates <- directional_rescue
  
  cat(
    "\nUsing directional rescue genes for heatmap.\n"
  )
  
} else {
  
  heatmap_candidates <- rescue_table
  
  warning(
    "No rescue genes passed thresholds. Using strongest signals."
  )
}


############################################################
# 53. Sort heatmap genes
############################################################

heatmap_candidates <- heatmap_candidates[
  order(
    heatmap_candidates$Combined_FDR,
    -abs(
      heatmap_candidates$Rescue_logFC
    ),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

heatmap_n <- min(
  TOP_HEATMAP_GENES,
  nrow(heatmap_candidates)
)

heatmap_genes <-
  heatmap_candidates$GeneSymbol[
    seq_len(heatmap_n)
  ]

heatmap_genes <- unique(
  heatmap_genes[
    heatmap_genes %in%
      rownames(expr_gene)
  ]
)


############################################################
# 54. Rescue heatmap
############################################################

if (length(heatmap_genes) >= 2) {
  
  heatmap_matrix <- expr_gene[
    heatmap_genes,
    ,
    drop = FALSE
  ]
  
  heatmap_z <- t(
    scale(
      t(
        heatmap_matrix
      )
    )
  )
  
  heatmap_z[
    !is.finite(heatmap_z)
  ] <- 0
  
  
  heatmap_annotation <- data.frame(
    
    Group =
      as.character(
        meta$group
      ),
    
    Donor =
      as.character(
        meta$donor
      ),
    
    stringsAsFactors = FALSE
  )
  
  rownames(
    heatmap_annotation
  ) <- rownames(meta)
  
  
  pheatmap(
    
    heatmap_z,
    
    annotation_col =
      heatmap_annotation,
    
    cluster_rows = TRUE,
    
    cluster_cols = TRUE,
    
    show_rownames = TRUE,
    
    show_colnames = TRUE,
    
    fontsize_row = 7,
    
    fontsize_col = 8,
    
    border_color = NA,
    
    main =
      "Top transcriptional rescue genes",
    
    filename =
      "figures/GSE44765_rescue/Top_rescue_genes_heatmap.pdf",
    
    width = 11,
    
    height = 10
  )
  
  
  pheatmap(
    
    heatmap_z,
    
    annotation_col =
      heatmap_annotation,
    
    cluster_rows = TRUE,
    
    cluster_cols = TRUE,
    
    show_rownames = TRUE,
    
    show_colnames = TRUE,
    
    fontsize_row = 7,
    
    fontsize_col = 8,
    
    border_color = NA,
    
    main =
      "Top transcriptional rescue genes",
    
    filename =
      "figures/GSE44765_rescue/Top_rescue_genes_heatmap.png",
    
    width = 11,
    
    height = 10
  )
}


############################################################
# 55. Select trajectory genes
############################################################

if (nrow(strict_rescue) > 0) {
  
  trajectory_candidates <-
    strict_rescue
  
} else {
  
  trajectory_candidates <-
    directional_rescue
}


############################################################
# 56. Trajectory plot
############################################################

if (
  nrow(
    trajectory_candidates
  ) > 0
) {
  
  trajectory_candidates <-
    trajectory_candidates[
      order(
        trajectory_candidates$Combined_FDR,
        -abs(
          trajectory_candidates$Rescue_logFC
        ),
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
  
  trajectory_n <- min(
    TOP_TRAJECTORY_GENES,
    nrow(
      trajectory_candidates
    )
  )
  
  trajectory_genes <-
    trajectory_candidates$GeneSymbol[
      seq_len(
        trajectory_n
      )
    ]
  
  trajectory_genes <- trajectory_genes[
    trajectory_genes %in%
      rownames(expr_gene)
  ]
  
  
  trajectory_samples <- rownames(meta)[
    meta$group %in%
      c(
        "Intact",
        "P3_2D",
        "P3_3D"
      )
  ]
  
  trajectory_meta <- meta[
    trajectory_samples,
    ,
    drop = FALSE
  ]
  
  
  trajectory_meta$group <- factor(
    
    as.character(
      trajectory_meta$group
    ),
    
    levels = c(
      "Intact",
      "P3_2D",
      "P3_3D"
    )
  )
  
  
  trajectory_long_list <- list()
  
  counter <- 1
  
  
  for (gene in trajectory_genes) {
    
    temp_df <- data.frame(
      
      Gene = gene,
      
      sample =
        trajectory_samples,
      
      expression =
        as.numeric(
          expr_gene[
            gene,
            trajectory_samples
          ]
        ),
      
      group =
        as.character(
          trajectory_meta[
            trajectory_samples,
            "group"
          ]
        ),
      
      donor =
        as.character(
          trajectory_meta[
            trajectory_samples,
            "donor"
          ]
        ),
      
      stringsAsFactors = FALSE
    )
    
    trajectory_long_list[[counter]] <-
      temp_df
    
    counter <- counter + 1
  }
  
  
  trajectory_long <- do.call(
    rbind,
    trajectory_long_list
  )
  
  
  trajectory_long$group <- factor(
    
    trajectory_long$group,
    
    levels = c(
      "Intact",
      "P3_2D",
      "P3_3D"
    )
  )
  
  
  p_trajectory <- ggplot(
    
    trajectory_long,
    
    aes(
      x = group,
      y = expression,
      group = donor,
      shape = donor
    )
    
  ) +
    
    geom_line(
      alpha = 0.65
    ) +
    
    geom_point(
      size = 2.4
    ) +
    
    facet_wrap(
      
      ~ Gene,
      
      scales = "free_y",
      
      ncol = 4
    ) +
    
    labs(
      
      title =
        "Top transcriptional rescue trajectories",
      
      subtitle =
        "Each line represents the same donor",
      
      x = NULL,
      
      y =
        "RMA log2 expression",
      
      shape =
        "Donor"
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        ),
      
      strip.text =
        element_text(
          face = "bold"
        )
    )
  
  
  ggsave(
    
    "figures/GSE44765_rescue/Top_rescue_gene_trajectories.pdf",
    
    p_trajectory,
    
    width = 12,
    
    height = 9
  )
  
  
  ggsave(
    
    "figures/GSE44765_rescue/Top_rescue_gene_trajectories.png",
    
    p_trajectory,
    
    width = 12,
    
    height = 9,
    
    dpi = 300
  )
}


############################################################
# 57. Global movement toward Intact
############################################################

intact_samples <- rownames(meta)[
  meta$group == "Intact"
]

p3_2d_samples <- rownames(meta)[
  meta$group == "P3_2D"
]

p3_3d_samples <- rownames(meta)[
  meta$group == "P3_3D"
]


############################################################
# 58. Condition centroids
############################################################

centroid_intact <- rowMeans(
  
  expr_gene[
    ,
    intact_samples,
    drop = FALSE
  ]
)

centroid_2d <- rowMeans(
  
  expr_gene[
    ,
    p3_2d_samples,
    drop = FALSE
  ]
)

centroid_3d <- rowMeans(
  
  expr_gene[
    ,
    p3_3d_samples,
    drop = FALSE
  ]
)


############################################################
# 59. Distance from Intact
############################################################

distance_2d_to_intact <-
  
  abs(
    centroid_2d -
      centroid_intact
  )

distance_3d_to_intact <-
  
  abs(
    centroid_3d -
      centroid_intact
  )


############################################################
# 60. Global rescue statistics
############################################################

global_rescue_fraction <-
  
  mean(
    distance_3d_to_intact <
      distance_2d_to_intact
  )

median_distance_2d <-
  
  median(
    distance_2d_to_intact
  )

median_distance_3d <-
  
  median(
    distance_3d_to_intact
  )


############################################################
# 61. Save global summary
############################################################

global_summary <- data.frame(
  
  Metric = c(
    
    "Fraction_genes_3D_closer_to_Intact",
    
    "Median_absolute_distance_P3_2D_to_Intact",
    
    "Median_absolute_distance_P3_3D_to_Intact"
  ),
  
  Value = c(
    
    global_rescue_fraction,
    
    median_distance_2d,
    
    median_distance_3d
  )
)

write.csv(
  
  global_summary,
  
  "results/GSE44765_rescue/GSE44765_global_rescue_summary.csv",
  
  row.names = FALSE,
  
  quote = FALSE
)


############################################################
# 62. Save design
############################################################

write.csv(
  
  design,
  
  "results/GSE44765_rescue/GSE44765_design_matrix.csv",
  
  quote = FALSE
)

write.csv(
  
  contrast_matrix,
  
  "results/GSE44765_rescue/GSE44765_contrast_matrix.csv",
  
  quote = FALSE
)


############################################################
# 63. Save workspace
############################################################

save(
  
  expr_gene,
  meta,
  design,
  contrast_matrix,
  fit2,
  tt_loss,
  tt_rescue,
  tt_final,
  rescue_table,
  strict_rescue,
  directional_rescue,
  global_summary,
  
  file =
    "results/GSE44765_rescue/GSE44765_rescue_workspace.RData"
)


############################################################
# 64. Session information
############################################################

sink(
  
  "results/GSE44765_rescue/sessionInfo_rescue_analysis.txt"
)

sessionInfo()

sink()


############################################################
# 65. Human-readable summary
############################################################

summary_file <-
  "results/GSE44765_rescue/ANALYSIS_SUMMARY.txt"

sink(summary_file)

cat(
  "GSE44765 transcriptional rescue analysis\n"
)

cat(
  "=========================================\n\n"
)

cat(
  "FDR cutoff:",
  FDR_CUTOFF,
  "\n"
)

cat(
  "Absolute log2FC cutoff:",
  LOGFC_CUTOFF,
  "\n\n"
)

cat(
  "Unique genes analyzed:",
  nrow(expr_gene),
  "\n\n"
)

cat(
  "Genes altered between Intact and P3 2D:",
  sum(
    rescue_table$Altered_in_2D,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Genes significantly changed between P3 3D and P3 2D:",
  sum(
    rescue_table$Changed_in_3D,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Directional rescue genes:",
  nrow(directional_rescue),
  "\n"
)

cat(
  "Strict rescue genes:",
  nrow(strict_rescue),
  "\n\n"
)

cat(
  "Fraction of all genes where P3 3D is closer to Intact:",
  round(
    global_rescue_fraction,
    4
  ),
  "\n"
)

cat(
  "Median distance P3 2D -> Intact:",
  round(
    median_distance_2d,
    4
  ),
  "\n"
)

cat(
  "Median distance P3 3D -> Intact:",
  round(
    median_distance_3d,
    4
  ),
  "\n"
)

sink()


############################################################
# 66. Final console report
############################################################

cat("\n")
cat("========================================================\n")
cat(" ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("========================================================\n\n")

cat(
  "Genes analyzed:",
  nrow(expr_gene),
  "\n"
)

cat(
  "Altered in P3 2D vs Intact:",
  sum(
    rescue_table$Altered_in_2D,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significantly changed by P3 3D:",
  sum(
    rescue_table$Changed_in_3D,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Directional rescue genes:",
  nrow(directional_rescue),
  "\n"
)

cat(
  "STRICT rescue genes:",
  nrow(strict_rescue),
  "\n\n"
)

cat(
  "Fraction of all genes closer to Intact after 3D:",
  round(
    100 *
      global_rescue_fraction,
    1
  ),
  "%\n"
)

cat(
  "Median distance P3 2D -> Intact:",
  round(
    median_distance_2d,
    3
  ),
  "\n"
)

cat(
  "Median distance P3 3D -> Intact:",
  round(
    median_distance_3d,
    3
  ),
  "\n\n"
)

cat("Main figures:\n")

cat(
  "figures/GSE44765_rescue/Volcano_Intact_vs_P3_2D.png\n"
)

cat(
  "figures/GSE44765_rescue/Volcano_P3_3D_vs_P3_2D.png\n"
)

cat(
  "figures/GSE44765_rescue/Rescue_landscape.png\n"
)

cat(
  "figures/GSE44765_rescue/Top_rescue_genes_heatmap.png\n"
)

cat(
  "figures/GSE44765_rescue/Top_rescue_gene_trajectories.png\n"
)

cat("\nMain tables:\n")

cat(
  "results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv\n"
)

cat(
  "results/GSE44765_rescue/GSE44765_directional_rescue_genes.csv\n"
)

cat(
  "results/GSE44765_rescue/GSE44765_rescue_analysis_all_genes.csv\n"
)

cat(
  "results/GSE44765_rescue/ANALYSIS_SUMMARY.txt\n"
)

cat("\n")
cat("========================================================\n")