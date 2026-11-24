############################################################
# GSE44765 pathway-level rescue analysis
# OFFLINE / NO-MSigdbr / ReactomePA-compatible VERSION
#
# Purpose:
#   Avoid msigdbr / Zenodo downloads.
#
# Required inputs from Script 02:
#   results/GSE44765_rescue/GSE44765_gene_expression_RMA.rds
#   results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv
#   results/GSE44765_rescue/DE_Intact_vs_P3_2D_all_genes.csv
#   results/GSE44765_rescue/DE_P3_3D_vs_P3_2D_all_genes.csv
#
# Analyses:
#   1. GO Biological Process ORA
#   2. Reactome ORA
#   3. GO BP GSEA
#   4. Reactome GSEA
#   5. Pathway-level rescue in GO and Reactome
#   6. GO pathway sample-level trajectory / heatmap
#
# This version DOES NOT use msigdbr and DOES NOT download
# gene sets from Zenodo during analysis.
############################################################


############################################################
# 0. Basic settings
############################################################

graphics.off()
options(stringsAsFactors = FALSE)
set.seed(12345)

cat("\n")
cat("========================================================\n")
cat(" GSE44765 pathway rescue analysis - OFFLINE VERSION\n")
cat("========================================================\n\n")
cat("Working directory:\n")
cat(getwd(), "\n\n")


############################################################
# 1. Fixed thresholds
############################################################

FDR_CUTOFF <- 0.05
TOP_ORA_TERMS <- 15
TOP_PATHWAY_LABELS <- 12
TOP_TRAJECTORY_PATHWAYS <- 12
MIN_GO_GENES_FOR_SCORE <- 5


############################################################
# 2. Output folders
############################################################

dir.create(
  "results/GSE44765_pathway",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/GSE44765_pathway",
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Install CRAN packages if missing
############################################################

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "stringr"
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
# 4. Install Bioconductor manager if missing
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


############################################################
# 5. Install Bioconductor packages if missing
############################################################

bioc_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "AnnotationDbi"
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
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(AnnotationDbi)
})

cat("Packages loaded successfully.\n\n")


############################################################
# 7. Input files
############################################################

expr_file <-
  "results/GSE44765_rescue/GSE44765_gene_expression_RMA.rds"

strict_file <-
  "results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv"

loss_file <-
  "results/GSE44765_rescue/DE_Intact_vs_P3_2D_all_genes.csv"

rescue_file <-
  "results/GSE44765_rescue/DE_P3_3D_vs_P3_2D_all_genes.csv"

required_files <- c(
  expr_file,
  strict_file,
  loss_file,
  rescue_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  cat("\nMissing files:\n")
  print(missing_files)
  stop(
    "\nPlease run Scripts 01 and 02 successfully first.\n"
  )
}


############################################################
# 8. Load data
############################################################

expr_gene <- readRDS(expr_file)

strict_rescue <- read.csv(
  strict_file,
  stringsAsFactors = FALSE
)

tt_loss <- read.csv(
  loss_file,
  stringsAsFactors = FALSE
)

tt_rescue <- read.csv(
  rescue_file,
  stringsAsFactors = FALSE
)

cat(
  "Expression matrix:",
  nrow(expr_gene),
  "genes x",
  ncol(expr_gene),
  "samples\n"
)

cat(
  "Strict rescue genes:",
  nrow(strict_rescue),
  "\n\n"
)


############################################################
# 9. Check required columns
############################################################

if (!all(
  c("GeneSymbol", "Rescue_type") %in%
  colnames(strict_rescue)
)) {
  stop(
    "ERROR: strict rescue table lacks GeneSymbol/Rescue_type."
  )
}

if (!all(
  c("GeneSymbol", "t") %in%
  colnames(tt_loss)
)) {
  stop(
    "ERROR: Intact vs P3_2D table lacks GeneSymbol/t."
  )
}

if (!all(
  c("GeneSymbol", "t") %in%
  colnames(tt_rescue)
)) {
  stop(
    "ERROR: P3_3D vs P3_2D table lacks GeneSymbol/t."
  )
}


############################################################
# 10. Sample metadata
############################################################

meta <- data.frame(
  sample = paste0("DP", 1:18),
  donor = rep(c("D5", "D6", "D7"), times = 6),
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
  levels = c("D5", "D6", "D7")
)

if (!all(
  colnames(expr_gene) %in% rownames(meta)
)) {
  stop(
    "ERROR: expression samples do not match DP1-DP18."
  )
}

meta <- meta[
  colnames(expr_gene),
  ,
  drop = FALSE
]


############################################################
# 11. Gene universe and strict rescue directions
############################################################

gene_universe <- unique(rownames(expr_gene))
gene_universe <- gene_universe[
  !is.na(gene_universe) &
    gene_universe != ""
]

rescue_up <- unique(
  strict_rescue$GeneSymbol[
    strict_rescue$Rescue_type ==
      "Restored_upward"
  ]
)

rescue_down <- unique(
  strict_rescue$GeneSymbol[
    strict_rescue$Rescue_type ==
      "Restored_downward"
  ]
)

rescue_up <- rescue_up[
  !is.na(rescue_up) &
    rescue_up != ""
]

rescue_down <- rescue_down[
  !is.na(rescue_down) &
    rescue_down != ""
]

cat(
  "Restored upward genes:",
  length(rescue_up),
  "\n"
)

cat(
  "Restored downward genes:",
  length(rescue_down),
  "\n\n"
)


############################################################
# 12. Helper functions
############################################################

save_csv <- function(x, filename) {
  write.csv(
    x,
    file = filename,
    row.names = FALSE,
    quote = FALSE
  )
}

gene_ratio_numeric <- function(x) {
  sapply(
    strsplit(as.character(x), "/"),
    function(z) {
      if (length(z) != 2) {
        return(NA_real_)
      }
      as.numeric(z[1]) / as.numeric(z[2])
    }
  )
}

plot_enrichment <- function(
    enrichment_df,
    title_text,
    output_prefix,
    top_n = TOP_ORA_TERMS
) {
  
  if (
    is.null(enrichment_df) ||
    nrow(enrichment_df) == 0 ||
    !"p.adjust" %in% colnames(enrichment_df)
  ) {
    significant_df <- data.frame()
  } else {
    significant_df <- enrichment_df[
      !is.na(enrichment_df$p.adjust) &
        enrichment_df$p.adjust < FDR_CUTOFF,
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(significant_df) == 0) {
    
    p <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = "No pathways reached FDR < 0.05",
        size = 5
      ) +
      xlim(-1, 1) +
      ylim(-1, 1) +
      labs(title = title_text) +
      theme_void()
    
  } else {
    
    significant_df <- significant_df[
      order(significant_df$p.adjust),
      ,
      drop = FALSE
    ]
    
    plot_n <- min(
      top_n,
      nrow(significant_df)
    )
    
    plot_df <- significant_df[
      seq_len(plot_n),
      ,
      drop = FALSE
    ]
    
    plot_df$GeneRatioNumeric <-
      gene_ratio_numeric(
        plot_df$GeneRatio
      )
    
    plot_df$minus_log10_FDR <-
      -log10(
        pmax(
          plot_df$p.adjust,
          1e-300
        )
      )
    
    plot_df$Description_plot <-
      stringr::str_wrap(
        plot_df$Description,
        width = 42
      )
    
    plot_df$Description_plot <- factor(
      plot_df$Description_plot,
      levels = rev(
        plot_df$Description_plot
      )
    )
    
    p <- ggplot(
      plot_df,
      aes(
        x = GeneRatioNumeric,
        y = Description_plot,
        size = Count,
        color = minus_log10_FDR
      )
    ) +
      geom_point() +
      labs(
        title = title_text,
        x = "Gene ratio",
        y = NULL,
        size = "Gene count",
        color = "-log10(FDR)"
      ) +
      theme_classic(base_size = 12)
  }
  
  ggsave(
    filename = paste0(
      "figures/GSE44765_pathway/",
      output_prefix,
      ".pdf"
    ),
    plot = p,
    width = 9,
    height = 7
  )
  
  ggsave(
    filename = paste0(
      "figures/GSE44765_pathway/",
      output_prefix,
      ".png"
    ),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  significant_df
}

symbol_to_entrez_table <- function(symbols) {
  
  symbols <- unique(
    symbols[
      !is.na(symbols) &
        symbols != ""
    ]
  )
  
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c(
      "SYMBOL",
      "ENTREZID"
    )
  )
  
  mapping <- mapping[
    !is.na(mapping$ENTREZID) &
      mapping$ENTREZID != "",
    ,
    drop = FALSE
  ]
  
  unique(mapping)
}

prepare_entrez_rank <- function(de_table) {
  
  temp <- de_table[
    !is.na(de_table$GeneSymbol) &
      de_table$GeneSymbol != "" &
      !is.na(de_table$t) &
      is.finite(de_table$t),
    c(
      "GeneSymbol",
      "t"
    ),
    drop = FALSE
  ]
  
  mapping <- symbol_to_entrez_table(
    temp$GeneSymbol
  )
  
  temp <- merge(
    temp,
    mapping,
    by.x = "GeneSymbol",
    by.y = "SYMBOL"
  )
  
  temp$abs_t <- abs(temp$t)
  
  temp <- temp[
    order(
      temp$ENTREZID,
      -temp$abs_t
    ),
    ,
    drop = FALSE
  ]
  
  temp <- temp[
    !duplicated(temp$ENTREZID),
    ,
    drop = FALSE
  ]
  
  rank_vector <- temp$t
  
  names(rank_vector) <-
    temp$ENTREZID
  
  rank_vector <- sort(
    rank_vector,
    decreasing = TRUE
  )
  
  rank_vector
}

flatten_gsea_result <- function(df) {
  
  out <- df
  
  list_cols <- vapply(
    out,
    is.list,
    logical(1)
  )
  
  if (any(list_cols)) {
    for (nm in names(out)[list_cols]) {
      out[[nm]] <- vapply(
        out[[nm]],
        function(z) {
          paste(
            z,
            collapse = ";"
          )
        },
        character(1)
      )
    }
  }
  
  out
}


############################################################
# 13. GO Biological Process ORA
############################################################

cat(
  "Running GO Biological Process ORA...\n"
)

go_up <- clusterProfiler::enrichGO(
  gene = rescue_up,
  universe = gene_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

go_down <- clusterProfiler::enrichGO(
  gene = rescue_down,
  universe = gene_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

go_up_df <- as.data.frame(go_up)
go_down_df <- as.data.frame(go_down)

save_csv(
  go_up_df,
  "results/GSE44765_pathway/GO_BP_restored_upward.csv"
)

save_csv(
  go_down_df,
  "results/GSE44765_pathway/GO_BP_restored_downward.csv"
)

go_up_sig <- plot_enrichment(
  go_up_df,
  "GO Biological Process: restored upward",
  "GO_BP_restored_upward"
)

go_down_sig <- plot_enrichment(
  go_down_df,
  "GO Biological Process: restored downward",
  "GO_BP_restored_downward"
)


############################################################
# 14. Reactome ORA
############################################################

cat(
  "Running Reactome ORA...\n"
)

rescue_up_map <-
  symbol_to_entrez_table(
    rescue_up
  )

rescue_down_map <-
  symbol_to_entrez_table(
    rescue_down
  )

universe_map <-
  symbol_to_entrez_table(
    gene_universe
  )

rescue_up_entrez <-
  unique(
    rescue_up_map$ENTREZID
  )

rescue_down_entrez <-
  unique(
    rescue_down_map$ENTREZID
  )

universe_entrez <-
  unique(
    universe_map$ENTREZID
  )

reactome_up <- ReactomePA::enrichPathway(
  gene = rescue_up_entrez,
  organism = "human",
  universe = universe_entrez,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  readable = TRUE
)

reactome_down <- ReactomePA::enrichPathway(
  gene = rescue_down_entrez,
  organism = "human",
  universe = universe_entrez,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  readable = TRUE
)

reactome_up_df <-
  as.data.frame(
    reactome_up
  )

reactome_down_df <-
  as.data.frame(
    reactome_down
  )

save_csv(
  reactome_up_df,
  "results/GSE44765_pathway/Reactome_restored_upward.csv"
)

save_csv(
  reactome_down_df,
  "results/GSE44765_pathway/Reactome_restored_downward.csv"
)

reactome_up_sig <- plot_enrichment(
  reactome_up_df,
  "Reactome: restored upward",
  "Reactome_restored_upward"
)

reactome_down_sig <- plot_enrichment(
  reactome_down_df,
  "Reactome: restored downward",
  "Reactome_restored_downward"
)


############################################################
# 15. Prepare ranked Entrez gene lists for GSEA
############################################################

cat(
  "Preparing ranked Entrez gene lists...\n"
)

rank_loss_entrez <-
  prepare_entrez_rank(
    tt_loss
  )

rank_rescue_entrez <-
  prepare_entrez_rank(
    tt_rescue
  )

cat(
  "Ranked genes for Intact vs P3_2D:",
  length(rank_loss_entrez),
  "\n"
)

cat(
  "Ranked genes for P3_3D vs P3_2D:",
  length(rank_rescue_entrez),
  "\n\n"
)



############################################################
# 16A. Version-compatible Reactome GSEA helper
#
# Different ReactomePA versions expose slightly different
# arguments. We therefore inspect the installed function
# and pass only arguments supported by that version.
############################################################

run_reactome_gsea <- function(gene_list) {
  
  available_args <- names(
    formals(
      ReactomePA::gsePathway
    )
  )
  
  args <- list(
    geneList = gene_list,
    organism = "human",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  # Some newer versions support this argument; older ones do not.
  if ("by" %in% available_args) {
    args$by <- "fgsea"
  }
  
  # Only add eps when supported by the installed version.
  if ("eps" %in% available_args) {
    args$eps <- 0
  }
  
  do.call(
    ReactomePA::gsePathway,
    args
  )
}


############################################################
# 16. GO BP GSEA
############################################################

cat(
  "Running GO BP GSEA: Intact vs P3_2D...\n"
)

gsea_go_loss <- clusterProfiler::gseGO(
  geneList = rank_loss_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  by = "fgsea"
)

cat(
  "Running GO BP GSEA: P3_3D vs P3_2D...\n"
)

gsea_go_rescue <- clusterProfiler::gseGO(
  geneList = rank_rescue_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  by = "fgsea"
)

go_loss_gsea_df <-
  as.data.frame(
    gsea_go_loss
  )

go_rescue_gsea_df <-
  as.data.frame(
    gsea_go_rescue
  )

save_csv(
  flatten_gsea_result(
    go_loss_gsea_df
  ),
  "results/GSE44765_pathway/GSEA_GO_BP_Intact_vs_P3_2D.csv"
)

save_csv(
  flatten_gsea_result(
    go_rescue_gsea_df
  ),
  "results/GSE44765_pathway/GSEA_GO_BP_P3_3D_vs_P3_2D.csv"
)


############################################################
# 17. Reactome GSEA
############################################################

cat(
  "Running Reactome GSEA: Intact vs P3_2D...
"
)

gsea_reactome_loss <- run_reactome_gsea(
  rank_loss_entrez
)

cat(
  "Running Reactome GSEA: P3_3D vs P3_2D...
"
)

gsea_reactome_rescue <- run_reactome_gsea(
  rank_rescue_entrez
)

reactome_loss_gsea_df <-
  as.data.frame(
    gsea_reactome_loss
  )

reactome_rescue_gsea_df <-
  as.data.frame(
    gsea_reactome_rescue
  )

save_csv(
  flatten_gsea_result(
    reactome_loss_gsea_df
  ),
  "results/GSE44765_pathway/GSEA_Reactome_Intact_vs_P3_2D.csv"
)

save_csv(
  flatten_gsea_result(
    reactome_rescue_gsea_df
  ),
  "results/GSE44765_pathway/GSEA_Reactome_P3_3D_vs_P3_2D.csv"
)


############################################################
# 18. Build pathway rescue table
############################################################

build_pathway_rescue_table <- function(
    loss_df,
    rescue_df
) {
  
  required_columns <- c(
    "ID",
    "Description",
    "NES",
    "pvalue",
    "p.adjust"
  )
  
  if (!all(
    required_columns %in%
    colnames(loss_df)
  )) {
    stop(
      "ERROR: loss GSEA table has unexpected columns."
    )
  }
  
  if (!all(
    required_columns %in%
    colnames(rescue_df)
  )) {
    stop(
      "ERROR: rescue GSEA table has unexpected columns."
    )
  }
  
  x1 <- loss_df[
    ,
    required_columns,
    drop = FALSE
  ]
  
  colnames(x1) <- c(
    "PathwayID",
    "Description",
    "Loss_NES",
    "Loss_P",
    "Loss_FDR"
  )
  
  x2 <- rescue_df[
    ,
    required_columns,
    drop = FALSE
  ]
  
  colnames(x2) <- c(
    "PathwayID",
    "Description_rescue",
    "Rescue_NES",
    "Rescue_P",
    "Rescue_FDR"
  )
  
  merged <- merge(
    x1,
    x2,
    by = "PathwayID",
    all = TRUE
  )
  
  merged$Description_final <- ifelse(
    !is.na(merged$Description),
    merged$Description,
    merged$Description_rescue
  )
  
  merged$Same_direction <-
    sign(
      merged$Loss_NES
    ) ==
    sign(
      merged$Rescue_NES
    )
  
  merged$Strict_pathway_rescue <-
    !is.na(
      merged$Loss_FDR
    ) &
    !is.na(
      merged$Rescue_FDR
    ) &
    merged$Loss_FDR <
    FDR_CUTOFF &
    merged$Rescue_FDR <
    FDR_CUTOFF &
    merged$Same_direction
  
  merged$Rescue_type <-
    "No_strict_rescue"
  
  merged$Rescue_type[
    !is.na(
      merged$Strict_pathway_rescue
    ) &
      merged$Strict_pathway_rescue &
      merged$Loss_NES > 0
  ] <- "Restored_upward"
  
  merged$Rescue_type[
    !is.na(
      merged$Strict_pathway_rescue
    ) &
      merged$Strict_pathway_rescue &
      merged$Loss_NES < 0
  ] <- "Restored_downward"
  
  merged$Combined_FDR <- pmax(
    merged$Loss_FDR,
    merged$Rescue_FDR,
    na.rm = FALSE
  )
  
  merged$Rescue_strength <- pmin(
    abs(
      merged$Loss_NES
    ),
    abs(
      merged$Rescue_NES
    )
  )
  
  merged <- merged[
    order(
      merged$Combined_FDR,
      -merged$Rescue_strength,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  merged
}


############################################################
# 19. GO / Reactome rescue tables
############################################################

go_pathway_rescue <-
  build_pathway_rescue_table(
    go_loss_gsea_df,
    go_rescue_gsea_df
  )

reactome_pathway_rescue <-
  build_pathway_rescue_table(
    reactome_loss_gsea_df,
    reactome_rescue_gsea_df
  )

save_csv(
  go_pathway_rescue,
  "results/GSE44765_pathway/GO_BP_pathway_rescue_table.csv"
)

save_csv(
  reactome_pathway_rescue,
  "results/GSE44765_pathway/Reactome_pathway_rescue_table.csv"
)

go_strict <- go_pathway_rescue[
  !is.na(
    go_pathway_rescue$Strict_pathway_rescue
  ) &
    go_pathway_rescue$Strict_pathway_rescue,
  ,
  drop = FALSE
]

reactome_strict <- reactome_pathway_rescue[
  !is.na(
    reactome_pathway_rescue$Strict_pathway_rescue
  ) &
    reactome_pathway_rescue$Strict_pathway_rescue,
  ,
  drop = FALSE
]

save_csv(
  go_strict,
  "results/GSE44765_pathway/GO_BP_STRICT_rescued_pathways.csv"
)

save_csv(
  reactome_strict,
  "results/GSE44765_pathway/Reactome_STRICT_rescued_pathways.csv"
)


############################################################
# 20. Pathway rescue scatter helper
############################################################

plot_pathway_rescue <- function(
    pathway_table,
    title_text,
    output_prefix
) {
  
  plot_df <- pathway_table[
    is.finite(
      pathway_table$Loss_NES
    ) &
      is.finite(
        pathway_table$Rescue_NES
      ),
    ,
    drop = FALSE
  ]
  
  plot_df$Plot_class <-
    "Other"
  
  plot_df$Plot_class[
    !is.na(
      plot_df$Strict_pathway_rescue
    ) &
      plot_df$Strict_pathway_rescue
  ] <- "Strict pathway rescue"
  
  label_df <- plot_df[
    plot_df$Plot_class ==
      "Strict pathway rescue",
    ,
    drop = FALSE
  ]
  
  label_df <- label_df[
    order(
      label_df$Combined_FDR,
      -label_df$Rescue_strength,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  label_n <- min(
    TOP_PATHWAY_LABELS,
    nrow(label_df)
  )
  
  if (label_n > 0) {
    label_df <- label_df[
      seq_len(label_n),
      ,
      drop = FALSE
    ]
  } else {
    label_df <- plot_df[
      0,
      ,
      drop = FALSE
    ]
  }
  
  label_df$Label <-
    stringr::str_wrap(
      label_df$Description_final,
      width = 34
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = Loss_NES,
      y = Rescue_NES,
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
      size = 2.1,
      alpha = 0.70
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        x = Loss_NES,
        y = Rescue_NES,
        label = Label
      ),
      inherit.aes = FALSE,
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      title = title_text,
      subtitle =
        "Same-sign significant enrichment in both contrasts indicates pathway-level rescue",
      x =
        "NES: Intact - P3 2D",
      y =
        "NES: P3 3D - P3 2D",
      color =
        "Classification"
    ) +
    theme_classic(
      base_size = 12
    )
  
  ggsave(
    filename = paste0(
      "figures/GSE44765_pathway/",
      output_prefix,
      ".pdf"
    ),
    plot = p,
    width = 9,
    height = 7
  )
  
  ggsave(
    filename = paste0(
      "figures/GSE44765_pathway/",
      output_prefix,
      ".png"
    ),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  p
}


############################################################
# 21. Save GO / Reactome rescue landscapes
############################################################

p_go_rescue <- plot_pathway_rescue(
  go_pathway_rescue,
  "GO BP pathway-level transcriptional rescue",
  "GSEA_GO_BP_pathway_rescue_landscape"
)

p_reactome_rescue <- plot_pathway_rescue(
  reactome_pathway_rescue,
  "Reactome pathway-level transcriptional rescue",
  "GSEA_Reactome_pathway_rescue_landscape"
)


############################################################
# 22. GO pathway trajectory / heatmap
#
# Uses full GO pathway membership from org.Hs.eg.db.
# This section is descriptive visualization.
############################################################

relevant_samples <- rownames(meta)[
  meta$group %in%
    c(
      "Intact",
      "P3_2D",
      "P3_3D"
    )
]

trajectory_meta <- meta[
  relevant_samples,
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

expr_trajectory <- expr_gene[
  ,
  relevant_samples,
  drop = FALSE
]

gene_z <- t(
  scale(
    t(
      expr_trajectory
    )
  )
)

gene_z[
  !is.finite(gene_z)
] <- 0


############################################################
# 23. Select top strict GO pathways
############################################################

if (nrow(go_strict) > 0) {
  
  go_strict <- go_strict[
    order(
      go_strict$Combined_FDR,
      -go_strict$Rescue_strength,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  n_top_go <- min(
    TOP_TRAJECTORY_PATHWAYS,
    nrow(go_strict)
  )
  
  top_go_ids <- go_strict$PathwayID[
    seq_len(n_top_go)
  ]
  
  top_go_desc <- go_strict$Description_final[
    seq_len(n_top_go)
  ]
  
  
  ##########################################################
  # 24. Retrieve full gene membership for each GO term
  ##########################################################
  
  pathway_score_rows <- list()
  score_counter <- 1
  
  for (i in seq_along(top_go_ids)) {
    
    go_id <- top_go_ids[i]
    go_desc <- top_go_desc[i]
    
    go_members <- tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = go_id,
        keytype = "GOALL",
        columns = c(
          "SYMBOL"
        )
      ),
      error = function(e) {
        data.frame()
      }
    )
    
    if (
      is.null(go_members) ||
      nrow(go_members) == 0 ||
      !"SYMBOL" %in%
      colnames(go_members)
    ) {
      next
    }
    
    pathway_genes <- unique(
      go_members$SYMBOL
    )
    
    pathway_genes <- pathway_genes[
      !is.na(pathway_genes) &
        pathway_genes != ""
    ]
    
    pathway_genes <- intersect(
      pathway_genes,
      rownames(gene_z)
    )
    
    if (
      length(pathway_genes) <
      MIN_GO_GENES_FOR_SCORE
    ) {
      next
    }
    
    pathway_scores <- colMeans(
      gene_z[
        pathway_genes,
        ,
        drop = FALSE
      ]
    )
    
    pathway_score_rows[[score_counter]] <-
      pathway_scores
    
    names(pathway_score_rows)[score_counter] <-
      go_desc
    
    score_counter <-
      score_counter + 1
  }
  
  
  ##########################################################
  # 25. Build pathway score matrix
  ##########################################################
  
  if (length(pathway_score_rows) > 0) {
    
    pathway_score_matrix <- do.call(
      rbind,
      pathway_score_rows
    )
    
    pathway_score_matrix <-
      pathway_score_matrix[
        ,
        relevant_samples,
        drop = FALSE
      ]
    
    pathway_score_long_list <- list()
    long_counter <- 1
    
    for (
      pathway_name in
      rownames(pathway_score_matrix)
    ) {
      
      temp <- data.frame(
        Pathway = pathway_name,
        Sample =
          colnames(
            pathway_score_matrix
          ),
        Score =
          as.numeric(
            pathway_score_matrix[
              pathway_name,
              ,
              drop = TRUE
            ]
          ),
        stringsAsFactors = FALSE
      )
      
      temp$Group <- as.character(
        trajectory_meta[
          temp$Sample,
          "group"
        ]
      )
      
      temp$Donor <- as.character(
        trajectory_meta[
          temp$Sample,
          "donor"
        ]
      )
      
      pathway_score_long_list[[long_counter]] <-
        temp
      
      long_counter <-
        long_counter + 1
    }
    
    pathway_score_long <- do.call(
      rbind,
      pathway_score_long_list
    )
    
    pathway_score_long$Group <- factor(
      pathway_score_long$Group,
      levels = c(
        "Intact",
        "P3_2D",
        "P3_3D"
      )
    )
    
    pathway_score_long$Pathway_label <-
      stringr::str_wrap(
        pathway_score_long$Pathway,
        width = 32
      )
    
    save_csv(
      pathway_score_long,
      "results/GSE44765_pathway/Top_rescued_GO_pathway_sample_scores.csv"
    )
    
    
    ########################################################
    # 26. Pathway trajectory plot
    ########################################################
    
    p_pathway_trajectory <- ggplot(
      pathway_score_long,
      aes(
        x = Group,
        y = Score,
        group = Donor,
        shape = Donor
      )
    ) +
      geom_line(
        alpha = 0.65
      ) +
      geom_point(
        size = 2.5
      ) +
      facet_wrap(
        ~ Pathway_label,
        scales = "free_y",
        ncol = 3
      ) +
      labs(
        title =
          "Top rescued GO BP pathway trajectories",
        subtitle =
          "Mean gene-wise Z score; each line represents the same donor",
        x = NULL,
        y =
          "Pathway expression score",
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
      filename =
        "figures/GSE44765_pathway/Top_rescued_GO_pathway_trajectories.pdf",
      plot =
        p_pathway_trajectory,
      width = 12,
      height = 11
    )
    
    ggsave(
      filename =
        "figures/GSE44765_pathway/Top_rescued_GO_pathway_trajectories.png",
      plot =
        p_pathway_trajectory,
      width = 12,
      height = 11,
      dpi = 300
    )
    
    
    ########################################################
    # 27. Pathway heatmap
    ########################################################
    
    heatmap_annotation <- data.frame(
      Group =
        as.character(
          trajectory_meta$group
        ),
      Donor =
        as.character(
          trajectory_meta$donor
        ),
      stringsAsFactors = FALSE
    )
    
    rownames(
      heatmap_annotation
    ) <- rownames(
      trajectory_meta
    )
    
    rownames(
      pathway_score_matrix
    ) <- stringr::str_wrap(
      rownames(
        pathway_score_matrix
      ),
      width = 36
    )
    
    pheatmap(
      pathway_score_matrix,
      annotation_col =
        heatmap_annotation,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      fontsize_row = 8,
      main =
        "Top rescued GO BP pathway scores",
      filename =
        "figures/GSE44765_pathway/Top_rescued_GO_pathway_heatmap.pdf",
      width = 10,
      height = 8
    )
    
    pheatmap(
      pathway_score_matrix,
      annotation_col =
        heatmap_annotation,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      fontsize_row = 8,
      main =
        "Top rescued GO BP pathway scores",
      filename =
        "figures/GSE44765_pathway/Top_rescued_GO_pathway_heatmap.png",
      width = 10,
      height = 8
    )
  }
}


############################################################
# 28. Count significant ORA pathways
############################################################

count_significant <- function(df) {
  
  if (
    is.null(df) ||
    nrow(df) == 0 ||
    !"p.adjust" %in%
    colnames(df)
  ) {
    return(0)
  }
  
  sum(
    !is.na(
      df$p.adjust
    ) &
      df$p.adjust <
      FDR_CUTOFF
  )
}

n_go_up <-
  count_significant(
    go_up_df
  )

n_go_down <-
  count_significant(
    go_down_df
  )

n_reactome_up <-
  count_significant(
    reactome_up_df
  )

n_reactome_down <-
  count_significant(
    reactome_down_df
  )


############################################################
# 29. Numeric summary
############################################################

pathway_summary <- data.frame(
  Metric = c(
    "Strict_rescue_genes_total",
    "Strict_rescue_genes_upward",
    "Strict_rescue_genes_downward",
    "GO_BP_ORA_significant_upward",
    "GO_BP_ORA_significant_downward",
    "Reactome_ORA_significant_upward",
    "Reactome_ORA_significant_downward",
    "GO_BP_strict_rescued_pathways_GSEA",
    "Reactome_strict_rescued_pathways_GSEA"
  ),
  Value = c(
    nrow(strict_rescue),
    length(rescue_up),
    length(rescue_down),
    n_go_up,
    n_go_down,
    n_reactome_up,
    n_reactome_down,
    nrow(go_strict),
    nrow(reactome_strict)
  )
)

save_csv(
  pathway_summary,
  "results/GSE44765_pathway/GSE44765_pathway_summary.csv"
)


############################################################
# 30. Human-readable summary
############################################################

summary_file <-
  "results/GSE44765_pathway/ANALYSIS_SUMMARY_pathways.txt"

sink(summary_file)

cat(
  "GSE44765 pathway-level rescue analysis\n"
)

cat(
  "======================================\n\n"
)

cat(
  "Version: offline / no-msigdbr / ReactomePA-version-compatible\n\n"
)

cat(
  "Strict rescue genes:",
  nrow(strict_rescue),
  "\n"
)

cat(
  "Restored upward genes:",
  length(rescue_up),
  "\n"
)

cat(
  "Restored downward genes:",
  length(rescue_down),
  "\n\n"
)

cat(
  "Significant GO BP ORA pathways - upward:",
  n_go_up,
  "\n"
)

cat(
  "Significant GO BP ORA pathways - downward:",
  n_go_down,
  "\n"
)

cat(
  "Significant Reactome ORA pathways - upward:",
  n_reactome_up,
  "\n"
)

cat(
  "Significant Reactome ORA pathways - downward:",
  n_reactome_down,
  "\n\n"
)

cat(
  "Strict GO BP pathway rescue by GSEA:",
  nrow(go_strict),
  "\n"
)

cat(
  "Strict Reactome pathway rescue by GSEA:",
  nrow(reactome_strict),
  "\n\n"
)

if (nrow(go_strict) > 0) {
  
  cat(
    "Top strict GO BP rescued pathways:\n"
  )
  
  top_n <- min(
    15,
    nrow(go_strict)
  )
  
  for (i in seq_len(top_n)) {
    cat(
      i,
      ". ",
      go_strict$Description_final[i],
      " | Loss NES=",
      round(
        go_strict$Loss_NES[i],
        3
      ),
      " | Rescue NES=",
      round(
        go_strict$Rescue_NES[i],
        3
      ),
      " | Combined FDR=",
      signif(
        go_strict$Combined_FDR[i],
        3
      ),
      "\n",
      sep = ""
    )
  }
  
  cat("\n")
}

if (nrow(reactome_strict) > 0) {
  
  cat(
    "Top strict Reactome rescued pathways:\n"
  )
  
  top_n <- min(
    15,
    nrow(
      reactome_strict
    )
  )
  
  for (i in seq_len(top_n)) {
    cat(
      i,
      ". ",
      reactome_strict$Description_final[i],
      " | Loss NES=",
      round(
        reactome_strict$Loss_NES[i],
        3
      ),
      " | Rescue NES=",
      round(
        reactome_strict$Rescue_NES[i],
        3
      ),
      " | Combined FDR=",
      signif(
        reactome_strict$Combined_FDR[i],
        3
      ),
      "\n",
      sep = ""
    )
  }
}

sink()


############################################################
# 31. Save session info
############################################################

sink(
  "results/GSE44765_pathway/sessionInfo_pathway_analysis.txt"
)

sessionInfo()

sink()


############################################################
# 32. Save workspace
############################################################

save(
  strict_rescue,
  rescue_up,
  rescue_down,
  go_up_df,
  go_down_df,
  reactome_up_df,
  reactome_down_df,
  go_loss_gsea_df,
  go_rescue_gsea_df,
  reactome_loss_gsea_df,
  reactome_rescue_gsea_df,
  go_pathway_rescue,
  reactome_pathway_rescue,
  go_strict,
  reactome_strict,
  pathway_summary,
  file =
    "results/GSE44765_pathway/GSE44765_pathway_workspace.RData"
)


############################################################
# 33. Final report
############################################################

cat("\n")
cat("========================================================\n")
cat(" PATHWAY ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("========================================================\n\n")

cat(
  "Strict rescue genes:",
  nrow(strict_rescue),
  "\n"
)

cat(
  "Restored upward genes:",
  length(rescue_up),
  "\n"
)

cat(
  "Restored downward genes:",
  length(rescue_down),
  "\n\n"
)

cat(
  "GO BP ORA significant pathways (up/down):",
  n_go_up,
  "/",
  n_go_down,
  "\n"
)

cat(
  "Reactome ORA significant pathways (up/down):",
  n_reactome_up,
  "/",
  n_reactome_down,
  "\n\n"
)

cat(
  "STRICT GO BP rescued pathways:",
  nrow(go_strict),
  "\n"
)

cat(
  "STRICT Reactome rescued pathways:",
  nrow(reactome_strict),
  "\n\n"
)

cat(
  "Main figures:\n"
)

cat(
  "figures/GSE44765_pathway/GO_BP_restored_upward.png\n"
)

cat(
  "figures/GSE44765_pathway/GO_BP_restored_downward.png\n"
)

cat(
  "figures/GSE44765_pathway/Reactome_restored_upward.png\n"
)

cat(
  "figures/GSE44765_pathway/Reactome_restored_downward.png\n"
)

cat(
  "figures/GSE44765_pathway/GSEA_GO_BP_pathway_rescue_landscape.png\n"
)

cat(
  "figures/GSE44765_pathway/GSEA_Reactome_pathway_rescue_landscape.png\n"
)

if (
  file.exists(
    "figures/GSE44765_pathway/Top_rescued_GO_pathway_trajectories.png"
  )
) {
  cat(
    "figures/GSE44765_pathway/Top_rescued_GO_pathway_trajectories.png\n"
  )
}

if (
  file.exists(
    "figures/GSE44765_pathway/Top_rescued_GO_pathway_heatmap.png"
  )
) {
  cat(
    "figures/GSE44765_pathway/Top_rescued_GO_pathway_heatmap.png\n"
  )
}

cat("\nMain result tables:\n")

cat(
  "results/GSE44765_pathway/GO_BP_STRICT_rescued_pathways.csv\n"
)

cat(
  "results/GSE44765_pathway/Reactome_STRICT_rescued_pathways.csv\n"
)

cat(
  "results/GSE44765_pathway/GSE44765_pathway_summary.csv\n"
)

cat(
  "results/GSE44765_pathway/ANALYSIS_SUMMARY_pathways.txt\n"
)

cat("\n")
cat(
  "This version did not use msigdbr or download MSigDB from Zenodo.\n"
)
cat(
  "GO and Reactome were analyzed using installed Bioconductor annotation/resources.\n"
)
cat("========================================================\n")
