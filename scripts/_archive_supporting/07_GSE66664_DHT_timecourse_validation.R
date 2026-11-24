############################################################
# 07_GSE66664_DHT_timecourse_validation.R
#
# COMPLETE ONE-CLICK VERSION
#
# Purpose
# -------
# Analyze the 140-array GSE66664 BAB/BAN DHT time course
# using the PREDEFINED, FROZEN GSE44765 strict rescue signature.
#
# This script answers four questions:
#
# 1) Baseline replication:
#    At 0 h, is non-balding BAN more rescue-like than balding BAB?
#
# 2) Dynamic rescue-score response:
#    How does the predefined rescue score change after DHT
#    in BAB and BAN at 1 nM and 10 nM?
#
# 3) Differential DHT response:
#    Does BAN respond differently from BAB?
#    Interaction = (BAN_t - BAN_0) - (BAB_t - BAB_0)
#
# 4) Cross-study temporal concordance:
#    At each time point, how strongly does the DHT response
#    align with the GSE44765 P3_3D - P3_2D rescue effect?
#
# IMPORTANT INTERPRETATION
# ------------------------
# GSE66664 contains repeated experimental measurements of
# BAB and BAN immortalized dermal papilla cell-line models.
# The 140 arrays are NOT 140 independent patients.
#
# Input
# -----
# data_raw/GSE66664/
#   GSE66664_non-normalized_data.txt.gz
#   OR a downloaded duplicate such as:
#   GSE66664_non-normalized_data.txt(1).gz
#
# GPL10558 annotation is searched automatically in:
#   data_raw/GSE66664/
#   data_raw/GSE66663/
#
# Frozen discovery signature:
#   results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv
# or any matching file under results/
#
# Main outputs
# ------------
# results/GSE66664_DHT/
# figures/GSE66664_DHT/
############################################################


############################################################
# 0. Basic settings
############################################################

graphics.off()

options(
  stringsAsFactors = FALSE
)

set.seed(
  66664
)

cat("\n")
cat("============================================================\n")
cat(" GSE66664 DHT TIME-COURSE VALIDATION STARTED\n")
cat("============================================================\n\n")

cat(
  "Working directory:\n",
  getwd(),
  "\n\n"
)


############################################################
# 1. Project paths
############################################################

raw_dir <-
  "data_raw/GSE66664"

annotation_search_dirs <- c(
  "data_raw/GSE66664",
  "data_raw/GSE66663"
)

results_dir <-
  "results/GSE66664_DHT"

figures_dir <-
  "figures/GSE66664_DHT"

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
# 2. Locate frozen GSE44765 signature
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
  
  discovery_signature_file <-
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
        "Cannot find GSE44765_STRICT_rescue_genes.csv under results/.\n",
        "Run the GSE44765 rescue analysis first.\n"
      )
    )
  }
  
  discovery_signature_file <-
    signature_candidates[1]
}

cat(
  "Frozen discovery signature:\n",
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

if (
  !requireNamespace(
    "limma",
    quietly = TRUE
  )
) {
  
  BiocManager::install(
    "limma",
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  
  library(
    limma
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
# 4. Locate GSE66664 non-normalized matrix
#
# This intentionally accepts duplicate-download names such as:
# GSE66664_non-normalized_data.txt(1).gz
############################################################

matrix_candidates <- list.files(
  path = raw_dir,
  pattern = "^GSE66664_non-normalized_data.*\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (
  length(
    matrix_candidates
  ) == 0
) {
  
  stop(
    paste0(
      "\nERROR:\n",
      "No GSE66664 non-normalized matrix found in:\n",
      raw_dir,
      "\n\nExpected a name beginning with:\n",
      "GSE66664_non-normalized_data\n"
    )
  )
}

candidate_sizes <- file.info(
  matrix_candidates
)$size

non_normalized_file <-
  matrix_candidates[
    which.max(
      candidate_sizes
    )
  ]

if (
  length(
    matrix_candidates
  ) > 1
) {
  
  cat(
    "Multiple candidate matrices found.\n",
    "Using the largest file:\n",
    non_normalized_file,
    "\n\n"
  )
  
} else {
  
  cat(
    "Non-normalized matrix:\n",
    non_normalized_file,
    "\n\n"
  )
}


############################################################
# 5. Read matrix
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

if (
  !"ILMN_ID" %in%
  colnames(
    raw_table
  )
) {
  
  stop(
    "ERROR: ILMN_ID column is missing."
  )
}

cat(
  "Matrix rows:",
  nrow(
    raw_table
  ),
  "\n"
)

cat(
  "Matrix columns:",
  ncol(
    raw_table
  ),
  "\n\n"
)


############################################################
# 6. Detect the 140 signal columns
############################################################

signal_columns <- grep(
  "^AVG_Signal\\.",
  colnames(
    raw_table
  ),
  value = TRUE
)

detection_columns <- grep(
  "^Detection\\.",
  colnames(
    raw_table
  ),
  value = TRUE
)

cat(
  "AVG_Signal columns:",
  length(
    signal_columns
  ),
  "\n"
)

cat(
  "Detection columns:",
  length(
    detection_columns
  ),
  "\n\n"
)

if (
  length(
    signal_columns
  ) != 140
) {
  
  stop(
    paste0(
      "ERROR: expected 140 AVG_Signal columns, found ",
      length(
        signal_columns
      ),
      "."
    )
  )
}

if (
  length(
    detection_columns
  ) != 140
) {
  
  stop(
    paste0(
      "ERROR: expected 140 Detection columns, found ",
      length(
        detection_columns
      ),
      "."
    )
  )
}


############################################################
# 7. Parse sample names directly from matrix header
############################################################

signal_suffix <- sub(
  "^AVG_Signal\\.",
  "",
  signal_columns
)

sample_pattern <-
  "^(BAB|BAN)_(1nM|10nM) DHT_([^_]+)_rep([A-Za-z0-9]+)$"

parsed_match <- regexec(
  sample_pattern,
  signal_suffix
)

parsed_values <- regmatches(
  signal_suffix,
  parsed_match
)

parse_ok <- lengths(
  parsed_values
) == 5

if (
  !all(
    parse_ok
  )
) {
  
  stop(
    paste0(
      "ERROR: failed to parse sample names:\n",
      paste(
        signal_suffix[
          !parse_ok
        ],
        collapse = "\n"
      )
    )
  )
}

CellLine <- vapply(
  parsed_values,
  function(x) {
    x[2]
  },
  FUN.VALUE = character(1)
)

Dose <- vapply(
  parsed_values,
  function(x) {
    x[3]
  },
  FUN.VALUE = character(1)
)

Time <- vapply(
  parsed_values,
  function(x) {
    x[4]
  },
  FUN.VALUE = character(1)
)

Replicate <- vapply(
  parsed_values,
  function(x) {
    x[5]
  },
  FUN.VALUE = character(1)
)


############################################################
# 8. Time conversion and canonical ordering
############################################################

time_levels <- c(
  "0h",
  "15min",
  "30min",
  "1h",
  "3h",
  "6h",
  "12h",
  "16h",
  "20h",
  "24h",
  "36h",
  "48h"
)

time_to_minutes <- function(
    x
) {
  
  if (
    grepl(
      "min$",
      x
    )
  ) {
    
    return(
      as.numeric(
        sub(
          "min$",
          "",
          x
        )
      )
    )
  }
  
  if (
    grepl(
      "h$",
      x
    )
  ) {
    
    return(
      60 *
        as.numeric(
          sub(
            "h$",
            "",
            x
          )
        )
    )
  }
  
  return(
    NA_real_
  )
}

TimeMinutes <- vapply(
  Time,
  time_to_minutes,
  FUN.VALUE = numeric(1)
)

if (
  any(
    !Time %in%
    time_levels
  )
) {
  
  stop(
    paste0(
      "ERROR: unexpected time labels detected:\n",
      paste(
        unique(
          Time[
            !Time %in%
              time_levels
          ]
        ),
        collapse = ", "
      )
    )
  )
}


############################################################
# 9. Build metadata
############################################################

meta <- data.frame(
  
  SignalColumn =
    signal_columns,
  
  Sample =
    signal_suffix,
  
  CellLine =
    CellLine,
  
  Dose =
    Dose,
  
  Time =
    Time,
  
  TimeMinutes =
    TimeMinutes,
  
  Replicate =
    Replicate,
  
  stringsAsFactors =
    FALSE
)

meta$CellLine <- factor(
  meta$CellLine,
  levels = c(
    "BAB",
    "BAN"
  )
)

meta$Dose <- factor(
  meta$Dose,
  levels = c(
    "1nM",
    "10nM"
  )
)

meta$Time <- factor(
  meta$Time,
  levels =
    time_levels,
  ordered = TRUE
)

meta <- meta[
  order(
    meta$CellLine,
    meta$Dose,
    meta$TimeMinutes,
    meta$Replicate
  ),
  ,
  drop = FALSE
]

rownames(
  meta
) <- meta$Sample

write.csv(
  meta,
  file.path(
    results_dir,
    "GSE66664_sample_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 10. Verify all signal columns have matching Detection column
############################################################

expected_detection_columns <- paste0(
  "Detection.",
  meta$Sample
)

if (
  !all(
    expected_detection_columns %in%
    colnames(
      raw_table
    )
  )
) {
  
  stop(
    "ERROR: at least one signal column lacks its matching Detection column."
  )
}


############################################################
# 11. Summarize experimental-cell replicate counts
############################################################

replicate_count <- aggregate(
  Replicate ~
    CellLine +
    Dose +
    Time,
  data = meta,
  FUN = length
)

colnames(
  replicate_count
)[
  colnames(
    replicate_count
  ) == "Replicate"
] <- "N"

replicate_count$TimeMinutes <- vapply(
  as.character(
    replicate_count$Time
  ),
  time_to_minutes,
  FUN.VALUE = numeric(1)
)

replicate_count <- replicate_count[
  order(
    replicate_count$CellLine,
    replicate_count$Dose,
    replicate_count$TimeMinutes
  ),
  ,
  drop = FALSE
]

write.csv(
  replicate_count,
  file.path(
    results_dir,
    "GSE66664_condition_replicate_counts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

incomplete_cells <- replicate_count[
  replicate_count$N < 3,
  ,
  drop = FALSE
]

cat(
  "Experimental cells:",
  nrow(
    replicate_count
  ),
  "\n"
)

cat(
  "Cells with fewer than 3 arrays:",
  nrow(
    incomplete_cells
  ),
  "\n\n"
)

if (
  nrow(
    incomplete_cells
  ) > 0
) {
  
  print(
    incomplete_cells
  )
  
  cat("\n")
}


############################################################
# 12. Build raw-signal matrix in metadata order
############################################################

raw_signal <- as.matrix(
  raw_table[
    ,
    meta$SignalColumn,
    drop = FALSE
  ]
)

storage.mode(
  raw_signal
) <- "numeric"

rownames(
  raw_signal
) <- raw_table$ILMN_ID

colnames(
  raw_signal
) <- meta$Sample

if (
  any(
    !is.finite(
      raw_signal
    )
  )
) {
  
  stop(
    "ERROR: non-finite raw signal values detected."
  )
}

if (
  any(
    raw_signal <= 0
  )
) {
  
  stop(
    "ERROR: non-positive raw signal values detected."
  )
}


############################################################
# 13. Detection-score matrix and QC
############################################################

detection_matrix <- as.matrix(
  raw_table[
    ,
    expected_detection_columns,
    drop = FALSE
  ]
)

storage.mode(
  detection_matrix
) <- "numeric"

rownames(
  detection_matrix
) <- raw_table$ILMN_ID

colnames(
  detection_matrix
) <- meta$Sample

detection_qc <- data.frame(
  
  Sample =
    meta$Sample,
  
  CellLine =
    as.character(
      meta$CellLine
    ),
  
  Dose =
    as.character(
      meta$Dose
    ),
  
  Time =
    as.character(
      meta$Time
    ),
  
  Replicate =
    meta$Replicate,
  
  MedianRawSignal =
    apply(
      raw_signal,
      2,
      median,
      na.rm = TRUE
    ),
  
  MeanRawSignal =
    colMeans(
      raw_signal,
      na.rm = TRUE
    ),
  
  FractionDetectionScore_ge_0_95 =
    colMeans(
      detection_matrix >= 0.95,
      na.rm = TRUE
    ),
  
  stringsAsFactors =
    FALSE
)

write.csv(
  detection_qc,
  file.path(
    results_dir,
    "GSE66664_raw_QC.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 14. Log2 + quantile normalization
#
# This follows the processing strategy used for GSE66663
# and the supplied GEO non-normalized signal matrices.
############################################################

log2_signal <- log2(
  raw_signal
)

expr_probe <- limma::normalizeBetweenArrays(
  log2_signal,
  method = "quantile"
)

saveRDS(
  expr_probe,
  file.path(
    results_dir,
    "GSE66664_probe_level_log2_quantile_expression.rds"
  )
)


############################################################
# 15. Read [Probes] section from GPL10558 annotation
############################################################

read_gpl10558_probes <- function(
    filename
) {
  
  con <- gzfile(
    filename,
    open = "rt"
  )
  
  on.exit(
    close(
      con
    ),
    add = TRUE
  )
  
  all_lines <- readLines(
    con,
    warn = FALSE
  )
  
  probe_section <- which(
    trimws(
      all_lines
    ) == "[Probes]"
  )
  
  control_section <- which(
    trimws(
      all_lines
    ) == "[Controls]"
  )
  
  if (
    length(
      probe_section
    ) == 0 ||
    length(
      control_section
    ) == 0
  ) {
    
    stop(
      paste0(
        "ERROR: cannot locate [Probes]/[Controls] in ",
        basename(
          filename
        )
      )
    )
  }
  
  probe_start <-
    probe_section[1] + 1
  
  eligible_controls <- control_section[
    control_section >
      probe_section[1]
  ]
  
  control_start <-
    eligible_controls[1]
  
  probe_end <-
    control_start - 1
  
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
  
  unlink(
    temp_file
  )
  
  required_columns <- c(
    "Probe_Id",
    "Symbol"
  )
  
  if (
    !all(
      required_columns %in%
      colnames(
        probe_table
      )
    )
  ) {
    
    stop(
      paste0(
        "ERROR: annotation ",
        basename(
          filename
        ),
        " lacks Probe_Id/Symbol."
      )
    )
  }
  
  probe_table
}


############################################################
# 16. Locate annotation files
############################################################

annotation_candidates <- character(0)

for (
  current_dir in annotation_search_dirs
) {
  
  if (
    dir.exists(
      current_dir
    )
  ) {
    
    current_files <- list.files(
      path = current_dir,
      pattern = "^GPL10558.*\\.txt\\.gz$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    annotation_candidates <- c(
      annotation_candidates,
      current_files
    )
  }
}

annotation_candidates <- unique(
  annotation_candidates
)

if (
  length(
    annotation_candidates
  ) == 0
) {
  
  stop(
    paste0(
      "\nERROR:\n",
      "No GPL10558 annotation file found in:\n",
      paste(
        annotation_search_dirs,
        collapse = "\n"
      ),
      "\n"
    )
  )
}

cat(
  "GPL10558 annotation candidates:",
  length(
    annotation_candidates
  ),
  "\n"
)


############################################################
# 17. Select annotation by exact probe overlap
############################################################

annotation_tables <- list()

annotation_match_rows <- list()

for (
  i in seq_along(
    annotation_candidates
  )
) {
  
  current_file <-
    annotation_candidates[i]
  
  current_annotation <-
    read_gpl10558_probes(
      current_file
    )
  
  n_overlap <- length(
    intersect(
      rownames(
        expr_probe
      ),
      current_annotation$Probe_Id
    )
  )
  
  annotation_tables[
    i
  ] <- list(
    current_annotation
  )
  
  annotation_match_rows[
    i
  ] <- list(
    data.frame(
      
      CandidateIndex =
        i,
      
      AnnotationFile =
        basename(
          current_file
        ),
      
      FullPath =
        current_file,
      
      AnnotationProbeRows =
        nrow(
          current_annotation
        ),
      
      ExpressionProbeRows =
        nrow(
          expr_probe
        ),
      
      MatchingProbeIDs =
        n_overlap,
      
      MatchFraction =
        n_overlap /
        nrow(
          expr_probe
        ),
      
      stringsAsFactors =
        FALSE
    )
  )
}

annotation_match_table <- do.call(
  rbind,
  annotation_match_rows
)

annotation_match_table <- annotation_match_table[
  order(
    -annotation_match_table$MatchingProbeIDs,
    -annotation_match_table$MatchFraction
  ),
  ,
  drop = FALSE
]

write.csv(
  annotation_match_table,
  file.path(
    results_dir,
    "GSE66664_annotation_probe_match.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

best_annotation_index <-
  annotation_match_table$CandidateIndex[1]

annotation <- getElement(
  annotation_tables,
  best_annotation_index
)

selected_annotation_file <-
  annotation_match_table$FullPath[1]

cat("\nAnnotation match:\n")
print(
  annotation_match_table
)

cat(
  "\nSelected annotation:\n",
  selected_annotation_file,
  "\n\n"
)


############################################################
# 18. Probe-to-gene mapping
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
  !is.na(
    annotation_small$Symbol
  ) &
    annotation_small$Symbol != "",
  ,
  drop = FALSE
]

annotation_small <- annotation_small[
  !duplicated(
    annotation_small$Probe_Id
  ),
  ,
  drop = FALSE
]

probe_to_symbol <-
  annotation_small$Symbol

names(
  probe_to_symbol
) <- annotation_small$Probe_Id

common_probe_ids <- intersect(
  rownames(
    expr_probe
  ),
  names(
    probe_to_symbol
  )
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
  nrow(
    expr_probe
  ),
  "\n"
)

cat(
  "Mapped probes with symbols:",
  nrow(
    expr_mapped
  ),
  "\n"
)


############################################################
# 19. Collapse probes to gene level
############################################################

expr_gene <- limma::avereps(
  expr_mapped,
  ID = gene_symbol
)

expr_gene <- expr_gene[
  !is.na(
    rownames(
      expr_gene
    )
  ) &
    rownames(
      expr_gene
    ) != "",
  ,
  drop = FALSE
]

cat(
  "Gene-level features:",
  nrow(
    expr_gene
  ),
  "\n\n"
)

saveRDS(
  expr_gene,
  file.path(
    results_dir,
    "GSE66664_gene_level_expression.rds"
  )
)


############################################################
# 20. PCA
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
  sum(
    pca$sdev^2
  )

pca_df <- data.frame(
  
  Sample =
    rownames(
      pca$x
    ),
  
  PC1 =
    pca$x[, 1],
  
  PC2 =
    pca$x[, 2],
  
  stringsAsFactors =
    FALSE
)

pca_df$CellLine <-
  meta[
    pca_df$Sample,
    "CellLine"
  ]

pca_df$Dose <-
  meta[
    pca_df$Sample,
    "Dose"
  ]

pca_df$Time <-
  meta[
    pca_df$Sample,
    "Time"
  ]

p_pca <- ggplot(
  
  pca_df,
  
  aes(
    x = PC1,
    y = PC2,
    color = CellLine,
    shape = Dose
  )
  
) +
  
  geom_point(
    size = 2.7,
    alpha = 0.75
  ) +
  
  labs(
    
    title =
      "GSE66664: BAB/BAN DHT time course",
    
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
    
    color =
      "Cell line",
    
    shape =
      "DHT series"
  ) +
  
  theme_classic(
    base_size = 13
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66664_PCA.png"
  ),
  p_pca,
  width = 8.5,
  height = 6.8,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66664_PCA.pdf"
  ),
  p_pca,
  width = 8.5,
  height = 6.8
)


############################################################
# 21. Load frozen discovery signature
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

signature_overlap <- signature[
  signature$GeneSymbol %in%
    rownames(
      expr_gene
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

signature_coverage <- nrow(
  signature_overlap
) /
  nrow(
    signature
  )

write.csv(
  signature_overlap,
  file.path(
    results_dir,
    "GSE44765_signature_overlap_in_GSE66664.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Discovery signature genes:",
  nrow(
    signature
  ),
  "\n"
)

cat(
  "Signature genes represented in GSE66664:",
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
  "Restored-upward overlap:",
  length(
    up_genes
  ),
  "\n"
)

cat(
  "Restored-downward overlap:",
  length(
    down_genes
  ),
  "\n\n"
)


############################################################
# 22. Predefined rescue score across all 140 arrays
#
# Score:
#   mean Z(restored-upward) - mean Z(restored-downward)
#
# Z scaling is gene-wise within GSE66664.
# This score is cohort-relative and intended for within-study
# time-course comparison, not as an absolute clinical score.
############################################################

signature_expr <- expr_gene[
  signature_overlap$GeneSymbol,
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

rescue_score <-
  up_score -
  down_score

score_df <- data.frame(
  
  Sample =
    colnames(
      expr_gene
    ),
  
  CellLine =
    as.character(
      meta[
        colnames(
          expr_gene
        ),
        "CellLine"
      ]
    ),
  
  Dose =
    as.character(
      meta[
        colnames(
          expr_gene
        ),
        "Dose"
      ]
    ),
  
  Time =
    as.character(
      meta[
        colnames(
          expr_gene
        ),
        "Time"
      ]
    ),
  
  TimeMinutes =
    meta[
      colnames(
        expr_gene
      ),
      "TimeMinutes"
    ],
  
  Replicate =
    meta[
      colnames(
        expr_gene
      ),
      "Replicate"
    ],
  
  Up_score =
    as.numeric(
      up_score[
        colnames(
          expr_gene
        )
      ]
    ),
  
  Down_score =
    as.numeric(
      down_score[
        colnames(
          expr_gene
        )
      ]
    ),
  
  Rescue_signature_score =
    as.numeric(
      rescue_score[
        colnames(
          expr_gene
        )
      ]
    ),
  
  stringsAsFactors =
    FALSE
)

score_df$CellLine <- factor(
  score_df$CellLine,
  levels = c(
    "BAB",
    "BAN"
  )
)

score_df$Dose <- factor(
  score_df$Dose,
  levels = c(
    "1nM",
    "10nM"
  )
)

score_df$Time <- factor(
  score_df$Time,
  levels =
    time_levels,
  ordered = TRUE
)

write.csv(
  score_df,
  file.path(
    results_dir,
    "GSE66664_predefined_rescue_scores_by_sample.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 23. Score summary by CellLine x Dose x Time
############################################################

score_split <- split(
  score_df$Rescue_signature_score,
  interaction(
    score_df$CellLine,
    score_df$Dose,
    score_df$Time,
    drop = TRUE
  )
)

score_summary_rows <- list()

score_counter <- 1

for (
  current_cell in c(
    "BAB",
    "BAN"
  )
) {
  
  for (
    current_dose in c(
      "1nM",
      "10nM"
    )
  ) {
    
    for (
      current_time in time_levels
    ) {
      
      values <- score_df$Rescue_signature_score[
        as.character(
          score_df$CellLine
        ) == current_cell &
          as.character(
            score_df$Dose
          ) == current_dose &
          as.character(
            score_df$Time
          ) == current_time
      ]
      
      if (
        length(
          values
        ) == 0
      ) {
        
        next
      }
      
      score_summary_rows[
        score_counter
      ] <- list(
        data.frame(
          
          CellLine =
            current_cell,
          
          Dose =
            current_dose,
          
          Time =
            current_time,
          
          TimeMinutes =
            time_to_minutes(
              current_time
            ),
          
          N =
            length(
              values
            ),
          
          Mean =
            mean(
              values
            ),
          
          SD =
            if (
              length(
                values
              ) > 1
            ) {
              sd(
                values
              )
            } else {
              NA_real_
            },
          
          Median =
            median(
              values
            ),
          
          stringsAsFactors =
            FALSE
        )
      )
      
      score_counter <-
        score_counter + 1
    }
  }
}

score_summary <- do.call(
  rbind,
  score_summary_rows
)

score_summary$Time <- factor(
  score_summary$Time,
  levels =
    time_levels,
  ordered = TRUE
)

write.csv(
  score_summary,
  file.path(
    results_dir,
    "GSE66664_predefined_rescue_score_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 24. Baseline means and change-from-baseline score
############################################################

baseline_means <- score_summary[
  as.character(
    score_summary$Time
  ) == "0h",
  c(
    "CellLine",
    "Dose",
    "Mean"
  ),
  drop = FALSE
]

colnames(
  baseline_means
)[
  colnames(
    baseline_means
  ) == "Mean"
] <- "BaselineMean"

score_df_with_baseline <- merge(
  score_df,
  baseline_means,
  by = c(
    "CellLine",
    "Dose"
  ),
  all.x = TRUE,
  sort = FALSE
)

score_df_with_baseline$DeltaFromBaseline <-
  score_df_with_baseline$Rescue_signature_score -
  score_df_with_baseline$BaselineMean

write.csv(
  score_df_with_baseline,
  file.path(
    results_dir,
    "GSE66664_rescue_score_change_from_baseline_by_sample.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 25. Mean delta-from-baseline by condition
############################################################

delta_summary_rows <- list()

delta_counter <- 1

for (
  current_cell in c(
    "BAB",
    "BAN"
  )
) {
  
  for (
    current_dose in c(
      "1nM",
      "10nM"
    )
  ) {
    
    for (
      current_time in time_levels
    ) {
      
      values <- score_df_with_baseline$DeltaFromBaseline[
        as.character(
          score_df_with_baseline$CellLine
        ) == current_cell &
          as.character(
            score_df_with_baseline$Dose
          ) == current_dose &
          as.character(
            score_df_with_baseline$Time
          ) == current_time
      ]
      
      if (
        length(
          values
        ) == 0
      ) {
        
        next
      }
      
      delta_summary_rows[
        delta_counter
      ] <- list(
        data.frame(
          
          CellLine =
            current_cell,
          
          Dose =
            current_dose,
          
          Time =
            current_time,
          
          TimeMinutes =
            time_to_minutes(
              current_time
            ),
          
          N =
            length(
              values
            ),
          
          MeanDelta =
            mean(
              values
            ),
          
          SDDelta =
            if (
              length(
                values
              ) > 1
            ) {
              sd(
                values
              )
            } else {
              NA_real_
            },
          
          stringsAsFactors =
            FALSE
        )
      )
      
      delta_counter <-
        delta_counter + 1
    }
  }
}

delta_summary <- do.call(
  rbind,
  delta_summary_rows
)

delta_summary$Time <- factor(
  delta_summary$Time,
  levels =
    time_levels,
  ordered = TRUE
)

write.csv(
  delta_summary,
  file.path(
    results_dir,
    "GSE66664_rescue_score_delta_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 26. BAN-minus-BAB interaction at score level
############################################################

interaction_score_rows <- list()

interaction_score_counter <- 1

for (
  current_dose in c(
    "1nM",
    "10nM"
  )
) {
  
  for (
    current_time in time_levels
  ) {
    
    bab_row <- delta_summary[
      delta_summary$CellLine == "BAB" &
        delta_summary$Dose == current_dose &
        as.character(
          delta_summary$Time
        ) == current_time,
      ,
      drop = FALSE
    ]
    
    ban_row <- delta_summary[
      delta_summary$CellLine == "BAN" &
        delta_summary$Dose == current_dose &
        as.character(
          delta_summary$Time
        ) == current_time,
      ,
      drop = FALSE
    ]
    
    if (
      nrow(
        bab_row
      ) == 1 &&
      nrow(
        ban_row
      ) == 1
    ) {
      
      interaction_score_rows[
        interaction_score_counter
      ] <- list(
        data.frame(
          
          Dose =
            current_dose,
          
          Time =
            current_time,
          
          TimeMinutes =
            time_to_minutes(
              current_time
            ),
          
          BAB_Delta =
            bab_row$MeanDelta,
          
          BAN_Delta =
            ban_row$MeanDelta,
          
          Interaction_BANminusBAB =
            ban_row$MeanDelta -
            bab_row$MeanDelta,
          
          stringsAsFactors =
            FALSE
        )
      )
      
      interaction_score_counter <-
        interaction_score_counter + 1
    }
  }
}

interaction_score <- do.call(
  rbind,
  interaction_score_rows
)

interaction_score$Time <- factor(
  interaction_score$Time,
  levels =
    time_levels,
  ordered = TRUE
)

write.csv(
  interaction_score,
  file.path(
    results_dir,
    "GSE66664_score_interaction_BANminusBAB.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 27. Rescue-score trajectory figure
############################################################

p_score_trajectory <- ggplot(
  
  score_summary,
  
  aes(
    x = Time,
    y = Mean,
    group = CellLine,
    color = CellLine
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  
  geom_line(
    linewidth = 0.85
  ) +
  
  geom_point(
    size = 2.4
  ) +
  
  geom_errorbar(
    
    aes(
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    
    width = 0.18,
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ Dose,
    ncol = 1,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "GSE66664 predefined rescue-state trajectory after DHT",
    
    subtitle =
      "Higher score = more GSE44765 rescue-like direction",
    
    x =
      "Time after DHT",
    
    y =
      "Predefined rescue signature score",
    
    color =
      "Cell line"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66664_rescue_score_timecourse.png"
  ),
  p_score_trajectory,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66664_rescue_score_timecourse.pdf"
  ),
  p_score_trajectory,
  width = 10,
  height = 8
)


############################################################
# 28. Score-level interaction figure
############################################################

p_score_interaction <- ggplot(
  
  interaction_score,
  
  aes(
    x = Time,
    y = Interaction_BANminusBAB,
    group = Dose,
    color = Dose
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2.6
  ) +
  
  labs(
    
    title =
      "Differential DHT response of the predefined rescue score",
    
    subtitle =
      "(BAN change from baseline) - (BAB change from baseline)",
    
    x =
      "Time after DHT",
    
    y =
      "BAN-minus-BAB rescue-score response",
    
    color =
      "DHT dose"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66664_rescue_score_interaction.png"
  ),
  p_score_interaction,
  width = 10,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66664_rescue_score_interaction.pdf"
  ),
  p_score_interaction,
  width = 10,
  height = 6.5
)


############################################################
# 29. Build limma design using all 48 experimental cells
############################################################

meta$Condition <- paste(
  as.character(
    meta$CellLine
  ),
  as.character(
    meta$Dose
  ),
  as.character(
    meta$Time
  ),
  sep = "_"
)

meta$Condition <- factor(
  meta$Condition
)

design <- model.matrix(
  ~ 0 + Condition,
  data = meta
)

colnames(
  design
) <- sub(
  "^Condition",
  "",
  colnames(
    design
  )
)

if (
  qr(
    design
  )$rank !=
  ncol(
    design
  )
) {
  
  stop(
    "ERROR: design matrix is not full rank."
  )
}

write.csv(
  design,
  file.path(
    results_dir,
    "GSE66664_design_matrix.csv"
  ),
  quote = FALSE
)


############################################################
# 30. Helper for contrast vectors
############################################################

make_contrast_vector <- function(
    named_weights,
    design_columns
) {
  
  v <- rep(
    0,
    length(
      design_columns
    )
  )
  
  names(
    v
  ) <- design_columns
  
  missing_terms <- setdiff(
    names(
      named_weights
    ),
    design_columns
  )
  
  if (
    length(
      missing_terms
    ) > 0
  ) {
    
    stop(
      paste0(
        "ERROR: contrast term(s) missing from design:\n",
        paste(
          missing_terms,
          collapse = "\n"
        )
      )
    )
  }
  
  v[
    names(
      named_weights
    )
  ] <- named_weights
  
  v
}


############################################################
# 31. Build all prespecified contrasts
############################################################

contrast_list <- list()

contrast_meta_rows <- list()

contrast_counter <- 1

design_columns <- colnames(
  design
)

add_contrast <- function(
    contrast_name,
    named_weights,
    contrast_type,
    cell_line_label,
    dose_label,
    time_label
) {
  
  v <- make_contrast_vector(
    named_weights,
    design_columns
  )
  
  contrast_list[
    contrast_name
  ] <<- list(
    v
  )
  
  contrast_meta_rows[
    contrast_counter
  ] <<- list(
    data.frame(
      
      Contrast =
        contrast_name,
      
      ContrastType =
        contrast_type,
      
      CellLine =
        cell_line_label,
      
      Dose =
        dose_label,
      
      Time =
        time_label,
      
      TimeMinutes =
        time_to_minutes(
          time_label
        ),
      
      stringsAsFactors =
        FALSE
    )
  )
  
  contrast_counter <<-
    contrast_counter + 1
}


############################################################
# 32. Baseline BAN - BAB contrasts
############################################################

for (
  current_dose in c(
    "1nM",
    "10nM"
  )
) {
  
  contrast_name <- paste0(
    "BAN_vs_BAB_",
    current_dose,
    "_0h"
  )
  
  weights <- c(
    1,
    -1
  )
  
  names(
    weights
  ) <- c(
    paste(
      "BAN",
      current_dose,
      "0h",
      sep = "_"
    ),
    paste(
      "BAB",
      current_dose,
      "0h",
      sep = "_"
    )
  )
  
  add_contrast(
    contrast_name =
      contrast_name,
    named_weights =
      weights,
    contrast_type =
      "BaselineCrossSectional",
    cell_line_label =
      "BANvsBAB",
    dose_label =
      current_dose,
    time_label =
      "0h"
  )
}

pooled_baseline_weights <- c(
  0.5,
  0.5,
  -0.5,
  -0.5
)

names(
  pooled_baseline_weights
) <- c(
  "BAN_1nM_0h",
  "BAN_10nM_0h",
  "BAB_1nM_0h",
  "BAB_10nM_0h"
)

add_contrast(
  contrast_name =
    "BAN_vs_BAB_pooled_0h",
  named_weights =
    pooled_baseline_weights,
  contrast_type =
    "BaselineCrossSectionalPooled",
  cell_line_label =
    "BANvsBAB",
  dose_label =
    "Pooled",
  time_label =
    "0h"
)


############################################################
# 33. Time-specific DHT response, interaction,
#     and cross-sectional contrasts
############################################################

post_times <- time_levels[
  time_levels != "0h"
]

for (
  current_dose in c(
    "1nM",
    "10nM"
  )
) {
  
  for (
    current_time in post_times
  ) {
    
    bab_t <- paste(
      "BAB",
      current_dose,
      current_time,
      sep = "_"
    )
    
    bab_0 <- paste(
      "BAB",
      current_dose,
      "0h",
      sep = "_"
    )
    
    ban_t <- paste(
      "BAN",
      current_dose,
      current_time,
      sep = "_"
    )
    
    ban_0 <- paste(
      "BAN",
      current_dose,
      "0h",
      sep = "_"
    )
    
    bab_name <- paste0(
      "BAB_",
      current_dose,
      "_",
      current_time,
      "_vs_0h"
    )
    
    bab_weights <- c(
      1,
      -1
    )
    
    names(
      bab_weights
    ) <- c(
      bab_t,
      bab_0
    )
    
    add_contrast(
      contrast_name =
        bab_name,
      named_weights =
        bab_weights,
      contrast_type =
        "DHTResponse",
      cell_line_label =
        "BAB",
      dose_label =
        current_dose,
      time_label =
        current_time
    )
    
    ban_name <- paste0(
      "BAN_",
      current_dose,
      "_",
      current_time,
      "_vs_0h"
    )
    
    ban_weights <- c(
      1,
      -1
    )
    
    names(
      ban_weights
    ) <- c(
      ban_t,
      ban_0
    )
    
    add_contrast(
      contrast_name =
        ban_name,
      named_weights =
        ban_weights,
      contrast_type =
        "DHTResponse",
      cell_line_label =
        "BAN",
      dose_label =
        current_dose,
      time_label =
        current_time
    )
    
    interaction_name <- paste0(
      "Interaction_BANminusBAB_",
      current_dose,
      "_",
      current_time
    )
    
    interaction_weights <- c(
      1,
      -1,
      -1,
      1
    )
    
    names(
      interaction_weights
    ) <- c(
      ban_t,
      ban_0,
      bab_t,
      bab_0
    )
    
    add_contrast(
      contrast_name =
        interaction_name,
      named_weights =
        interaction_weights,
      contrast_type =
        "Interaction",
      cell_line_label =
        "BANminusBAB",
      dose_label =
        current_dose,
      time_label =
        current_time
    )
    
    cross_name <- paste0(
      "BAN_vs_BAB_",
      current_dose,
      "_",
      current_time
    )
    
    cross_weights <- c(
      1,
      -1
    )
    
    names(
      cross_weights
    ) <- c(
      ban_t,
      bab_t
    )
    
    add_contrast(
      contrast_name =
        cross_name,
      named_weights =
        cross_weights,
      contrast_type =
        "CrossSectional",
      cell_line_label =
        "BANvsBAB",
      dose_label =
        current_dose,
      time_label =
        current_time
    )
  }
}

contrast_matrix <- do.call(
  cbind,
  contrast_list
)

contrast_meta <- do.call(
  rbind,
  contrast_meta_rows
)

rownames(
  contrast_meta
) <- NULL

write.csv(
  contrast_meta,
  file.path(
    results_dir,
    "GSE66664_contrast_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  contrast_matrix,
  file.path(
    results_dir,
    "GSE66664_contrast_matrix.csv"
  ),
  quote = FALSE
)

cat(
  "Prespecified limma contrasts:",
  ncol(
    contrast_matrix
  ),
  "\n\n"
)


############################################################
# 34. Fit limma model
############################################################

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


############################################################
# 35. Save full all-gene result for pooled baseline
############################################################

baseline_de <- limma::topTable(
  fit2,
  coef =
    "BAN_vs_BAB_pooled_0h",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

baseline_de$GeneSymbol <- rownames(
  baseline_de
)

baseline_de <- baseline_de[
  ,
  c(
    "GeneSymbol",
    setdiff(
      colnames(
        baseline_de
      ),
      "GeneSymbol"
    )
  ),
  drop = FALSE
]

write.csv(
  baseline_de,
  file.path(
    results_dir,
    "GSE66664_BAN_vs_BAB_pooled_0h_all_genes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 36. Extract signature-gene effects for every contrast
############################################################

signature_effect_rows <- list()

concordance_rows <- list()

effect_counter <- 1

for (
  contrast_name in colnames(
    contrast_matrix
  )
) {
  
  current_tt <- limma::topTable(
    fit2,
    coef =
      contrast_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )
  
  match_index <- match(
    signature_overlap$GeneSymbol,
    rownames(
      current_tt
    )
  )
  
  current_effects <- data.frame(
    
    GeneSymbol =
      signature_overlap$GeneSymbol,
    
    Rescue_type =
      signature_overlap$Rescue_type,
    
    Discovery_Rescue_logFC =
      signature_overlap$Rescue_logFC,
    
    External_logFC =
      current_tt$logFC[
        match_index
      ],
    
    External_t =
      current_tt$t[
        match_index
      ],
    
    External_P =
      current_tt$P.Value[
        match_index
      ],
    
    External_FDR =
      current_tt$adj.P.Val[
        match_index
      ],
    
    Contrast =
      contrast_name,
    
    stringsAsFactors =
      FALSE
  )
  
  current_meta <- contrast_meta[
    contrast_meta$Contrast ==
      contrast_name,
    ,
    drop = FALSE
  ]
  
  current_effects$ContrastType <-
    current_meta$ContrastType
  
  current_effects$CellLine <-
    current_meta$CellLine
  
  current_effects$Dose <-
    current_meta$Dose
  
  current_effects$Time <-
    current_meta$Time
  
  current_effects$TimeMinutes <-
    current_meta$TimeMinutes
  
  current_effects$ConcordantDirection <- sign(
    current_effects$Discovery_Rescue_logFC
  ) ==
    sign(
      current_effects$External_logFC
    )
  
  current_effects$ExpectedSignedExternalLogFC <- sign(
    current_effects$Discovery_Rescue_logFC
  ) *
    current_effects$External_logFC
  
  valid <- is.finite(
    current_effects$Discovery_Rescue_logFC
  ) &
    is.finite(
      current_effects$External_logFC
    )
  
  rho <- suppressWarnings(
    cor(
      current_effects$Discovery_Rescue_logFC[
        valid
      ],
      current_effects$External_logFC[
        valid
      ],
      method = "spearman"
    )
  )
  
  directional_fraction <- mean(
    current_effects$ConcordantDirection[
      valid
    ],
    na.rm = TRUE
  )
  
  median_expected <- median(
    current_effects$ExpectedSignedExternalLogFC[
      valid
    ],
    na.rm = TRUE
  )
  
  mean_expected <- mean(
    current_effects$ExpectedSignedExternalLogFC[
      valid
    ],
    na.rm = TRUE
  )
  
  signature_effect_rows[
    effect_counter
  ] <- list(
    current_effects
  )
  
  concordance_rows[
    effect_counter
  ] <- list(
    data.frame(
      
      Contrast =
        contrast_name,
      
      ContrastType =
        current_meta$ContrastType,
      
      CellLine =
        current_meta$CellLine,
      
      Dose =
        current_meta$Dose,
      
      Time =
        current_meta$Time,
      
      TimeMinutes =
        current_meta$TimeMinutes,
      
      N_signature_genes =
        sum(
          valid
        ),
      
      Spearman_rho =
        rho,
      
      Directional_concordance =
        directional_fraction,
      
      Median_expected_signed_external_logFC =
        median_expected,
      
      Mean_expected_signed_external_logFC =
        mean_expected,
      
      stringsAsFactors =
        FALSE
    )
  )
  
  effect_counter <-
    effect_counter + 1
}

signature_effects_long <- do.call(
  rbind,
  signature_effect_rows
)

concordance_summary <- do.call(
  rbind,
  concordance_rows
)

rownames(
  signature_effects_long
) <- NULL

rownames(
  concordance_summary
) <- NULL

write.csv(
  signature_effects_long,
  file.path(
    results_dir,
    "GSE66664_signature_gene_effects_all_contrasts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  concordance_summary,
  file.path(
    results_dir,
    "GSE66664_cross_study_concordance_all_contrasts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 37. CAMERA for the two frozen signature arms
############################################################

signature_index <- list(
  
  Restored_upward =
    which(
      rownames(
        expr_gene
      ) %in%
        up_genes
    ),
  
  Restored_downward =
    which(
      rownames(
        expr_gene
      ) %in%
        down_genes
    )
)

camera_rows <- list()

camera_counter <- 1

for (
  j in seq_len(
    ncol(
      contrast_matrix
    )
  )
) {
  
  contrast_name <- colnames(
    contrast_matrix
  )[j]
  
  current_camera <- limma::camera(
    
    y =
      expr_gene,
    
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
  
  current_camera$GeneSet <- rownames(
    current_camera
  )
  
  current_camera$Contrast <-
    contrast_name
  
  current_meta <- contrast_meta[
    contrast_meta$Contrast ==
      contrast_name,
    ,
    drop = FALSE
  ]
  
  current_camera$ContrastType <-
    current_meta$ContrastType
  
  current_camera$CellLine <-
    current_meta$CellLine
  
  current_camera$Dose <-
    current_meta$Dose
  
  current_camera$Time <-
    current_meta$Time
  
  current_camera$TimeMinutes <-
    current_meta$TimeMinutes
  
  current_camera <- current_camera[
    ,
    c(
      "Contrast",
      "ContrastType",
      "CellLine",
      "Dose",
      "Time",
      "TimeMinutes",
      "GeneSet",
      setdiff(
        colnames(
          current_camera
        ),
        c(
          "Contrast",
          "ContrastType",
          "CellLine",
          "Dose",
          "Time",
          "TimeMinutes",
          "GeneSet"
        )
      )
    ),
    drop = FALSE
  ]
  
  camera_rows[
    camera_counter
  ] <- list(
    current_camera
  )
  
  camera_counter <-
    camera_counter + 1
}

camera_results <- do.call(
  rbind,
  camera_rows
)

rownames(
  camera_results
) <- NULL

camera_results$WithinContrastFDR <-
  camera_results$FDR

camera_results$GlobalFDR <- p.adjust(
  camera_results$PValue,
  method = "BH"
)

write.csv(
  camera_results,
  file.path(
    results_dir,
    "GSE66664_CAMERA_all_contrasts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 38. Cross-study concordance time course:
#     DHT response in BAB and BAN
############################################################

dht_concordance <- concordance_summary[
  concordance_summary$ContrastType ==
    "DHTResponse",
  ,
  drop = FALSE
]

dht_concordance$Time <- factor(
  dht_concordance$Time,
  levels =
    time_levels,
  ordered = TRUE
)

p_dht_concordance <- ggplot(
  
  dht_concordance,
  
  aes(
    x = Time,
    y = Spearman_rho,
    group = CellLine,
    color = CellLine
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2.4
  ) +
  
  facet_wrap(
    ~ Dose,
    ncol = 1
  ) +
  
  labs(
    
    title =
      "Temporal concordance of DHT response with GSE44765 rescue",
    
    subtitle =
      "Positive rho = DHT response moves in the GSE44765 rescue direction",
    
    x =
      "Time after DHT",
    
    y =
      "Spearman rho",
    
    color =
      "Cell line"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66664_DHT_response_rescue_concordance_timecourse.png"
  ),
  p_dht_concordance,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66664_DHT_response_rescue_concordance_timecourse.pdf"
  ),
  p_dht_concordance,
  width = 10,
  height = 8
)


############################################################
# 39. Interaction concordance time course
#
# Positive rho:
#   BAN response is more rescue-like than BAB response.
#
# Negative rho:
#   BAB response is more rescue-like than BAN response.
############################################################

interaction_concordance <- concordance_summary[
  concordance_summary$ContrastType ==
    "Interaction",
  ,
  drop = FALSE
]

interaction_concordance$Time <- factor(
  interaction_concordance$Time,
  levels =
    time_levels,
  ordered = TRUE
)

p_interaction_concordance <- ggplot(
  
  interaction_concordance,
  
  aes(
    x = Time,
    y = Spearman_rho,
    group = Dose,
    color = Dose
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2.6
  ) +
  
  labs(
    
    title =
      "Differential BAB/BAN DHT-response concordance",
    
    subtitle =
      "Effect = (BAN_t - BAN_0) - (BAB_t - BAB_0)",
    
    x =
      "Time after DHT",
    
    y =
      "Spearman rho with GSE44765 rescue effect",
    
    color =
      "DHT dose"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66664_interaction_rescue_concordance_timecourse.png"
  ),
  p_interaction_concordance,
  width = 10,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66664_interaction_rescue_concordance_timecourse.pdf"
  ),
  p_interaction_concordance,
  width = 10,
  height = 6.5
)


############################################################
# 40. Cross-sectional BAN-vs-BAB concordance over time
############################################################

cross_sectional_concordance <- concordance_summary[
  concordance_summary$ContrastType %in% c(
    "BaselineCrossSectional",
    "CrossSectional"
  ),
  ,
  drop = FALSE
]

cross_sectional_concordance$Time <- factor(
  cross_sectional_concordance$Time,
  levels =
    time_levels,
  ordered = TRUE
)

p_cross_sectional <- ggplot(
  
  cross_sectional_concordance,
  
  aes(
    x = Time,
    y = Spearman_rho,
    group = Dose,
    color = Dose
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2.6
  ) +
  
  labs(
    
    title =
      "BAN-versus-BAB alignment with the GSE44765 rescue axis",
    
    subtitle =
      "Positive rho = BAN is more aligned with the rescue direction",
    
    x =
      "Time after DHT",
    
    y =
      "Spearman rho: GSE44765 rescue vs BAN-BAB",
    
    color =
      "DHT dose"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66664_BAN_vs_BAB_rescue_concordance_timecourse.png"
  ),
  p_cross_sectional,
  width = 10,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66664_BAN_vs_BAB_rescue_concordance_timecourse.pdf"
  ),
  p_cross_sectional,
  width = 10,
  height = 6.5
)


############################################################
# 41. Baseline score summary
############################################################

baseline_score_df <- score_df[
  as.character(
    score_df$Time
  ) == "0h",
  ,
  drop = FALSE
]

baseline_score_summary_rows <- list()

baseline_score_counter <- 1

for (
  current_dose in c(
    "1nM",
    "10nM"
  )
) {
  
  bab_values <- baseline_score_df$Rescue_signature_score[
    as.character(
      baseline_score_df$CellLine
    ) == "BAB" &
      as.character(
        baseline_score_df$Dose
      ) == current_dose
  ]
  
  ban_values <- baseline_score_df$Rescue_signature_score[
    as.character(
      baseline_score_df$CellLine
    ) == "BAN" &
      as.character(
        baseline_score_df$Dose
      ) == current_dose
  ]
  
  baseline_score_summary_rows[
    baseline_score_counter
  ] <- list(
    data.frame(
      
      Dose =
        current_dose,
      
      BAB_N =
        length(
          bab_values
        ),
      
      BAN_N =
        length(
          ban_values
        ),
      
      BAB_Mean =
        mean(
          bab_values
        ),
      
      BAN_Mean =
        mean(
          ban_values
        ),
      
      BAN_minus_BAB =
        mean(
          ban_values
        ) -
        mean(
          bab_values
        ),
      
      stringsAsFactors =
        FALSE
    )
  )
  
  baseline_score_counter <-
    baseline_score_counter + 1
}

all_bab_0h <- baseline_score_df$Rescue_signature_score[
  as.character(
    baseline_score_df$CellLine
  ) == "BAB"
]

all_ban_0h <- baseline_score_df$Rescue_signature_score[
  as.character(
    baseline_score_df$CellLine
  ) == "BAN"
]

baseline_score_summary_rows[
  baseline_score_counter
] <- list(
  data.frame(
    
    Dose =
      "PooledDescriptive",
    
    BAB_N =
      length(
        all_bab_0h
      ),
    
    BAN_N =
      length(
        all_ban_0h
      ),
    
    BAB_Mean =
      mean(
        all_bab_0h
      ),
    
    BAN_Mean =
      mean(
        all_ban_0h
      ),
    
    BAN_minus_BAB =
      mean(
        all_ban_0h
      ) -
      mean(
        all_bab_0h
      ),
    
    stringsAsFactors =
      FALSE
  )
)

baseline_score_summary <- do.call(
  rbind,
  baseline_score_summary_rows
)

write.csv(
  baseline_score_summary,
  file.path(
    results_dir,
    "GSE66664_baseline_rescue_score_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 42. Baseline signature score plot
############################################################

p_baseline_score <- ggplot(
  
  baseline_score_df,
  
  aes(
    x = CellLine,
    y = Rescue_signature_score,
    color = CellLine,
    shape = Dose,
    label = paste0(
      Dose,
      "_rep",
      Replicate
    )
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_point(
    size = 3.5,
    position = position_jitter(
      width = 0.05,
      height = 0
    )
  ) +
  
  ggrepel::geom_text_repel(
    size = 2.8,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  labs(
    
    title =
      "GSE66664 baseline predefined rescue score",
    
    subtitle =
      "0 h samples from the two dose-series experiments",
    
    x =
      NULL,
    
    y =
      "Predefined rescue signature score",
    
    color =
      "Cell line",
    
    shape =
      "DHT series"
  ) +
  
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    figures_dir,
    "GSE66664_baseline_rescue_score.png"
  ),
  p_baseline_score,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "GSE66664_baseline_rescue_score.pdf"
  ),
  p_baseline_score,
  width = 7.5,
  height = 6.5
)


############################################################
# 43. Optional same-model replication against GSE66663
#
# This is NOT independent patient validation.
# It is only a consistency check within the BAB/BAN model.
############################################################

gse66663_candidates <- list.files(
  path = "results",
  pattern = "^GSE66663_BAN_vs_BAB_all_genes\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

same_model_replication <- NULL

if (
  length(
    gse66663_candidates
  ) > 0
) {
  
  gse66663_de <- read.csv(
    gse66663_candidates[1],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (
    all(
      c(
        "GeneSymbol",
        "logFC"
      ) %in%
      colnames(
        gse66663_de
      )
    )
  ) {
    
    gse66664_small <- baseline_de[
      ,
      c(
        "GeneSymbol",
        "logFC"
      ),
      drop = FALSE
    ]
    
    colnames(
      gse66664_small
    )[2] <-
      "GSE66664_Baseline_logFC"
    
    gse66663_small <- gse66663_de[
      ,
      c(
        "GeneSymbol",
        "logFC"
      ),
      drop = FALSE
    ]
    
    colnames(
      gse66663_small
    )[2] <-
      "GSE66663_logFC"
    
    same_model_merged <- merge(
      gse66663_small,
      gse66664_small,
      by = "GeneSymbol"
    )
    
    same_model_valid <- is.finite(
      same_model_merged$GSE66663_logFC
    ) &
      is.finite(
        same_model_merged$GSE66664_Baseline_logFC
      )
    
    same_model_replication <- data.frame(
      
      NCommonGenes =
        sum(
          same_model_valid
        ),
      
      SpearmanRho =
        suppressWarnings(
          cor(
            same_model_merged$GSE66663_logFC[
              same_model_valid
            ],
            same_model_merged$GSE66664_Baseline_logFC[
              same_model_valid
            ],
            method = "spearman"
          )
        ),
      
      DirectionalConcordance =
        mean(
          sign(
            same_model_merged$GSE66663_logFC[
              same_model_valid
            ]
          ) ==
            sign(
              same_model_merged$GSE66664_Baseline_logFC[
                same_model_valid
              ]
            ),
          na.rm = TRUE
        ),
      
      stringsAsFactors =
        FALSE
    )
    
    write.csv(
      same_model_replication,
      file.path(
        results_dir,
        "GSE66663_vs_GSE66664_baseline_same_model_replication.csv"
      ),
      row.names = FALSE,
      quote = FALSE
    )
  }
}


############################################################
# 44. Extract key baseline concordance rows
############################################################

baseline_concordance <- concordance_summary[
  concordance_summary$ContrastType %in% c(
    "BaselineCrossSectional",
    "BaselineCrossSectionalPooled"
  ),
  ,
  drop = FALSE
]

write.csv(
  baseline_concordance,
  file.path(
    results_dir,
    "GSE66664_baseline_cross_study_concordance.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 45. Significant CAMERA summary
############################################################

camera_significant <- camera_results[
  camera_results$GlobalFDR < 0.05,
  ,
  drop = FALSE
]

write.csv(
  camera_significant,
  file.path(
    results_dir,
    "GSE66664_CAMERA_significant_globalFDR_0.05.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


############################################################
# 46. Human-readable summary
############################################################

summary_file <- file.path(
  results_dir,
  "ANALYSIS_SUMMARY_GSE66664.txt"
)

sink(
  summary_file
)

cat(
  "GSE66664 DHT time-course validation\n"
)

cat(
  "========================================\n\n"
)

cat(
  "Model:\n"
)

cat(
  "BAB = immortalized balding dermal papilla cell line\n"
)

cat(
  "BAN = immortalized non-balding dermal papilla cell line\n"
)

cat(
  "Arrays are experimental replicates/time-course measurements,\n"
)

cat(
  "NOT independent patients.\n\n"
)

cat(
  "Total arrays:",
  nrow(
    meta
  ),
  "\n"
)

cat(
  "Experimental cells:",
  nrow(
    replicate_count
  ),
  "\n"
)

cat(
  "Time points:",
  paste(
    time_levels,
    collapse = ", "
  ),
  "\n"
)

cat(
  "DHT series: 1nM, 10nM\n\n"
)

cat(
  "Cells with fewer than 3 arrays:\n"
)

if (
  nrow(
    incomplete_cells
  ) == 0
) {
  
  cat(
    "None\n"
  )
  
} else {
  
  print(
    incomplete_cells
  )
}

cat("\n")

cat(
  "Supplied expression probes:",
  nrow(
    raw_signal
  ),
  "\n"
)

cat(
  "Selected annotation:",
  selected_annotation_file,
  "\n"
)

cat(
  "Best annotation probe matches:",
  annotation_match_table$MatchingProbeIDs[1],
  "/",
  nrow(
    raw_signal
  ),
  "\n"
)

cat(
  "Mapped gene-level features:",
  nrow(
    expr_gene
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
  "Signature genes represented:",
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
  "Restored-upward overlap:",
  length(
    up_genes
  ),
  "\n"
)

cat(
  "Restored-downward overlap:",
  length(
    down_genes
  ),
  "\n\n"
)

cat(
  "Baseline rescue score summary:\n"
)

print(
  baseline_score_summary
)

cat("\n")

cat(
  "Baseline cross-study concordance with GSE44765 rescue effect:\n"
)

print(
  baseline_concordance
)

cat("\n")

cat(
  "Score-level differential DHT response:\n"
)

print(
  interaction_score
)

cat("\n")

cat(
  "DHT-response concordance with GSE44765 rescue effect:\n"
)

print(
  dht_concordance
)

cat("\n")

cat(
  "Interaction concordance:\n"
)

print(
  interaction_concordance
)

cat("\n")

cat(
  "Cross-sectional BAN-vs-BAB concordance across time:\n"
)

print(
  cross_sectional_concordance
)

cat("\n")

cat(
  "CAMERA tests significant after global BH FDR < 0.05:",
  nrow(
    camera_significant
  ),
  "/",
  nrow(
    camera_results
  ),
  "\n\n"
)

if (
  !is.null(
    same_model_replication
  )
) {
  
  cat(
    "Optional GSE66663 vs GSE66664 0h same-model replication:\n"
  )
  
  print(
    same_model_replication
  )
  
  cat("\n")
}

cat(
  "Interpretation guide:\n"
)

cat(
  "- Baseline BAN-BAB positive concordance supports replication of the\n"
)

cat(
  "  non-balding/rescue-like relationship seen in GSE66663.\n"
)

cat(
  "- For DHTResponse contrasts, positive rho means DHT moves that cell line\n"
)

cat(
  "  toward the GSE44765 rescue direction; negative rho means away from it.\n"
)

cat(
  "- For Interaction contrasts, positive rho means the BAN DHT response is\n"
)

cat(
  "  more rescue-like than the BAB response relative to their own baselines.\n"
)

cat(
  "- These are cell-line-model dynamics, not patient-level clinical effects.\n"
)

sink()


############################################################
# 47. Save workspace
############################################################

save(
  
  raw_table,
  meta,
  replicate_count,
  incomplete_cells,
  raw_signal,
  detection_matrix,
  expr_probe,
  annotation_match_table,
  annotation,
  expr_gene,
  pca,
  signature,
  signature_overlap,
  score_df,
  score_summary,
  delta_summary,
  interaction_score,
  design,
  contrast_matrix,
  contrast_meta,
  fit2,
  signature_effects_long,
  concordance_summary,
  camera_results,
  baseline_de,
  baseline_score_summary,
  same_model_replication,
  
  file = file.path(
    results_dir,
    "GSE66664_DHT_timecourse_workspace.RData"
  )
)


############################################################
# 48. Session information
############################################################

sink(
  file.path(
    results_dir,
    "sessionInfo_GSE66664.txt"
  )
)

sessionInfo()

sink()


############################################################
# 49. Final console report
############################################################

cat("\n")
cat("============================================================\n")
cat(" GSE66664 DHT TIME-COURSE VALIDATION COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "Arrays analyzed:",
  nrow(
    meta
  ),
  "\n"
)

cat(
  "Expression probes:",
  nrow(
    raw_signal
  ),
  "\n"
)

cat(
  "Gene-level features:",
  nrow(
    expr_gene
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

cat("\nIncomplete replicate cells:\n")
print(
  incomplete_cells
)

cat("\nBaseline score summary:\n")
print(
  baseline_score_summary
)

cat("\nBaseline GSE44765 concordance:\n")
print(
  baseline_concordance
)

cat("\nScore interaction trajectory:\n")
print(
  interaction_score
)

cat("\nMain figures:\n")

cat(
  file.path(
    figures_dir,
    "GSE66664_PCA.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66664_baseline_rescue_score.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66664_rescue_score_timecourse.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66664_rescue_score_interaction.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66664_DHT_response_rescue_concordance_timecourse.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66664_interaction_rescue_concordance_timecourse.png"
  ),
  "\n"
)

cat(
  file.path(
    figures_dir,
    "GSE66664_BAN_vs_BAB_rescue_concordance_timecourse.png"
  ),
  "\n"
)

cat("\nSummary file:\n")
cat(
  summary_file,
  "\n"
)

cat("\nIMPORTANT:\n")
cat(
  "Do not treat the 140 arrays as 140 independent patients.\n"
)

cat("============================================================\n")
