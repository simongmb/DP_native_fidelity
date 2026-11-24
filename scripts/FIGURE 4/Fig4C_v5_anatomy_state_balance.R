# ============================================================
# Figure 4C v5
# Ordered sample view across the two frozen rescue programmes
#
# Input:
#   subexpr = 1553 genes x 110 samples
#
# Frozen programmes:
#   Native-aligned rescue      = 81 genes
#   Alternative-state rescue    = 99 genes
#
# Score construction:
#   1. gene-wise z-score across samples
#   2. align by restored-up (+1) / restored-down (-1)
#   3. average aligned z-score within each frozen program
#
# Sample ordering:
#   Balance index = Native score - Alternative score
#
# Plot:
#   3-row ordered heatmap
#     Native-aligned programme score
#     Alternative-state programme score
#     Balance index
#
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "tidyr",
  "ggplot2",
  "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "\nMissing packages:\n",
      paste(missing_pkgs, collapse = ", "),
      "\n\nInstall with:\n",
      "install.packages(c(",
      paste0('"', missing_pkgs, '"', collapse = ", "),
      "))\n"
    )
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scales)
})


# ------------------------------------------------------------
# 1. Locate project root
# ------------------------------------------------------------

get_script_dir <- function() {

  args <- commandArgs(
    trailingOnly = FALSE
  )

  file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
  )

  if (length(file_arg) > 0) {
    return(
      dirname(
        normalizePath(
          sub(
            "^--file=",
            "",
            file_arg[1]
          )
        )
      )
    )
  }

  getwd()
}


start_dir <- get_script_dir()

candidate_roots <- unique(
  c(
    start_dir,
    getwd(),
    dirname(start_dir),
    dirname(dirname(start_dir)),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )
)

is_valid_root <- function(x) {

  frozen_dir <- file.path(
    x,
    "results",
    "FROZEN_PROGRAMS"
  )

  frozen_ok <- all(
    file.exists(
      c(
        file.path(
          frozen_dir,
          "Native_identity_restoration_81_FROZEN.csv"
        ),
        file.path(
          frozen_dir,
          "Alternative_state_rescue_99_FROZEN.csv"
        )
      )
    )
  )

  audit_ok <- length(
    list.files(
      x,
      pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
  ) > 0

  frozen_ok && audit_ok
}

valid_roots <- candidate_roots[
  vapply(
    candidate_roots,
    is_valid_root,
    logical(1)
  )
]

if (length(valid_roots) == 0) {
  stop("\nCannot identify project root.\n")
}

project_root <- valid_roots[1]

cat(
  "\nProject root:\n",
  project_root,
  "\n"
)


# ------------------------------------------------------------
# 2. Load clean-session subexpr
# ------------------------------------------------------------

loader <- file.path(
  project_root,
  "scripts",
  "load_subexpr_clean.R"
)

if (!file.exists(loader)) {
  stop(
    paste0(
      "
Cannot find clean-session loader:
",
      loader,
      "
"
    )
  )
}

source(loader, local = FALSE)

if (!exists("subexpr", inherits = TRUE)) {
  stop(
    "
QC STOP: load_subexpr_clean.R did not create subexpr.
"
  )
}

expr_mat <- get(
  "subexpr",
  inherits = TRUE
)

if (
  !is.matrix(expr_mat) &&
  !is.data.frame(expr_mat)
) {
  stop(
    "
subexpr must be a matrix/data.frame.
"
  )
}

expr_mat <- as.matrix(
  expr_mat
)

storage.mode(expr_mat) <- "numeric"

if (
  is.null(rownames(expr_mat)) ||
  is.null(colnames(expr_mat))
) {
  stop(
    "
subexpr must have gene rownames and sample colnames.
"
  )
}

rownames(expr_mat) <- toupper(
  stringr::str_trim(
    rownames(expr_mat)
  )
)

if (nrow(expr_mat) != 1553L) {
  stop(
    paste0(
      "
QC STOP: expected 1553 genes in clean-session subexpr, found ",
      nrow(expr_mat),
      ".
"
    )
  )
}

if (ncol(expr_mat) != 110L) {
  stop(
    paste0(
      "
QC STOP: expected 110 profiles in clean-session subexpr, found ",
      ncol(expr_mat),
      ".
"
    )
  )
}

sample_ids <- colnames(
  expr_mat
)

cat(
  "
Clean-session subexpr dimensions:",
  nrow(expr_mat),
  "x",
  ncol(expr_mat),
  "
"
)


# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------

clean_gene <- function(x) {

  x <- toupper(
    stringr::str_trim(
      as.character(x)
    )
  )

  x[
    x == ""
  ] <- NA_character_

  unique(
    x[
      !is.na(x)
    ]
  )
}


read_frozen <- function(path) {

  x <- readr::read_csv(
    path,
    show_col_types = FALSE
  )

  if (
    !"GeneSymbol" %in%
      names(x)
  ) {
    stop(
      paste0(
        "\nGeneSymbol missing in:\n",
        path,
        "\n"
      )
    )
  }

  clean_gene(
    x$GeneSymbol
  )
}


pick_col <- function(
  df,
  candidates
) {

  hit <- candidates[
    candidates %in%
      names(df)
  ]

  if (
    length(hit) == 0
  ) {
    return(
      NA_character_
    )
  }

  hit[1]
}


# ------------------------------------------------------------
# 4. Frozen programs
# ------------------------------------------------------------

frozen_dir <- file.path(
  project_root,
  "results",
  "FROZEN_PROGRAMS"
)

native81 <- read_frozen(
  file.path(
    frozen_dir,
    "Native_identity_restoration_81_FROZEN.csv"
  )
)

alternative99 <- read_frozen(
  file.path(
    frozen_dir,
    "Alternative_state_rescue_99_FROZEN.csv"
  )
)

if (
  length(native81) != 81
) {
  stop(
    "\nNative frozen program != 81 genes.\n"
  )
}

if (
  length(alternative99) != 99
) {
  stop(
    "\nAlternative frozen program != 99 genes.\n"
  )
}

if (
  !all(
    native81 %in%
      rownames(expr_mat)
  )
) {
  stop(
    "\nsubexpr does not contain all 81 native genes.\n"
  )
}

if (
  !all(
    alternative99 %in%
      rownames(expr_mat)
  )
) {
  stop(
    "\nsubexpr does not contain all 99 alternative genes.\n"
  )
}


# ------------------------------------------------------------
# 5. Audit table + restoration direction
# ------------------------------------------------------------

audit_files <- list.files(
  project_root,
  pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (
  length(audit_files) == 0
) {
  stop(
    "\nCannot find 02_gene_level_native_fidelity_audit.csv\n"
  )
}

if (
  length(audit_files) > 1
) {

  in_results <- grepl(
    "results",
    audit_files,
    ignore.case = TRUE
  )

  if (
    any(
      in_results
    )
  ) {
    audit_files <- audit_files[
      in_results
    ]
  }
}

audit_file <- audit_files[1]

audit0 <- readr::read_csv(
  audit_file,
  show_col_types = FALSE
)

gene_col <- pick_col(
  audit0,
  c(
    "GeneSymbol",
    "gene_symbol",
    "Gene",
    "gene",
    "symbol",
    "SYMBOL"
  )
)

if (
  is.na(
    gene_col
  )
) {
  stop(
    "\nCannot identify gene column in audit table.\n"
  )
}

direction_candidates <- c(
  "restoration_direction",
  "Restoration_direction",
  "restored_direction",
  "Restored_direction",
  "rescue_direction",
  "Rescue_direction",
  "strict_rescue_direction",
  "Strict_rescue_direction",
  "direction_restored",
  "restored_up_down",
  "Restored_up_down",
  "rescue_up_down",
  "direction"
)

up_flag_candidates <- c(
  "restored_up",
  "Restored_up",
  "is_restored_up",
  "is_restored_upward",
  "restored_upward"
)

signed_effect_candidates <- c(
  "rescue_logFC",
  "Rescue_logFC",
  "rescue_effect",
  "Rescue_effect",
  "restoration_logFC",
  "Restoration_logFC",
  "treatment_logFC",
  "Treatment_logFC",
  "perturbation_logFC"
)

direction_col <- pick_col(
  audit0,
  direction_candidates
)

up_flag_col <- pick_col(
  audit0,
  up_flag_candidates
)

signed_effect_col <- pick_col(
  audit0,
  signed_effect_candidates
)

direction_sign <- rep(
  NA_real_,
  nrow(audit0)
)

direction_source <- NA_character_


if (
  !is.na(
    direction_col
  )
) {

  raw_dir <- tolower(
    stringr::str_trim(
      as.character(
        audit0[[
          direction_col
        ]]
      )
    )
  )

  direction_sign[
    grepl(
      "up|upward|increase|positive|\\+1",
      raw_dir
    )
  ] <- 1

  direction_sign[
    grepl(
      "down|downward|decrease|negative|-1",
      raw_dir
    )
  ] <- -1

  direction_source <- paste0(
    "categorical column: ",
    direction_col
  )
}


if (
  all(
    is.na(
      direction_sign
    )
  ) &&
  !is.na(
    up_flag_col
  )
) {

  flag <- audit0[[
    up_flag_col
  ]]

  if (
    is.logical(
      flag
    )
  ) {

    direction_sign <- ifelse(
      flag,
      1,
      -1
    )

    direction_source <- paste0(
      "logical column: ",
      up_flag_col
    )
  }
}


if (
  all(
    is.na(
      direction_sign
    )
  ) &&
  !is.na(
    signed_effect_col
  )
) {

  effect <- suppressWarnings(
    as.numeric(
      audit0[[
        signed_effect_col
      ]]
    )
  )

  direction_sign <- sign(
    effect
  )

  direction_sign[
    !is.finite(
      direction_sign
    ) |
    direction_sign == 0
  ] <- NA_real_

  direction_source <- paste0(
    "signed numeric column: ",
    signed_effect_col
  )
}


if (
  all(
    is.na(
      direction_sign
    )
  )
) {
  stop(
    "\nCannot recover restored-up/restored-down direction.\n"
  )
}

direction_tbl <- tibble::tibble(

  GeneSymbol =
    toupper(
      stringr::str_trim(
        as.character(
          audit0[[
            gene_col
          ]]
        )
      )
    ),

  Direction_sign =
    direction_sign
) %>%

  dplyr::filter(
    !is.na(
      GeneSymbol
    ),
    is.finite(
      Direction_sign
    ),
    Direction_sign %in%
      c(
        -1,
        1
      )
  ) %>%

  dplyr::distinct(
    GeneSymbol,
    .keep_all = TRUE
  )

native_direction <- tibble::tibble(
  GeneSymbol =
    native81
) %>%
  dplyr::left_join(
    direction_tbl,
    by = "GeneSymbol"
  )

alt_direction <- tibble::tibble(
  GeneSymbol =
    alternative99
) %>%
  dplyr::left_join(
    direction_tbl,
    by = "GeneSymbol"
  )

if (
  any(
    !is.finite(
      native_direction$Direction_sign
    )
  ) ||
  any(
    !is.finite(
      alt_direction$Direction_sign
    )
  )
) {
  stop(
    "\nDirection missing for frozen genes.\n"
  )
}

native_up <- sum(
  native_direction$Direction_sign == 1
)

native_down <- sum(
  native_direction$Direction_sign == -1
)

alt_up <- sum(
  alt_direction$Direction_sign == 1
)

alt_down <- sum(
  alt_direction$Direction_sign == -1
)

if (
  native_up != 55 ||
  native_down != 26 ||
  alt_up != 30 ||
  alt_down != 69
) {
  stop(
    "\nDirection QC does not reproduce frozen Figure 3D counts.\n"
  )
}


# ------------------------------------------------------------
# 6. Gene-wise z scoring
# ------------------------------------------------------------

z_by_gene <- function(mat) {

  z <- t(
    apply(
      mat,
      1,
      function(v) {

        s <- stats::sd(
          v,
          na.rm = TRUE
        )

        if (
          !is.finite(
            s
          ) ||
          s == 0
        ) {

          return(
            rep(
              0,
              length(v)
            )
          )
        }

        (
          v -
          mean(
            v,
            na.rm = TRUE
          )
        ) / s
      }
    )
  )

  rownames(z) <- rownames(
    mat
  )

  colnames(z) <- colnames(
    mat
  )

  z
}


native_z <- z_by_gene(
  expr_mat[
    native81,
    ,
    drop = FALSE
  ]
)

alt_z <- z_by_gene(
  expr_mat[
    alternative99,
    ,
    drop = FALSE
  ]
)


# ------------------------------------------------------------
# 7. Direction alignment
# ------------------------------------------------------------

native_sign <- native_direction$Direction_sign[
  match(
    native81,
    native_direction$GeneSymbol
  )
]

alt_sign <- alt_direction$Direction_sign[
  match(
    alternative99,
    alt_direction$GeneSymbol
  )
]

native_aligned <- sweep(
  native_z,
  1,
  native_sign,
  "*"
)

alt_aligned <- sweep(
  alt_z,
  1,
  alt_sign,
  "*"
)


# ------------------------------------------------------------
# 8. Sample scores
# ------------------------------------------------------------

native_score <- colMeans(
  native_aligned,
  na.rm = TRUE
)

alternative_score <- colMeans(
  alt_aligned,
  na.rm = TRUE
)

sample_ids <- colnames(
  expr_mat
)

score_df <- tibble::tibble(

  Sample =
    sample_ids,

  Native_score =
    as.numeric(
      native_score[
        sample_ids
      ]
    ),

  Alternative_score =
    as.numeric(
      alternative_score[
        sample_ids
      ]
    )
) %>%

  dplyr::mutate(

    Balance_index =
      Native_score -
      Alternative_score
  )


# ------------------------------------------------------------
# 9. Order samples by anatomy-referenced state balance
# ------------------------------------------------------------

score_df <- score_df %>%

  dplyr::arrange(
    Balance_index
  ) %>%

  dplyr::mutate(

    Sample_order =
      dplyr::row_number(),

    Sample_factor =
      factor(
        Sample,
        levels =
          Sample
      )
  )


# ------------------------------------------------------------
# 10. Standardize the 3 displayed rows for common color scale
#
# Each row is z-scaled across samples.
# This is ONLY for visualization.
# The original values remain exported separately.
# ------------------------------------------------------------

row_z <- function(x) {

  s <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(
      s
    ) ||
    s == 0
  ) {

    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (
    x -
    mean(
      x,
      na.rm = TRUE
    )
  ) / s
}


heat_df <- score_df %>%

  dplyr::transmute(

    Sample =
      Sample,

    Sample_order =
      Sample_order,

    `Native-aligned programme` =
      row_z(
        Native_score
      ),

    `Alternative-state programme` =
      row_z(
        Alternative_score
      ),

    `Native-aligned - Alternative-state balance` =
      row_z(
        Balance_index
      )
  ) %>%

  tidyr::pivot_longer(

    cols =
      c(
        `Native-aligned programme`,
        `Alternative-state programme`,
        `Native-aligned - Alternative-state balance`
      ),

    names_to =
      "Metric",

    values_to =
      "Z"
  )


heat_df$Metric <- factor(

  heat_df$Metric,

  levels = c(
    "Native-aligned programme",
    "Alternative-state programme",
    "Native-aligned - Alternative-state balance"
  )
)


# ------------------------------------------------------------
# 11. Optional sample-group strip
#
# If affected_gsm / aff_cols are available,
# export membership but do not encode it in the main heatmap
# unless the identity is unambiguous.
# ------------------------------------------------------------

affected_ids <- character(0)

if (
  exists(
    "affected_gsm",
    envir = .GlobalEnv
  )
) {

  ag <- get(
    "affected_gsm",
    envir = .GlobalEnv
  )

  if (
    is.atomic(
      ag
    ) &&
    is.null(
      dim(
        ag
      )
    )
  ) {

    affected_ids <- union(
      affected_ids,
      as.character(
        ag
      )
    )
  }
}

if (
  exists(
    "aff_cols",
    envir = .GlobalEnv
  )
) {

  ac <- get(
    "aff_cols",
    envir = .GlobalEnv
  )

  if (
    is.atomic(
      ac
    ) &&
    is.null(
      dim(
        ac
      )
    )
  ) {

    affected_ids <- union(
      affected_ids,
      as.character(
        ac
      )
    )
  }
}

score_df <- score_df %>%

  dplyr::mutate(

    Affected_set =
      Sample %in%
      affected_ids
  )


# ------------------------------------------------------------
# 12. Correlations for reporting
# ------------------------------------------------------------

cor_native_alt <- suppressWarnings(
  stats::cor.test(
    score_df$Native_score,
    score_df$Alternative_score,
    method = "spearman",
    exact = FALSE
  )
)

rho_native_alt <- unname(
  cor_native_alt$estimate
)

p_native_alt <- cor_native_alt$p.value


# ------------------------------------------------------------
# 13. Plot
# ------------------------------------------------------------

COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"

max_abs_z <- max(
  abs(
    heat_df$Z
  ),
  na.rm = TRUE
)

fill_limit <- min(
  max(
    2,
    max_abs_z
  ),
  3
)


p <- ggplot(
  heat_df,
  aes(
    x = Sample_order,
    y = Metric,
    fill = Z
  )
) +

  geom_tile(
    height = 0.78,
    width = 1.0
  ) +

  scale_fill_gradient2(
    low = "#5A7F91",
    mid = "#F5F5F2",
    high = "#C97857",
    midpoint = 0,
    limits = c(
      -fill_limit,
      fill_limit
    ),
    oob = scales::squish,
    name = "Row z-score"
  ) +

  scale_x_continuous(
    expand = c(
      0,
      0
    ),
    breaks = NULL
  ) +

  labs(

    x =
      "Samples ordered by Native-aligned - Alternative-state balance",

    y =
      NULL,

    title =
      "Graded anatomy-referenced state balance across samples",

    subtitle =
      paste0(
        "Ordered by state balance; ",
        "programme correlation rho = ",
        sprintf(
          "%.2f",
          rho_native_alt
        )
      ),

    caption =
      paste0(
        "n = ",
        nrow(
          score_df
        ),
        " samples. ",
        "Rows are z-standardized independently for visualization."
      )
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(

    plot.title =
      element_text(
        size = 10.4,
        face = "bold",
        colour = COL_TEXT,
        margin = margin(
          b = 4
        )
      ),

    plot.subtitle =
      element_text(
        size = 8.0,
        colour = COL_MUTED,
        margin = margin(
          b = 8
        )
      ),

    plot.caption =
      element_text(
        size = 7.2,
        colour = "#95A4AA",
        hjust = 0,
        margin = margin(
          t = 7
        )
      ),

    axis.title.x =
      element_text(
        size = 9.5,
        colour = COL_TEXT,
        margin = margin(
          t = 7
        )
      ),

    axis.text.y =
      element_text(
        size = 9.0,
        colour = COL_TEXT
      ),

    axis.line =
      element_blank(),

    axis.ticks =
      element_blank(),

    legend.position =
      "right",

    legend.title =
      element_text(
        size = 8.3,
        colour = COL_TEXT
      ),

    legend.text =
      element_text(
        size = 7.8,
        colour = COL_MUTED
      ),

    plot.margin =
      margin(
        8,
        10,
        8,
        10,
        unit = "pt"
      )
  )


# ------------------------------------------------------------
# 14. Output
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 4",
  "Panel_C_orderedStateBalance_v5"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  filename = file.path(
    out_dir,
    "Fig4_PanelC_orderedStateBalance_v5.png"
  ),
  plot = p,
  width = 7.0,
  height = 3.6,
  units = "in",
  dpi = 600,
  bg = "white"
)


pdf(
  file = file.path(
    out_dir,
    "Fig4_PanelC_orderedStateBalance_v5.pdf"
  ),
  width = 7.0,
  height = 3.6,
  useDingbats = FALSE
)

print(
  p
)

dev.off()


ggsave(
  filename = file.path(
    out_dir,
    "Fig4_PanelC_orderedStateBalance_v5.tiff"
  ),
  plot = p,
  width = 7.0,
  height = 3.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 15. Export data
# ------------------------------------------------------------

readr::write_csv(
  score_df,
  file.path(
    out_dir,
    "Fig4C_ordered_sample_scores.csv"
  )
)

readr::write_csv(
  heat_df,
  file.path(
    out_dir,
    "Fig4C_heatmap_plot_data.csv"
  )
)


stats_df <- tibble::tibble(

  N_samples =
    nrow(
      score_df
    ),

  Spearman_rho_native_vs_alternative =
    rho_native_alt,

  Spearman_P =
    p_native_alt,

  Native_program_n =
    81,

  Alternative_program_n =
    99,

  Native_restored_up =
    native_up,

  Native_restored_down =
    native_down,

  Alternative_restored_up =
    alt_up,

  Alternative_restored_down =
    alt_down,

  Affected_set_samples_detected =
    sum(
      score_df$Affected_set
    )
)

readr::write_csv(
  stats_df,
  file.path(
    out_dir,
    "Fig4C_statistics.csv"
  )
)


qc_df <- tibble::tibble(

  Item = c(
    "subexpr genes",
    "subexpr samples",
    "Native genes represented",
    "Alternative genes represented",
    "Native restored-up",
    "Native restored-down",
    "Alternative restored-up",
    "Alternative restored-down"
  ),

  Value = c(
    nrow(
      expr_mat
    ),
    ncol(
      expr_mat
    ),
    length(
      native81
    ),
    length(
      alternative99
    ),
    native_up,
    native_down,
    alt_up,
    alt_down
  ),

  Expected = c(
    "1553",
    "110",
    "81",
    "99",
    "55",
    "26",
    "30",
    "69"
  )
)

readr::write_csv(
  qc_df,
  file.path(
    out_dir,
    "Fig4C_QC.csv"
  )
)


readme <- c(

  "Figure 4C ordered anatomy-referenced state balance v5",

  "",

  "Input:",
  paste0(
    "subexpr = ",
    nrow(
      expr_mat
    ),
    " genes x ",
    ncol(
      expr_mat
    ),
    " samples"
  ),

  "",

  "Frozen programmes:",
  "Native-aligned rescue = 81 genes",
  "Alternative-state rescue = 99 genes",

  "",

  "Sample order:",
  "Balance index = Native score - Alternative score",

  "",

  "Displayed rows:",
  "1. Native-aligned programme score",
  "2. Alternative-state programme score",
  "3. Native-aligned - Alternative-state balance",

  "",

  "For visualization only, each displayed row is z-standardized across samples.",

  "",

  "Display refinement in v5: shortened title/subtitle to prevent clipping; no numerical or analytical changes.",

  "",

  paste0(
    "Spearman rho between Native and Alternative scores = ",
    sprintf(
      "%.4f",
      rho_native_alt
    )
  ),

  paste0(
    "Direction source: ",
    direction_source
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig4C_orderedStateBalance_v5.txt"
  )
)


cat(
  "\n============================================\n",
  "FIGURE 4C COMPLETE\n",
  "============================================\n\n",
  sep = ""
)

print(
  stats_df,
  n = Inf
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)
