# ============================================================
# Figure 5A v6
# Non-programme gene associations with the anatomy-referenced state balance
#
# PURPOSE
# ------------------------------------------------------------
# Test whether the Native-aligned - Alternative-state balance
# extends beyond the 180 genes used to define the two frozen
# programmes.
#
# Key design:
#   - Reuse the frozen Figure 4A sample-level balance defined from:
#       81 Native-aligned rescue genes
#       99 Alternative-state rescue genes
#
#   - Then test ONLY non-programme genes from the frozen
#     1,553-gene analyzable background.
#
# Therefore, every tested gene is excluded from the defining programmes
# and there is no direct self-inclusion in the state-balance target.
#
# Expected held-out background:
#   1553 - 81 - 99 = 1373 genes
#
# For each non-programme gene:
#   Spearman correlation(expression, Native - Alternative balance)
#   BH FDR across all non-programme genes
#
# Main figure:
#   rho vs -log10(FDR), with strongest positive/negative
#   non-programme genes labelled.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "ggplot2",
  "ggrepel"
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
  library(ggplot2)
  library(ggrepel)
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

  all(
    file.exists(
      c(
        file.path(
          frozen_dir,
          "Analyzable_background_1553_FROZEN.csv"
        ),
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
  ) &&
    length(
      list.files(
        x,
        pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    ) > 0
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

  x[x == ""] <- NA_character_

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

  if (!"GeneSymbol" %in% names(x)) {
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

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


# ------------------------------------------------------------
# 4. Frozen programmes + background
# ------------------------------------------------------------

frozen_dir <- file.path(
  project_root,
  "results",
  "FROZEN_PROGRAMS"
)

background1553 <- read_frozen(
  file.path(
    frozen_dir,
    "Analyzable_background_1553_FROZEN.csv"
  )
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

if (length(background1553) != 1553) {
  stop(
    "\nFrozen background != 1553 genes.\n"
  )
}

if (length(native81) != 81) {
  stop(
    "\nNative frozen program != 81 genes.\n"
  )
}

if (length(alternative99) != 99) {
  stop(
    "\nAlternative frozen program != 99 genes.\n"
  )
}

if (
  length(
    intersect(
      native81,
      alternative99
    )
  ) != 0
) {
  stop(
    "\nFrozen 81/99 programs overlap unexpectedly.\n"
  )
}

if (
  !all(
    native81 %in%
      rownames(expr_mat)
  ) ||
  !all(
    alternative99 %in%
      rownames(expr_mat)
  )
) {
  stop(
    "\nsubexpr does not contain all frozen program genes.\n"
  )
}


# ------------------------------------------------------------
# 5. Non-programme background
# ------------------------------------------------------------

program180 <- union(
  native81,
  alternative99
)

heldout_genes <- setdiff(
  background1553,
  program180
)

heldout_present <- intersect(
  heldout_genes,
  rownames(expr_mat)
)

cat(
  "\nHeld-out background:\n",
  "Expected = 1373\n",
  "Defined = ",
  length(heldout_genes),
  "\n",
  "Present in subexpr = ",
  length(heldout_present),
  "\n",
  sep = ""
)

if (length(heldout_genes) != 1373) {
  stop(
    "\nHeld-out background is not exactly 1,373 genes.\n"
  )
}

if (length(heldout_present) < 0.95 * 1373) {
  stop(
    "\nToo many held-out background genes are missing from subexpr.\n"
  )
}


# ------------------------------------------------------------
# 6. Audit + restoration direction
# ------------------------------------------------------------

audit_files <- list.files(
  project_root,
  pattern = "^02_gene_level_native_fidelity_audit\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(audit_files) == 0) {
  stop(
    "\nCannot find 02_gene_level_native_fidelity_audit.csv\n"
  )
}

if (length(audit_files) > 1) {

  in_results <- grepl(
    "results",
    audit_files,
    ignore.case = TRUE
  )

  if (any(in_results)) {
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

if (is.na(gene_col)) {
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


if (!is.na(direction_col)) {

  raw_dir <- tolower(
    stringr::str_trim(
      as.character(
        audit0[[direction_col]]
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
  all(is.na(direction_sign)) &&
  !is.na(up_flag_col)
) {

  flag <- audit0[[up_flag_col]]

  if (is.logical(flag)) {

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
  all(is.na(direction_sign)) &&
  !is.na(signed_effect_col)
) {

  effect <- suppressWarnings(
    as.numeric(
      audit0[[signed_effect_col]]
    )
  )

  direction_sign <- sign(
    effect
  )

  direction_sign[
    !is.finite(direction_sign) |
    direction_sign == 0
  ] <- NA_real_

  direction_source <- paste0(
    "signed numeric column: ",
    signed_effect_col
  )
}

if (all(is.na(direction_sign))) {
  stop(
    "\nCannot recover restoration direction.\n"
  )
}

direction_tbl <- tibble::tibble(
  GeneSymbol =
    toupper(
      stringr::str_trim(
        as.character(
          audit0[[gene_col]]
        )
      )
    ),
  Direction_sign =
    direction_sign
) %>%
  dplyr::filter(
    !is.na(GeneSymbol),
    is.finite(Direction_sign),
    Direction_sign %in%
      c(-1, 1)
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
  any(!is.finite(native_direction$Direction_sign)) ||
  any(!is.finite(alt_direction$Direction_sign))
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
# 7. Recover the frozen Figure 4A state-balance target
# ------------------------------------------------------------
# Use the already-exported Figure 4A sample scores as the only
# target source. Do not silently redefine the state balance here.

fig4_score_files <- list.files(
  project_root,
  pattern = "^Fig4A_sample_program_scores\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(fig4_score_files) == 0) {
  stop(
    paste0(
      "\nQC STOP\n",
      "Cannot find Fig4A_sample_program_scores.csv.\n",
      "Figure 5A will NOT reconstruct the target independently.\n"
    )
  )
}

score_file_rank <- vapply(
  fig4_score_files,
  function(path) {
    s <- 0
    if (grepl("sampleProgrammes_Nature_wide_v7", path, ignore.case = TRUE)) s <- s + 100
    if (grepl("FIGURE 4", path, ignore.case = TRUE)) s <- s + 10
    s
  },
  numeric(1)
)

fig4_score_file <- fig4_score_files[which.max(score_file_rank)]

fig4_scores <- readr::read_csv(
  fig4_score_file,
  show_col_types = FALSE
)

required_score_cols <- c(
  "Sample",
  "Native_identity_score",
  "Alternative_state_score"
)

if (!all(required_score_cols %in% names(fig4_scores))) {
  stop(
    paste0(
      "\nQC STOP\n",
      "Fig4A_sample_program_scores.csv does not contain the required columns:\n",
      paste(required_score_cols, collapse = ", "),
      "\n"
    )
  )
}

if (!setequal(as.character(fig4_scores$Sample), sample_ids)) {
  stop(
    "\nQC STOP\nFigure 4A sample IDs do not match the current expression matrix.\n"
  )
}

fig4_scores <- fig4_scores[
  match(sample_ids, as.character(fig4_scores$Sample)),
  ,
  drop = FALSE
]

balance_score <- as.numeric(fig4_scores$Native_identity_score) -
  as.numeric(fig4_scores$Alternative_state_score)

names(balance_score) <- sample_ids

if (
  length(balance_score) != length(sample_ids) ||
  any(!is.finite(balance_score))
) {
  stop(
    "\nQC STOP\nInvalid frozen Figure 4A state-balance target.\n"
  )
}

target_source <- normalizePath(fig4_score_file)

cat(
  "\nFrozen Figure 4A target loaded from:\n",
  target_source,
  "\n"
)


# ------------------------------------------------------------
# 8. Non-programme gene associations
# ------------------------------------------------------------

association_rows <- lapply(
  heldout_present,
  function(g) {

    values <- as.numeric(
      expr_mat[
        g,
        sample_ids
      ]
    )

    ct <- suppressWarnings(
      stats::cor.test(
        values,
        balance_score,
        method = "spearman",
        exact = FALSE
      )
    )

    tibble::tibble(
      GeneSymbol =
        g,
      Spearman_rho =
        unname(
          ct$estimate
        ),
      P_value =
        ct$p.value
    )
  }
)

assoc_df <- dplyr::bind_rows(
  association_rows
) %>%

  dplyr::mutate(
    FDR =
      p.adjust(
        P_value,
        method = "BH"
      ),
    NegLog10FDR =
      -log10(
        pmax(
          FDR,
          1e-300
        )
      ),
    Direction =
      dplyr::case_when(
        Spearman_rho > 0 ~
          "Native-aligned associated",
        Spearman_rho < 0 ~
          "Alternative-aligned associated",
        TRUE ~
          "Neutral"
      )
  ) %>%

  dplyr::arrange(
    Spearman_rho
  )


# ------------------------------------------------------------
# 9. Summary / top labels
# ------------------------------------------------------------

n_tested <- nrow(
  assoc_df
)

n_fdr <- sum(
  assoc_df$FDR < 0.05,
  na.rm = TRUE
)

n_pos_fdr <- sum(
  assoc_df$FDR < 0.05 &
  assoc_df$Spearman_rho > 0,
  na.rm = TRUE
)

n_neg_fdr <- sum(
  assoc_df$FDR < 0.05 &
  assoc_df$Spearman_rho < 0,
  na.rm = TRUE
)

TOP_N <- 8

top_pos <- assoc_df %>%
  dplyr::filter(
    Spearman_rho > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      abs(
        Spearman_rho
      )
    )
  ) %>%
  dplyr::slice_head(
    n = TOP_N
  )

top_neg <- assoc_df %>%
  dplyr::filter(
    Spearman_rho < 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      abs(
        Spearman_rho
      )
    )
  ) %>%
  dplyr::slice_head(
    n = TOP_N
  )

label_df <- dplyr::bind_rows(
  top_pos,
  top_neg
) %>%
  dplyr::distinct(
    GeneSymbol,
    .keep_all = TRUE
  )


cat(
  "\nNon-programme gene association analysis:\n",
  "Genes tested = ",
  n_tested,
  "\n",
  "FDR < 0.05 = ",
  n_fdr,
  "\n",
  "Positive significant = ",
  n_pos_fdr,
  "\n",
  "Negative significant = ",
  n_neg_fdr,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Plot — Nature-style wide / short layout
# ------------------------------------------------------------

COL_POS <- "#3D8C82"
COL_NEG <- "#C97857"
COL_NS <- "#B9C2C6"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"
COL_ZERO <- "#A8B5BA"

assoc_df <- assoc_df %>%
  dplyr::mutate(
    Plot_class =
      dplyr::case_when(
        FDR < 0.05 &
        Spearman_rho > 0 ~
          "Native-aligned associated",
        FDR < 0.05 &
        Spearman_rho < 0 ~
          "Alternative-aligned associated",
        TRUE ~
          "Not significant"
      )
  )

plot_cols <- c(
  "Native-aligned associated" = COL_POS,
  "Alternative-aligned associated" = COL_NEG,
  "Not significant" = COL_NS
)

p <- ggplot(
  assoc_df,
  aes(
    x = Spearman_rho,
    y = NegLog10FDR,
    colour = Plot_class
  )
) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.28,
    colour = COL_ZERO
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "22",
    linewidth = 0.28,
    colour = COL_ZERO
  ) +

  geom_point(
    size = 1.05,
    alpha = 0.68
  ) +

  # Fixed, deterministic label lanes for the 8 strongest genes on each arm.
  # This avoids any remaining overlap from automatic repulsion.
  geom_segment(
    data = top_neg %>%
      dplyr::arrange(dplyr::desc(NegLog10FDR)) %>%
      dplyr::mutate(
        label_y = seq(38.5, 14.0, length.out = dplyr::n()),
        label_x = -0.93
      ),
    aes(
      x = Spearman_rho,
      y = NegLog10FDR,
      xend = label_x + 0.025,
      yend = label_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.18,
    colour = "#AEB9BE",
    alpha = 0.75
  ) +

  geom_text(
    data = top_neg %>%
      dplyr::arrange(dplyr::desc(NegLog10FDR)) %>%
      dplyr::mutate(
        label_y = seq(38.5, 14.0, length.out = dplyr::n()),
        label_x = -0.93
      ),
    aes(
      x = label_x,
      y = label_y,
      label = GeneSymbol
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 1.48,
    family = "Helvetica",
    colour = COL_TEXT
  ) +

  geom_segment(
    data = top_pos %>%
      dplyr::arrange(dplyr::desc(NegLog10FDR)) %>%
      dplyr::mutate(
        label_y = seq(38.5, 14.0, length.out = dplyr::n()),
        label_x = 0.93
      ),
    aes(
      x = Spearman_rho,
      y = NegLog10FDR,
      xend = label_x - 0.025,
      yend = label_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.18,
    colour = "#AEB9BE",
    alpha = 0.75
  ) +

  geom_text(
    data = top_pos %>%
      dplyr::arrange(dplyr::desc(NegLog10FDR)) %>%
      dplyr::mutate(
        label_y = seq(38.5, 14.0, length.out = dplyr::n()),
        label_x = 0.93
      ),
    aes(
      x = label_x,
      y = label_y,
      label = GeneSymbol
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 1.48,
    family = "Helvetica",
    colour = COL_TEXT
  ) +

  scale_colour_manual(
    values = plot_cols,
    name = NULL
  ) +

  scale_x_continuous(
    limits = c(-1.02, 1.02),
    breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    limits = c(0, 41),
    expand = c(0, 0)
  ) +

  labs(
    x = "Spearman rho with anatomy-referenced state balance",
    y = expression(-log[10]("BH FDR")),
    title = "Non-programme genes associate with the anatomy-referenced state balance",
    subtitle = paste0(
      n_fdr,
      " of ",
      n_tested,
      " non-programme genes reach BH FDR < 0.05"
    ),
    caption = paste0(
      "Only non-programme genes were tested; dashed line marks BH FDR = 0.05."
    )
  ) +

  theme_classic(
    base_family = "Helvetica",
    base_size = 7
  ) +

  theme(
    plot.title = element_text(
      size = 5.7,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      margin = margin(b = 1.0)
    ),

    plot.subtitle = element_text(
      size = 4.5,
      colour = COL_MUTED,
      hjust = 0.5,
      margin = margin(b = 2.5)
    ),

    plot.caption = element_text(
      size = 4.1,
      colour = "#95A4AA",
      hjust = 0.5,
      margin = margin(t = 2.5)
    ),

    axis.title.x = element_text(
      size = 5.0,
      colour = COL_TEXT,
      margin = margin(t = 3)
    ),

    axis.title.y = element_text(
      size = 5.0,
      colour = COL_TEXT,
      margin = margin(r = 3)
    ),

    axis.text = element_text(
      size = 4.6,
      colour = COL_MUTED
    ),

    axis.line = element_line(
      linewidth = 0.30,
      colour = COL_TEXT
    ),

    axis.ticks = element_line(
      linewidth = 0.26,
      colour = COL_TEXT
    ),

    axis.ticks.length = unit(
      1.0,
      "mm"
    ),

    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",

    legend.text = element_text(
      size = 4.5,
      colour = COL_MUTED
    ),

    legend.key.width = unit(
      3.0,
      "mm"
    ),

    legend.key.height = unit(
      2.3,
      "mm"
    ),

    legend.spacing.x = unit(
      1.6,
      "mm"
    ),

    legend.margin = margin(
      t = 0,
      r = 0,
      b = 1.5,
      l = 0
    ),

    plot.margin = margin(
      4,
      6,
      4,
      6,
      unit = "pt"
    )
  )


# ------------------------------------------------------------
# 11. Output
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 5",
  "Panel_A_nonProgrammeAssociation_v6"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    out_dir,
    "Fig5_PanelA_nonProgrammeAssociation_v6.png"
  ),
  plot = p,
  width = 88,
  height = 48,
  units = "mm",
  dpi = 600,
  bg = "white"
)

pdf(
  file = file.path(
    out_dir,
    "Fig5_PanelA_nonProgrammeAssociation_v6.pdf"
  ),
  width = 88 / 25.4,
  height = 48 / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)

print(p)
dev.off()

ggsave(
  filename = file.path(
    out_dir,
    "Fig5_PanelA_nonProgrammeAssociation_v6.tiff"
  ),
  plot = p,
  width = 88,
  height = 48,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 12. Export
# ------------------------------------------------------------

readr::write_csv(
  assoc_df,
  file.path(
    out_dir,
    "Fig5A_all_nonprogramme_gene_associations.csv"
  )
)

readr::write_csv(
  label_df,
  file.path(
    out_dir,
    "Fig5A_labeled_nonprogramme_genes.csv"
  )
)

stats_df <- tibble::tibble(
  Samples =
    length(
      sample_ids
    ),
  Frozen_background =
    length(
      background1553
    ),
  Program_genes_excluded =
    length(
      program180
    ),
  Heldout_genes_defined =
    length(
      heldout_genes
    ),
  Heldout_genes_tested =
    n_tested,
  FDR_lt_0_05 =
    n_fdr,
  Positive_FDR_lt_0_05 =
    n_pos_fdr,
  Negative_FDR_lt_0_05 =
    n_neg_fdr
)

readr::write_csv(
  stats_df,
  file.path(
    out_dir,
    "Fig5A_statistics.csv"
  )
)

qc_df <- tibble::tibble(
  Item = c(
    "subexpr genes",
    "subexpr samples",
    "Frozen background",
    "Native program",
    "Alternative program",
    "Program overlap",
    "Held-out background defined",
    "Held-out background tested",
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
      background1553
    ),
    length(
      native81
    ),
    length(
      alternative99
    ),
    length(
      intersect(
        native81,
        alternative99
      )
    ),
    length(
      heldout_genes
    ),
    n_tested,
    native_up,
    native_down,
    alt_up,
    alt_down
  ),
  Expected = c(
    "1553",
    "110",
    "1553",
    "81",
    "99",
    "0",
    "1373",
    ">=1304",
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
    "Fig5A_QC.csv"
  )
)

readme <- c(
  "Figure 5A — non-programme gene association analysis v6",
  "",
  "Purpose:",
  "Test whether the anatomy-referenced state balance is reflected in genes outside the 180 defining programme genes.",
  "",
  paste0(
    "Frozen background = ",
    length(
      background1553
    )
  ),
  paste0(
    "Programme genes excluded = ",
    length(
      program180
    )
  ),
  paste0(
    "Non-programme genes defined = ",
    length(
      heldout_genes
    )
  ),
  paste0(
    "Non-programme genes tested = ",
    n_tested
  ),
  paste0(
    "BH FDR < 0.05 = ",
    n_fdr
  ),
  paste0(
    "Positive significant = ",
    n_pos_fdr
  ),
  paste0(
    "Negative significant = ",
    n_neg_fdr
  ),
  "",
  "No tested gene contributes to the Native-aligned or Alternative-state programme scores.",
  paste0(
    "State-balance target source = ",
    target_source
  ),
  paste0(
    "Direction source used to define frozen programmes = ",
    direction_source
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig5A_nonProgrammeAssociation_v6.txt"
  )
)

cat(
  "\n============================================\n",
  "FIGURE 5A NON-PROGRAMME ASSOCIATION ANALYSIS COMPLETE\n",
  "============================================\n\n",
  sep = ""
)

print(
  stats_df,
  n = Inf
)

cat(
  "\nTop positive non-programme genes:\n"
)

print(
  top_pos %>%
    dplyr::select(
      GeneSymbol,
      Spearman_rho,
      FDR
    ),
  n = Inf
)

cat(
  "\nTop negative non-programme genes:\n"
)

print(
  top_neg %>%
    dplyr::select(
      GeneSymbol,
      Spearman_rho,
      FDR
    ),
  n = Inf
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)
