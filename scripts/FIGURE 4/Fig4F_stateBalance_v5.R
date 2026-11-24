# ============================================================
# Figure 4F v5 — terminology-aligned clean-session reproducible
# Global concordance of all 180 frozen genes with the
# anatomy-referenced state balance
#
# Purpose:
# Complement Figure 4E by showing ALL 180 leave-one-gene-out
# associations rather than only the strongest genes.
#
# Method:
#   - same frozen 81 / 99 programmes
#   - same subexpr matrix
#   - same direction alignment
#   - for each gene, remove it from its own programme score
#   - correlate direction-aligned gene expression with the
#     leave-one-out Native - Alternative balance
#
# Figure:
#   Distribution of leave-one-out Spearman rho by programme
#
# Statistics:
#   Wilcoxon rank-sum test
#   Cliff's delta
#   sign-concordance proportions
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
  "ggplot2"
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

loader_file <- file.path(
  project_root,
  "scripts",
  "load_subexpr_clean.R"
)

if (!file.exists(loader_file)) {
  stop(
    paste0(
      "\nCannot find clean-session loader:\n",
      loader_file,
      "\n"
    )
  )
}

source(
  loader_file,
  local = FALSE
)

expr_mat <- subexpr

if (
  !is.matrix(expr_mat) &&
  !is.data.frame(expr_mat)
) {
  stop("\nsubexpr must be a matrix or data.frame.\n")
}

expr_mat <- as.matrix(expr_mat)
storage.mode(expr_mat) <- "numeric"

if (
  is.null(rownames(expr_mat)) ||
  is.null(colnames(expr_mat))
) {
  stop("\nsubexpr must have gene rownames and sample colnames.\n")
}

rownames(expr_mat) <- toupper(
  stringr::str_trim(
    rownames(expr_mat)
  )
)

if (nrow(expr_mat) != 1553L) {
  stop(
    paste0(
      "\nQC STOP: expected 1553 genes in subexpr, found ",
      nrow(expr_mat),
      ".\n"
    )
  )
}

if (ncol(expr_mat) != 110L) {
  stop(
    paste0(
      "\nQC STOP: expected 110 profiles in subexpr, found ",
      ncol(expr_mat),
      ".\n"
    )
  )
}

cat(
  "\nsubexpr dimensions:",
  nrow(expr_mat),
  "x",
  ncol(expr_mat),
  "\n"
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

format_p <- function(p) {

  if (!is.finite(p)) {
    return("NA")
  }

  if (p == 0 || p < 2.2e-16) {
    return("< 2.2e-16")
  }

  if (p < 0.001) {
    return(
      format(
        p,
        scientific = TRUE,
        digits = 2
      )
    )
  }

  sprintf(
    "%.3f",
    p
  )
}


# ------------------------------------------------------------
# 4. Frozen programmes
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

if (length(native81) != 81) {
  stop("\nNative frozen program != 81 genes.\n")
}

if (length(alternative99) != 99) {
  stop("\nAlternative frozen program != 99 genes.\n")
}

if (
  !all(
    native81 %in%
      rownames(expr_mat)
  )
) {
  stop("\nsubexpr does not contain all 81 native genes.\n")
}

if (
  !all(
    alternative99 %in%
      rownames(expr_mat)
  )
) {
  stop("\nsubexpr does not contain all 99 alternative genes.\n")
}


# ------------------------------------------------------------
# 5. Audit + direction
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
  stop("\nCannot identify gene column in audit table.\n")
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
    "\nCannot recover restored-up/restored-down direction.\n"
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


# ------------------------------------------------------------
# 6. Frozen direction QC
# ------------------------------------------------------------

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
# 7. Gene-wise z scoring + direction alignment
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
          !is.finite(s) ||
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

  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)

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

native_full <- colMeans(
  native_aligned,
  na.rm = TRUE
)

alt_full <- colMeans(
  alt_aligned,
  na.rm = TRUE
)

sample_ids <- colnames(
  expr_mat
)


# ------------------------------------------------------------
# 8. Leave-one-gene-out associations
# ------------------------------------------------------------

driver_rows <- list()
idx <- 1


for (g in native81) {

  native_loo <- colMeans(
    native_aligned[
      setdiff(
        native81,
        g
      ),
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  balance_loo <- native_loo -
    alt_full

  gene_aligned <- as.numeric(
    native_aligned[
      g,
      sample_ids
    ]
  )

  ct <- suppressWarnings(
    stats::cor.test(
      gene_aligned,
      balance_loo,
      method = "spearman",
      exact = FALSE
    )
  )

  driver_rows[[idx]] <- tibble::tibble(
    GeneSymbol = g,
    Program = "Native-aligned programme",
    Spearman_rho = unname(
      ct$estimate
    ),
    P_value = ct$p.value
  )

  idx <- idx + 1
}


for (g in alternative99) {

  alt_loo <- colMeans(
    alt_aligned[
      setdiff(
        alternative99,
        g
      ),
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  balance_loo <- native_full -
    alt_loo

  gene_aligned <- as.numeric(
    alt_aligned[
      g,
      sample_ids
    ]
  )

  ct <- suppressWarnings(
    stats::cor.test(
      gene_aligned,
      balance_loo,
      method = "spearman",
      exact = FALSE
    )
  )

  driver_rows[[idx]] <- tibble::tibble(
    GeneSymbol = g,
    Program = "Alternative-state programme",
    Spearman_rho = unname(
      ct$estimate
    ),
    P_value = ct$p.value
  )

  idx <- idx + 1
}


drivers <- dplyr::bind_rows(
  driver_rows
) %>%

  dplyr::mutate(
    FDR =
      p.adjust(
        P_value,
        method = "BH"
      )
  )


# ------------------------------------------------------------
# 9. Programme-level statistics
# ------------------------------------------------------------

native_rho <- drivers %>%
  dplyr::filter(
    Program ==
      "Native-aligned programme"
  ) %>%
  dplyr::pull(
    Spearman_rho
  )

alt_rho <- drivers %>%
  dplyr::filter(
    Program ==
      "Alternative-state programme"
  ) %>%
  dplyr::pull(
    Spearman_rho
  )


wilcox_res <- suppressWarnings(
  stats::wilcox.test(
    native_rho,
    alt_rho,
    alternative = "two.sided",
    exact = FALSE
  )
)

wilcox_p <- wilcox_res$p.value


# Cliff's delta:
# positive => native rho values larger than alternative rho values
delta_matrix <- outer(
  native_rho,
  alt_rho,
  FUN = "-"
)

cliffs_delta <- mean(
  sign(
    delta_matrix
  )
)

abs_delta <- abs(
  cliffs_delta
)

delta_magnitude <- dplyr::case_when(
  abs_delta < 0.147 ~ "negligible",
  abs_delta < 0.330 ~ "small",
  abs_delta < 0.474 ~ "medium",
  TRUE ~ "large"
)


native_positive_fraction <- mean(
  native_rho > 0
)

alt_negative_fraction <- mean(
  alt_rho < 0
)

native_fdr_positive_fraction <- mean(
  drivers$Program ==
    "Native-aligned programme" &
  drivers$Spearman_rho > 0 &
  drivers$FDR < 0.05
) / mean(
  drivers$Program ==
    "Native-aligned programme"
)

alt_fdr_negative_fraction <- mean(
  drivers$Program ==
    "Alternative-state programme" &
  drivers$Spearman_rho < 0 &
  drivers$FDR < 0.05
) / mean(
  drivers$Program ==
    "Alternative-state programme"
)


summary_df <- drivers %>%

  dplyr::group_by(
    Program
  ) %>%

  dplyr::summarise(
    n =
      dplyr::n(),
    Median_rho =
      median(
        Spearman_rho
      ),
    Q1 =
      as.numeric(
        quantile(
          Spearman_rho,
          0.25
        )
      ),
    Q3 =
      as.numeric(
        quantile(
          Spearman_rho,
          0.75
        )
      ),
    FDR_lt_0_05 =
      sum(
        FDR < 0.05
      ),
    Positive_rho =
      sum(
        Spearman_rho > 0
      ),
    Negative_rho =
      sum(
        Spearman_rho < 0
      ),
    .groups =
      "drop"
  )


cat(
  "\nProgramme-level LOO rho summary:\n"
)

print(
  summary_df,
  n = Inf
)

cat(
  "\nWilcoxon P = ",
  format_p(
    wilcox_p
  ),
  "\n",
  "Cliff's delta = ",
  sprintf(
    "%.4f",
    cliffs_delta
  ),
  " (",
  delta_magnitude,
  ")\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Plot — horizontal wide / short layout
# ------------------------------------------------------------

COL_NATIVE <- "#3D8C82"
COL_ALT <- "#C97857"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"
COL_ZERO <- "#A8B5BA"

program_cols <- c(
  "Native-aligned programme" = COL_NATIVE,
  "Alternative-state programme" = COL_ALT
)

plot_df <- drivers %>%
  dplyr::mutate(
    Program = factor(
      Program,
      levels = c(
        "Native-aligned programme",
        "Alternative-state programme"
      )
    )
  )

p_label <- if (
  wilcox_p == 0 ||
  wilcox_p < 2.2e-16
) {
  "Wilcoxon P < 2.2e-16"
} else {
  paste0(
    "Wilcoxon P = ",
    format_p(wilcox_p)
  )
}

stat_label <- paste0(
  p_label,
  "; Cliff's delta = ",
  sprintf("%.2f", cliffs_delta),
  " (",
  delta_magnitude,
  ")"
)

p <- ggplot(
  plot_df,
  aes(
    y = Program,
    x = Spearman_rho,
    fill = Program,
    colour = Program
  )
) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.30,
    colour = COL_ZERO,
    linetype = "22"
  ) +

  geom_violin(
    width = 0.78,
    trim = TRUE,
    alpha = 0.15,
    linewidth = 0.45
  ) +

  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.30,
    linewidth = 0.52
  ) +

  geom_jitter(
    height = 0.07,
    width = 0,
    size = 0.95,
    alpha = 0.48,
    stroke = 0
  ) +

  scale_fill_manual(
    values = program_cols,
    guide = "none"
  ) +

  scale_colour_manual(
    values = program_cols,
    guide = "none"
  ) +

  scale_y_discrete(
    labels = c(
      "Native-aligned programme" =
        "Native-aligned programme",
      "Alternative-state programme" =
        "Alternative-state programme"
    )
  ) +

  labs(
    x = "Leave-one-out Spearman rho",
    y = NULL,
    title = "Frozen programmes show opposing state-balance associations",
    subtitle = stat_label,
    caption = paste0(
      "All 180 frozen genes shown; each gene excluded from its own programme score before testing."
    )
  ) +

  theme_classic(
    base_family = "Helvetica",
    base_size = 7
  ) +

  theme(
    plot.title = element_text(
      size = 5.1,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      margin = margin(b = 3)
    ),

    plot.subtitle = element_text(
      size = 4.3,
      colour = COL_MUTED,
      hjust = 0.5,
      margin = margin(b = 5)
    ),

    plot.caption = element_text(
      size = 4.1,
      colour = "#95A4AA",
      hjust = 0.5,
      margin = margin(t = 3.2)
    ),

    axis.title.x = element_text(
      size = 5.0,
      colour = COL_TEXT,
      margin = margin(t = 3)
    ),

    axis.text.x = element_text(
      size = 4.6,
      colour = COL_MUTED
    ),

    axis.text.y = element_text(
      size = 5.0,
      colour = COL_TEXT,
      margin = margin(r = 4)
    ),

    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),

    axis.line.x = element_line(
      linewidth = 0.30,
      colour = COL_TEXT
    ),

    axis.ticks.x = element_line(
      linewidth = 0.26,
      colour = COL_TEXT
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
  "FIGURE 4",
  "Panel_F_globalConcordance_stateBalance_v5"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  filename = file.path(
    out_dir,
    "Fig4_PanelF_globalConcordance_stateBalance_v5.png"
  ),
  plot = p,
  width = 88,
  height = 50,
  units = "mm",
  dpi = 600,
  bg = "white"
)


pdf(
  file = file.path(
    out_dir,
    "Fig4_PanelF_globalConcordance_stateBalance_v5.pdf"
  ),
  width = 88 / 25.4,
  height = 50 / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)

print(
  p
)

dev.off()


ggsave(
  filename = file.path(
    out_dir,
    "Fig4_PanelF_globalConcordance_stateBalance_v5.tiff"
  ),
  plot = p,
  width = 88,
  height = 50,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


readr::write_csv(
  drivers,
  file.path(
    out_dir,
    "Fig4F_all_LOO_gene_associations.csv"
  )
)

readr::write_csv(
  summary_df,
  file.path(
    out_dir,
    "Fig4F_program_summary.csv"
  )
)


stats_df <- tibble::tibble(
  Test = "Wilcoxon rank-sum",
  Native_n = length(
    native_rho
  ),
  Alternative_n = length(
    alt_rho
  ),
  Native_median_rho = median(
    native_rho
  ),
  Alternative_median_rho = median(
    alt_rho
  ),
  Wilcoxon_W = unname(
    wilcox_res$statistic
  ),
  Wilcoxon_P = wilcox_p,
  Cliffs_delta = cliffs_delta,
  Delta_magnitude = delta_magnitude,
  Native_fraction_rho_positive = native_positive_fraction,
  Alternative_fraction_rho_negative = alt_negative_fraction
)

readr::write_csv(
  stats_df,
  file.path(
    out_dir,
    "Fig4F_statistics.csv"
  )
)


qc_df <- tibble::tibble(
  Item = c(
    "subexpr genes",
    "subexpr samples",
    "Native frozen genes",
    "Alternative frozen genes",
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
    "Fig4F_QC.csv"
  )
)


readme <- c(
  "Figure 4F global leave-one-out concordance v5",
  "",
  "Purpose:",
  "Show the full distribution of leave-one-out gene associations for all 180 frozen genes.",
  "",
  "Method:",
  "Each gene was removed from its own programme score before calculating Native-aligned - Alternative-state balance.",
  "Direction-aligned gene expression was then correlated with that leave-one-out balance.",
  "",
  paste0(
    "Native median rho = ",
    sprintf(
      "%.4f",
      median(
        native_rho
      )
    )
  ),
  paste0(
    "Alternative median rho = ",
    sprintf(
      "%.4f",
      median(
        alt_rho
      )
    )
  ),
  paste0(
    "Wilcoxon P = ",
    format_p(
      wilcox_p
    )
  ),
  paste0(
    "Cliff's delta = ",
    sprintf(
      "%.4f",
      cliffs_delta
    ),
    " (",
    delta_magnitude,
    ")"
  ),
  paste0(
    "Native fraction rho > 0 = ",
    sprintf(
      "%.4f",
      native_positive_fraction
    )
  ),
  paste0(
    "Alternative fraction rho < 0 = ",
    sprintf(
      "%.4f",
      alt_negative_fraction
    )
  ),
  "",
  paste0(
    "Direction source: ",
    direction_source
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig4F_globalConcordance_stateBalance_v5.txt"
  )
)


cat(
  "\n============================================\n",
  "FIGURE 4F v5 COMPLETE\n",
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
