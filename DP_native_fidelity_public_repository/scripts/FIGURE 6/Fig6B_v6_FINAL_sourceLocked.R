# ============================================================
# Figure 6B v5 — DONOR-BLOCKED + CLEAN-SESSION SUBEXPR
# Marginal association versus multivariate model prominence
#
# Question:
# Are the recurrent Figure 6A model features simply the genes with the
# strongest univariate association with the anatomy-referenced state balance,
# or do some genes show disproportionate multivariate prominence?
#
# X axis:
#   Spearman rho between each non-programme gene and anatomy-referenced state balance
#
# Y axis:
#   Mean model prominence from Figure 6A v7
#
# Inputs:
#   Fig6A_all_heldout_prominence_donorBlocked_v5.csv
#   Fig5C_sample_predictions.csv
#   subexpr
#
# Interpretation:
#   - strong |rho| + high prominence:
#       concordant marginal and multivariate signal
#   - modest |rho| + high prominence:
#       multivariate-specific model contribution
#
# This analysis is descriptive and does not imply causality.
# ============================================================

required_pkgs <- c(
  "ggplot2",
  "readr",
  "dplyr",
  "tibble",
  "stringr",
  "scales",
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
      "\n\nInstall with:\ninstall.packages(c(",
      paste0('"', missing_pkgs, '"', collapse = ", "),
      "))\n"
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(scales)
  library(ggrepel)
})

EXPECTED_N_SAMPLES <- 110
EXPECTED_HELDOUT_N <- 1373
TOP_LABEL_N <- 12

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(
      dirname(
        normalizePath(
          sub("^--file=", "", file_arg[1])
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
    dirname(dirname(dirname(start_dir))),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )
)

candidate_roots <- candidate_roots[
  dir.exists(candidate_roots)
]

score_root <- function(x) {
  s <- 0

  if (
    length(
      list.files(
        x,
        pattern = "^Fig6A_all_nonprogramme_prominence_donorBlocked_v7\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    ) > 0
  ) {
    s <- s + 10
  }

  if (
    length(
      list.files(
        x,
        pattern = "^Fig5C_sample_predictions\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    ) > 0
  ) {
    s <- s + 5
  }

  if (dir.exists(file.path(x, "figures"))) {
    s <- s + 2
  }

  s
}

root_scores <- vapply(
  candidate_roots,
  score_root,
  numeric(1)
)

project_root <- candidate_roots[
  which.max(root_scores)
]

cat(
  "\nProject root:\n",
  project_root,
  "\n"
)

clean_gene <- function(x) {
  x <- toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
  x[x == ""] <- NA_character_
  x
}

safe_read_csv <- function(path) {
  tryCatch(
    readr::read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    ),
    error = function(e) NULL
  )
}

# ------------------------------------------------------------
# Load clean-session subexpr
# ------------------------------------------------------------
loader_file <- file.path(
  project_root,
  "scripts",
  "load_subexpr_clean.R"
)

if (!file.exists(loader_file)) {
  stop(
    paste0(
      "\nQC STOP\nCannot find clean-session loader:\n",
      loader_file,
      "\n"
    )
  )
}

source(loader_file, local = FALSE)

expr <- as.matrix(subexpr)
storage.mode(expr) <- "numeric"

if (nrow(expr) != 1553L) {
  stop(
    paste0(
      "\nQC STOP\nExpected 1553 genes in subexpr, found ",
      nrow(expr),
      ".\n"
    )
  )
}

if (ncol(expr) != EXPECTED_N_SAMPLES) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      EXPECTED_N_SAMPLES,
      " samples in subexpr, found ",
      ncol(expr),
      ".\n"
    )
  )
}

if (is.null(rownames(expr))) {
  stop("\nQC STOP\nsubexpr has no gene row names.\n")
}

rownames(expr) <- clean_gene(rownames(expr))


# ------------------------------------------------------------
# Recover Figure 6A prominence table
# ------------------------------------------------------------
prom_files <- list.files(
  project_root,
  pattern = "^Fig6A_all_nonprogramme_prominence_donorBlocked_v7\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(prom_files) == 0) {
  stop(
    "\nQC STOP\nCannot find Fig6A_all_nonprogramme_prominence_donorBlocked_v7.csv\n"
  )
}

prom_scores <- vapply(
  prom_files,
  function(path) {
    s <- 0
    if (grepl("donorBlocked|donor_blocked", path, ignore.case = TRUE)) s <- s + 100
    if (grepl("FIGURE 6|Fig6", path, ignore.case = TRUE)) s <- s + 10
    s
  },
  numeric(1)
)

prom_file <- prom_files[
  which.max(prom_scores)
]

prom <- safe_read_csv(
  prom_file
)

required_prom_cols <- c(
  "GeneSymbol",
  "Mean_prominence",
  "Top10_frequency",
  "Sign_consistency",
  "Median_standardized_effect",
  "Direction"
)

if (
  is.null(prom) ||
  !all(
    required_prom_cols %in%
      names(prom)
  )
) {
  stop(
    "\nQC STOP\nInvalid Figure 6A prominence table.\n"
  )
}

prom <- prom %>%
  dplyr::mutate(
    GeneSymbol =
      clean_gene(
        GeneSymbol
      )
  )

if (
  nrow(prom) !=
    EXPECTED_HELDOUT_N
) {
  stop(
    paste0(
      "\nQC STOP\nExpected ",
      EXPECTED_HELDOUT_N,
      " non-programme genes in Figure 6A table, found ",
      nrow(prom),
      ".\n"
    )
  )
}

if (
  length(
    setdiff(
      prom$GeneSymbol,
      rownames(expr)
    )
  ) > 0
) {
  stop(
    "\nQC STOP\nSome Figure 6A non-programme genes are missing from subexpr.\n"
  )
}

# ------------------------------------------------------------
# Recover frozen observed balance from Figure 5C
# ------------------------------------------------------------
pred_files <- list.files(
  project_root,
  pattern = "^Fig5C_sample_predictions\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(pred_files) == 0) {
  stop(
    "\nQC STOP\nCannot find Fig5C_sample_predictions.csv\n"
  )
}

score_pred_file <- function(path) {
  s <- 0

  if (grepl("stateBalanceReconstruction_v12|Panel_C_stateBalanceReconstruction_v12", path, ignore.case = TRUE)) {
    s <- s + 200
  }

  if (grepl("donor_blocked_v8|Panel_C_donor", path, ignore.case = TRUE)) {
    s <- s + 100
  }

  if (grepl("Panel_C_v3", path, ignore.case = TRUE)) {
    s <- s + 30
  }

  if (grepl("Panel_C_v2", path, ignore.case = TRUE)) {
    s <- s + 20
  }

  if (grepl("FIGURE 5", path, ignore.case = TRUE)) {
    s <- s + 10
  }

  s
}

pred_scores <- vapply(
  pred_files,
  score_pred_file,
  numeric(1)
)

pred_file <- pred_files[
  which.max(pred_scores)
]

pred_tab <- safe_read_csv(
  pred_file
)

if (
  is.null(pred_tab) ||
  !"Observed_balance" %in%
    names(pred_tab) ||
  nrow(pred_tab) !=
    EXPECTED_N_SAMPLES
) {
  stop(
    "\nQC STOP\nInvalid Figure 5C prediction table.\n"
  )
}

balance <- as.numeric(
  pred_tab$Observed_balance
)

expr_samples <- colnames(expr)

if (
  "Sample" %in%
    names(pred_tab) &&
  !is.null(expr_samples) &&
  setequal(
    as.character(
      pred_tab$Sample
    ),
    expr_samples
  )
) {
  balance <- balance[
    match(
      expr_samples,
      as.character(
        pred_tab$Sample
      )
    )
  ]

  alignment_mode <-
    "Exact sample-name alignment"
} else {
  alignment_mode <-
    "Figure 5C row order retained"
}

if (
  length(balance) !=
    EXPECTED_N_SAMPLES ||
  any(!is.finite(balance))
) {
  stop(
    "\nQC STOP\nInvalid anatomy-referenced state-balance vector.\n"
  )
}

# ------------------------------------------------------------
# Compute marginal Spearman association for all non-programme genes
# ------------------------------------------------------------
cat(
  "\nComputing marginal Spearman associations for",
  EXPECTED_HELDOUT_N,
  "non-programme genes...\n"
)

assoc_rows <- lapply(
  prom$GeneSymbol,
  function(g) {
    xg <- as.numeric(
      expr[
        g,
        ,
        drop = TRUE
      ]
    )

    ct <- suppressWarnings(
      cor.test(
        xg,
        balance,
        method = "spearman",
        exact = FALSE
      )
    )

    tibble::tibble(
      GeneSymbol = g,
      Spearman_rho =
        unname(
          ct$estimate
        ),
      P_value =
        ct$p.value
    )
  }
)

assoc <- dplyr::bind_rows(
  assoc_rows
) %>%
  dplyr::mutate(
    FDR =
      p.adjust(
        P_value,
        method = "BH"
      ),
    Abs_rho =
      abs(
        Spearman_rho
      )
  )

# ------------------------------------------------------------
# Join with Figure 6A prominence
# ------------------------------------------------------------
joint <- prom %>%
  dplyr::left_join(
    assoc,
    by = "GeneSymbol"
  )

if (
  any(
    !is.finite(
      joint$Spearman_rho
    )
  )
) {
  stop(
    "\nQC STOP\nNon-finite Spearman rho detected after merge.\n"
  )
}

# ------------------------------------------------------------
# Quantify marginal-vs-multivariate relationship
# ------------------------------------------------------------
rho_prominence <- suppressWarnings(
  cor(
    joint$Abs_rho,
    joint$Mean_prominence,
    method = "spearman"
  )
)

rho_signed <- suppressWarnings(
  cor(
    joint$Spearman_rho,
    joint$Median_standardized_effect,
    method = "spearman"
  )
)

# Candidate multivariate-specific genes:
# high prominence, modest marginal association.
prom_q90 <- unname(
  quantile(
    joint$Mean_prominence,
    0.90,
    na.rm = TRUE
  )
)

rho_q50 <- unname(
  quantile(
    joint$Abs_rho,
    0.50,
    na.rm = TRUE
  )
)

joint <- joint %>%
  dplyr::mutate(
    Multivariate_specific =
      Mean_prominence >=
        prom_q90 &
      Abs_rho <=
        rho_q50
  )

# ------------------------------------------------------------
# Select labels
# ------------------------------------------------------------
label_top_prom <- joint %>%
  dplyr::arrange(
    dplyr::desc(
      Mean_prominence
    )
  ) %>%
  dplyr::slice_head(
    n = 8
  )

label_specific <- joint %>%
  dplyr::filter(
    Multivariate_specific
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Mean_prominence
    ),
    Abs_rho
  ) %>%
  dplyr::slice_head(
    n = 6
  )

label_df <- dplyr::bind_rows(
  label_top_prom,
  label_specific
) %>%
  dplyr::distinct(
    GeneSymbol,
    .keep_all = TRUE
  ) %>%
  dplyr::slice_head(
    n = TOP_LABEL_N
  )

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 6",
  "Panel_B_marginalVsMultivariate_donorBlocked_v5"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Export tables
# ------------------------------------------------------------
readr::write_csv(
  joint,
  file.path(
    out_dir,
    "Fig6B_all_nonprogramme_marginal_vs_multivariate.csv"
  )
)

readr::write_csv(
  label_df,
  file.path(
    out_dir,
    "Fig6B_labeled_genes.csv"
  )
)

readr::write_csv(
  joint %>%
    dplyr::filter(
      Multivariate_specific
    ) %>%
    dplyr::arrange(
      dplyr::desc(
        Mean_prominence
      )
    ),
  file.path(
    out_dir,
    "Fig6B_multivariate_specific_candidates.csv"
  )
)

# ------------------------------------------------------------
# Main plot -- Nature-style Figure 6B
# ------------------------------------------------------------
# Scientific calculations above are unchanged.
# Only the visual layout is revised for the final 2 x 2 Figure 6 assembly.

COL_NATIVE <- "#3D8C82"
COL_ALT <- "#C97857"
COL_OTHER <- "#B9C5CA"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"

p <- ggplot(
  joint,
  aes(
    x = Spearman_rho,
    y = Mean_prominence
  )
) +
  geom_hline(
    yintercept = prom_q90,
    linetype = "dashed",
    linewidth = 0.28,
    colour = "#AAB7BC"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.28,
    colour = "#AAB7BC"
  ) +
  geom_point(
    aes(colour = Direction),
    size = 1.15,
    alpha = 0.48
  ) +
  geom_point(
    data = label_df,
    aes(
      x = Spearman_rho,
      y = Mean_prominence,
      colour = Direction
    ),
    size = 1.75,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = GeneSymbol),
    size = 1.70,
    family = "Helvetica",
    colour = COL_TEXT,
    box.padding = 0.22,
    point.padding = 0.14,
    segment.size = 0.20,
    segment.colour = "#9AABB2",
    max.overlaps = Inf,
    min.segment.length = 0,
    force = 1.3,
    force_pull = 0.8,
    direction = "both",
    seed = 20261030,
    max.time = 2
  ) +
  scale_colour_manual(
    values = c(
      "Native-aligned" = COL_NATIVE,
      "Alternative-aligned" = COL_ALT,
      "Neutral" = COL_OTHER
    ),
    breaks = c(
      "Alternative-aligned",
      "Native-aligned"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    limits = c(0, 1.03),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_x_continuous(
    breaks = c(-0.2, 0, 0.2, 0.4),
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  labs(
    x = "Marginal Spearman correlation with anatomy-referenced state balance",
    y = "Mean multivariate model prominence",
    title = "Multivariate prominence broadly tracks marginal association",
    subtitle = paste0(
      "Non-programme genes: n = ",
      nrow(joint),
      "; Spearman |rho| vs prominence = ",
      sprintf("%.2f", rho_prominence)
    )
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 1.8,
        alpha = 1
      )
    )
  ) +
  theme_classic(
    base_family = "Helvetica",
    base_size = 5.0
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      size = 6.0,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      lineheight = 1.0,
      margin = margin(b = 1.0, unit = "pt")
    ),
    plot.subtitle = element_text(
      size = 4.45,
      colour = COL_MUTED,
      hjust = 0.5,
      margin = margin(b = 1.0, unit = "pt")
    ),

    # Compact legend directly below subtitle.
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.margin = margin(0, 0, 1.0, 0, unit = "pt"),
    legend.spacing.x = grid::unit(5.0, "pt"),
    legend.key.width = grid::unit(6.0, "pt"),
    legend.key.height = grid::unit(5.0, "pt"),
    legend.text = element_text(
      size = 4.15,
      colour = COL_MUTED
    ),

    axis.title.x = element_text(
      size = 4.75,
      colour = COL_TEXT,
      margin = margin(t = 1.8, unit = "pt")
    ),
    axis.title.y = element_text(
      size = 4.75,
      colour = COL_TEXT,
      margin = margin(r = 1.8, unit = "pt")
    ),
    axis.text = element_text(
      size = 4.35,
      colour = COL_MUTED
    ),
    axis.line = element_line(
      linewidth = 0.30,
      colour = COL_TEXT
    ),
    axis.ticks = element_line(
      linewidth = 0.24,
      colour = COL_TEXT
    ),
    axis.ticks.length = grid::unit(1.2, "pt"),

    plot.margin = margin(
      t = 2.0,
      r = 2.8,
      b = 2.0,
      l = 2.5,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# Nature-style export
# ------------------------------------------------------------
# Match Figure 6A so A/B form a strict top-row pair.
FIG_W_MM <- 88
FIG_H_MM <- 62

out_dir_nature <- out_dir

dir.create(
  out_dir_nature,
  recursive = TRUE,
  showWarnings = FALSE
)

png_file <- file.path(
  out_dir_nature,
  "Fig6_PanelB_marginalVsMultivariate_donorBlocked_v5.png"
)

pdf_file <- file.path(
  out_dir_nature,
  "Fig6_PanelB_marginalVsMultivariate_donorBlocked_v5.pdf"
)

tiff_file <- file.path(
  out_dir_nature,
  "Fig6_PanelB_marginalVsMultivariate_donorBlocked_v5.tiff"
)

ggsave(
  filename = png_file,
  plot = p,
  width = FIG_W_MM,
  height = FIG_H_MM,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# Base PDF device for macOS compatibility.
grDevices::pdf(
  file = pdf_file,
  width = FIG_W_MM / 25.4,
  height = FIG_H_MM / 25.4,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE,
  paper = "special"
)
print(p)
grDevices::dev.off()

if (!file.exists(pdf_file) ||
    is.na(file.info(pdf_file)$size) ||
    file.info(pdf_file)$size <= 0) {
  stop(
    paste0(
      "PDF export failed: ",
      pdf_file
    )
  )
}

ggsave(
  filename = tiff_file,
  plot = p,
  width = FIG_W_MM,
  height = FIG_H_MM,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

cat(
  "\nNature-style Figure 6B written to:\n",
  out_dir_nature,
  "\n"
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------
qc <- tibble::tibble(
  Metric = c(
    "Non-programme genes",
    "Spearman corr |rho| vs mean prominence",
    "Spearman corr signed rho vs median standardized effect",
    "90th percentile prominence threshold",
    "Median |rho| threshold",
    "Multivariate-specific candidate genes",
    "Genes with FDR < 0.05 marginal association"
  ),
  Value = c(
    nrow(joint),
    rho_prominence,
    rho_signed,
    prom_q90,
    rho_q50,
    sum(
      joint$Multivariate_specific
    ),
    sum(
      joint$FDR <
        0.05
    )
  )
)

readr::write_csv(
  qc,
  file.path(
    out_dir,
    "Fig6B_QC_summary.csv"
  )
)

# ------------------------------------------------------------
# README
# ------------------------------------------------------------
readme <- c(
  "Figure 6B v5 — donor-blocked complete archive",
  "",
  "Marginal association versus donor-blocked multivariate model prominence.",
  "",
  "Purpose:",
  "Determine whether Figure 6A recurrent model features are simply the strongest univariate correlates of the anatomy-referenced state balance.",
  "",
  "X axis:",
  "Spearman rho between each non-programme gene and the frozen anatomy-referenced state balance.",
  "",
  "Y axis:",
  "Mean model prominence from Figure 6A v7 donor-blocked analysis.",
  "",
  "Interpretation:",
  "The analysis quantifies concordance between marginal association and donor-blocked multivariate prominence; the threshold-defined candidate subset is reported descriptively.",
  "",
  paste0(
    "Spearman correlation between |rho| and prominence = ",
    sprintf(
      "%.4f",
      rho_prominence
    )
  ),
  paste0(
    "Spearman correlation between signed rho and median standardized effect = ",
    sprintf(
      "%.4f",
      rho_signed
    )
  ),
  paste0(
    "Multivariate-specific candidates = ",
    sum(
      joint$Multivariate_specific
    )
  ),
  "",
  "Limitations:",
  "This is a descriptive decomposition of marginal association versus multivariate model prominence.",
  "It does not establish causal or mechanistic roles.",
  "",
  paste0("Figure 6A prominence source: ", normalizePath(prom_file)),
  paste0("Figure 5C target source: ", normalizePath(pred_file)),
  paste0("Sample alignment: ", alignment_mode)
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig6_PanelB_donorBlocked_v5.txt"
  )
)

cat(
  "\n============================================\n"
)

cat(
  "FIGURE 6B v5 DONOR-BLOCKED COMPLETE\n"
)

cat(
  "============================================\n\n"
)

print(
  qc,
  n = Inf
)

cat(
  "\nLabeled genes:\n"
)

print(
  label_df %>%
    dplyr::select(
      GeneSymbol,
      Spearman_rho,
      FDR,
      Mean_prominence,
      Top10_frequency,
      Sign_consistency,
      Direction,
      Multivariate_specific
    ),
  n = Inf
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)

cat(
  "\nMain figure:\n",
  file.path(
    out_dir,
    "Fig6_PanelB_marginalVsMultivariate_v1.pdf"
  ),
  "\n"
)
