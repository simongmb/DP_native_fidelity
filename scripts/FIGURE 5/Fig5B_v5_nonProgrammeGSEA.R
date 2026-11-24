# ============================================================
# Figure 5B v5 — OFFLINE terminology-aligned
# Non-programme gene preranked GO Biological Process GSEA
#
# PURPOSE
# ------------------------------------------------------------
# Same biological question as Fig5B v1, but WITHOUT msigdbr
# and WITHOUT any internet/Zenodo download.
#
# Uses local Bioconductor annotations:
#   org.Hs.eg.db
#   GO.db
#   AnnotationDbi
#   fgsea
#
# Input:
#   subexpr = 1553 genes x 110 samples
#
# Ranking:
#   1373 non-programme genes
#   Spearman rho with anatomy-referenced state balance
#
# Interpretation:
#   NES > 0 -> Native-aligned end enriched BP
#   NES < 0 -> Alternative-aligned end enriched BP
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
  "fgsea",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "GO.db"
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
      "\n\nInstall CRAN packages with install.packages().\n",
      "Install Bioconductor packages with BiocManager::install().\n"
    )
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
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

wrap_term <- function(x, width = 52) {
  paste(
    strwrap(
      x,
      width = width
    ),
    collapse = "\n"
  )
}


# ------------------------------------------------------------
# 4. Frozen programmes + non-programme background
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
  stop("\nFrozen background != 1553 genes.\n")
}

if (length(native81) != 81) {
  stop("\nNative frozen program != 81 genes.\n")
}

if (length(alternative99) != 99) {
  stop("\nAlternative frozen program != 99 genes.\n")
}

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

if (length(heldout_genes) != 1373) {
  stop(
    "\nHeld-out background is not exactly 1,373 genes.\n"
  )
}

if (length(heldout_present) < 0.95 * 1373) {
  stop(
    "\nToo many held-out genes are missing from subexpr.\n"
  )
}


# ------------------------------------------------------------
# 5. Audit + restoration direction
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
# 6. Programme scoring
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

native_score <- colMeans(
  native_aligned,
  na.rm = TRUE
)

alt_score <- colMeans(
  alt_aligned,
  na.rm = TRUE
)

balance_score <- native_score -
  alt_score

names(balance_score) <- sample_ids


# ------------------------------------------------------------
# 7. Non-programme gene ranking
# ------------------------------------------------------------

rank_rows <- lapply(
  heldout_present,
  function(g) {

    values <- as.numeric(
      expr_mat[
        g,
        sample_ids
      ]
    )

    rho <- suppressWarnings(
      stats::cor(
        values,
        balance_score,
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )

    tibble::tibble(
      GeneSymbol =
        g,
      Spearman_rho =
        rho
    )
  }
)

rank_df <- dplyr::bind_rows(
  rank_rows
) %>%
  dplyr::filter(
    is.finite(
      Spearman_rho
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Spearman_rho
    )
  )

ranks <- rank_df$Spearman_rho

names(ranks) <- rank_df$GeneSymbol

ranks <- sort(
  ranks,
  decreasing = TRUE
)

cat(
  "\nNon-programme ranked genes:",
  length(ranks),
  "\n"
)


# ------------------------------------------------------------
# 8. Build GO-BP pathways OFFLINE
#
# SYMBOL -> GO mapping via org.Hs.eg.db
# Keep ontology == BP.
# ------------------------------------------------------------

cat(
  "\nBuilding local GO Biological Process pathways...\n"
)

go_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = names(ranks),
  keytype = "SYMBOL",
  columns = c(
    "SYMBOL",
    "GO",
    "ONTOLOGY"
  )
)

go_map <- as.data.frame(
  go_map,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(
    !is.na(SYMBOL),
    !is.na(GO),
    ONTOLOGY == "BP",
    SYMBOL %in%
      names(ranks)
  ) %>%
  dplyr::distinct(
    SYMBOL,
    GO
  )

if (nrow(go_map) == 0) {
  stop(
    "\nNo local GO-BP annotations recovered for ranked genes.\n"
  )
}


# GO term names
go_ids <- unique(
  go_map$GO
)

go_terms <- AnnotationDbi::select(
  GO.db,
  keys = go_ids,
  keytype = "GOID",
  columns = c(
    "GOID",
    "TERM",
    "ONTOLOGY"
  )
)

go_terms <- as.data.frame(
  go_terms,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(
    ONTOLOGY == "BP"
  ) %>%
  dplyr::distinct(
    GOID,
    .keep_all = TRUE
  )


pathway_tbl <- go_map %>%
  dplyr::inner_join(
    go_terms %>%
      dplyr::select(
        GOID,
        TERM
      ),
    by = c(
      "GO" =
        "GOID"
    )
  )


pathways <- split(
  pathway_tbl$SYMBOL,
  pathway_tbl$GO
)

pathways <- lapply(
  pathways,
  unique
)


# ------------------------------------------------------------
# 9. Preranked GSEA
# ------------------------------------------------------------

set.seed(
  20261028
)

fg <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = ranks,
  minSize = 10,
  maxSize = 500,
  eps = 0
)

fg_df <- as.data.frame(
  fg,
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(
    go_terms %>%
      dplyr::select(
        GOID,
        TERM
      ),
    by = c(
      "pathway" =
        "GOID"
    )
  ) %>%
  dplyr::mutate(
    Direction =
      dplyr::case_when(
        NES > 0 ~
          "Native-aligned enriched",
        NES < 0 ~
          "Alternative-aligned enriched",
        TRUE ~
          "Neutral"
      )
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(
      abs(
        NES
      )
    )
  )


# ------------------------------------------------------------
# 10. Significant pathway selection
# ------------------------------------------------------------

sig_df <- fg_df %>%
  dplyr::filter(
    is.finite(
      NES
    ),
    is.finite(
      padj
    ),
    padj < 0.05,
    !is.na(
      TERM
    )
  )

n_sig <- nrow(
  sig_df
)

n_native_sig <- sum(
  sig_df$NES > 0
)

n_alt_sig <- sum(
  sig_df$NES < 0
)

TOP_PER_DIRECTION <- 8


top_native <- sig_df %>%
  dplyr::filter(
    NES > 0
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(
      NES
    )
  ) %>%
  dplyr::slice_head(
    n = TOP_PER_DIRECTION
  )


top_alt <- sig_df %>%
  dplyr::filter(
    NES < 0
  ) %>%
  dplyr::arrange(
    padj,
    NES
  ) %>%
  dplyr::slice_head(
    n = TOP_PER_DIRECTION
  )


plot_df <- dplyr::bind_rows(
  top_alt,
  top_native
) %>%
  dplyr::distinct(
    pathway,
    .keep_all = TRUE
  ) %>%
  dplyr::arrange(
    NES
  ) %>%
  dplyr::mutate(
    Term_wrapped =
      vapply(
        TERM,
        wrap_term,
        character(1)
      ),
    Term_plot =
      factor(
        Term_wrapped,
        levels =
          Term_wrapped
      ),
    NegLog10FDR =
      -log10(
        pmax(
          padj,
          1e-300
        )
      )
  )


if (nrow(plot_df) == 0) {
  stop(
    paste0(
      "\nNo local GO Biological Process pathways reach ",
      "BH FDR < 0.05.\n",
      "All GSEA results will still be written before stop only ",
      "if you modify the script; no main figure is forced.\n"
    )
  )
}


cat(
  "\nGSEA summary:\n",
  "Significant GO-BP pathways = ",
  n_sig,
  "\n",
  "Native-aligned = ",
  n_native_sig,
  "\n",
  "Alternative-aligned = ",
  n_alt_sig,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Plot — Nature-style wide / compact layout
# ------------------------------------------------------------

COL_NATIVE <- "#3D8C82"
COL_ALT <- "#C97857"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"
COL_ZERO <- "#A8B5BA"

direction_cols <- c(
  "Native-aligned enriched" = COL_NATIVE,
  "Alternative-aligned enriched" = COL_ALT
)

p <- ggplot(
  plot_df,
  aes(
    x = NES,
    y = Term_plot,
    colour = Direction
  )
) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.28,
    colour = COL_ZERO
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = Term_plot
    ),
    linewidth = 0.52,
    alpha = 0.76
  ) +

  geom_point(
    aes(
      size = NegLog10FDR
    ),
    alpha = 0.96
  ) +

  scale_colour_manual(
    values = direction_cols,
    name = NULL
  ) +

  scale_size_continuous(
    range = c(
      1.7,
      3.1
    ),
    name =
      expression(
        -log[10]("FDR")
      )
  ) +

  scale_x_continuous(
    breaks = c(
      -2,
      -1,
      0,
      1,
      2
    ),
    expand = expansion(
      mult = c(
        0.05,
        0.06
      )
    )
  ) +

  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    title = "Non-programme genes show opposing biological-process associations",
    subtitle = paste0(
      "Preranked GO-BP GSEA across ",
      length(ranks),
      " non-programme genes; positive NES = Native-aligned end"
    ),
    caption = "Local GO-BP annotations; displayed terms meet BH FDR < 0.05."
  ) +

  theme_classic(
    base_family = "Helvetica",
    base_size = 7
  ) +

  theme(
    plot.title = element_text(
      size = 5.6,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      margin = margin(
        b = 1.0
      )
    ),

    plot.subtitle = element_text(
      size = 4.4,
      colour = COL_MUTED,
      hjust = 0.5,
      margin = margin(
        b = 2.5
      )
    ),

    plot.caption = element_text(
      size = 4.0,
      colour = "#95A4AA",
      hjust = 0.5,
      margin = margin(
        t = 2.5
      )
    ),

    axis.title.x = element_text(
      size = 5.0,
      colour = COL_TEXT,
      margin = margin(
        t = 3
      )
    ),

    axis.text.x = element_text(
      size = 4.6,
      colour = COL_MUTED
    ),

    axis.text.y = element_text(
      size = 4.15,
      colour = COL_TEXT,
      lineheight = 1.02,
      margin = margin(
        r = 8
      )
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

    axis.ticks.length = unit(
      1.0,
      "mm"
    ),

    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical",

    legend.title = element_text(
      size = 4.5,
      colour = COL_TEXT
    ),

    legend.text = element_text(
      size = 4.35,
      colour = COL_MUTED
    ),

    legend.key.height = unit(
      2.6,
      "mm"
    ),

    legend.key.width = unit(
      3.0,
      "mm"
    ),

    legend.spacing.y = unit(
      0.8,
      "mm"
    ),

    legend.margin = margin(
      l = 3
    ),

    plot.margin = margin(
      4,
      6,
      4,
      10,
      unit = "pt"
    )
  )


# ------------------------------------------------------------
# 12. Output
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 5",
  "Panel_B_nonProgrammeGSEA_Nature_wide_v5"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  filename = file.path(
    out_dir,
    "Fig5_PanelB_nonProgrammeGSEA_Nature_wide_v5.png"
  ),
  plot = p,
  width = 94,
  height = 56,
  units = "mm",
  dpi = 600,
  bg = "white"
)


pdf(
  file = file.path(
    out_dir,
    "Fig5_PanelB_nonProgrammeGSEA_Nature_wide_v5.pdf"
  ),
  width = 94 / 25.4,
  height = 56 / 25.4,
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
    "Fig5_PanelB_nonProgrammeGSEA_Nature_wide_v5.tiff"
  ),
  plot = p,
  width = 94,
  height = 56,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 13. Export
# ------------------------------------------------------------

readr::write_csv(
  rank_df,
  file.path(
    out_dir,
    "Fig5B_nonprogramme_gene_ranks.csv"
  )
)

readr::write_csv(
  fg_df,
  file.path(
    out_dir,
    "Fig5B_all_GO_BP_GSEA_results.csv"
  )
)

readr::write_csv(
  plot_df,
  file.path(
    out_dir,
    "Fig5B_plot_data.csv"
  )
)

stats_df <- tibble::tibble(
  Ranked_genes =
    length(
      ranks
    ),
  GO_BP_gene_sets_mapped =
    length(
      pathways
    ),
  Significant_pathways_FDR_lt_0_05 =
    n_sig,
  Native_axis_significant =
    n_native_sig,
  Alternative_axis_significant =
    n_alt_sig
)

readr::write_csv(
  stats_df,
  file.path(
    out_dir,
    "Fig5B_statistics.csv"
  )
)

qc_df <- tibble::tibble(
  Item = c(
    "subexpr genes",
    "subexpr samples",
    "Frozen background",
    "Native program",
    "Alternative program",
    "Held-out genes defined",
    "Held-out genes ranked",
    "GO-BP annotated ranked genes",
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
      heldout_genes
    ),
    length(
      ranks
    ),
    length(
      unique(
        go_map$SYMBOL
      )
    ),
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
    "1373",
    ">=1304",
    "report",
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
    "Fig5B_QC.csv"
  )
)

readme <- c(
  "Figure 5B — OFFLINE non-programme gene GO-BP preranked GSEA v5",
  "",
  "No msigdbr and no internet connection are used.",
  "",
  "Ranking:",
  "Non-programme genes only.",
  "Ranking metric = Spearman rho with anatomy-referenced state balance.",
  "",
  "Gene sets:",
  "Local GO Biological Process annotations from org.Hs.eg.db + GO.db.",
  "",
  "Interpretation:",
  "NES > 0 = Native-aligned end enriched",
  "NES < 0 = Alternative-aligned end enriched",
  "",
  paste0(
    "Ranked genes = ",
    length(
      ranks
    )
  ),
  paste0(
    "GO-BP pathways mapped = ",
    length(
      pathways
    )
  ),
  paste0(
    "Significant GO-BP pathways = ",
    n_sig
  ),
  paste0(
    "Native-aligned significant = ",
    n_native_sig
  ),
  paste0(
    "Alternative-aligned significant = ",
    n_alt_sig
  ),
  "",
  paste0(
    "Direction source for programme scoring = ",
    direction_source
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig5B_nonProgrammeGSEA_Nature_wide_v5.txt"
  )
)


cat(
  "\n============================================\n",
  "FIGURE 5B OFFLINE GSEA COMPLETE\n",
  "============================================\n\n",
  sep = ""
)

print(
  stats_df,
  n = Inf
)

cat(
  "\nTop Native-aligned GO-BP terms:\n"
)

print(
  top_native %>%
    dplyr::select(
      TERM,
      NES,
      padj
    ),
  n = Inf
)

cat(
  "\nTop Alternative-aligned GO-BP terms:\n"
)

print(
  top_alt %>%
    dplyr::select(
      TERM,
      NES,
      padj
    ),
  n = Inf
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)
