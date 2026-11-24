# ============================================================
# Supplementary Figure S4B v5
# Matched random-programme permutation test
#
# Question:
# Is the observed inverse association between the frozen
# 81-gene Native-aligned rescue programme and 99-gene Alternative-state rescue programme
# stronger than expected for random programmes of the same size
# and same restored-up/restored-down composition?
#
# Observed:
#   Native-aligned rescue programme      = 81 genes = 55 up + 26 down
#   Alternative-state rescue programme = 99 genes = 30 up + 69 down
#
# Null model:
#   From the frozen 1,553-gene analyzable background,
#   repeatedly draw two DISJOINT random programmes while matching
#   the exact up/down composition above.
#
# For each permutation:
#   1. gene-wise z-score across samples
#   2. direction alignment
#   3. sample-level programme scores
#   4. Spearman rho between the two random programme scores
#
# Permutations = 10,000
#
# ============================================================

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

cat("\nProject root:\n", project_root, "\n")


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
    candidates %in% names(df)
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


# ------------------------------------------------------------
# 4. Frozen gene sets
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
  stop("\nNative-aligned rescue frozen programme != 81 genes.\n")
}

if (length(alternative99) != 99) {
  stop("\nAlternative-state rescue frozen programme != 99 genes.\n")
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
  stop("\nCannot find 02_gene_level_native_fidelity_audit.csv\n")
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
    Direction_sign %in% c(-1, 1)
  ) %>%

  dplyr::distinct(
    GeneSymbol,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 6. Eligible background
#
# Must be:
#   - frozen background
#   - present in subexpr
#   - have valid direction sign
# ------------------------------------------------------------

eligible_tbl <- tibble::tibble(
  GeneSymbol =
    background1553
) %>%

  dplyr::left_join(
    direction_tbl,
    by = "GeneSymbol"
  ) %>%

  dplyr::filter(
    GeneSymbol %in%
      rownames(expr_mat),
    is.finite(
      Direction_sign
    )
  )

eligible_up <- eligible_tbl %>%
  dplyr::filter(
    Direction_sign == 1
  ) %>%
  dplyr::pull(
    GeneSymbol
  )

eligible_down <- eligible_tbl %>%
  dplyr::filter(
    Direction_sign == -1
  ) %>%
  dplyr::pull(
    GeneSymbol
  )

cat(
  "\nEligible background:\n",
  "Total = ",
  nrow(eligible_tbl),
  "\n",
  "Up = ",
  length(eligible_up),
  "\n",
  "Down = ",
  length(eligible_down),
  "\n",
  sep = ""
)

if (
  length(eligible_up) <
    (55 + 30) ||
  length(eligible_down) <
    (26 + 69)
) {
  stop(
    "\nInsufficient eligible genes for matched disjoint random programmes.\n"
  )
}


# ------------------------------------------------------------
# 7. Frozen direction QC
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
# 8. Gene-wise z scores for all eligible genes
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

eligible_genes <- eligible_tbl$GeneSymbol

z_all <- z_by_gene(
  expr_mat[
    eligible_genes,
    ,
    drop = FALSE
  ]
)

sign_all <- eligible_tbl$Direction_sign[
  match(
    rownames(z_all),
    eligible_tbl$GeneSymbol
  )
]

aligned_all <- sweep(
  z_all,
  1,
  sign_all,
  "*"
)


# ------------------------------------------------------------
# 9. Observed full-programme rho
# ------------------------------------------------------------

native_score_obs <- colMeans(
  aligned_all[
    native81,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)

alt_score_obs <- colMeans(
  aligned_all[
    alternative99,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)

observed_rho <- suppressWarnings(
  stats::cor(
    native_score_obs,
    alt_score_obs,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
)

cat(
  "\nObserved rho =",
  sprintf(
    "%.4f",
    observed_rho
  ),
  "\n"
)


# ------------------------------------------------------------
# 10. Matched random-programme permutations
# ------------------------------------------------------------

N_PERM <- 10000

set.seed(
  20261028
)

perm_rho <- numeric(
  N_PERM
)


for (
  i in seq_len(
    N_PERM
  )
) {

  # Native random set:
  # 55 up + 26 down
  native_up_rand <- sample(
    eligible_up,
    size = 55,
    replace = FALSE
  )

  native_down_rand <- sample(
    eligible_down,
    size = 26,
    replace = FALSE
  )

  native_rand <- c(
    native_up_rand,
    native_down_rand
  )


  # Alternative must be disjoint from native random set.
  remaining_up <- setdiff(
    eligible_up,
    native_up_rand
  )

  remaining_down <- setdiff(
    eligible_down,
    native_down_rand
  )

  alt_up_rand <- sample(
    remaining_up,
    size = 30,
    replace = FALSE
  )

  alt_down_rand <- sample(
    remaining_down,
    size = 69,
    replace = FALSE
  )

  alt_rand <- c(
    alt_up_rand,
    alt_down_rand
  )


  native_score_rand <- colMeans(
    aligned_all[
      native_rand,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  alt_score_rand <- colMeans(
    aligned_all[
      alt_rand,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )


  perm_rho[i] <- suppressWarnings(
    stats::cor(
      native_score_rand,
      alt_score_rand,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
}


if (
  any(!is.finite(perm_rho))
) {
  stop(
    "\nOne or more permutation rho values are non-finite.\n"
  )
}


perm_df <- tibble::tibble(
  Iteration =
    seq_len(
      N_PERM
    ),
  Spearman_rho =
    perm_rho
)


# ------------------------------------------------------------
# 11. Empirical P
#
# One-sided:
# probability under matched null of rho <= observed rho.
#
# +1 correction prevents zero P.
# ------------------------------------------------------------

empirical_p <- (
  sum(
    perm_rho <= observed_rho
  ) + 1
) / (
  N_PERM + 1
)

null_median <- median(
  perm_rho
)

null_q025 <- as.numeric(
  quantile(
    perm_rho,
    0.025
  )
)

null_q975 <- as.numeric(
  quantile(
    perm_rho,
    0.975
  )
)

null_min <- min(
  perm_rho
)

null_max <- max(
  perm_rho
)

cat(
  "\nPermutation summary:\n",
  "Null median rho = ",
  sprintf("%.4f", null_median),
  "\n",
  "Null 95% interval = ",
  sprintf("%.4f", null_q025),
  " to ",
  sprintf("%.4f", null_q975),
  "\n",
  "Observed rho = ",
  sprintf("%.4f", observed_rho),
  "\n",
  "Empirical P = ",
  format(
    empirical_p,
    scientific = TRUE,
    digits = 3
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Plot
# ------------------------------------------------------------

COL_NULL <- "#A7B7BD"
COL_DENSITY <- "#6F8C95"
COL_OBS <- "#C97857"
COL_ZERO <- "#A8B5BA"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"

subtitle_text <- paste0(
  "Matched null median rho = ",
  sprintf(
    "%.2f",
    null_median
  ),
  "; observed rho = ",
  sprintf(
    "%.2f",
    observed_rho
  ),
  "; empirical P = ",
  format(
    empirical_p,
    scientific = TRUE,
    digits = 2
  )
)


p <- ggplot(
  perm_df,
  aes(
    x = Spearman_rho
  )
) +

  geom_histogram(
    bins = 45,
    fill = COL_NULL,
    colour = "white",
    linewidth = 0.18,
    alpha = 0.92
  ) +

  geom_density(
    aes(
      y = after_stat(
        density *
        N_PERM *
        diff(
          range(
            perm_rho
          )
        ) /
        45
      )
    ),
    colour = COL_DENSITY,
    linewidth = 0.48,
    adjust = 1.0
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "22",
    linewidth = 0.28,
    colour = COL_ZERO
  ) +

  geom_vline(
    xintercept = observed_rho,
    linewidth = 0.62,
    colour = COL_OBS
  ) +

  annotate(
    "text",
    x = observed_rho,
    y = Inf,
    label = paste0(
      "Observed\nrho = ",
      sprintf(
        "%.2f",
        observed_rho
      )
    ),
    hjust = -0.06,
    vjust = 1.20,
    size = 4.8 / ggplot2::.pt,
    colour = COL_OBS
  ) +

  labs(
    x =
      "Spearman rho under matched random programmes",

    y =
      "Number of permutations",

    title =
      "Observed programme opposition exceeds matched random expectation",

    subtitle =
      subtitle_text,

    caption =
      paste0(
        N_PERM,
        " permutations. Random programmes matched the exact sizes and\n",
        "restored-up/restored-down composition of the frozen 81- and 99-gene programmes."
      )
  ) +

  theme_classic(
    base_size = 6
  ) +

  theme(
    plot.title =
      element_text(
        size = 6.2,
        face = "bold",
        colour = COL_TEXT,
        hjust = 0.5,
        margin = margin(
          b = 2.5
        )
      ),

    plot.title.position = "plot",

    plot.subtitle =
      element_text(
        size = 4.8,
        colour = COL_MUTED,
        hjust = 0.5,
        margin = margin(
          b = 4
        )
      ),

    plot.caption =
      element_text(
        size = 4.0,
        colour = "#81939B",
        hjust = 0.5,
        lineheight = 0.95,
        margin = margin(
          t = 3
        )
      ),

    axis.title =
      element_text(
        size = 5.6,
        colour = COL_TEXT
      ),

    axis.text =
      element_text(
        size = 4.8,
        colour = COL_MUTED
      ),

    plot.margin =
      margin(
        5,
        5,
        9,
        5,
        unit = "pt"
      )
  )


# ------------------------------------------------------------
# 13. Output
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "figures",
  "SUPPLEMENTARY FIGURES",
  "S4_Panel_B_matchedNull_v4"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    out_dir,
    "S4_Panel_B_matchedNull_v4.png"
  ),
  plot = p,
  width = 88,
  height = 62,
  units = "mm",
  dpi = 600,
  bg = "white"
)

pdf(
  file = file.path(
    out_dir,
    "S4_Panel_B_matchedNull_v4.pdf"
  ),
  width = 88 / 25.4,
  height = 62 / 25.4,
  useDingbats = FALSE
)

print(p)
dev.off()

ggsave(
  filename = file.path(
    out_dir,
    "S4_Panel_B_matchedNull_v4.tiff"
  ),
  plot = p,
  width = 88,
  height = 62,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


readr::write_csv(
  perm_df,
  file.path(
    out_dir,
    "S4B_matched_null_rho.csv"
  )
)

summary_df <- tibble::tibble(
  Metric = c(
    "Observed rho",
    "Null median rho",
    "Null 2.5% quantile",
    "Null 97.5% quantile",
    "Null minimum rho",
    "Null maximum rho",
    "Empirical one-sided P",
    "Permutations",
    "Eligible background genes",
    "Eligible restored-up genes",
    "Eligible restored-down genes"
  ),
  Value = c(
    observed_rho,
    null_median,
    null_q025,
    null_q975,
    null_min,
    null_max,
    empirical_p,
    N_PERM,
    nrow(eligible_tbl),
    length(eligible_up),
    length(eligible_down)
  )
)

readr::write_csv(
  summary_df,
  file.path(
    out_dir,
    "S4B_matched_null_summary.csv"
  )
)

qc_df <- tibble::tibble(
  Item = c(
    "Frozen background genes",
    "subexpr genes",
    "subexpr samples",
    "Eligible background genes",
    "Native-aligned rescue genes",
    "Alternative-state rescue genes",
    "Native restored-up",
    "Native restored-down",
    "Alternative restored-up",
    "Alternative restored-down",
    "Permutations"
  ),
  Value = c(
    length(background1553),
    nrow(expr_mat),
    ncol(expr_mat),
    nrow(eligible_tbl),
    length(native81),
    length(alternative99),
    native_up,
    native_down,
    alt_up,
    alt_down,
    N_PERM
  ),
  Expected = c(
    "1553",
    "1553",
    "110",
    "report",
    "81",
    "99",
    "55",
    "26",
    "30",
    "69",
    "10000"
  )
)

readr::write_csv(
  qc_df,
  file.path(
    out_dir,
    "S4B_QC.csv"
  )
)

readme <- c(
  "Supplementary Figure S4B matched random-programme permutation test v5",
  "",
  "Question:",
  "Is the observed inverse association stronger than expected for random programmes?",
  "",
  "Null model:",
  "Two disjoint random programmes sampled from the frozen analyzable background.",
  "Programme sizes and restored-up/restored-down composition are matched exactly.",
  "",
  "Native-aligned rescue null programme: 81 genes = 55 up + 26 down",
  "Alternative-state rescue null programme: 99 genes = 30 up + 69 down",
  "",
  paste0(
    "Observed rho = ",
    sprintf("%.4f", observed_rho)
  ),
  paste0(
    "Null median rho = ",
    sprintf("%.4f", null_median)
  ),
  paste0(
    "Null 95% interval = ",
    sprintf("%.4f", null_q025),
    " to ",
    sprintf("%.4f", null_q975)
  ),
  paste0(
    "Empirical one-sided P = ",
    format(
      empirical_p,
      scientific = TRUE,
      digits = 4
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
    "README_S4B_matchedNull_v5.txt"
  )
)

cat(
  "\n============================================\n",
  "FIGURE 4D COMPLETE\n",
  "============================================\n\n",
  sep = ""
)

print(
  summary_df,
  n = Inf
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)
