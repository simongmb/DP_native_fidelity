# ============================================================
# Supplementary Figure S4A v6 — clean-session reproducible
# Bootstrap robustness of the inverse two-programme relationship
#
# Goal:
# Test whether the strong inverse association observed in Fig. 4A
# is robust to subsampling genes from the two frozen programmes.
#
# Fixed input:
#   subexpr = 1553 genes x 110 samples
#
# Frozen programs:
#   Native-aligned rescue      = 81 genes
#   Alternative-state rescue   = 99 genes
#
# Bootstrap:
#   1000 iterations
#   sample 70% of genes from each programme without replacement
#   recalculate direction-aligned sample scores
#   calculate Spearman rho across samples
#
# Main output:
#   distribution of bootstrap Spearman rho
#   observed rho from full programmes shown as reference
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
  stop(
    "\nCannot identify project root.\n"
  )
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


format_p <- function(p) {

  if (is.na(p)) {
    return("NA")
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


if (
  length(native81) != 81
) {

  stop(
    "\nNative-aligned rescue frozen programme != 81 genes.\n"
  )
}


if (
  length(alternative99) != 99
) {

  stop(
    "\nAlternative-state rescue frozen programme != 99 genes.\n"
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
    "\nFrozen programmes unexpectedly overlap.\n"
  )
}


# ------------------------------------------------------------
# 5. Coverage QC
# ------------------------------------------------------------

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
# 6. Audit table + restoration direction
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


# categorical
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


# logical flag
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


# signed numeric effect
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
# 7. Gene-wise z scores
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


native_z_all <- z_by_gene(
  expr_mat[
    native81,
    ,
    drop = FALSE
  ]
)


alt_z_all <- z_by_gene(
  expr_mat[
    alternative99,
    ,
    drop = FALSE
  ]
)


native_sign_all <- native_direction$Direction_sign[
  match(
    native81,
    native_direction$GeneSymbol
  )
]


alt_sign_all <- alt_direction$Direction_sign[
  match(
    alternative99,
    alt_direction$GeneSymbol
  )
]


native_aligned_all <- sweep(
  native_z_all,
  1,
  native_sign_all,
  "*"
)


alt_aligned_all <- sweep(
  alt_z_all,
  1,
  alt_sign_all,
  "*"
)


# ------------------------------------------------------------
# 8. Observed full-programme relationship
# ------------------------------------------------------------

native_score_full <- colMeans(
  native_aligned_all,
  na.rm = TRUE
)


alt_score_full <- colMeans(
  alt_aligned_all,
  na.rm = TRUE
)


observed_cor <- suppressWarnings(
  stats::cor.test(
    native_score_full,
    alt_score_full,
    method = "spearman",
    exact = FALSE
  )
)


observed_rho <- unname(
  observed_cor$estimate
)


observed_p <- observed_cor$p.value


cat(
  "\nObserved full-programme Spearman rho =",
  sprintf(
    "%.4f",
    observed_rho
  ),
  "\n"
)


# ------------------------------------------------------------
# 9. Bootstrap settings
# ------------------------------------------------------------

N_BOOT <- 1000

FRACTION <- 0.70

N_NATIVE_SAMPLE <- floor(
  length(
    native81
  ) *
    FRACTION
)


N_ALT_SAMPLE <- floor(
  length(
    alternative99
  ) *
    FRACTION
)


if (
  N_NATIVE_SAMPLE < 2 ||
  N_ALT_SAMPLE < 2
) {

  stop(
    "\nBootstrap gene subset sizes are too small.\n"
  )
}


cat(
  "\nBootstrap settings:\n",
  "Iterations = ",
  N_BOOT,
  "\n",
  "Native genes per iteration = ",
  N_NATIVE_SAMPLE,
  "/81\n",
  "Alternative genes per iteration = ",
  N_ALT_SAMPLE,
  "/99\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Bootstrap
# ------------------------------------------------------------

set.seed(
  20261028
)


bootstrap_rho <- numeric(
  N_BOOT
)


for (
  b in seq_len(
    N_BOOT
  )
) {

  native_idx <- sample(
    seq_along(
      native81
    ),
    size = N_NATIVE_SAMPLE,
    replace = FALSE
  )


  alt_idx <- sample(
    seq_along(
      alternative99
    ),
    size = N_ALT_SAMPLE,
    replace = FALSE
  )


  native_score_b <- colMeans(
    native_aligned_all[
      native_idx,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )


  alt_score_b <- colMeans(
    alt_aligned_all[
      alt_idx,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )


  bootstrap_rho[b] <- suppressWarnings(
    stats::cor(
      native_score_b,
      alt_score_b,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
}


if (
  any(
    !is.finite(
      bootstrap_rho
    )
  )
) {

  stop(
    "\nOne or more bootstrap correlations are non-finite.\n"
  )
}


boot_df <- tibble::tibble(

  Iteration =
    seq_len(
      N_BOOT
    ),

  Spearman_rho =
    bootstrap_rho
)


# ------------------------------------------------------------
# 11. Bootstrap summary
# ------------------------------------------------------------

rho_median <- median(
  bootstrap_rho
)


rho_q025 <- as.numeric(
  quantile(
    bootstrap_rho,
    probs = 0.025,
    type = 7
  )
)


rho_q975 <- as.numeric(
  quantile(
    bootstrap_rho,
    probs = 0.975,
    type = 7
  )
)


fraction_negative <- mean(
  bootstrap_rho < 0
)


fraction_strong_negative <- mean(
  bootstrap_rho <= -0.5
)


fraction_le_observed <- mean(
  bootstrap_rho <= observed_rho
)


summary_df <- tibble::tibble(

  Metric = c(
    "Observed full-programme rho",
    "Bootstrap median rho",
    "Bootstrap 2.5% quantile",
    "Bootstrap 97.5% quantile",
    "Fraction rho < 0",
    "Fraction rho <= -0.5",
    "Fraction rho <= observed rho",
    "Bootstrap iterations",
    "Native genes sampled per iteration",
    "Alternative genes sampled per iteration"
  ),

  Value = c(
    observed_rho,
    rho_median,
    rho_q025,
    rho_q975,
    fraction_negative,
    fraction_strong_negative,
    fraction_le_observed,
    N_BOOT,
    N_NATIVE_SAMPLE,
    N_ALT_SAMPLE
  )
)


cat(
  "\nBootstrap summary:\n"
)


print(
  summary_df,
  n = Inf
)


# ------------------------------------------------------------
# 12. Plot
# ------------------------------------------------------------

COL_FILL <- "#91AAA9"
COL_LINE <- "#3D8C82"
COL_OBS <- "#C97857"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"


subtitle_text <- paste0(

  "Median rho = ",
  sprintf(
    "%.2f",
    rho_median
  ),

  "; 95% bootstrap interval ",
  sprintf(
    "%.2f",
    rho_q025
  ),

  " to ",
  sprintf(
    "%.2f",
    rho_q975
  ),

  "; ",
  sprintf(
    "%.1f",
    100 *
      fraction_negative
  ),

  "% of iterations < 0"
)


p <- ggplot(
  boot_df,
  aes(
    x = Spearman_rho
  )
) +

  geom_histogram(
    bins = 35,
    fill = COL_FILL,
    colour = "white",
    linewidth = 0.18,
    alpha = 0.90
  ) +

  geom_density(
    aes(
      y = after_stat(
        density *
        N_BOOT *
        diff(
          range(
            bootstrap_rho
          )
        ) /
        35
      )
    ),
    colour = COL_LINE,
    linewidth = 0.48,
    adjust = 1.0
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "22",
    linewidth = 0.28,
    colour = "#A8B5BA"
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
      "Full programmes\nrho = ",
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
      "Bootstrap Spearman rho",

    y =
      "Number of iterations",

    title =
      "Inverse programme association remains robust to gene subsampling",

    subtitle =
      subtitle_text,

    caption =
      paste0(
        N_BOOT,
        " iterations; each iteration sampled ",
        N_NATIVE_SAMPLE,
        "/81 Native-aligned rescue genes and ",
        N_ALT_SAMPLE,
        "/99 Alternative-state rescue genes\nwithout replacement."
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
        lineheight = 1.0,
        margin = margin(
          t = 4
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
  "S4_Panel_A_bootstrap_v5"
)


dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  filename = file.path(
    out_dir,
    "S4_Panel_A_bootstrap_v5.png"
  ),
  plot = p,
  width = 88,
  height = 66,
  units = "mm",
  dpi = 600,
  bg = "white"
)


pdf(
  file = file.path(
    out_dir,
    "S4_Panel_A_bootstrap_v5.pdf"
  ),
  width = 88 / 25.4,
  height = 66 / 25.4,
  useDingbats = FALSE
)

print(
  p
)

dev.off()


ggsave(
  filename = file.path(
    out_dir,
    "S4_Panel_A_bootstrap_v5.tiff"
  ),
  plot = p,
  width = 88,
  height = 66,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


readr::write_csv(
  boot_df,
  file.path(
    out_dir,
    "S4A_bootstrap_rho_values.csv"
  )
)


readr::write_csv(
  summary_df,
  file.path(
    out_dir,
    "S4A_bootstrap_summary.csv"
  )
)


qc_df <- tibble::tibble(

  Item = c(
    "subexpr genes",
    "subexpr samples",
    "Native-aligned rescue genes",
    "Alternative-state rescue genes",
    "Native restored-up",
    "Native restored-down",
    "Alternative restored-up",
    "Alternative restored-down",
    "Bootstrap iterations",
    "Native genes sampled",
    "Alternative genes sampled"
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
    alt_down,
    N_BOOT,
    N_NATIVE_SAMPLE,
    N_ALT_SAMPLE
  ),

  Expected = c(
    "1553",
    "110",
    "81",
    "99",
    "55",
    "26",
    "30",
    "69",
    "1000",
    as.character(
      N_NATIVE_SAMPLE
    ),
    as.character(
      N_ALT_SAMPLE
    )
  )
)


readr::write_csv(
  qc_df,
  file.path(
    out_dir,
    "S4A_QC.csv"
  )
)


readme <- c(

  "Supplementary Figure S4A bootstrap robustness v6",

  "",

  "Purpose:",
  paste0(
    "Test whether the inverse sample-level association between ",
    "the frozen 81-gene and 99-gene programmes is robust to ",
    "subsampling genes."
  ),

  "",

  paste0(
    "Observed full-programme Spearman rho = ",
    sprintf(
      "%.4f",
      observed_rho
    )
  ),

  paste0(
    "Observed full-programme P = ",
    format_p(
      observed_p
    )
  ),

  "",

  paste0(
    "Bootstrap iterations = ",
    N_BOOT
  ),

  paste0(
    "Native genes sampled each iteration = ",
    N_NATIVE_SAMPLE,
    "/81"
  ),

  paste0(
    "Alternative genes sampled each iteration = ",
    N_ALT_SAMPLE,
    "/99"
  ),

  "",

  paste0(
    "Bootstrap median rho = ",
    sprintf(
      "%.4f",
      rho_median
    )
  ),

  paste0(
    "95% bootstrap interval = ",
    sprintf(
      "%.4f",
      rho_q025
    ),
    " to ",
    sprintf(
      "%.4f",
      rho_q975
    )
  ),

  paste0(
    "Fraction of iterations with rho < 0 = ",
    sprintf(
      "%.4f",
      fraction_negative
    )
  ),

  paste0(
    "Fraction of iterations with rho <= -0.5 = ",
    sprintf(
      "%.4f",
      fraction_strong_negative
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
    "README_S4A_bootstrap_v6.txt"
  )
)


cat(
  "\n============================================\n",
  "FIGURE 4B BOOTSTRAP COMPLETE\n",
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
