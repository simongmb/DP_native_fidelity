# ============================================================
# Figure 4A v7
# Sample-level relationship between the two frozen programmes
#
# Terminology-aligned revision:
# - uses fixed input: subexpr
# - frozen 81 / 99 programs unchanged
# - direction-aligned scoring unchanged
# - removes useless "Other samples" legend
# - adds Spearman rho and P value
# - updates display terminology to Native-aligned / Alternative-state programmes
# - keeps all analysis, scoring, donor-aware inference and data unchanged
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

candidate_roots <- unique(c(
  start_dir,
  getwd(),
  dirname(start_dir),
  dirname(dirname(start_dir)),
  dirname(getwd()),
  dirname(dirname(getwd()))
))

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
      "\nCannot find clean-session loader:\n",
      loader,
      "\n"
    )
  )
}

source(loader, local = FALSE)

if (!exists("subexpr", inherits = TRUE)) {
  stop("\nQC STOP: load_subexpr_clean.R did not create subexpr.\n")
}

expr_mat <- get("subexpr", inherits = TRUE)

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
      "\nQC STOP: expected 1553 genes in clean-session subexpr, found ",
      nrow(expr_mat),
      ".\n"
    )
  )
}

if (ncol(expr_mat) != 110L) {
  stop(
    paste0(
      "\nQC STOP: expected 110 profiles in clean-session subexpr, found ",
      ncol(expr_mat),
      ".\n"
    )
  )
}

cat(
  "\nClean-session subexpr dimensions:",
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

pick_col <- function(df, candidates) {
  hit <- candidates[
    candidates %in% names(df)
  ]

  if (length(hit) == 0) {
    return(NA_character_)
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

  sprintf("%.3f", p)
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

if (length(native81) != 81) {
  stop("\nNative frozen program != 81 genes.\n")
}

if (length(alternative99) != 99) {
  stop("\nAlternative frozen program != 99 genes.\n")
}


# ------------------------------------------------------------
# 5. Coverage QC
# ------------------------------------------------------------

native_present <- intersect(
  native81,
  rownames(expr_mat)
)

alt_present <- intersect(
  alternative99,
  rownames(expr_mat)
)

if (length(native_present) != 81) {
  stop("\nsubexpr no longer contains all 81 native genes.\n")
}

if (length(alt_present) != 99) {
  stop("\nsubexpr no longer contains all 99 alternative genes.\n")
}


# ------------------------------------------------------------
# 6. Audit table and direction recovery
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

  direction_sign <- sign(effect)

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
  GeneSymbol = toupper(
    stringr::str_trim(
      as.character(
        audit0[[gene_col]]
      )
    )
  ),
  Direction_sign = direction_sign
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
# 7. Frozen direction QC
# ------------------------------------------------------------

native_direction <- tibble::tibble(
  GeneSymbol = native81
) %>%
  dplyr::left_join(
    direction_tbl,
    by = "GeneSymbol"
  )

alt_direction <- tibble::tibble(
  GeneSymbol = alternative99
) %>%
  dplyr::left_join(
    direction_tbl,
    by = "GeneSymbol"
  )

if (
  any(!is.finite(native_direction$Direction_sign)) ||
  any(!is.finite(alt_direction$Direction_sign))
) {
  stop("\nDirection missing for frozen genes.\n")
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
# 8. Gene-wise z score
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


# ------------------------------------------------------------
# 9. Direction alignment
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
# 10. Program scores
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
  Sample = sample_ids,
  Native_identity_score = as.numeric(
    native_score[
      sample_ids
    ]
  ),
  Alternative_state_score = as.numeric(
    alternative_score[
      sample_ids
    ]
  )
)


# ------------------------------------------------------------
# 11. Spearman association + donor-aware inference
# ------------------------------------------------------------

cor_res <- suppressWarnings(
  stats::cor.test(
    score_df$Native_identity_score,
    score_df$Alternative_state_score,
    method = "spearman",
    exact = FALSE
  )
)

rho <- unname(
  cor_res$estimate
)

# Nominal sample-level P is retained only for audit/QC.
# It is NOT used as the primary inferential P value because the
# 110 anatomical profiles include multiple profiles per donor.
cor_p_nominal <- cor_res$p.value

# Recover donor IDs from sample names using the same convention as
# the donor-blocked Figure 5/6 analyses.
extract_donor <- function(ids) {
  out <- stringr::str_match(
    ids,
    "^.*_([A-J])_[^_]+$"
  )[, 2]

  # Fallback for IDs formatted as COMPARTMENT_DONOR_...
  if (anyNA(out)) {
    parts <- strsplit(ids, "_", fixed = TRUE)
    if (all(lengths(parts) >= 2L)) {
      candidate <- vapply(parts, `[`, character(1), 2)
      if (length(unique(candidate)) == 10L) {
        out <- candidate
      }
    }
  }

  out
}

donor_ids <- extract_donor(score_df$Sample)

if (
  anyNA(donor_ids) ||
  length(unique(donor_ids)) != 10L
) {
  stop(
    "\nQC STOP: could not recover exactly 10 donor IDs for donor-aware permutation.\n"
  )
}

permute_within_donor <- function(x, donor) {
  out <- x
  for (d in unique(donor)) {
    idx <- which(donor == d)
    out[idx] <- sample(
      x[idx],
      length(idx),
      replace = FALSE
    )
  }
  out
}

N_PERM_DONOR <- 5000L
DONOR_PERM_SEED <- 20261102L

set.seed(DONOR_PERM_SEED)

extreme_count <- 0L

for (b in seq_len(N_PERM_DONOR)) {
  alt_perm <- permute_within_donor(
    score_df$Alternative_state_score,
    donor_ids
  )

  rho_perm <- suppressWarnings(
    cor(
      score_df$Native_identity_score,
      alt_perm,
      method = "spearman"
    )
  )

  if (
    is.finite(rho_perm) &&
    abs(rho_perm) >= abs(rho)
  ) {
    extreme_count <- extreme_count + 1L
  }
}

donor_perm_p <- (
  extreme_count + 1L
) / (
  N_PERM_DONOR + 1L
)

stat_label <- paste0(
  "Spearman rho = ",
  sprintf("%.2f", rho),
  "   within-donor P = ",
  format_p(donor_perm_p)
)

cat(
  "\n",
  stat_label,
  "\n",
  "Exact rho = ",
  sprintf("%.6f", rho),
  "\nNominal sample-level P (audit only) = ",
  format_p(cor_p_nominal),
  "\nWithin-donor permutation P = ",
  format_p(donor_perm_p),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Quadrant counts
# ------------------------------------------------------------

score_df <- score_df %>%
  dplyr::mutate(
    State = dplyr::case_when(
      Native_identity_score >= 0 &
      Alternative_state_score >= 0 ~
        "Dual-positive",

      Native_identity_score >= 0 &
      Alternative_state_score < 0 ~
        "Native-dominant",

      Native_identity_score < 0 &
      Alternative_state_score >= 0 ~
        "Alternative-state dominant",

      TRUE ~
        "Dual-low"
    )
  )

state_summary <- score_df %>%
  dplyr::count(
    State,
    name = "n"
  ) %>%
  dplyr::mutate(
    Fraction = n / sum(n)
  )


# ------------------------------------------------------------
# 13. Plot — Nature-style wide / short layout
# ------------------------------------------------------------

COL_POINT <- "#758D96"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"
COL_NATIVE <- "#3D8C82"
COL_ALT <- "#C97857"
COL_ZERO <- "#D5DEE2"

p <- ggplot(
  score_df,
  aes(
    x = Native_identity_score,
    y = Alternative_state_score
  )
) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.28,
    colour = COL_ZERO
  ) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.28,
    colour = COL_ZERO
  ) +

  geom_point(
    size = 1.65,
    alpha = 0.72,
    colour = COL_POINT
  ) +

  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = stat_label,
    hjust = 1.08,
    vjust = 1.35,
    size = 1.72,
    family = "Helvetica",
    colour = COL_MUTED
  ) +

  labs(
    x = "Native-aligned programme score",
    y = "Alternative-state programme score",
    title = "Native-aligned and Alternative-state scores are inversely associated",
    caption = paste0(
      "Input: subexpr (",
      nrow(expr_mat),
      " genes x ",
      ncol(expr_mat),
      " samples); gene-wise values z-standardized and direction-aligned."
    )
  ) +

  theme_classic(
    base_family = "Helvetica",
    base_size = 7
  ) +

  theme(
    plot.title = element_text(
      size = 5.3,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      margin = margin(b = 1.0)
    ),

    plot.caption = element_text(
      size = 4.2,
      colour = "#8FA0A7",
      hjust = 0.5,
      margin = margin(t = 3)
    ),

    axis.title.x = element_text(
      size = 5.4,
      colour = COL_TEXT,
      margin = margin(t = 3)
    ),

    axis.title.y = element_text(
      size = 5.4,
      colour = COL_TEXT,
      margin = margin(r = 3)
    ),

    axis.text = element_text(
      size = 4.8,
      colour = COL_MUTED
    ),

    axis.line = element_line(
      linewidth = 0.32,
      colour = COL_TEXT
    ),

    axis.ticks = element_line(
      linewidth = 0.28,
      colour = COL_TEXT
    ),

    axis.ticks.length = unit(
      1.2,
      "mm"
    ),

    plot.margin = margin(
      t = 3,
      r = 5,
      b = 4,
      l = 5,
      unit = "pt"
    )
  )


# ------------------------------------------------------------
# 13. Output
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 4",
  "Panel_A_sampleProgrammes_Nature_wide_v7"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    out_dir,
    "Fig4_PanelA_Nature_wide_v7.png"
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
    "Fig4_PanelA_Nature_wide_v7.pdf"
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
    "Fig4_PanelA_Nature_wide_v7.tiff"
  ),
  plot = p,
  width = 88,
  height = 48,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

readr::write_csv(
  score_df,
  file.path(
    out_dir,
    "Fig4A_sample_program_scores.csv"
  )
)


stats_df <- tibble::tibble(
  N_samples = nrow(score_df),
  Spearman_rho = rho,
  Spearman_P_nominal = cor_p_nominal,
  Within_donor_permutation_P = donor_perm_p,
  Within_donor_permutations = N_PERM_DONOR,
  Native_program_n = 81,
  Alternative_program_n = 99,
  Native_restored_up = native_up,
  Native_restored_down = native_down,
  Alternative_restored_up = alt_up,
  Alternative_restored_down = alt_down
)

readr::write_csv(
  stats_df,
  file.path(
    out_dir,
    "Fig4A_statistics.csv"
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
    nrow(expr_mat),
    ncol(expr_mat),
    length(native_present),
    length(alt_present),
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
    "Fig4A_QC.csv"
  )
)

readme <- c(
  "Figure 4A sample-level frozen-programme relationship v6",
  "",
  paste0(
    "Input: subexpr = ",
    nrow(expr_mat),
    " genes x ",
    ncol(expr_mat),
    " samples"
  ),
  "",
  "Frozen programmes:",
  "Native-aligned rescue = 81 genes",
  "Alternative-state rescue = 99 genes",
  "",
  paste0(
    "Direction source: ",
    direction_source
  ),
  "",
  "Scoring:",
  "1. gene-wise z-score across samples",
  "2. restored-up genes multiplied by +1",
  "3. restored-down genes multiplied by -1",
  "4. mean aligned z-score per frozen programme",
  "",
  paste0(
    "Spearman rho = ",
    sprintf("%.4f", rho)
  ),
  paste0(
    "Nominal sample-level Spearman P (audit only) = ",
    format_p(cor_p_nominal)
  ),
  paste0(
    "Within-donor permutation P = ",
    format_p(donor_perm_p)
  ),
  paste0(
    "Within-donor permutations = ",
    N_PERM_DONOR
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig4A_sampleProgrammes_v7.txt"
  )
)

cat(
  "\n============================================\n",
  "FIGURE 4A v7 COMPLETE\n",
  "============================================\n\n",
  sep = ""
)

print(
  stats_df,
  n = Inf
)

cat(
  "\nState summary:\n"
)

print(
  state_summary,
  n = Inf
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)
