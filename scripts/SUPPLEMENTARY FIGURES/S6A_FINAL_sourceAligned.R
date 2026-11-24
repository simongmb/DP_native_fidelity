# ============================================================
# Supplementary Figure S6A v6 — compact descriptive version
#
# Purpose:
# Descriptive comparison of the full 1,373-gene non-programme model and the
# nested top-25 reconstruction across the same 20 repeated
# donor-blocked nested-CV runs.
#
# IMPORTANT:
# - No Wilcoxon test is performed.
# - No repeat-level inferential P value is generated.
# - Paired lines are retained only to visualize stability and
#   within-repeat consistency.
#
# Visual revision:
# - Slightly narrower panel
# - Larger fonts
# - Reduced inter-panel spacing
# - Reduced outer margins
# - Slightly larger repeat points and median diamonds
#
# Outputs:
#   Supplementary_Fig6A_descriptive.png
#   Supplementary_Fig6A_descriptive.pdf
#   Supplementary_Fig6A_descriptive.tiff
#   S6A_repeat_descriptive_summary.csv
#   S6A_within_repeat_differences_descriptive.csv
#   S6A_plot_source_data.csv
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr"
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
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})


# ------------------------------------------------------------
# 1. Constants
# ------------------------------------------------------------

EXPECTED_REPEATS <- 20L

FULL_MODEL_NAME <-
  "All 1,373 non-programme genes"

TOP25_MODEL_NAME <-
  "Nested top-25 genes"


# ------------------------------------------------------------
# 2. Locate project root
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
  
  if (length(file_arg) > 0L) {
    
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


find_project_root <- function(
    start = get_script_dir()) {
  
  start <- normalizePath(
    start,
    winslash = "/",
    mustWork = TRUE
  )
  
  candidates <- unique(
    c(
      start,
      getwd(),
      dirname(start),
      dirname(
        dirname(start)
      ),
      dirname(
        dirname(
          dirname(start)
        )
      ),
      dirname(
        getwd()
      ),
      dirname(
        dirname(
          getwd()
        )
      )
    )
  )
  
  candidates <- candidates[
    dir.exists(candidates)
  ]
  
  for (x in candidates) {
    
    if (
      dir.exists(
        file.path(
          x,
          "figures"
        )
      ) &&
      dir.exists(
        file.path(
          x,
          "scripts"
        )
      )
    ) {
      
      return(x)
    }
  }
  
  stop(
    paste0(
      "\nCannot identify project root ",
      "containing figures/ and scripts/.\n"
    )
  )
}


project_root <- find_project_root()

cat(
  "\nProject root:\n",
  project_root,
  "\n"
)


# ------------------------------------------------------------
# 3. Locate Figure 6D repeat-level metrics
# ------------------------------------------------------------

candidate_files <- list.files(
  project_root,
  pattern =
    "^Fig6D_v[0-9]+_repeat_metrics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(candidate_files) == 0L) {
  
  stop(
    paste0(
      "\nQC STOP\n",
      "Cannot find ",
      "a finalized Fig6D repeat-metrics file.\n",
      "Run the final Figure 6 ",
      "top-25 script first.\n"
    )
  )
}


score_file <- function(path) {
  
  s <- 0
  
  if (
    grepl(
      "Panel_D_nestedTop25_donorBlocked",
      path,
      ignore.case = TRUE
    )
  ) {
    
    s <- s + 100
  }
  
  if (
    grepl(
      "FIGURE 6",
      path,
      ignore.case = TRUE
    )
  ) {
    s <- s + 20
  }

  ver <- suppressWarnings(
    as.integer(
      stringr::str_match(
        basename(path),
        "^Fig6D_v([0-9]+)_repeat_metrics\\.csv$"
      )[, 2]
    )
  )
  if (is.finite(ver)) s <- s + ver

  s
}


scores <- vapply(
  candidate_files,
  score_file,
  numeric(1)
)

metrics_file <- candidate_files[
  which.max(scores)
]

cat(
  "\nRepeat metrics source:\n",
  metrics_file,
  "\n"
)


# ------------------------------------------------------------
# 4. Read repeat-level data
# ------------------------------------------------------------

dat <- readr::read_csv(
  metrics_file,
  show_col_types = FALSE,
  progress = FALSE
)


required_cols <- c(
  "Model",
  "Repeat",
  "Spearman_rho",
  "CV_R2"
)

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0L) {
  
  stop(
    paste0(
      "\nQC STOP\n",
      "Missing required columns:\n",
      paste(
        missing_cols,
        collapse = ", "
      ),
      "\n"
    )
  )
}


dat <- dat %>%
  dplyr::filter(
    Model %in%
      c(
        FULL_MODEL_NAME,
        TOP25_MODEL_NAME
      )
  ) %>%
  dplyr::mutate(
    
    Repeat =
      as.integer(
        Repeat
      ),
    
    Model = factor(
      Model,
      levels = c(
        FULL_MODEL_NAME,
        TOP25_MODEL_NAME
      ),
      labels = c(
        "Full 1,373-gene\nnon-programme model",
        "Nested top-25\nreconstruction"
      )
    )
  )


# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------

qc_counts <- dat %>%
  dplyr::count(
    Model,
    name = "n"
  )

if (
  nrow(qc_counts) != 2L ||
  any(
    qc_counts$n !=
    EXPECTED_REPEATS
  )
) {
  
  stop(
    paste0(
      "\nQC STOP\n",
      "Expected 20 repeat-level rows ",
      "per model.\n\n",
      paste(
        capture.output(
          print(
            qc_counts
          )
        ),
        collapse = "\n"
      ),
      "\n"
    )
  )
}


if (
  dplyr::n_distinct(
    dat$Repeat
  ) !=
  EXPECTED_REPEATS
) {
  
  stop(
    paste0(
      "\nQC STOP\n",
      "Expected exactly 20 ",
      "repeat identifiers.\n"
    )
  )
}


if (
  any(
    !is.finite(
      dat$Spearman_rho
    )
  ) ||
  any(
    !is.finite(
      dat$CV_R2
    )
  )
) {
  
  stop(
    paste0(
      "\nQC STOP\n",
      "Non-finite repeat-level ",
      "performance values detected.\n"
    )
  )
}


cat(
  "\nQC passed:\n"
)

print(
  qc_counts,
  n = Inf
)


# ------------------------------------------------------------
# 6. Descriptive summaries only
# ------------------------------------------------------------

summary_tbl <- dat %>%
  dplyr::group_by(
    Model
  ) %>%
  dplyr::summarise(
    
    n_repeats =
      dplyr::n(),
    
    median_Spearman_rho =
      median(
        Spearman_rho,
        na.rm = TRUE
      ),
    
    min_Spearman_rho =
      min(
        Spearman_rho,
        na.rm = TRUE
      ),
    
    max_Spearman_rho =
      max(
        Spearman_rho,
        na.rm = TRUE
      ),
    
    median_CV_R2 =
      median(
        CV_R2,
        na.rm = TRUE
      ),
    
    min_CV_R2 =
      min(
        CV_R2,
        na.rm = TRUE
      ),
    
    max_CV_R2 =
      max(
        CV_R2,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 7. Within-repeat descriptive differences
# ------------------------------------------------------------

wide <- dat %>%
  dplyr::select(
    Model,
    Repeat,
    Spearman_rho,
    CV_R2
  ) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = c(
      Spearman_rho,
      CV_R2
    )
  )


diff_tbl <- tibble::tibble(
  
  Repeat =
    wide$Repeat,
  
  Spearman_difference_top25_minus_full =
    wide$`Spearman_rho_Nested top-25\nreconstruction` -
    wide$`Spearman_rho_Full 1,373-gene\nnon-programme model`,
  
  CV_R2_difference_top25_minus_full =
    wide$`CV_R2_Nested top-25\nreconstruction` -
    wide$`CV_R2_Full 1,373-gene\nnon-programme model`
)


# ------------------------------------------------------------
# 8. Long-format plotting data
# ------------------------------------------------------------

plot_dat <- dat %>%
  dplyr::select(
    Model,
    Repeat,
    Spearman_rho,
    CV_R2
  ) %>%
  tidyr::pivot_longer(
    
    cols = c(
      Spearman_rho,
      CV_R2
    ),
    
    names_to =
      "Metric",
    
    values_to =
      "Value"
  ) %>%
  dplyr::mutate(
    
    Metric =
      dplyr::recode(
        Metric,
        Spearman_rho =
          "Spearman rho",
        CV_R2 =
          "CV R2"
      ),
    
    Metric = factor(
      Metric,
      levels = c(
        "Spearman rho",
        "CV R2"
      )
    )
  )


median_dat <- plot_dat %>%
  dplyr::group_by(
    Metric,
    Model
  ) %>%
  dplyr::summarise(
    
    Value =
      median(
        Value,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 9. Colours
# ------------------------------------------------------------

COL_TEXT <-
  "#20343D"

COL_MUTED <-
  "#718993"

COL_LINE <-
  "#91A5AD"

COL_POINT <-
  "#456B78"


# ------------------------------------------------------------
# 10. Plot
# ------------------------------------------------------------

p <- ggplot(
  plot_dat,
  aes(
    x = Model,
    y = Value,
    group =
      interaction(
        Repeat,
        Metric
      )
  )
) +
  
  # paired repeat-level trajectories
  geom_line(
    linewidth = 0.34,
    alpha = 0.30,
    colour = COL_LINE
  ) +
  
  # individual repeat-level estimates
  geom_point(
    size = 1.45,
    alpha = 0.75,
    colour = COL_POINT
  ) +
  
  # median diamonds
  geom_point(
    data = median_dat,
    aes(
      x = Model,
      y = Value
    ),
    inherit.aes = FALSE,
    shape = 23,
    size = 3.0,
    stroke = 0.50,
    fill = "white",
    colour = COL_TEXT
  ) +
  
  facet_wrap(
    ~ Metric,
    nrow = 1,
    scales = "free_y"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    
    title =
      paste0(
        "Repeat-level stability of nested top-25 reconstruction"
      ),
    
    subtitle =
      paste0(
        EXPECTED_REPEATS,
        " repeated donor-blocked nested-CV runs; ",
        "lines connect the same repeat"
      )
  ) +
  
  theme_classic(
    base_family = "Helvetica",
    base_size = 7
  ) +
  
  theme(
    
    plot.title.position =
      "plot",
    
    plot.title =
      element_text(
        face = "bold",
        size = 8.2,
        colour = COL_TEXT,
        hjust = 0.5,
        lineheight = 1.0,
        margin = margin(
          b = 1.2,
          unit = "pt"
        )
      ),
    
    plot.subtitle =
      element_text(
        size = 6.0,
        colour = COL_MUTED,
        hjust = 0.5,
        lineheight = 1.0,
        margin = margin(
          b = 2.0,
          unit = "pt"
        )
      ),
    
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold",
        size = 6.8,
        colour = COL_TEXT,
        margin = margin(
          b = 1.5,
          unit = "pt"
        )
      ),
    
    axis.title =
      element_blank(),
    
    axis.text.x =
      element_text(
        size = 5.4,
        colour = COL_TEXT,
        lineheight = 0.90,
        margin = margin(
          t = 2.5,
          unit = "pt"
        )
      ),
    
    axis.text.y =
      element_text(
        size = 5.7,
        colour = COL_MUTED
      ),
    
    axis.line =
      element_line(
        linewidth = 0.32,
        colour = COL_TEXT
      ),
    
    axis.ticks.x =
      element_blank(),
    
    axis.ticks.y =
      element_line(
        linewidth = 0.25,
        colour = COL_TEXT
      ),
    
    axis.ticks.length =
      grid::unit(
        1.2,
        "pt"
      ),
    
    # reduce space between Spearman and R2 panels
    panel.spacing.x =
      grid::unit(
        4,
        "pt"
      ),
    
    # tighter outer margins
    plot.margin =
      margin(
        t = 2.5,
        r = 2.5,
        b = 7.0,
        l = 2.5,
        unit = "pt"
      )
  )


# ------------------------------------------------------------
# 11. Output directory
# ------------------------------------------------------------

out_dir <- file.path(
  project_root,
  "figures",
  "SUPPLEMENTARY FIGURES",
  "S6A_descriptive_noWilcoxon_compact"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 12. Save source data
# ------------------------------------------------------------

readr::write_csv(
  summary_tbl,
  file.path(
    out_dir,
    "S6A_repeat_descriptive_summary.csv"
  )
)

readr::write_csv(
  diff_tbl,
  file.path(
    out_dir,
    "S6A_within_repeat_differences_descriptive.csv"
  )
)

readr::write_csv(
  plot_dat,
  file.path(
    out_dir,
    "S6A_plot_source_data.csv"
  )
)


# ------------------------------------------------------------
# 13. Figure dimensions
# ------------------------------------------------------------
#
# Previous version:
#   178 x 70 mm
#
# Compact version:
#   158 x 74 mm
#
# Slightly narrower + slightly taller to reduce empty appearance
# while retaining readable typography.
# ------------------------------------------------------------

FIG_W_MM <- 166
FIG_H_MM <- 74


png_file <- file.path(
  out_dir,
  "Supplementary_Fig6A_descriptive.png"
)

pdf_file <- file.path(
  out_dir,
  "Supplementary_Fig6A_descriptive.pdf"
)

tiff_file <- file.path(
  out_dir,
  "Supplementary_Fig6A_descriptive.tiff"
)


# ------------------------------------------------------------
# 14. PNG export
# ------------------------------------------------------------

ggsave(
  filename = png_file,
  plot = p,
  width = FIG_W_MM,
  height = FIG_H_MM,
  units = "mm",
  dpi = 600,
  bg = "white"
)


if (
  !file.exists(png_file) ||
  is.na(
    file.info(
      png_file
    )$size
  ) ||
  file.info(
    png_file
  )$size <= 0
) {
  
  stop(
    paste0(
      "PNG export failed: ",
      png_file
    )
  )
}


# ------------------------------------------------------------
# 15. PDF export
# ------------------------------------------------------------

grDevices::pdf(
  file = pdf_file,
  width =
    FIG_W_MM / 25.4,
  height =
    FIG_H_MM / 25.4,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE,
  paper = "special"
)

print(p)

grDevices::dev.off()


if (
  !file.exists(pdf_file) ||
  is.na(
    file.info(
      pdf_file
    )$size
  ) ||
  file.info(
    pdf_file
  )$size <= 0
) {
  
  stop(
    paste0(
      "PDF export failed: ",
      pdf_file
    )
  )
}


# ------------------------------------------------------------
# 16. TIFF export
# ------------------------------------------------------------

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


if (
  !file.exists(tiff_file) ||
  is.na(
    file.info(
      tiff_file
    )$size
  ) ||
  file.info(
    tiff_file
  )$size <= 0
) {
  
  stop(
    paste0(
      "TIFF export failed: ",
      tiff_file
    )
  )
}


# ------------------------------------------------------------
# 17. README
# ------------------------------------------------------------

readme <- c(
  
  "Supplementary Figure S6A v6",
  "",
  
  "Purpose:",
  paste0(
    "Descriptive comparison of the full 1,373-gene non-programme model ",
    "and nested top-25 reconstruction across the same ",
    "20 repeated donor-blocked nested-CV runs."
  ),
  "",
  
  "Inferential policy:",
  paste0(
    "The 20 repeated runs reuse the same donors and profiles ",
    "and are therefore not independent experimental replicates."
  ),
  
  paste0(
    "Repeat-level differences are shown only to illustrate ",
    "stability and consistency."
  ),
  
  paste0(
    "No Wilcoxon test or other repeat-level inferential ",
    "P value is calculated."
  ),
  "",
  
  "Visual revision:",
  "Figure width reduced from 178 mm to 158 mm.",
  "Figure height increased from 70 mm to 74 mm.",
  "Fonts enlarged.",
  "Inter-panel spacing reduced.",
  "Outer margins reduced.",
  "Repeat points and median diamonds slightly enlarged.",
  "",
  
  "Outputs:",
  "S6A_repeat_descriptive_summary.csv",
  "S6A_within_repeat_differences_descriptive.csv",
  "S6A_plot_source_data.csv",
  "Supplementary_Fig6A_descriptive.png",
  "Supplementary_Fig6A_descriptive.pdf",
  "Supplementary_Fig6A_descriptive.tiff"
)


writeLines(
  readme,
  file.path(
    out_dir,
    "README_S6A_v6.txt"
  )
)


# ------------------------------------------------------------
# 18. Console summary
# ------------------------------------------------------------

cat(
  "\n",
  "============================================\n",
  sep = ""
)

cat(
  "SUPPLEMENTARY FIGURE S6A v6 COMPLETE\n"
)

cat(
  "============================================\n\n"
)


cat(
  "Repeat-level descriptive summary:\n"
)

print(
  summary_tbl,
  n = Inf
)


cat(
  "\nNo Wilcoxon test was performed.\n"
)

cat(
  paste0(
    "Repeat-level comparisons are descriptive ",
    "only and are used to assess stability ",
    "and consistency.\n"
  )
)


cat(
  "\nOutput files:\n"
)

cat(
  "PNG:  ",
  png_file,
  "\n",
  sep = ""
)

cat(
  "PDF:  ",
  pdf_file,
  "\n",
  sep = ""
)

cat(
  "TIFF: ",
  tiff_file,
  "\n",
  sep = ""
)


cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)