
# ============================================================
# Supplementary Figure S2A v8 — Anatomical-reference robustness
# Rescue-fidelity decoupling across follicular compartments
#
# Purpose:
#   Demonstrate that the weak association between Rescue fraction
#   and native DP fidelity is not driven by one particular
#   anatomical reference comparison.
#
# Changes from v1:
#   1. Shorter subtitle to prevent clipping.
#   2. X-axis extended to +0.03 so rho = 0 is clearly inside
#      the plotting region instead of near the right boundary.
#   3. No changes to data, statistics, ordering, or interpretation.
#
# Data source:
#   02_gene_level_native_fidelity_audit.csv
#
# Included contrasts:
#   DP vs DC
#   DP vs CTS
#   DP vs HM1
#   DP vs HM2
#   DP vs IRS1
#   DP vs ORS1
#   DP vs ORS2
#   Composite native fidelity
#
# NOTE:
#   native_fidelity_raw_DP_minus_nonDP is intentionally excluded
#   because it is not direction-aligned and therefore does not
#   represent the same fidelity definition used in the main analysis.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
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

  return(
    getwd()
  )
}


start_dir <- get_script_dir()


candidate_roots <- unique(
  c(
    start_dir,
    getwd(),
    dirname(start_dir),
    dirname(
      dirname(start_dir)
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


candidate_roots <- candidate_roots[
  dir.exists(
    candidate_roots
  )
]


score_root <- function(x) {

  score <- 0

  if (
    dir.exists(
      file.path(
        x,
        "results"
      )
    )
  ) {
    score <- score + 2
  }

  if (
    dir.exists(
      file.path(
        x,
        "figures"
      )
    )
  ) {
    score <- score + 1
  }

  score
}


scores <- sapply(
  candidate_roots,
  score_root
)


project_root <- candidate_roots[
  which.max(
    scores
  )
]


cat(
  "Project root:",
  project_root,
  "\n"
)


# ------------------------------------------------------------
# 2. Locate native-fidelity audit table
# ------------------------------------------------------------
target_name <- "02_gene_level_native_fidelity_audit.csv"


search_dirs <- unique(
  c(
    file.path(
      project_root,
      "results"
    ),
    project_root
  )
)


search_dirs <- search_dirs[
  dir.exists(
    search_dirs
  )
]


matches <- character(0)


for (d in search_dirs) {

  m <- list.files(
    d,
    pattern = paste0(
      "^",
      target_name,
      "$"
    ),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )

  matches <- c(
    matches,
    m
  )
}


matches <- unique(
  matches
)


if (length(matches) == 0) {

  stop(
    paste0(
      "\nCannot find ",
      target_name,
      "\nSearched under:\n",
      paste(
        search_dirs,
        collapse = "\n"
      )
    )
  )
}


if (length(matches) > 1) {

  in_results <- grepl(
    "results",
    matches,
    ignore.case = TRUE
  )

  if (any(in_results)) {

    matches <- matches[
      in_results
    ]
  }
}


audit_file <- matches[1]


cat(
  "Using audit file:",
  audit_file,
  "\n"
)


# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------
dat <- read_csv(
  audit_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 4. Detect Rescue fraction column
# ------------------------------------------------------------
rescue_candidates <- c(
  "Rescue_fraction",
  "rescue_fraction",
  "RescueFraction",
  "rescue.frac",
  "rescue"
)


rescue_col <- rescue_candidates[
  rescue_candidates %in% names(dat)
][1]


if (is.na(rescue_col)) {

  stop(
    "Could not detect Rescue_fraction column."
  )
}


cat(
  "Rescue column:",
  rescue_col,
  "\n"
)


# ------------------------------------------------------------
# 5. Define anatomical-reference columns
# ------------------------------------------------------------
contrast_map <- c(

  "aligned_DP_vs_DC" =
    "DP vs DC",

  "aligned_DP_vs_CTS" =
    "DP vs CTS",

  "aligned_DP_vs_HM1" =
    "DP vs HM1",

  "aligned_DP_vs_HM2" =
    "DP vs HM2",

  "aligned_DP_vs_IRS1" =
    "DP vs IRS1",

  "aligned_DP_vs_ORS1" =
    "DP vs ORS1",

  "aligned_DP_vs_ORS2" =
    "DP vs ORS2",

  "native_fidelity_direction_aligned" =
    "Composite native fidelity"
)


missing_cols <- setdiff(
  names(
    contrast_map
  ),
  names(
    dat
  )
)


if (length(missing_cols) > 0) {

  stop(
    paste0(
      "\nMissing required anatomical-reference columns:\n",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 6. Calculate Spearman correlation for each reference
# ------------------------------------------------------------
calc_corr <- function(
  colname,
  label
) {

  x <- suppressWarnings(
    as.numeric(
      dat[[rescue_col]]
    )
  )

  y <- suppressWarnings(
    as.numeric(
      dat[[colname]]
    )
  )


  ok <- is.finite(x) &
    is.finite(y)


  ct <- suppressWarnings(
    cor.test(
      x[ok],
      y[ok],
      method = "spearman",
      exact = FALSE
    )
  )


  tibble(
    fidelity_column = colname,
    reference = label,
    n = sum(ok),
    rho = unname(
      ct$estimate
    ),
    p_value = ct$p.value
  )
}


plot_dat <- bind_rows(
  lapply(
    names(
      contrast_map
    ),
    function(colname) {

      calc_corr(
        colname,
        contrast_map[[colname]]
      )
    }
  )
)


# ------------------------------------------------------------
# 7. Ordering and labels
# ------------------------------------------------------------
display_order <- c(
  "Composite native fidelity",
  "DP vs DC",
  "DP vs CTS",
  "DP vs HM1",
  "DP vs HM2",
  "DP vs IRS1",
  "DP vs ORS1",
  "DP vs ORS2"
)


plot_dat <- plot_dat %>%
  dplyr::mutate(

    reference = factor(
      reference,
      levels = rev(
        display_order
      )
    ),

    Type = ifelse(
      as.character(reference) ==
        "Composite native fidelity",
      "Composite",
      "Individual anatomy"
    ),

    rho_label = sprintf(
      "%.3f",
      rho
    )
  )


cat(
  "\nFigure 2E statistics:\n"
)


print(
  plot_dat %>%
    dplyr::mutate(
      reference =
        as.character(reference)
    ) %>%
    dplyr::arrange(
      match(
        reference,
        display_order
      )
    )
)


# ------------------------------------------------------------
# 8. Colors
# ------------------------------------------------------------
type_cols <- c(
  "Composite" =
    "#3D8C82",

  "Individual anatomy" =
    "#8EA5AF"
)


# ------------------------------------------------------------
# 9. Plot
# ------------------------------------------------------------
# Fixed x limits:
# Left side keeps enough room for rho values around -0.12.
# Right side extends beyond zero to +0.03 so the zero reference
# line appears inside the panel rather than against the border.

plot_xlim <- c(
  -0.16,
  0.03
)


p <- ggplot(
  plot_dat,
  aes(
    x = rho,
    y = reference,
    colour = Type
  )
) +

  # No-association reference line
  geom_vline(
    xintercept = 0,
    linewidth = 0.30,
    colour = "#C9D5DA"
  ) +

  # Stem from zero to each rho
  geom_segment(
    aes(
      x = 0,
      xend = rho,
      yend = reference
    ),
    linewidth = 0.62,
    alpha = 0.42,
    show.legend = FALSE
  ) +

  # Correlation point
  geom_point(
    size = 2.05,
    alpha = 0.92
  ) +

  # Numeric rho label
  geom_text(
    aes(
      label = rho_label
    ),
    hjust = 1.18,
    colour = "#657983",
    size = 5.0 / ggplot2::.pt,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = type_cols,
    name = NULL
  ) +

  scale_x_continuous(
    limits = plot_xlim,
    breaks = c(
      -0.15,
      -0.10,
      -0.05,
      0.00
    ),
    labels = c(
      "-0.15",
      "-0.10",
      "-0.05",
      "0.00"
    ),
    name =
      "Spearman rho with Rescue fraction",
    expand = c(
      0,
      0
    )
  ) +

  labs(
    y = NULL,
    title =
      "Weak rescue-fidelity association across references"
  ) +

  theme_classic(
    base_size = 6
  ) +

  theme(

    # --------------------------------------------------------
    # Title
    # --------------------------------------------------------
    plot.title = element_text(
      size = 6.2,
      face = "bold",
      colour = "#273C4B",
      hjust = 0.5,
      margin = margin(
        b = 2.8
      )
    ),
    plot.title.position = "plot",

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------
    legend.position = "top",

    legend.justification = "center",

    legend.direction = "horizontal",

    legend.text = element_text(
      size = 5.4,
      colour = "#526A74"
    ),

    legend.key.width = unit(
      2.2,
      "mm"
    ),

    legend.key.height = unit(
      2.0,
      "mm"
    ),

    legend.spacing.x = unit(
      1.4,
      "mm"
    ),

    legend.margin = margin(
      t = 0,
      r = 0,
      b = 1,
      l = 0
    ),

    # --------------------------------------------------------
    # X axis
    # --------------------------------------------------------
    axis.title.x = element_text(
      size = 6.1,
      colour = "#20343D",
      margin = margin(
        t = 4
      )
    ),

    axis.text.x = element_text(
      size = 5.3,
      colour = "#718993",
      margin = margin(
        t = 1.5
      )
    ),

    # --------------------------------------------------------
    # Y labels
    # --------------------------------------------------------
    axis.text.y = element_text(
      size = 5.4,
      colour = "#20343D",
      margin = margin(
        r = 3
      )
    ),

    axis.title.y =
      element_blank(),

    axis.line.y =
      element_blank(),

    axis.ticks.y =
      element_blank(),

    # --------------------------------------------------------
    # X-axis line / ticks
    # --------------------------------------------------------
    axis.line.x = element_line(
      linewidth = 0.30,
      colour = "#20343D"
    ),

    axis.ticks.x = element_line(
      linewidth = 0.28,
      colour = "#20343D"
    ),

    axis.ticks.length = unit(
      2,
      "mm"
    ),

    # --------------------------------------------------------
    # Margins
    # --------------------------------------------------------
    plot.margin = margin(
      t = 4,
      r = 6,
      b = 5,
      l = 5,
      unit = "pt"
    ),

    panel.grid =
      element_blank()
  )


# ------------------------------------------------------------
# 10. Output
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "SUPPLEMENTARY FIGURES",
  "S2_Panel_A_v8"
)


dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# PNG
# ------------------------------------------------------------
ggsave(
  filename = file.path(
    out_dir,
    "S2_Panel_A_v8.png"
  ),
  plot = p,
  width = 74,
  height = 58,
  units = "mm",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# PDF
# ------------------------------------------------------------
pdf(
  file = file.path(
    out_dir,
    "S2_Panel_A_v8.pdf"
  ),
  width = 74 / 25.4,
  height = 58 / 25.4,
  useDingbats = FALSE
)

print(p)

dev.off()


# ------------------------------------------------------------
# TIFF
# ------------------------------------------------------------
ggsave(
  filename = file.path(
    out_dir,
    "S2_Panel_A_v8.tiff"
  ),
  plot = p,
  width = 74,
  height = 58,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 11. Export statistics
# ------------------------------------------------------------
write_csv(
  plot_dat %>%
    dplyr::mutate(
      reference =
        as.character(reference)
    ) %>%
    dplyr::select(
      reference,
      fidelity_column,
      n,
      rho,
      p_value,
      Type
    ),
  file.path(
    out_dir,
    "S2_Panel_A_v8_statistics.csv"
  )
)


# ------------------------------------------------------------
# 12. README / provenance
# ------------------------------------------------------------
readme <- c(

  "Supplementary Figure S2A v8 — Anatomical-reference robustness",

  "",

  "Purpose:",

  paste(
    "Test whether the weak association between Rescue fraction and",
    "native DP fidelity is robust to the choice of anatomical reference."
  ),

  "",

  paste(
    "Input:",
    normalizePath(audit_file)
  ),

  "",

  "Included direction-aligned references:",

  paste(
    display_order,
    collapse = "; "
  ),

  "",

  "Excluded:",

  paste(
    "native_fidelity_raw_DP_minus_nonDP because it is not",
    "direction-aligned and therefore does not represent the",
    "same fidelity definition used in the main analysis."
  ),

  "",

  "Changes from v1:",

  "- Reformatted as Supplementary Figure S2A at 74 x 58 mm final size.",

  "- Scientific values, ordering, x-axis limits and correlations unchanged.",

  "- rho = 0 now appears clearly inside the plotting region.",

  "- Data, correlations, ordering, colors, and interpretation unchanged.",

  "",

  "Interpretation:",

  paste(
    "Direction-aligned anatomical comparisons consistently show",
    "only weak associations with Rescue fraction, supporting the",
    "separability of rescue magnitude and native DP fidelity."
  )
)


writeLines(
  readme,
  file.path(
    out_dir,
    "README_S2_Panel_A_v8.txt"
  )
)


cat(
  "\nDone.\n"
)

cat(
  "Output directory:\n"
)

cat(
  out_dir,
  "\n"
)
