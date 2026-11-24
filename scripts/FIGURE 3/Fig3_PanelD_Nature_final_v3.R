
# ============================================================
# Figure 3D v3 — Programme direction composition
#
# Changes from v1:
#   1. Shorten subtitle to prevent clipping.
#   2. Make stack order match legend order:
#        left  = Restored up
#        right = Restored down
#   3. Keep counts, percentages, colors, and statistics unchanged.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(grid)
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

  return(getwd())
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

candidate_roots <- candidate_roots[
  dir.exists(candidate_roots)
]

score_root <- function(x) {

  score <- 0

  if (dir.exists(file.path(x, "results"))) {
    score <- score + 2
  }

  if (dir.exists(file.path(x, "figures"))) {
    score <- score + 1
  }

  score
}

scores <- sapply(
  candidate_roots,
  score_root
)

project_root <- candidate_roots[
  which.max(scores)
]

cat(
  "Project root:",
  project_root,
  "\n"
)


# ------------------------------------------------------------
# 2. Frozen composition data
# ------------------------------------------------------------
plot_dat <- tibble(
  program = c(
    "Native-aligned rescue",
    "Native-aligned rescue",
    "Alternative-state rescue",
    "Alternative-state rescue"
  ),

  direction = c(
    "Restored up",
    "Restored down",
    "Restored up",
    "Restored down"
  ),

  n = c(
    55,
    26,
    30,
    69
  )
) %>%
  group_by(program) %>%
  mutate(
    total = sum(n),
    fraction = n / total,
    percent = 100 * fraction,
    label = paste0(
      n,
      " (",
      round(percent),
      "%)"
    )
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 3. QC
# ------------------------------------------------------------
qc <- plot_dat %>%
  dplyr::group_by(program) %>%
  dplyr::summarise(
    total = sum(n),
    up = sum(
      n[
        direction == "Restored up"
      ]
    ),
    down = sum(
      n[
        direction == "Restored down"
      ]
    ),
    .groups = "drop"
  )

print(qc)

native_total <- qc$total[
  qc$program ==
    "Native-aligned rescue"
]

alt_total <- qc$total[
  qc$program ==
    "Alternative-state rescue"
]

if (
  length(native_total) != 1 ||
    native_total != 81
) {
  stop(
    "Native-aligned rescue programme does not sum to 81."
  )
}

if (
  length(alt_total) != 1 ||
    alt_total != 99
) {
  stop(
    "Alternative-state rescue programme does not sum to 99."
  )
}


# ------------------------------------------------------------
# 4. Factor ordering
# ------------------------------------------------------------
plot_dat <- plot_dat %>%
  mutate(
    program = factor(
      program,
      levels = rev(
        c(
          "Native-aligned rescue",
          "Alternative-state rescue"
        )
      )
    ),

    direction = factor(
      direction,
      levels = c(
        "Restored up",
        "Restored down"
      )
    )
  )


# ------------------------------------------------------------
# 5. Colors
# ------------------------------------------------------------
COL_UP <- "#C97857"
COL_DOWN <- "#3D8C82"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"

direction_cols <- c(
  "Restored up" = COL_UP,
  "Restored down" = COL_DOWN
)


# ------------------------------------------------------------
# 6. Plot — Nature-style compact panel
# ------------------------------------------------------------
p <- ggplot(
  plot_dat,
  aes(
    x = fraction,
    y = program,
    fill = direction
  )
) +

  geom_col(
    width = 0.46,
    colour = "white",
    linewidth = 0.30,
    position = position_stack(
      reverse = TRUE
    )
  ) +

  geom_text(
    aes(
      label = label
    ),
    position = position_fill(
      vjust = 0.5,
      reverse = TRUE
    ),
    size = 2.05,
    family = "Helvetica",
    colour = "white",
    fontface = "bold"
  ) +

  scale_fill_manual(
    values = direction_cols,
    breaks = c(
      "Restored up",
      "Restored down"
    ),
    name = NULL
  ) +

  scale_x_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    limits = c(
      0,
      1
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    x = "Programme composition",
    y = NULL,
    subtitle = "Programme direction composition"
  ) +

  theme_classic(
    base_family = "Helvetica",
    base_size = 7
  ) +

  theme(
    plot.subtitle = element_text(
      size = 6.2,
      colour = COL_MUTED,
      hjust = 0,
      margin = margin(
        b = 4
      )
    ),

    axis.title.x = element_text(
      size = 5.6,
      colour = COL_TEXT,
      margin = margin(
        t = 5
      )
    ),

    axis.text.x = element_text(
      size = 5.1,
      colour = COL_MUTED
    ),

    axis.text.y = element_text(
      size = 5.5,
      colour = COL_TEXT,
      margin = margin(
        r = 5
      )
    ),

    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),

    axis.line.x = element_line(
      linewidth = 0.32,
      colour = COL_TEXT
    ),

    axis.ticks.x = element_line(
      linewidth = 0.28,
      colour = COL_TEXT
    ),

    axis.ticks.length = unit(
      1.3,
      "mm"
    ),

    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",

    legend.text = element_text(
      size = 5.2,
      colour = "#526A74"
    ),

    legend.key.width = unit(
      3.5,
      "mm"
    ),

    legend.key.height = unit(
      2.5,
      "mm"
    ),

    legend.spacing.x = unit(
      2.0,
      "mm"
    ),

    legend.margin = margin(
      t = 0,
      r = 0,
      b = 2,
      l = 0
    ),

    plot.margin = margin(
      t = 4,
      r = 5,
      b = 5,
      l = 5,
      unit = "pt"
    ),

    panel.grid = element_blank()
  )


# ------------------------------------------------------------
# 7. Output
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 3",
  "Panel_D_Nature_final_v3"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# PNG
ggsave(
  filename = file.path(
    out_dir,
    "Fig3_PanelD_Nature_final_v3.png"
  ),
  plot = p,
  width = 88,
  height = 43,
  units = "mm",
  dpi = 600,
  bg = "white"
)


# PDF
pdf(
  file = file.path(
    out_dir,
    "Fig3_PanelD_Nature_final_v3.pdf"
  ),
  width = 88 / 25.4,
  height = 43 / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)

print(p)

dev.off()


# TIFF
ggsave(
  filename = file.path(
    out_dir,
    "Fig3_PanelD_Nature_final_v3.tiff"
  ),
  plot = p,
  width = 88,
  height = 43,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 8. Export source data
# ------------------------------------------------------------
write_csv(
  plot_dat %>%
    mutate(
      program = as.character(program),
      direction = as.character(direction)
    ),
  file.path(
    out_dir,
    "Fig3_PanelD_Nature_final_v3_composition.csv"
  )
)


# ------------------------------------------------------------
# 9. README
# ------------------------------------------------------------
readme <- c(

  "Figure 3D v3 — Nature final | Programme direction composition",

  "",

  "Frozen counts:",

  "Native-aligned rescue programme:",
  "  Restored up = 55",
  "  Restored down = 26",
  "  Total = 81",

  "",

  "Alternative-state rescue programme:",
  "  Restored up = 30",
  "  Restored down = 69",
  "  Total = 99",

  "",

  "Nature-style layout changes:",

  "- Compact 88 x 43 mm final panel.",
  "- Helvetica typography and reduced final-size text.",
  "- Compact left-aligned legend.",
  "- Thinner bars, axes and dividers for Nature-style balance.",
  "- No counts, percentages, colors, or scientific interpretation changed."
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig3_PanelD_Nature_final_v3.txt"
  )
)


cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
