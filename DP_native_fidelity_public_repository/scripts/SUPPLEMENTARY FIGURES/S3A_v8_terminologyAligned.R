
# ============================================================
# Supplementary Figure S3A v8 — Contractile-state overlap within the
# Alternative-state rescue programme
#
# Changes from v2:
#   1. Soften the title to avoid implying statistically
#      significant enrichment.
#   2. Restore legend order:
#      Native-aligned rescue first,
#      Alternative-state rescue second.
#   3. Slightly tighten the caption wording.
#
# All statistics, positions, representative genes,
# and interpretation are otherwise unchanged.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(grid)
  library(patchwork)
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
  if (dir.exists(file.path(x, "results"))) score <- score + 2
  if (dir.exists(file.path(x, "figures"))) score <- score + 1
  score
}

scores <- sapply(candidate_roots, score_root)

project_root <- candidate_roots[
  which.max(scores)
]

cat(
  "Project root:",
  project_root,
  "\n"
)


# ------------------------------------------------------------
# 2. Frozen overlap statistics
# ------------------------------------------------------------
overlap_df <- tibble(
  reference = c(
    "DS contractility core",
    "DS contractility core",
    "vSMC smooth-muscle\ncontraction core",
    "vSMC smooth-muscle\ncontraction core"
  ),

  program = c(
    "Native-aligned rescue",
    "Alternative-state rescue",
    "Native-aligned rescue",
    "Alternative-state rescue"
  ),

  hits = c(
    2,
    6,
    0,
    2
  ),

  program_size = c(
    81,
    99,
    81,
    99
  ),

  p_value = c(
    0.30,
    0.30,
    0.50,
    0.50
  )
) %>%
  mutate(
    overlap_fraction = hits / program_size,
    overlap_percent = 100 * overlap_fraction,
    hit_label = paste0(
      hits,
      "/",
      program_size
    )
  )


# ------------------------------------------------------------
# 3. Programme order and manual y positions
# ------------------------------------------------------------
overlap_df <- overlap_df %>%
  mutate(
    program = factor(
      program,
      levels = c(
        "Native-aligned rescue",
        "Alternative-state rescue"
      )
    ),

    ref_base = case_when(
      reference == "DS contractility core" ~ 2,
      TRUE ~ 1
    ),

    y_pos = case_when(
      program == "Native-aligned rescue" ~ ref_base + 0.10,
      TRUE ~ ref_base - 0.10
    )
  )


# ------------------------------------------------------------
# 4. Representative contractile genes
# ------------------------------------------------------------
contractile_genes <- tibble(
  gene = c(
    "MYL9",
    "LMOD1",
    "FLNA",
    "MYLK",
    "FERMT2",
    "VCL"
  ),

  in_DS_core = TRUE,

  in_vSMC_core = c(
    TRUE,
    FALSE,
    FALSE,
    TRUE,
    FALSE,
    FALSE
  )
)


# ------------------------------------------------------------
# 5. Colors
# ------------------------------------------------------------
COL_NATIVE <- "#3D8C82"
COL_ALT <- "#C97857"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"

program_cols <- c(
  "Native-aligned rescue" =
    COL_NATIVE,

  "Alternative-state rescue" =
    COL_ALT
)


# ------------------------------------------------------------
# 6. Main overlap plot
# ------------------------------------------------------------
p_main <- ggplot(
  overlap_df,
  aes(
    x = overlap_percent,
    y = y_pos,
    colour = program
  )
) +

  geom_segment(
    aes(
      x = 0,
      xend = overlap_percent,
      y = y_pos,
      yend = y_pos
    ),
    linewidth = 0.58,
    alpha = 0.34,
    show.legend = FALSE
  ) +

  geom_point(
    size = 2.35
  ) +

  geom_text(
    aes(
      label = hit_label
    ),
    hjust = -0.38,
    size = 5.0 / ggplot2::.pt,
    colour = COL_TEXT,
    show.legend = FALSE
  ) +

  annotate(
    "text",
    x = 8.15,
    y = 2,
    label = "P = 0.30",
    size = 5.0 / ggplot2::.pt,
    colour = COL_MUTED
  ) +

  annotate(
    "text",
    x = 8.15,
    y = 1,
    label = "P = 0.50",
    size = 5.0 / ggplot2::.pt,
    colour = COL_MUTED
  ) +

  scale_colour_manual(
    values = program_cols,
    breaks = c(
      "Native-aligned rescue",
      "Alternative-state rescue"
    ),
    name = NULL
  ) +

  scale_y_continuous(
    breaks = c(
      2,
      1
    ),
    labels = c(
      "DS contractility core",
      "vSMC smooth-muscle\ncontraction core"
    ),
    limits = c(
      0.68,
      2.32
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_x_continuous(
    limits = c(
      0,
      9.2
    ),
    breaks = c(
      0,
      2,
      4,
      6,
      8
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = c(
      0,
      0
    )
  ) +

  labs(
    x = "Programme overlap",
    y = NULL,
    title =
      "Contractile overlap in the Alternative-state rescue programme"
  ) +

  theme_classic(
    base_size = 6
  ) +

  theme(
    plot.title = element_text(
      size = 5.9,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      margin = margin(
        b = 4
      )
    ),
    plot.title.position = "plot",

    axis.title.x = element_text(
      size = 5.8,
      colour = COL_TEXT,
      margin = margin(
        t = 4
      )
    ),

    axis.text.x = element_text(
      size = 5.0,
      colour = COL_MUTED
    ),

    axis.text.y = element_text(
      size = 5.2,
      colour = COL_TEXT,
      margin = margin(
        r = 4
      ),
      lineheight = 0.95
    ),

    axis.line.y = element_blank(),

    axis.ticks.y = element_blank(),

    axis.line.x = element_line(
      linewidth = 0.30,
      colour = COL_TEXT
    ),

    axis.ticks.x = element_line(
      linewidth = 0.28,
      colour = COL_TEXT
    ),

    legend.position = "top",

    legend.justification = "center",

    legend.text = element_text(
      size = 5.0,
      colour = "#526A74"
    ),

    legend.key.width = unit(
      2.4,
      "mm"
    ),

    legend.key.height = unit(
      1.8,
      "mm"
    ),

    plot.margin = margin(
      t = 4,
      r = 6,
      b = 3,
      l = 5,
      unit = "pt"
    ),

    panel.grid = element_blank()
  )


# ------------------------------------------------------------
# 7. Representative-gene matrix
# ------------------------------------------------------------
gene_long <- contractile_genes %>%
  pivot_longer(
    cols = c(
      in_DS_core,
      in_vSMC_core
    ),
    names_to = "context",
    values_to = "present"
  ) %>%
  mutate(
    context = recode(
      context,
      in_DS_core =
        "DS contractility",

      in_vSMC_core =
        "vSMC contraction"
    ),

    gene = factor(
      gene,
      levels = rev(
        c(
          "MYL9",
          "LMOD1",
          "FLNA",
          "MYLK",
          "FERMT2",
          "VCL"
        )
      )
    )
  )


p_gene <- ggplot(
  gene_long,
  aes(
    x = context,
    y = gene
  )
) +

  geom_tile(
    aes(
      fill = present
    ),
    colour = "white",
    linewidth = 0.55,
    width = 0.92,
    height = 0.92
  ) +

  scale_fill_manual(
    values = c(
      "TRUE" = COL_ALT,
      "FALSE" = "#EEF2F3"
    ),
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = NULL,
    title =
      "Representative contractile genes in the 99-gene programme"
  ) +

  theme_classic(
    base_size = 6
  ) +

  theme(
    plot.title = element_text(
      size = 5.5,
      face = "bold",
      colour = COL_TEXT,
      hjust = 0.5,
      margin = margin(
        b = 3
      )
    ),
    plot.title.position = "plot",

    axis.text.x = element_text(
      size = 5.0,
      colour = COL_TEXT
    ),

    axis.text.y = element_text(
      size = 5.0,
      face = "italic",
      colour = COL_TEXT,
      margin = margin(
        r = 3
      )
    ),

    axis.line = element_blank(),

    axis.ticks = element_blank(),

    plot.margin = margin(
      t = 6,
      r = 5,
      b = 2,
      l = 5,
      unit = "pt"
    ),

    panel.grid = element_blank()
  )


# ------------------------------------------------------------
# 8. Combine
# ------------------------------------------------------------
p <- p_main /
  p_gene +
  plot_layout(
    heights = c(
      1.22,
      1
    )
  )


# ------------------------------------------------------------
# 9. Output
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "SUPPLEMENTARY FIGURES",
  "S3_Panel_A_v8"
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
    "S3_Panel_A_v8.png"
  ),
  plot = p,
  width = 80,
  height = 72,
  units = "mm",
  dpi = 600,
  bg = "white"
)


# PDF
pdf(
  file = file.path(
    out_dir,
    "S3_Panel_A_v8.pdf"
  ),
  width = 80 / 25.4,
  height = 72 / 25.4,
  useDingbats = FALSE
)

print(p)

dev.off()


# TIFF
ggsave(
  filename = file.path(
    out_dir,
    "S3_Panel_A_v8.tiff"
  ),
  plot = p,
  width = 80,
  height = 72,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 10. Export source data
# ------------------------------------------------------------
write_csv(
  overlap_df,
  file.path(
    out_dir,
    "S3_Panel_A_v8_overlap_summary.csv"
  )
)

write_csv(
  contractile_genes,
  file.path(
    out_dir,
    "S3_Panel_A_v8_contractile_genes.csv"
  )
)


# ------------------------------------------------------------
# 11. README
# ------------------------------------------------------------
readme <- c(

  "Supplementary Figure S3A v8 — Contractile-state overlap within the Alternative-state rescue programme",

  "",

  "Changes from v2:",

  paste(
    "The subtitle was softened from a comparative-frequency claim",
    "to a descriptive overlap statement because the programme-wide",
    "differences are not statistically significant."
  ),

  paste(
    "Legend order was restored to Native-aligned rescue first",
    "and Alternative-state rescue second."
  ),

  paste(
    "The caption was shortened to:",
    "'Contextual overlap; no significant programme-wide enrichment'."
  ),

  "",

  "Statistics unchanged:",

  "DS contractility core: 2/81 vs 6/99; P approximately 0.30",
  "vSMC contraction core: 0/81 vs 2/99; P approximately 0.50",

  "",

  "Interpretation:",

  paste(
    "The Alternative-state rescue programme contains a recognizable",
    "contractile component, but this is contextual mechanistic support",
    "rather than significant enrichment of the entire 99-gene programme."
  )
)


writeLines(
  readme,
  file.path(
    out_dir,
    "README_S3_Panel_A_v8.txt"
  )
)


cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
