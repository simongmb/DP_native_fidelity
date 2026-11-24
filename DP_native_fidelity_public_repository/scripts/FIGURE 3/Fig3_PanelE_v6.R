# ============================================================
# Figure 3E v6 — Representative genes illustrate two
# anatomy-referenced rescue outcomes
#
# IMPORTANT:
#   Lines indicate programme membership only.
#   They do NOT represent gene-gene regulatory interactions.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(grid)
})

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

COL_NATIVE <- "#3D8C82"
COL_ALT <- "#C97857"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"
COL_NATIVE_FILL <- "#EAF4F2"
COL_ALT_FILL <- "#F9EEEA"
COL_LIGHT <- "#D8E1E5"

native_genes <- c(
  "LRP4",
  "NDP",
  "PTGDS",
  "ETV1",
  "ENPP2"
)

alternative_genes <- c(
  "MYLK",
  "MYL9",
  "LMOD1",
  "FLNA",
  "VCL",
  "FERMT2"
)

native_nodes <- tibble(
  gene = native_genes,
  x = c(
    0.16,
    0.34,
    0.16,
    0.34,
    0.25
  ),
  y = c(
    0.625,
    0.625,
    0.445,
    0.445,
    0.265
  ),
  module = "Native identity"
)

alternative_nodes <- tibble(
  gene = alternative_genes,
  x = c(
    0.66,
    0.84,
    0.66,
    0.84,
    0.66,
    0.84
  ),
  y = c(
    0.625,
    0.625,
    0.445,
    0.445,
    0.265,
    0.265
  ),
  module = "Alternative state"
)

node_df <- dplyr::bind_rows(
  native_nodes,
  alternative_nodes
)

module_centers <- tibble(
  module = c(
    "Native identity",
    "Alternative state"
  ),
  x = c(
    0.25,
    0.75
  ),
  y = c(
    0.735,
    0.735
  )
)

edge_df <- node_df %>%
  dplyr::left_join(
    module_centers,
    by = "module",
    suffix = c(
      "_gene",
      "_center"
    )
  )

source_df <- dplyr::bind_rows(
  tibble(
    program =
      "Native-aligned rescue",
    biological_component =
      "Alignment with native DP identity",
    gene =
      native_genes
  ),
  tibble(
    program =
      "Alternative-state rescue",
    biological_component =
      "Alternative low-fidelity rescue pattern",
    gene =
      alternative_genes
  )
)

p <- ggplot() +

  annotate(
    "rect",
    xmin = 0.025, xmax = 0.475,
    ymin = 0.17, ymax = 0.875,
    fill = "#F5FAF9",
    colour = "#D7E7E3",
    linewidth = 0.24
  ) +

  annotate(
    "rect",
    xmin = 0.525, xmax = 0.975,
    ymin = 0.17, ymax = 0.875,
    fill = "#FCF7F4",
    colour = "#EBDDD6",
    linewidth = 0.24
  ) +

  annotate(
    "segment",
    x = 0.50, xend = 0.50,
    y = 0.20, yend = 0.84,
    linewidth = 0.28,
    colour = COL_LIGHT
  ) +

  geom_segment(
    data = edge_df,
    aes(
      x = x_center,
      y = y_center - 0.035,
      xend = x_gene,
      yend = y_gene + 0.035,
      colour = module
    ),
    linewidth = 0.28,
    alpha = 0.32,
    show.legend = FALSE
  ) +

  geom_point(
    data = native_nodes,
    aes(x = x, y = y),
    shape = 21,
    size = 7.0,
    stroke = 0.55,
    colour = COL_NATIVE,
    fill = "white"
  ) +

  geom_text(
    data = native_nodes,
    aes(x = x, y = y, label = gene),
    size = 1.72,
    family = "Helvetica",
    fontface = "italic",
    colour = COL_TEXT
  ) +

  geom_point(
    data = alternative_nodes,
    aes(x = x, y = y),
    shape = 21,
    size = 7.0,
    stroke = 0.55,
    colour = COL_ALT,
    fill = "white"
  ) +

  geom_text(
    data = alternative_nodes,
    aes(x = x, y = y, label = gene),
    size = 1.72,
    family = "Helvetica",
    fontface = "italic",
    colour = COL_TEXT
  ) +

  annotate(
    "text",
    x = 0.250, y = 0.835,
    label = "Native-aligned rescue",
    size = 2.35,
    family = "Helvetica",
    colour = COL_NATIVE,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = 0.250, y = 0.780,
    label = "Native-aligned rescue programme",
    size = 1.75,
    family = "Helvetica",
    colour = COL_MUTED
  ) +

  annotate(
    "text",
    x = 0.750, y = 0.835,
    label = "Alternative-state rescue",
    size = 2.35,
    family = "Helvetica",
    colour = COL_ALT,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = 0.750, y = 0.780,
    label = "Alternative-state rescue programme",
    size = 1.75,
    family = "Helvetica",
    colour = COL_MUTED
  ) +

  annotate(
    "text",
    x = 0.250, y = 0.102,
    label = "Alignment with native\nDP identity",
    size = 1.78,
    family = "Helvetica",
    colour = COL_NATIVE,
    fontface = "bold",
    lineheight = 1.08
  ) +

  annotate(
    "text",
    x = 0.750, y = 0.102,
    label = "Alternative low-fidelity\nrescue pattern",
    size = 1.78,
    family = "Helvetica",
    colour = COL_ALT,
    fontface = "bold",
    lineheight = 1.08
  ) +

  labs(
    subtitle = "Selected genes illustrate distinct anatomy-referenced rescue outcomes",
    caption = paste(
      "Selected genes are illustrative members of the frozen programmes",
      "and were not used to define either programme."
    )
  ) +

  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0.045, 0.885),
    clip = "off"
  ) +

  scale_colour_manual(
    values = c(
      "Native identity" = COL_NATIVE,
      "Alternative state" = COL_ALT
    )
  ) +

  theme_void(
    base_family = "Helvetica"
  ) +

  theme(
    plot.subtitle = element_text(
      size = 6.1,
      colour = COL_MUTED,
      hjust = 0.5,
      margin = margin(b = 3)
    ),

    plot.caption = element_text(
      size = 4.45,
      colour = "#98A7AE",
      hjust = 0.5,
      margin = margin(t = 4)
    ),

    plot.margin = margin(
      t = 5,
      r = 6,
      b = 5,
      l = 6,
      unit = "pt"
    )
  )

out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 3",
  "Panel_E_Nature_final_v6"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    out_dir,
    "Fig3_PanelE_Nature_final_v6.png"
  ),
  plot = p,
  width = 88,
  height = 47,
  units = "mm",
  dpi = 600,
  bg = "white"
)

pdf(
  file = file.path(
    out_dir,
    "Fig3_PanelE_Nature_final_v6.pdf"
  ),
  width = 88 / 25.4,
  height = 47 / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)

print(p)

dev.off()

ggsave(
  filename = file.path(
    out_dir,
    "Fig3_PanelE_Nature_final_v6.tiff"
  ),
  plot = p,
  width = 88,
  height = 47,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write_csv(
  source_df,
  file.path(
    out_dir,
    "Fig3_PanelE_Nature_final_v6_representative_genes.csv"
  )
)

readme <- c(

  "Figure 3E — Nature final | Representative genes illustrate two anatomy-referenced rescue outcomes",

  "",

  "Native-aligned rescue:",

  paste(
    native_genes,
    collapse = ", "
  ),

  "",

  "Alternative-state rescue:",

  paste(
    alternative_genes,
    collapse = ", "
  ),

  "",

  "Important:",

  paste(
    "The displayed genes are illustrative members selected from the frozen",
    "81-gene Native-aligned rescue and 99-gene Alternative-state rescue programmes."
  ),

  "They were not used to define either programme.",

  paste(
    "Lines indicate programme membership only and do not represent",
    "experimentally validated gene-gene regulatory interactions."
  ),

  paste(
    "The panel illustrates contrasting anatomy-referenced rescue outcomes",
    "and does not establish biological fate or functional regenerative potency."
  )
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig3_PanelE_Nature_final_v6.txt"
  )
)

cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
