# ============================================================
# Figure 3B | Nature-style landscape layout v3
# Anatomical consistency of the 81-gene Native-aligned rescue programme
#
# Layout change only:
# - Genes on x-axis (81 columns)
# - Seven anatomical contrasts on y-axis
# - Wider / shorter panel for Figure 3 assembly
# - Per-compartment positive fractions shown at row ends
# - Scientific data, ordering logic, representative labels unchanged
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(grid)
})

# ------------------------------------------------------------
# 1. Locate project root
# ------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
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

candidate_roots <- candidate_roots[dir.exists(candidate_roots)]

score_root <- function(x) {
  score <- 0
  if (dir.exists(file.path(x, "results"))) score <- score + 2
  if (dir.exists(file.path(x, "figures"))) score <- score + 1
  score
}

scores <- sapply(candidate_roots, score_root)
project_root <- candidate_roots[which.max(scores)]

cat("Project root:", project_root, "\n")

# ------------------------------------------------------------
# 2. Locate audit CSV
# ------------------------------------------------------------
target_name <- "02_gene_level_native_fidelity_audit.csv"

search_dirs <- unique(c(
  file.path(project_root, "results"),
  project_root
))
search_dirs <- search_dirs[dir.exists(search_dirs)]

matches <- character(0)

for (d in search_dirs) {
  m <- list.files(
    d,
    pattern = paste0("^", target_name, "$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )
  matches <- c(matches, m)
}

matches <- unique(matches)

if (length(matches) == 0) {
  stop(
    paste0(
      "\nCannot find ", target_name,
      "\nSearched under:\n",
      paste(search_dirs, collapse = "\n")
    )
  )
}

if (length(matches) > 1) {
  in_results <- grepl("results", matches, ignore.case = TRUE)
  if (any(in_results)) matches <- matches[in_results]
}

audit_file <- matches[1]
cat("Using audit file:", audit_file, "\n")

# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------
dat <- read_csv(
  audit_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 4. Required columns
# ------------------------------------------------------------
required_cols <- c(
  "GeneSymbol",
  "anatomy_consistent_candidate",
  "aligned_DP_vs_DC",
  "aligned_DP_vs_CTS",
  "aligned_DP_vs_HM1",
  "aligned_DP_vs_HM2",
  "aligned_DP_vs_IRS1",
  "aligned_DP_vs_ORS1",
  "aligned_DP_vs_ORS2"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "\nMissing required columns:\n",
      paste(missing_cols, collapse = ", "),
      "\n\nDetected columns:\n",
      paste(names(dat), collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 5. Logical helper
# ------------------------------------------------------------
as_flag <- function(x) {
  if (is.logical(x)) return(x)

  x_chr <- toupper(trimws(as.character(x)))

  x_chr %in% c(
    "TRUE", "T", "1", "YES", "Y"
  )
}

dat <- dat %>%
  mutate(
    anatomy_flag = as_flag(anatomy_consistent_candidate)
  )

# ------------------------------------------------------------
# 6. Frozen 81-gene Native-aligned rescue programme
# ------------------------------------------------------------
native81 <- dat %>%
  filter(anatomy_flag)

if (nrow(native81) != 81) {
  stop(
    paste0(
      "Expected 81 anatomy-consistent genes, but found ",
      nrow(native81),
      "."
    )
  )
}

# ------------------------------------------------------------
# 7. Anatomical contrasts
# ------------------------------------------------------------
contrast_map <- c(
  "aligned_DP_vs_DC"   = "DP vs DC",
  "aligned_DP_vs_CTS"  = "DP vs CTS",
  "aligned_DP_vs_HM1"  = "DP vs HM1",
  "aligned_DP_vs_HM2"  = "DP vs HM2",
  "aligned_DP_vs_IRS1" = "DP vs IRS1",
  "aligned_DP_vs_ORS1" = "DP vs ORS1",
  "aligned_DP_vs_ORS2" = "DP vs ORS2"
)

# ------------------------------------------------------------
# 8. Long-format heatmap data
# ------------------------------------------------------------
heat_dat <- native81 %>%
  dplyr::select(
    GeneSymbol,
    all_of(names(contrast_map))
  ) %>%
  pivot_longer(
    cols = all_of(names(contrast_map)),
    names_to = "contrast_col",
    values_to = "aligned_value"
  ) %>%
  mutate(
    contrast = recode(
      contrast_col,
      !!!contrast_map
    )
  )

# ------------------------------------------------------------
# 9. Gene ordering
# ------------------------------------------------------------
row_summary <- heat_dat %>%
  group_by(GeneSymbol) %>%
  summarise(
    mean_aligned = mean(aligned_value, na.rm = TRUE),
    positive_fraction = mean(aligned_value > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(positive_fraction),
    desc(mean_aligned)
  )

gene_order <- row_summary$GeneSymbol

# ------------------------------------------------------------
# 10. Representative labels
# ------------------------------------------------------------
preferred_label_genes <- c(
  "LRP4",
  "ENPP2",
  "PTGDS",
  "ETV1",
  "NDP",
  "TBX3",
  "LGR5",
  "CLEC11A",
  "GATA6",
  "CYGB",
  "PRDM1",
  "BMPR2"
)

label_genes <- preferred_label_genes[
  preferred_label_genes %in% native81$GeneSymbol
]

gene_label_function <- function(x) {
  ifelse(x %in% label_genes, x, "")
}

# ------------------------------------------------------------
# 11. Per-compartment summaries
# ------------------------------------------------------------
column_summary <- heat_dat %>%
  group_by(contrast) %>%
  summarise(
    positive_fraction = mean(aligned_value > 0, na.rm = TRUE),
    median_aligned = median(aligned_value, na.rm = TRUE),
    n_nonmissing = sum(is.finite(aligned_value)),
    .groups = "drop"
  ) %>%
  mutate(
    pct_label = paste0(round(positive_fraction * 100), "%")
  )

overall_positive_fraction <- mean(
  heat_dat$aligned_value > 0,
  na.rm = TRUE
)

median_gene_consistency <- median(
  row_summary$positive_fraction,
  na.rm = TRUE
)

cat("\n81-gene Native-aligned rescue programme summary:\n")
cat(
  "Cell-level positive fraction =",
  round(overall_positive_fraction, 4),
  "\n"
)
cat(
  "Median gene-level positive fraction =",
  round(median_gene_consistency, 4),
  "\n"
)
cat("\nPer-compartment positive fractions:\n")
print(column_summary)

# ------------------------------------------------------------
# 12. Robust symmetric color scale
# ------------------------------------------------------------
abs_limit <- unname(
  quantile(
    abs(heat_dat$aligned_value),
    probs = 0.95,
    na.rm = TRUE
  )
)

if (!is.finite(abs_limit) || abs_limit <= 0) {
  abs_limit <- max(
    abs(heat_dat$aligned_value),
    na.rm = TRUE
  )
}

if (!is.finite(abs_limit) || abs_limit <= 0) {
  abs_limit <- 1
}

# ------------------------------------------------------------
# 13. Landscape orientation
# ------------------------------------------------------------
# Keep gene ordering identical to v3, but place genes along x.
heat_dat <- heat_dat %>%
  mutate(
    GeneSymbol = factor(
      GeneSymbol,
      levels = gene_order
    ),
    contrast = factor(
      contrast,
      levels = rev(unname(contrast_map))
    )
  )

column_summary_plot <- column_summary %>%
  mutate(
    contrast = factor(
      as.character(contrast),
      levels = rev(unname(contrast_map))
    )
  )


# Compact row labels: compartment + positive fraction.
row_label_map <- setNames(
  paste0(
    as.character(column_summary_plot$contrast),
    "  (",
    column_summary_plot$pct_label,
    ")"
  ),
  as.character(column_summary_plot$contrast)
)

# ------------------------------------------------------------
# 14. Plot
# ------------------------------------------------------------
p <- ggplot(
  heat_dat,
  aes(
    x = GeneSymbol,
    y = contrast,
    fill = aligned_value
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.07
  ) +

  scale_fill_gradient2(
    low = "#C97857",
    mid = "#F4F7F8",
    high = "#3D8C82",
    midpoint = 0,
    limits = c(-abs_limit, abs_limit),
    oob = scales::squish,
    name = "DP-aligned\ncontrast"
  ) +

  # No isolated gene labels: all 81 genes are shown and ordered,
  # but the axis is described collectively to avoid implying a second gene set.
  scale_x_discrete(
    labels = function(x) rep("", length(x)),
    expand = expansion(mult = c(0.002, 0.002))
  ) +

  # The positive fraction is incorporated directly into each row label,
  # eliminating the separate right-side percentage column.
  scale_y_discrete(
    labels = row_label_map,
    expand = expansion(mult = c(0.05, 0.05))
  ) +

  labs(
    x = "81 Native-aligned rescue genes (ordered by anatomical consistency)",
    y = NULL,
    subtitle = "Consistent DP-aligned anatomy across seven follicular compartments"
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(17, "mm"),
      barwidth = unit(2.3, "mm"),
      ticks = TRUE,
      frame.colour = NA
    )
  ) +

  theme_classic(
    base_family = "Helvetica",
    base_size = 7
  ) +

  theme(
    plot.subtitle = element_text(
      size = 6.2,
      colour = "#718993",
      hjust = 0,
      margin = margin(b = 4)
    ),

    axis.title.x = element_text(
      size = 5.4,
      colour = "#526A74",
      margin = margin(t = 5)
    ),

    axis.text.x = element_blank(),

    axis.text.y = element_text(
      size = 5.4,
      colour = "#20343D",
      margin = margin(r = 4)
    ),

    axis.line = element_blank(),
    axis.ticks = element_blank(),

    legend.position = "right",
    legend.justification = "center",

    legend.title = element_text(
      size = 5.1,
      colour = "#20343D",
      lineheight = 0.92,
      hjust = 0.5
    ),

    legend.text = element_text(
      size = 4.7,
      colour = "#718993"
    ),

    legend.margin = margin(l = 4),

    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5,
      unit = "pt"
    ),

    panel.grid = element_blank()
  )


# ------------------------------------------------------------
# 15. Output
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 3",
  "Panel_B_Nature_landscape_v3"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

WIDTH_MM  <- 88
HEIGHT_MM <- 46
DPI       <- 600

ggsave(
  filename = file.path(
    out_dir,
    "Fig3_PanelB_Nature_landscape_v3.png"
  ),
  plot = p,
  width = WIDTH_MM,
  height = HEIGHT_MM,
  units = "mm",
  dpi = DPI,
  bg = "white"
)

pdf(
  file = file.path(
    out_dir,
    "Fig3_PanelB_Nature_landscape_v3.pdf"
  ),
  width = WIDTH_MM / 25.4,
  height = HEIGHT_MM / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)

print(p)
dev.off()

ggsave(
  filename = file.path(
    out_dir,
    "Fig3_PanelB_Nature_landscape_v3.tiff"
  ),
  plot = p,
  width = WIDTH_MM,
  height = HEIGHT_MM,
  units = "mm",
  dpi = DPI,
  compression = "lzw",
  bg = "white"
)

# ------------------------------------------------------------
# 16. Export data / provenance
# ------------------------------------------------------------
write_csv(
  heat_dat %>%
    mutate(
      GeneSymbol = as.character(GeneSymbol),
      contrast = as.character(contrast)
    ),
  file.path(
    out_dir,
    "Fig3_PanelB_Nature_landscape_v3_heatmap_data.csv"
  )
)

write_csv(
  row_summary,
  file.path(
    out_dir,
    "Fig3_PanelB_Nature_landscape_v3_gene_consistency_summary.csv"
  )
)

write_csv(
  column_summary %>%
    mutate(
      contrast = as.character(contrast)
    ),
  file.path(
    out_dir,
    "Fig3_PanelB_Nature_landscape_v3_compartment_summary.csv"
  )
)

summary_df <- tibble(
  metric = c(
    "Native-aligned rescue programme size",
    "Cell-level positive fraction",
    "Median gene-level positive fraction",
    "Heatmap symmetric color limit",
    "Displayed gene labels"
  ),
  value = c(
    "81",
    sprintf("%.6f", overall_positive_fraction),
    sprintf("%.6f", median_gene_consistency),
    sprintf("%.6f", abs_limit),
    paste(label_genes, collapse = "; ")
  )
)

write_csv(
  summary_df,
  file.path(
    out_dir,
    "Fig3_PanelB_Nature_landscape_v3_summary.csv"
  )
)

readme <- c(
  "Figure 3B — Nature landscape layout",
  "",
  "Purpose:",
  "Show anatomical consistency of the frozen 81-gene Native-aligned rescue programme.",
  "",
  paste("Input:", normalizePath(audit_file)),
  "",
  "Layout change:",
  "- Genes are on the x-axis.",
  "- Seven anatomical contrasts are on the y-axis.",
  "- The panel is wider and shorter for the Figure 3 assembly.",
  "- Per-compartment positive fractions are integrated into the row labels.",
  "- Individual gene labels are removed; the x-axis describes all 81 ordered genes collectively.",
  "- No data, gene ordering logic, or scientific definition was changed.",
  "",
  paste0("Canvas: ", WIDTH_MM, " x ", HEIGHT_MM, " mm")
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig3_PanelB_Nature_landscape_v3.txt"
  )
)

cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
