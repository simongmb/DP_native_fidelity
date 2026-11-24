
# ============================================================
# Figure 2D | Nature-style final panel
# High rescue can accompany divergent native fidelity
#
# Purpose:
#   Show concrete examples of high-rescue genes with strongly
#   positive versus strongly negative native DP fidelity.
#
# Nature-style final layout:
# - Final panel size: 88 x 56 mm, aligned with Figure 2B/2C
# - Typography and strokes sized for final print dimensions
# - Compact legend and margins; gene labels and rescue values retained
# - Same genes, values, colors, and scientific logic
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(forcats)
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

  if (any(in_results)) {
    matches <- matches[in_results]
  }
}

audit_file <- matches[1]

cat("Using audit file:", audit_file, "\n")


# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------
dat0 <- read_csv(
  audit_file,
  show_col_types = FALSE
)

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


gene_col <- pick_col(
  dat0,
  c(
    "GeneSymbol",
    "gene_symbol",
    "Gene",
    "symbol"
  )
)

rescue_col <- pick_col(
  dat0,
  c(
    "Rescue_fraction",
    "rescue_fraction",
    "RescueFraction",
    "rescue.frac",
    "rescue"
  )
)

fidelity_col <- pick_col(
  dat0,
  c(
    "native_fidelity_direction_aligned",
    "native_fidelity",
    "Native_fidelity",
    "Native_Fidelity",
    "direction_aligned_native_fidelity",
    "donor_aware_native_fidelity",
    "native_DP_fidelity"
  )
)


if (is.na(gene_col) || is.na(rescue_col) || is.na(fidelity_col)) {
  stop(
    paste0(
      "\nRequired columns were not found.\n",
      "Detected columns:\n",
      paste(names(dat0), collapse = ", "),
      "\n\nExpected gene symbol, rescue fraction, and native-fidelity columns."
    )
  )
}


dat <- dat0 %>%
  transmute(
    GeneSymbol = as.character(.data[[gene_col]]),
    Rescue_fraction = suppressWarnings(
      as.numeric(.data[[rescue_col]])
    ),
    Native_DP_fidelity = suppressWarnings(
      as.numeric(.data[[fidelity_col]])
    )
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    is.finite(Rescue_fraction),
    is.finite(Native_DP_fidelity)
  )


# ------------------------------------------------------------
# 4. High-rescue threshold
# ------------------------------------------------------------
rescue_q3 <- unname(
  quantile(
    dat$Rescue_fraction,
    probs = 0.75,
    na.rm = TRUE,
    type = 7
  )
)

cat("Rescue Q3 =", sprintf("%.9f", rescue_q3), "\n")


# ------------------------------------------------------------
# 5. Fixed representative candidates
# ------------------------------------------------------------
upper_candidates <- c(
  "LRP4",
  "ENPP2",
  "PTGDS",
  "ETV1",
  "NDP"
)

lower_candidates <- c(
  "MYL9",
  "AXL",
  "MYLK",
  "APOE",
  "NOTCH3"
)


# ------------------------------------------------------------
# 6. Select exactly these genes when available
# ------------------------------------------------------------
upper_sel <- dat %>%
  filter(
    GeneSymbol %in% upper_candidates
  ) %>%
  mutate(
    Group = "High fidelity",
    order_within = match(
      GeneSymbol,
      upper_candidates
    )
  ) %>%
  arrange(order_within)


lower_sel <- dat %>%
  filter(
    GeneSymbol %in% lower_candidates
  ) %>%
  mutate(
    Group = "Low fidelity",
    order_within = match(
      GeneSymbol,
      lower_candidates
    )
  ) %>%
  arrange(order_within)


# If any predefined genes are missing, stop rather than silently change the set.
missing_upper <- setdiff(
  upper_candidates,
  upper_sel$GeneSymbol
)

missing_lower <- setdiff(
  lower_candidates,
  lower_sel$GeneSymbol
)

if (length(missing_upper) > 0 || length(missing_lower) > 0) {

  stop(
    paste0(
      "\nSome predefined representative genes were not found.\n",
      "Missing upper-fidelity genes: ",
      paste(missing_upper, collapse = ", "),
      "\nMissing lower-fidelity genes: ",
      paste(missing_lower, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 7. Build plotting table with group gap
# ------------------------------------------------------------
rep_dat <- bind_rows(
  upper_sel,
  lower_sel
) %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "High fidelity",
        "Low fidelity"
      )
    )
  )


# Manual y positions with a visible gap between groups
upper_positions <- c(10, 9, 8, 7, 6)
lower_positions <- c(4.7, 3.7, 2.7, 1.7, 0.7)

rep_dat <- rep_dat %>%
  mutate(
    y_pos = case_when(
      Group == "High fidelity" ~
        upper_positions[match(GeneSymbol, upper_candidates)],

      Group == "Low fidelity" ~
        lower_positions[match(GeneSymbol, lower_candidates)]
    ),
    rescue_label = paste0(
      "R = ",
      sprintf("%.2f", Rescue_fraction)
    )
  )


cat("\nRepresentative genes used:\n")
print(
  rep_dat %>%
    dplyr::select(
      GeneSymbol,
      Rescue_fraction,
      Native_DP_fidelity,
      Group
    )
)


# ------------------------------------------------------------
# 8. Colors
# ------------------------------------------------------------
group_cols <- c(
  "High fidelity" = "#3D8C82",
  "Low fidelity"  = "#C97857"
)


# ------------------------------------------------------------
# 9. X range and label offsets
# ------------------------------------------------------------
x_min <- min(
  c(
    rep_dat$Native_DP_fidelity,
    0
  )
)

x_max <- max(
  c(
    rep_dat$Native_DP_fidelity,
    0
  )
)

x_span <- x_max - x_min

left_label_offset <- x_span * 0.075
right_label_offset <- x_span * 0.075


# ------------------------------------------------------------
# 10. Plot — sized for the final 88 x 56 mm Nature panel
# ------------------------------------------------------------
p <- ggplot(
  rep_dat,
  aes(
    x = Native_DP_fidelity,
    y = y_pos,
    colour = Group
  )
) +
  # Central zero line
  geom_vline(
    xintercept = 0,
    linewidth = 0.28,
    colour = "#C8D4D9"
  ) +
  # Subtle divider between the two gene groups
  geom_hline(
    yintercept = 5.35,
    linewidth = 0.22,
    colour = "#E2E8EB"
  ) +
  # Lollipop stems
  geom_segment(
    aes(
      x = 0,
      xend = Native_DP_fidelity,
      y = y_pos,
      yend = y_pos
    ),
    linewidth = 0.48,
    alpha = 0.48,
    show.legend = FALSE
  ) +
  # Fidelity points
  geom_point(
    size = 1.65,
    alpha = 0.92
  ) +
  # Rescue-fraction labels on the right of positive-fidelity genes
  geom_text(
    data = rep_dat %>% filter(Native_DP_fidelity >= 0),
    aes(
      x = Native_DP_fidelity + right_label_offset,
      y = y_pos,
      label = rescue_label
    ),
    inherit.aes = FALSE,
    colour = "#71858D",
    size = 1.78,
    hjust = 0,
    vjust = 0.5,
    family = "Helvetica"
  ) +
  # Rescue-fraction labels on the left of negative-fidelity genes
  geom_text(
    data = rep_dat %>% filter(Native_DP_fidelity < 0),
    aes(
      x = Native_DP_fidelity - left_label_offset,
      y = y_pos,
      label = rescue_label
    ),
    inherit.aes = FALSE,
    colour = "#71858D",
    size = 1.78,
    hjust = 1,
    vjust = 0.5,
    family = "Helvetica"
  ) +
  scale_y_continuous(
    breaks = rep_dat$y_pos,
    labels = rep_dat$GeneSymbol,
    limits = c(0.15, 10.55),
    expand = c(0, 0)
  ) +
  scale_colour_manual(
    values = group_cols,
    name = NULL
  ) +
  scale_x_continuous(
    name = "Native DP fidelity",
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  labs(
    y = NULL,
    subtitle = "High rescue can accompany divergent native fidelity"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(
    base_size = 6,
    base_family = "Helvetica"
  ) +
  theme(
    plot.subtitle = element_text(
      size = 6.3,
      colour = "#526A74",
      hjust = 0.5,
      margin = margin(b = 3, unit = "pt")
    ),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(
      size = 5.2,
      colour = "#526A74"
    ),
    legend.key.width = unit(2.2, "mm"),
    legend.key.height = unit(2.0, "mm"),
    legend.spacing.x = unit(1.4, "mm"),
    legend.margin = margin(t = 0, r = 0, b = 1, l = 0, unit = "pt"),
    axis.title.x = element_text(
      size = 6.1,
      colour = "#20343D",
      margin = margin(t = 2.5, unit = "pt")
    ),
    axis.text.x = element_text(
      size = 5.0,
      colour = "#718993",
      margin = margin(t = 1.5, unit = "pt")
    ),
    axis.text.y = element_text(
      size = 5.2,
      face = "italic",
      colour = "#20343D",
      margin = margin(r = 2.5, unit = "pt")
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.32,
      colour = "#20343D"
    ),
    axis.ticks.x = element_line(
      linewidth = 0.28,
      colour = "#20343D"
    ),
    axis.ticks.length = unit(1.0, "mm"),
    plot.margin = margin(
      t = 2.0,
      r = 8.0,
      b = 3.0,
      l = 3.0,
      unit = "pt"
    ),
    panel.grid = element_blank()
  )


# ------------------------------------------------------------
# 11. Output
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 2",
  "Panel_D_Nature_final"
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
    "Fig2_PanelD_Nature_final.png"
  ),
  plot = p,
  width = 88,
  height = 56,
  units = "mm",
  dpi = 600,
  bg = "white"
)


# PDF
pdf(
  file = file.path(
    out_dir,
    "Fig2_PanelD_Nature_final.pdf"
  ),
  width = 88 / 25.4,
  height = 56 / 25.4,
  useDingbats = FALSE
)

print(p)

dev.off()


# TIFF
ggsave(
  filename = file.path(
    out_dir,
    "Fig2_PanelD_Nature_final.tiff"
  ),
  plot = p,
  width = 88,
  height = 56,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 12. Export representative genes
# ------------------------------------------------------------
write_csv(
  rep_dat %>%
    dplyr::select(
      GeneSymbol,
      Rescue_fraction,
      Native_DP_fidelity,
      Group
    ),
  file.path(
    out_dir,
    "Fig2_PanelD_Nature_final_representative_genes.csv"
  )
)


# ------------------------------------------------------------
# 13. README / provenance
# ------------------------------------------------------------
readme <- c(
  "Figure 2D — Nature-style final panel; representative high-rescue genes",
  "",
  "Purpose:",
  paste(
    "Illustrate that high transcriptional rescue can coexist with",
    "strongly positive or strongly negative native DP fidelity."
  ),
  "",
  paste(
    "Input:",
    normalizePath(audit_file)
  ),
  paste(
    "Rescue Q3 =",
    sprintf("%.9f", rescue_q3)
  ),
  "",
  "High-fidelity genes:",
  paste(
    upper_candidates,
    collapse = ", "
  ),
  "",
  "Low-fidelity genes:",
  paste(
    lower_candidates,
    collapse = ", "
  ),
  "",
  "Visual changes for final layout:",
  "- Exported at 88 x 56 mm for direct Figure 2 assembly.",
  "- Helvetica typography and strokes sized for final dimensions.",
  "- Compact margins and legend; rescue labels retained.",
  "- Gene set, values, colors, and scientific logic are unchanged."
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig2_PanelD_Nature_final.txt"
  )
)


cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
