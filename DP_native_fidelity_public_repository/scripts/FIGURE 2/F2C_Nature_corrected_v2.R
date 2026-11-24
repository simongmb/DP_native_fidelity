# ============================================================
# Figure 2C | Nature-style corrected panel
# Rescue fraction vs direction-aligned native DP fidelity
#
# This version retains the approved Figure 2C visual structure:
#   - all analyzable genes in a light neutral background
#   - descriptive high-rescue / upper-fidelity group in teal
#   - descriptive high-rescue / lower-fidelity group in orange
#   - subtle rescue-Q3 and fidelity-Q1/Q3 reference lines
#   - Spearman statistic and compact legend above the plotting area
#
# The highlighted regions are descriptive quartile-defined groups only;
# they are NOT the frozen 81/99 programs used in later figures.
# Output is written to a separate folder and only named output files are
# replaced on rerun. No existing project folders are deleted.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(grid)
})

# ------------------------------------------------------------
# 1. Locate project root and input table
# ------------------------------------------------------------
TARGET_NAME <- "02_gene_level_native_fidelity_audit.csv"
PROJECT_FILE <- "HairRegeneration_AI.Rproj"

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/", mustWork = FALSE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

start_dir <- get_script_dir()
wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
candidate_roots <- unique(c(
  start_dir, wd,
  dirname(start_dir), dirname(dirname(start_dir)),
  dirname(dirname(dirname(start_dir))),
  dirname(wd), dirname(dirname(wd)), dirname(dirname(dirname(wd)))
))
candidate_roots <- candidate_roots[dir.exists(candidate_roots)]

score_root <- function(path) {
  score <- 0
  if (file.exists(file.path(path, PROJECT_FILE))) score <- score + 100
  if (dir.exists(file.path(path, "results"))) score <- score + 20
  if (dir.exists(file.path(path, "figures"))) score <- score + 10
  if (file.exists(file.path(path, "results", TARGET_NAME))) score <- score + 50
  score
}

root_scores <- vapply(candidate_roots, score_root, numeric(1))
if (length(root_scores) == 0 || max(root_scores) == 0) {
  stop(
    "Could not identify the HairRegeneration_AI project root. ",
    "Open the project in RStudio and run this script again."
  )
}
project_root <- candidate_roots[which.max(root_scores)]
cat("Project root: ", project_root, "\n", sep = "")

search_dirs <- unique(c(file.path(project_root, "results"), project_root))
search_dirs <- search_dirs[dir.exists(search_dirs)]
matches <- unique(unlist(lapply(search_dirs, function(d) {
  list.files(
    d, pattern = paste0("^", TARGET_NAME, "$"),
    recursive = TRUE, full.names = TRUE, ignore.case = FALSE
  )
}), use.names = FALSE))

if (length(matches) == 0) {
  stop("Cannot find ", TARGET_NAME, ". Searched under:\n",
       paste(search_dirs, collapse = "\n"))
}

preferred <- matches[grepl("DP_native_fidelity_reconstruction", matches,
                           ignore.case = TRUE)]
if (length(preferred) > 0) matches <- preferred
if (length(matches) > 1) {
  stop(
    "More than one candidate audit table was found. To avoid using the wrong ",
    "analysis, specify the intended input or remove duplicate audit tables:\n",
    paste(matches, collapse = "\n")
  )
}

audit_file <- normalizePath(matches[1], winslash = "/", mustWork = TRUE)
cat("Input audit table: ", audit_file, "\n", sep = "")

# ------------------------------------------------------------
# 2. Read data using the primary analysis columns
# ------------------------------------------------------------
dat0 <- read_csv(audit_file, show_col_types = FALSE)

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

gene_col <- pick_col(dat0, c("GeneSymbol", "gene_symbol", "Gene", "symbol"))
rescue_col <- pick_col(dat0, c(
  "Rescue_fraction", "rescue_fraction", "RescueFraction", "rescue.frac", "rescue"
))
fidelity_col <- pick_col(dat0, c(
  "native_fidelity_direction_aligned", "native_fidelity", "Native_fidelity",
  "Native_Fidelity", "direction_aligned_native_fidelity",
  "donor_aware_native_fidelity", "native_DP_fidelity"
))

if (is.na(rescue_col) || is.na(fidelity_col)) {
  stop(
    "Required rescue-fraction or native-fidelity column not found.\n",
    "Available columns:\n", paste(names(dat0), collapse = ", ")
  )
}

dat <- dat0 %>%
  transmute(
    GeneSymbol = if (!is.na(gene_col)) as.character(.data[[gene_col]]) else NA_character_,
    Rescue_fraction = suppressWarnings(as.numeric(.data[[rescue_col]])),
    Native_DP_fidelity = suppressWarnings(as.numeric(.data[[fidelity_col]]))
  ) %>%
  filter(is.finite(Rescue_fraction), is.finite(Native_DP_fidelity))

if (nrow(dat) < 3) {
  stop("Fewer than three complete genes remain; correlation cannot be estimated.")
}

# ------------------------------------------------------------
# 3. Primary statistic and descriptive quartile regions
# ------------------------------------------------------------
sp <- suppressWarnings(cor.test(
  dat$Rescue_fraction, dat$Native_DP_fidelity,
  method = "spearman", exact = FALSE
))
rho <- unname(sp$estimate)
pval <- sp$p.value
n_gene <- nrow(dat)

rescue_q3 <- unname(quantile(dat$Rescue_fraction, probs = 0.75,
                             na.rm = TRUE, type = 7))
fidelity_q1 <- unname(quantile(dat$Native_DP_fidelity, probs = 0.25,
                               na.rm = TRUE, type = 7))
fidelity_q3 <- unname(quantile(dat$Native_DP_fidelity, probs = 0.75,
                               na.rm = TRUE, type = 7))

dat <- dat %>%
  mutate(
    Region = case_when(
      Rescue_fraction >= rescue_q3 & Native_DP_fidelity >= fidelity_q3 ~
        "Upper-fidelity high-rescue",
      Rescue_fraction >= rescue_q3 & Native_DP_fidelity <= fidelity_q1 ~
        "Lower-fidelity high-rescue",
      TRUE ~ "Other"
    )
  )

n_upper <- sum(dat$Region == "Upper-fidelity high-rescue")
n_lower <- sum(dat$Region == "Lower-fidelity high-rescue")
n_other <- sum(dat$Region == "Other")

cat("Spearman rho = ", sprintf("%.6f", rho), "\n", sep = "")
cat("P = ", format(pval, scientific = TRUE), "\n", sep = "")
cat("n = ", format(n_gene, big.mark = ","), "\n", sep = "")
cat("Rescue Q3 = ", sprintf("%.6f", rescue_q3), "\n", sep = "")
cat("Fidelity Q1 = ", sprintf("%.6f", fidelity_q1), "\n", sep = "")
cat("Fidelity Q3 = ", sprintf("%.6f", fidelity_q3), "\n", sep = "")
cat("Highlighted groups: upper = ", n_upper, "; lower = ", n_lower,
    "; other = ", n_other, "\n", sep = "")

# Use a plain-text statistical subtitle. This avoids the malformed chained
# plotmath equality expression that caused source() to stop with "unexpected ==".
stat_expr <- sprintf(
  "Spearman rho = %.3f; P = %.1e; n = %s",
  rho, pval, format(n_gene, big.mark = ",", scientific = FALSE)
)

# ------------------------------------------------------------
# 4. Final panel size, colors, and plotting ranges
# ------------------------------------------------------------
# Match the 88 x 56 mm panel footprint used for Figure 2B.
PANEL_WIDTH_MM <- 88
PANEL_HEIGHT_MM <- 56
DPI <- 600
FONT <- "Helvetica"

INK <- "#263A45"
MUTED_TEXT <- "#718993"
POINT_OTHER <- "#B8C7CD"
POINT_UPPER <- "#3D8C82"
POINT_LOWER <- "#C97857"
ZERO_LINE <- "#D4DDE0"
THRESHOLD_LINE <- "#AAB9BF"

x_upper <- max(2, ceiling(max(dat$Rescue_fraction, na.rm = TRUE) / 0.5) * 0.5)
yr <- range(dat$Native_DP_fidelity, na.rm = TRUE)
ypad <- max(diff(yr) * 0.035, 0.08)

# Plot background points first, then the two highlighted groups.
dat_other <- dat %>% filter(Region == "Other")
dat_highlight <- dat %>%
  filter(Region != "Other") %>%
  mutate(Region = factor(
    Region,
    levels = c("Upper-fidelity high-rescue", "Lower-fidelity high-rescue")
  ))

region_cols <- c(
  "Upper-fidelity high-rescue" = POINT_UPPER,
  "Lower-fidelity high-rescue" = POINT_LOWER
)
legend_labels <- c(
  "Upper-fidelity high-rescue" = "High rescue + upper-fidelity quartile",
  "Lower-fidelity high-rescue" = "High rescue + lower-fidelity quartile"
)

# ------------------------------------------------------------
# 5. Plot — preserve the existing scientific visual encoding
# ------------------------------------------------------------
p <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = ZERO_LINE) +
  geom_vline(
    xintercept = rescue_q3, linewidth = 0.25,
    linetype = "22", colour = THRESHOLD_LINE
  ) +
  geom_hline(
    yintercept = c(fidelity_q1, fidelity_q3), linewidth = 0.25,
    linetype = "22", colour = THRESHOLD_LINE
  ) +
  geom_point(
    data = dat_other,
    aes(x = Rescue_fraction, y = Native_DP_fidelity),
    shape = 16, size = 0.78, alpha = 0.23, colour = POINT_OTHER
  ) +
  geom_point(
    data = dat_highlight,
    aes(x = Rescue_fraction, y = Native_DP_fidelity, colour = Region),
    shape = 16, size = 1.05, alpha = 0.88
  ) +
  scale_colour_manual(
    values = region_cols,
    breaks = c("Upper-fidelity high-rescue", "Lower-fidelity high-rescue"),
    labels = legend_labels,
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1.4, alpha = 1))
  ) +
  scale_x_continuous(
    breaks = seq(0, x_upper, by = 0.5),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = -2:2,
    labels = c("-2", "-1", "0", "1", "2"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = c(0, x_upper + 0.015 * x_upper),
    ylim = c(yr[1] - ypad, yr[2] + ypad),
    clip = "on"
  ) +
  labs(x = "Rescue fraction", y = "Native DP fidelity", subtitle = stat_expr) +
  theme_classic(base_size = 5.5, base_family = FONT) +
  theme(
    plot.subtitle = element_text(
      size = 5.2, colour = INK, hjust = 0.5,
      margin = margin(t = 0, b = 1.2, unit = "pt")
    ),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 4.9, colour = MUTED_TEXT),
    legend.key.width = unit(2.0, "mm"),
    legend.key.height = unit(1.8, "mm"),
    legend.spacing.x = unit(1.4, "mm"),
    legend.margin = margin(t = 0, r = 0, b = 1.0, l = 0, unit = "pt"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    axis.title.x = element_text(
      size = 5.8, colour = INK, margin = margin(t = 1.8, unit = "pt")
    ),
    axis.title.y = element_text(
      size = 5.8, colour = INK, margin = margin(r = 1.8, unit = "pt")
    ),
    axis.text.x = element_text(
      size = 5.0, colour = MUTED_TEXT, margin = margin(t = 0.8, unit = "pt")
    ),
    axis.text.y = element_text(
      size = 5.0, colour = MUTED_TEXT, margin = margin(r = 0.8, unit = "pt")
    ),
    axis.line = element_line(linewidth = 0.32, colour = INK),
    axis.ticks = element_line(linewidth = 0.28, colour = INK),
    axis.ticks.length = unit(0.8, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1.5, r = 2.0, b = 1.5, l = 2.0, unit = "pt")
  )

# ------------------------------------------------------------
# 6. Export to a separate folder; safe to rerun
# ------------------------------------------------------------
out_dir <- file.path(project_root, "figures", "FIGURE 2", "Panel_C_Nature_corrected")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
base_name <- "Fig2_PanelC_Nature_corrected"

# Vector PDF is intended for Illustrator assembly.
ggsave(
  filename = file.path(out_dir, paste0(base_name, ".pdf")),
  plot = p, device = grDevices::pdf,
  width = PANEL_WIDTH_MM, height = PANEL_HEIGHT_MM, units = "mm",
  family = FONT, useDingbats = FALSE, bg = "white"
)

ggsave(
  filename = file.path(out_dir, paste0(base_name, ".png")),
  plot = p, width = PANEL_WIDTH_MM, height = PANEL_HEIGHT_MM,
  units = "mm", dpi = DPI, bg = "white"
)

ggsave(
  filename = file.path(out_dir, paste0(base_name, ".tiff")),
  plot = p, width = PANEL_WIDTH_MM, height = PANEL_HEIGHT_MM,
  units = "mm", dpi = DPI, compression = "lzw", bg = "white"
)

write_csv(dat, file.path(out_dir, paste0(base_name, "_plot_data.csv")))

summary_df <- tibble(
  input_file = audit_file,
  rescue_column = rescue_col,
  fidelity_column = fidelity_col,
  n_total = n_gene,
  spearman_rho = rho,
  spearman_p = pval,
  rescue_Q3 = rescue_q3,
  fidelity_Q1 = fidelity_q1,
  fidelity_Q3 = fidelity_q3,
  n_upper_fidelity_high_rescue = n_upper,
  n_lower_fidelity_high_rescue = n_lower,
  n_other = n_other,
  panel_width_mm = PANEL_WIDTH_MM,
  panel_height_mm = PANEL_HEIGHT_MM,
  point_size_other_mm = 0.78,
  point_size_highlight_mm = 1.05,
  descriptive_quartile_highlights_shown = TRUE,
  frozen_81_99_programs_shown = FALSE
)
write_csv(summary_df, file.path(out_dir, paste0(base_name, "_summary.csv")))

readme <- c(
  "Figure 2C — Nature-style corrected panel",
  "",
  "Primary message: rescue fraction is weakly associated with native DP fidelity.",
  "Visual encoding retained from the existing Figure 2C: neutral background genes;",
  "teal upper-fidelity/high-rescue region; orange lower-fidelity/high-rescue region;",
  "and subtle descriptive quartile thresholds.",
  "The highlighted regions are descriptive and are not the frozen 81/99 programs.",
  "",
  paste("Input:", audit_file),
  paste("Rescue column:", rescue_col),
  paste("Native-fidelity column:", fidelity_col),
  paste("n =", format(n_gene, big.mark = ",")),
  paste("Spearman rho =", sprintf("%.6f", rho)),
  paste("Spearman P =", format(pval, scientific = TRUE)),
  paste("Rescue Q3 =", sprintf("%.6f", rescue_q3)),
  paste("Fidelity Q1 =", sprintf("%.6f", fidelity_q1)),
  paste("Fidelity Q3 =", sprintf("%.6f", fidelity_q3)),
  paste("Upper-fidelity/high-rescue n =", n_upper),
  paste("Lower-fidelity/high-rescue n =", n_lower),
  "",
  "Panel footprint: 88 x 56 mm.",
  "Helvetica typography; compact text sized for final placed size.",
  "PDF is vector output for Illustrator. Existing folders/files are not deleted.",
  paste("Output directory:", out_dir)
)
writeLines(readme, file.path(out_dir, paste0(base_name, "_README.txt")))

cat("\nFigure 2C export complete.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Panel size: ", PANEL_WIDTH_MM, " x ", PANEL_HEIGHT_MM, " mm\n", sep = "")
