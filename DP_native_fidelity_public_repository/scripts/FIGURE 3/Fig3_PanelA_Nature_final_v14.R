
# ============================================================
# Figure 3A Nature final v14 — Frozen definition of the 81 / 99 programmes
# Symmetric two-branch layout
#
# Purpose:
#   Show how the frozen 81-gene and 99-gene programmes were defined.
#
# Main changes from v2:
#   - Make the left and right branches vertically symmetric
#   - Keep final 81 / 99 programme boxes on the same baseline
#   - Move the six-gene exclusion note into the right branch
#     without pushing the final orange box downward
#   - Keep typography hierarchy consistent with Figure 2
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
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
  "Rescue_fraction",
  "native_fidelity_direction_aligned",
  "high_rescue",
  "high_native_fidelity",
  "low_native_fidelity",
  "anatomy_consistent_candidate",
  "non_native_candidate_pre_exclusion",
  "historical_proliferation_exclusion",
  "program"
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
# 5. Frozen thresholds
# ------------------------------------------------------------
rescue_q3 <- unname(
  quantile(
    dat$Rescue_fraction,
    probs = 0.75,
    na.rm = TRUE,
    type = 7
  )
)

fidelity_q1 <- unname(
  quantile(
    dat$native_fidelity_direction_aligned,
    probs = 0.25,
    na.rm = TRUE,
    type = 7
  )
)

fidelity_q3 <- unname(
  quantile(
    dat$native_fidelity_direction_aligned,
    probs = 0.75,
    na.rm = TRUE,
    type = 7
  )
)


# ------------------------------------------------------------
# 6. Logical conversion helper
# ------------------------------------------------------------
as_flag <- function(x) {
  if (is.logical(x)) return(x)

  x_chr <- toupper(trimws(as.character(x)))

  x_chr %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "Y"
  )
}

dat <- dat %>%
  mutate(
    high_rescue_flag =
      as_flag(high_rescue),

    high_fidelity_flag =
      as_flag(high_native_fidelity),

    low_fidelity_flag =
      as_flag(low_native_fidelity),

    anatomy_candidate_flag =
      as_flag(anatomy_consistent_candidate),

    low_candidate_flag =
      as_flag(non_native_candidate_pre_exclusion),

    proliferation_exclusion_flag =
      as_flag(historical_proliferation_exclusion)
  )


# ------------------------------------------------------------
# 7. Counts
# ------------------------------------------------------------
n_total <- nrow(dat)

n_high_rescue <- sum(
  dat$high_rescue_flag,
  na.rm = TRUE
)

n_native_candidate <- sum(
  dat$anatomy_candidate_flag,
  na.rm = TRUE
)

n_low_candidate <- sum(
  dat$low_candidate_flag,
  na.rm = TRUE
)

n_excluded <- sum(
  dat$low_candidate_flag &
    dat$proliferation_exclusion_flag,
  na.rm = TRUE
)

n_native_final <- n_native_candidate

n_alternative_final <- sum(
  dat$low_candidate_flag &
    !dat$proliferation_exclusion_flag,
  na.rm = TRUE
)



# ------------------------------------------------------------
# 7b. Reconstruct frozen programme rules for QC
# ------------------------------------------------------------
aligned_cols <- c(
  "aligned_DP_vs_DC",
  "aligned_DP_vs_CTS",
  "aligned_DP_vs_HM1",
  "aligned_DP_vs_HM2",
  "aligned_DP_vs_IRS1",
  "aligned_DP_vs_ORS1",
  "aligned_DP_vs_ORS2"
)

if (all(aligned_cols %in% names(dat))) {
  dat <- dat %>%
    mutate(
      anatomy_consistency_n_of_7_recalc =
        rowSums(across(all_of(aligned_cols), ~ .x > 0), na.rm = TRUE)
    )

  native_rule_recalc <-
    dat$Rescue_fraction > rescue_q3 &
    dat$native_fidelity_direction_aligned > fidelity_q3 &
    dat$anatomy_consistency_n_of_7_recalc >= 5

  low_rule_recalc <-
    dat$Rescue_fraction > rescue_q3 &
    dat$native_fidelity_direction_aligned < fidelity_q1 &
    dat$anatomy_consistency_n_of_7_recalc <= 2

  if (sum(native_rule_recalc, na.rm = TRUE) != 81) {
    stop("Frozen Native-aligned rescue rule does not reproduce n = 81.")
  }

  if (sum(low_rule_recalc, na.rm = TRUE) != 105) {
    stop("Frozen low-fidelity pre-exclusion rule does not reproduce n = 105.")
  }
}

# ------------------------------------------------------------
# 8. Frozen exclusion list
# ------------------------------------------------------------
expected_exclusions <- c(
  "ATAD2",
  "BIRC5",
  "CCNA2",
  "CENPA",
  "CIP2A",
  "MCM10"
)

observed_exclusions <- dat %>%
  filter(
    low_candidate_flag,
    proliferation_exclusion_flag
  ) %>%
  pull(GeneSymbol) %>%
  unique() %>%
  sort()


# ------------------------------------------------------------
# 9. QC hard stops
# ------------------------------------------------------------
if (n_total != 1553) {
  stop(
    paste0(
      "Expected 1,553 analyzable genes, but found ",
      n_total,
      "."
    )
  )
}

if (n_native_final != 81) {
  stop(
    paste0(
      "Expected 81 native-identity genes, but found ",
      n_native_final,
      "."
    )
  )
}

if (n_low_candidate != 105) {
  stop(
    paste0(
      "Expected 105 low-fidelity candidates, but found ",
      n_low_candidate,
      "."
    )
  )
}

if (n_excluded != 6) {
  stop(
    paste0(
      "Expected 6 proliferation exclusions, but found ",
      n_excluded,
      "."
    )
  )
}

if (n_alternative_final != 99) {
  stop(
    paste0(
      "Expected 99 alternative-state genes, but found ",
      n_alternative_final,
      "."
    )
  )
}

if (!setequal(observed_exclusions, expected_exclusions)) {
  stop(
    paste0(
      "Observed proliferation exclusions do not match frozen list.\n",
      "Observed: ",
      paste(observed_exclusions, collapse = ", "),
      "\nExpected: ",
      paste(expected_exclusions, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 10. Colors and final-size typography
# ------------------------------------------------------------
COL_TEAL <- "#3D8C82"
COL_ORANGE <- "#C97857"
COL_TEXT <- "#20343D"
COL_MUTED <- "#718993"
COL_BORDER <- "#BFCED4"
COL_GRAY_FILL <- "#F5F8F9"
COL_TEAL_FILL <- "#EDF6F4"
COL_ORANGE_FILL <- "#FBF1ED"
COL_EXCLUDE_FILL <- "#F7F4F2"
FONT_FAMILY <- "Helvetica"

# ggplot text size is in mm. These values correspond approximately to
# 5.0–6.5 pt at the FINAL exported panel size.
TXT_TITLE <- 2.25
TXT_BODY  <- 2.00
TXT_SMALL <- 1.78

# Final panel dimensions for a two-column Nature-style assembly.
WIDTH_MM  <- 88
HEIGHT_MM <- 72
DPI       <- 600

# ------------------------------------------------------------
# 11. Figure
# ------------------------------------------------------------
p <- ggplot() +

  # ==========================================================
  # TOP BOX
  # ==========================================================
  annotate(
    "rect",
    xmin = 0.235, xmax = 0.765,
    ymin = 0.890, ymax = 0.975,
    fill = COL_GRAY_FILL,
    colour = COL_BORDER,
    linewidth = 0.28
  ) +
  annotate(
    "text",
    x = 0.50, y = 0.945,
    label = "Analyzable strict-rescue genes",
    size = TXT_TITLE,
    family = FONT_FAMILY,
    colour = COL_TEXT,
    fontface = "bold",
    lineheight = 0.90
  ) +
  annotate(
    "text",
    x = 0.50, y = 0.908,
    label = paste0("n = ", format(n_total, big.mark = ",")),
    size = TXT_SMALL,
    family = FONT_FAMILY,
    colour = COL_MUTED
  ) +

  # ==========================================================
  # ARROW TO HIGH-RESCUE BOX
  # ==========================================================
  annotate(
    "segment",
    x = 0.50, xend = 0.50,
    y = 0.885, yend = 0.845,
    linewidth = 0.30,
    colour = "#93A6AE",
    arrow = arrow(length = unit(1.10, "mm"), type = "closed")
  ) +

  # ==========================================================
  # HIGH-RESCUE BOX
  # ==========================================================
  annotate(
    "rect",
    xmin = 0.22, xmax = 0.78,
    ymin = 0.705, ymax = 0.820,
    fill = "#FAFBFC",
    colour = "#C7D3D8",
    linewidth = 0.26
  ) +
  annotate(
    "text",
    x = 0.50, y = 0.782,
    label = "High rescue (upper quartile)",
    size = TXT_BODY,
    family = FONT_FAMILY,
    colour = COL_TEXT,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.50, y = 0.738,
    label = paste0("Rescue fraction > ", sprintf("%.3f", rescue_q3), "; n = ", n_high_rescue),
    size = TXT_SMALL,
    family = FONT_FAMILY,
    colour = COL_MUTED
  ) +

  # ==========================================================
  # SYMMETRIC BRANCH ARROWS
  # ==========================================================
  annotate(
    "segment",
    x = 0.43, xend = 0.22,
    y = 0.700, yend = 0.635,
    linewidth = 0.32,
    colour = COL_TEAL,
    arrow = arrow(length = unit(1.45, "mm"), type = "closed")
  ) +
  annotate(
    "segment",
    x = 0.57, xend = 0.78,
    y = 0.700, yend = 0.635,
    linewidth = 0.32,
    colour = COL_ORANGE,
    arrow = arrow(length = unit(1.45, "mm"), type = "closed")
  ) +

  # ==========================================================
  # LEFT RULE BOX
  # Equal top/bottom padding and fixed 0.050 vertical text spacing.
  # ==========================================================
  annotate(
    "rect",
    xmin = 0.000, xmax = 0.440,
    ymin = 0.410, ymax = 0.620,
    fill = COL_TEAL_FILL,
    colour = COL_TEAL,
    linewidth = 0.30
  ) +
  annotate(
    "text",
    x = 0.22, y = 0.565,
    label = "Native-aligned branch",
    size = TXT_TITLE,
    family = FONT_FAMILY,
    colour = COL_TEAL,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.22, y = 0.500,
    label = paste0("Native fidelity > ", sprintf("%.3f", fidelity_q3)),
    size = TXT_BODY,
    family = FONT_FAMILY,
    colour = COL_TEXT
  ) +
  annotate(
    "text",
    x = 0.22, y = 0.438,
    label = "Anatomy consistency >= 5/7",
    size = TXT_SMALL,
    family = FONT_FAMILY,
    colour = COL_MUTED
  ) +

  # ==========================================================
  # RIGHT RULE BOX
  # EXACTLY mirrors the left rule box geometry.
  # ==========================================================
  annotate(
    "rect",
    xmin = 0.560, xmax = 1.000,
    ymin = 0.410, ymax = 0.620,
    fill = COL_ORANGE_FILL,
    colour = COL_ORANGE,
    linewidth = 0.30
  ) +
  annotate(
    "text",
    x = 0.78, y = 0.565,
    label = "Low-fidelity branch",
    size = TXT_TITLE,
    family = FONT_FAMILY,
    colour = COL_ORANGE,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.78, y = 0.500,
    label = paste0("Native fidelity < ", sprintf("%.3f", fidelity_q1)),
    size = TXT_BODY,
    family = FONT_FAMILY,
    colour = COL_TEXT
  ) +
  annotate(
    "text",
    x = 0.78, y = 0.438,
    label = paste0("Anatomy consistency <= 2/7; n = ", n_low_candidate),
    size = TXT_SMALL,
    family = FONT_FAMILY,
    colour = COL_MUTED
  ) +

  # ==========================================================
  # LEFT ARROW TO FINAL PROGRAM
  # ==========================================================
  annotate(
    "segment",
    x = 0.22, xend = 0.22,
    y = 0.398, yend = 0.235,
    linewidth = 0.32,
    colour = COL_TEAL,
    arrow = arrow(length = unit(1.45, "mm"), type = "closed")
  ) +

  # ==========================================================
  # RIGHT EXCLUSION STEP
  # Exactly two text rows: title + one-line gene list.
  # Arrow shafts and arrowheads are kept fully outside neighbouring boxes.
  # ==========================================================
  annotate(
    "segment",
    x = 0.78, xend = 0.78,
    y = 0.398, yend = 0.365,
    linewidth = 0.30,
    colour = COL_ORANGE,
    arrow = arrow(length = unit(1.35, "mm"), type = "closed")
  ) +
  annotate(
    "rect",
    xmin = 0.505, xmax = 1.055,
    ymin = 0.225, ymax = 0.355,
    fill = COL_EXCLUDE_FILL,
    colour = "#C9BFB9",
    linewidth = 0.25
  ) +
  annotate(
    "text",
    x = 0.78, y = 0.315,
    label = "Exclude proliferation genes (n = 6)",
    size = TXT_SMALL,
    family = FONT_FAMILY,
    colour = COL_TEXT,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.78, y = 0.262,
    label = "ATAD2   BIRC5   CCNA2   CENPA   CIP2A   MCM10",
    size = TXT_SMALL,
    family = FONT_FAMILY,
    colour = COL_MUTED
  ) +
  annotate(
    "segment",
    x = 0.78, xend = 0.78,
    y = 0.216, yend = 0.187,
    linewidth = 0.32,
    colour = COL_ORANGE,
    arrow = arrow(length = unit(1.35, "mm"), type = "closed")
  ) +

  # ==========================================================
  # FINAL PROGRAM BOXES
  # Same width, height, baseline and identical internal line spacing.
  # ==========================================================
  annotate(
    "rect",
    xmin = 0.005, xmax = 0.405,
    ymin = 0.028, ymax = 0.175,
    fill = "#E4F2EF",
    colour = COL_TEAL,
    linewidth = 0.34
  ) +
  annotate(
    "text",
    x = 0.205, y = 0.138,
    label = "Native-aligned rescue",
    size = TXT_TITLE,
    family = FONT_FAMILY,
    colour = COL_TEAL,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.205, y = 0.100,
    label = "programme",
    size = TXT_BODY,
    family = FONT_FAMILY,
    colour = COL_TEXT
  ) +
  annotate(
    "text",
    x = 0.205, y = 0.061,
    label = paste0("n = ", n_native_final),
    size = TXT_BODY,
    family = FONT_FAMILY,
    colour = COL_MUTED,
    fontface = "bold"
  ) +

  annotate(
    "rect",
    xmin = 0.595, xmax = 0.995,
    ymin = 0.028, ymax = 0.175,
    fill = "#F8E9E3",
    colour = COL_ORANGE,
    linewidth = 0.34
  ) +
  annotate(
    "text",
    x = 0.795, y = 0.138,
    label = "Alternative-state rescue",
    size = TXT_TITLE,
    family = FONT_FAMILY,
    colour = COL_ORANGE,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.795, y = 0.100,
    label = "programme",
    size = TXT_BODY,
    family = FONT_FAMILY,
    colour = COL_TEXT
  ) +
  annotate(
    "text",
    x = 0.795, y = 0.061,
    label = paste0("n = ", n_alternative_final),
    size = TXT_BODY,
    family = FONT_FAMILY,
    colour = COL_MUTED,
    fontface = "bold"
  ) +

  # ==========================================================
  # FOOTNOTE
  # ==========================================================
  annotate(
    "text",
    x = 0.50, y = 0.008,
    label = "Programmes were frozen before downstream analyses",
    size = TXT_SMALL,
    family = FONT_FAMILY,
    colour = COL_MUTED
  ) +

  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  theme_void(base_family = FONT_FAMILY) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 4, r = 5, b = 4, l = 5, unit = "pt")
  )


# ------------------------------------------------------------
# 12. Output
# ------------------------------------------------------------
out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 3",
  "Panel_A_Nature_final_v14"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# PNG
ggsave(
  filename = file.path(out_dir, "Fig3_PanelA_Nature_final_v14.png"),
  plot = p,
  width = WIDTH_MM,
  height = HEIGHT_MM,
  units = "mm",
  dpi = DPI,
  bg = "white"
)

# PDF
pdf(
  file = file.path(out_dir, "Fig3_PanelA_Nature_final_v14.pdf"),
  width = WIDTH_MM / 25.4,
  height = HEIGHT_MM / 25.4,
  family = FONT_FAMILY,
  useDingbats = FALSE
)
print(p)
dev.off()

# TIFF
ggsave(
  filename = file.path(out_dir, "Fig3_PanelA_Nature_final_v14.tiff"),
  plot = p,
  width = WIDTH_MM,
  height = HEIGHT_MM,
  units = "mm",
  dpi = DPI,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 13. QC export
# ------------------------------------------------------------
qc_df <- tibble(
  metric = c(
    "Analyzable strict-rescue genes",
    "Rescue_fraction Q3",
    "Native fidelity Q1",
    "Native fidelity Q3",
    "High-rescue genes",
    "Native-aligned rescue final programme",
    "Low-fidelity pre-exclusion candidates",
    "Proliferation exclusions",
    "Alternative-state rescue final programme"
  ),

  value = c(
    as.character(n_total),
    sprintf("%.9f", rescue_q3),
    sprintf("%.9f", fidelity_q1),
    sprintf("%.9f", fidelity_q3),
    as.character(n_high_rescue),
    as.character(n_native_final),
    as.character(n_low_candidate),
    as.character(n_excluded),
    as.character(n_alternative_final)
  )
)

write_csv(
  qc_df,
  file.path(
    out_dir,
    "Fig3_PanelA_v3_QC.csv"
  )
)


# ------------------------------------------------------------
# 14. README
# ------------------------------------------------------------
readme <- c(

  "Figure 3A Nature final v14 — Frozen definition of the 81 / 99 programmes",

  "",

  paste(
    "Input:",
    normalizePath(audit_file)
  ),

  "",

  "Nature-style layout refinements:",

  "- Left and right branches are vertically symmetric and share common baselines.",
  "- All boxes use increased horizontal padding; long labels were shortened or split across lines.",
  "- The proliferation-exclusion box is separated from the final 99-gene box with no overlap.",
  "- Final 81- and 99-gene boxes use identical geometry and typography.",
  "- Helvetica typography and line weights are tuned for final-size Nature-style assembly.",
  "- Scientific definition updated to match the frozen reconstruction: strict >/< threshold rules and explicit >=5/7 vs <=2/7 anatomical-consistency criteria."
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_Fig3_PanelA_Nature_final_v14.txt"
  )
)


cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
