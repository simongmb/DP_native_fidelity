
# ============================================================
# Figure 7 v20 — DONOR-BLOCKED FINAL CLEAN-SESSION LOCK
# Nature-style conceptual model
# Central visual argument:
# anatomy-referenced state balance + broader non-programme encoding + recurrent compact model core
# ============================================================

suppressPackageStartupMessages(library(grid))

# -----------------------------
# Palette
# -----------------------------
COL_TEXT   <- "#20343D"
COL_MUTED  <- "#6F8189"
COL_RULE   <- "#AAB8BD"
COL_LIGHT  <- "#DCE4E7"
COL_NATIVE <- "#3D8C82"
COL_ALT    <- "#C97857"
COL_BG     <- "white"

# -----------------------------
# Helpers
# -----------------------------
tg <- function(x, y, label, size = 5.2, face = "plain",
               col = COL_TEXT, just = "centre") {
  textGrob(
    label,
    x = unit(x, "npc"),
    y = unit(y, "npc"),
    just = just,
    gp = gpar(
      fontsize = size,
      fontface = face,
      col = col,
      fontfamily = "Helvetica"
    )
  )
}

sg <- function(x1, y1, x2, y2, col = COL_RULE,
               lwd = 0.55, arrow_on = FALSE) {
  segmentsGrob(
    x0 = unit(x1, "npc"), y0 = unit(y1, "npc"),
    x1 = unit(x2, "npc"), y1 = unit(y2, "npc"),
    gp = gpar(col = col, lwd = lwd),
    arrow = if (arrow_on) {
      grid::arrow(length = unit(0.010, "npc"), type = "closed")
    } else NULL
  )
}

add <- function(lst, g) {
  lst[[length(lst) + 1]] <- g
  lst
}

# -----------------------------
# Figure body
# -----------------------------
g <- list()
g <- add(g, rectGrob(gp = gpar(fill = COL_BG, col = NA)))

# ============================================================
# 1. Programme anchors
# ============================================================
g <- add(g, tg(0.19, 0.915, "Native-aligned rescue programme",
               6.2, "bold", COL_NATIVE))
g <- add(g, tg(0.19, 0.884, "81 genes",
               5.0, "plain", COL_MUTED))

g <- add(g, tg(0.81, 0.915, "Alternative-state rescue programme",
               6.2, "bold", COL_ALT))
g <- add(g, tg(0.81, 0.884, "99 genes",
               5.0, "plain", COL_MUTED))

# ============================================================
# 2. Central continuum + distributed encoding as ONE visual
# ============================================================
g <- add(g, tg(0.50, 0.830,
               "Anatomy-referenced state balance",
               6.7, "bold"))

# wide continuum band
g <- add(g, rectGrob(
  x = unit(0.50, "npc"),
  y = unit(0.765, "npc"),
  width = unit(0.64, "npc"),
  height = unit(0.032, "npc"),
  gp = gpar(fill = "#F1F4F5", col = NA)
))
g <- add(g, rectGrob(
  x = unit(0.34, "npc"),
  y = unit(0.765, "npc"),
  width = unit(0.32, "npc"),
  height = unit(0.032, "npc"),
  gp = gpar(fill = COL_NATIVE, col = NA)
))
g <- add(g, rectGrob(
  x = unit(0.66, "npc"),
  y = unit(0.765, "npc"),
  width = unit(0.32, "npc"),
  height = unit(0.032, "npc"),
  gp = gpar(fill = COL_ALT, col = NA)
))

g <- add(g, tg(0.18, 0.725, "Native-aligned",
               5.0, "bold", COL_NATIVE))
g <- add(g, tg(0.82, 0.725, "Alternative-aligned",
               5.0, "bold", COL_ALT))

# feeders
g <- add(g, sg(0.27, 0.865, 0.34, 0.790, COL_NATIVE, 0.65))
g <- add(g, sg(0.73, 0.865, 0.66, 0.790, COL_ALT, 0.65))

# non-programme association title
g <- add(g, tg(0.50, 0.665,
               "Broader non-programme transcriptomic association",
               6.5, "bold"))
g <- add(g, tg(0.50, 0.635,
               "1,373 non-programme genes",
               5.0, "plain", COL_MUTED))

# connect continuum to distributed field
g <- add(g, sg(0.50, 0.700, 0.50, 0.655, COL_RULE, 0.55, TRUE))

# distributed architecture field:
# broad layered band with many short marks, but integrated with continuum
field_y <- 0.535
field_h <- 0.150
field_x1 <- 0.20
field_x2 <- 0.80

g <- add(g, roundrectGrob(
  x = unit(0.50, "npc"),
  y = unit(field_y, "npc"),
  width = unit(field_x2 - field_x1, "npc"),
  height = unit(field_h, "npc"),
  r = unit(0.008, "snpc"),
  gp = gpar(fill = "#FCFDFD", col = adjustcolor(COL_LIGHT, alpha.f = 0.75), lwd = 0.42)
))

g <- add(g, tg(0.50, 0.585,
               "Broader non-programme transcriptomic association",
               6.3, "bold"))

# deterministic micro-marks suggesting breadth only
xx <- seq(field_x1 + 0.022, field_x2 - 0.022, length.out = 56)
for (i in seq_along(xx)) {
  # Two vertically staggered rows to avoid barcode appearance
  y0 <- if (i %% 2 == 0) 0.512 else 0.488
  len <- 0.018 + 0.006 * (0.5 + 0.5 * sin(i * 1.22))
  col_i <- if (xx[i] < 0.42) {
    adjustcolor(COL_NATIVE, alpha.f = 0.38)
  } else if (xx[i] > 0.58) {
    adjustcolor(COL_ALT, alpha.f = 0.38)
  } else {
    adjustcolor(COL_RULE, alpha.f = 0.42)
  }
  g <- add(g, sg(xx[i], y0 - len/2, xx[i], y0 + len/2, col_i, 0.38))
}

g <- add(g, tg(0.50, 0.445,
               "Non-programme transcriptomic variation strongly reflects the prespecified state balance",
               5.0, "plain", COL_MUTED))

# ============================================================
# 3. Evidence labels embedded around the central field
# ============================================================
g <- add(g, sg(0.33, 0.468, 0.27, 0.402, COL_NATIVE, 0.50))
g <- add(g, sg(0.50, 0.458, 0.50, 0.402, COL_RULE, 0.50))
g <- add(g, sg(0.67, 0.468, 0.73, 0.402, COL_ALT, 0.50))

g <- add(g, tg(0.22, 0.372, "Recurrent model features",
               6.0, "bold", COL_NATIVE))
g <- add(g, tg(0.22, 0.344, "15 genes selected in all 200 outer folds",
               5.1, "plain", adjustcolor(COL_TEXT, alpha.f = 0.72)))

g <- add(g, tg(0.50, 0.372, "Marginal–multivariate concordance",
               6.0, "bold", COL_TEXT))
g <- add(g, tg(0.50, 0.344, "Absolute rho versus model prominence: rho = 0.73",
               5.1, "plain", adjustcolor(COL_TEXT, alpha.f = 0.72)))

g <- add(g, tg(0.78, 0.372, "Nested top-25 signature",
               6.0, "bold", COL_ALT))
g <- add(g, tg(0.78, 0.344, "Nested top-25 retains most state information",
               5.1, "plain", adjustcolor(COL_TEXT, alpha.f = 0.72)))

# ============================================================
# 4. Quantitative support
# ============================================================
g <- add(g, sg(0.30, 0.305, 0.70, 0.305, COL_LIGHT, 0.55))

g <- add(g, tg(0.35, 0.277, "Full model Spearman rho = 0.977",
               5.5, "bold"))
g <- add(g, tg(0.50, 0.277, "CV R² = 0.982",
               5.5, "bold"))
g <- add(g, tg(0.65, 0.277, "Fixed-prediction permutation P = 1 x 10^-4",
               4.9, "plain", COL_MUTED))

# ============================================================
# 5. Final compression test / conclusion
# ============================================================
g <- add(g, sg(0.50, 0.245, 0.50, 0.228, COL_RULE, 0.48, TRUE))

g <- add(g, roundrectGrob(
  x = unit(0.50, "npc"),
  y = unit(0.178, "npc"),
  width = unit(0.48, "npc"),
  height = unit(0.078, "npc"),
  r = unit(0.008, "snpc"),
  gp = gpar(fill = "#FAFBFB", col = COL_LIGHT, lwd = 0.50)
))

g <- add(g, tg(0.50, 0.196,
               "State-balance information is substantially compressible",
               6.5, "bold"))

g <- add(g, tg(0.50, 0.165,
               "Nested top-25: rho = 0.958, R² = 0.960   vs   full model: rho = 0.977, R² = 0.982",
               5.0, "plain", COL_MUTED))

fig <- gTree(children = do.call(gList, g))

# -----------------------------
# Output
# -----------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  }
  getwd()
}

start_dir <- get_script_dir()

candidate_roots <- unique(c(
  start_dir,
  getwd(),
  dirname(start_dir),
  dirname(dirname(start_dir)),
  dirname(getwd())
))
candidate_roots <- candidate_roots[dir.exists(candidate_roots)]

score_root <- function(x) {
  s <- 0
  if (dir.exists(file.path(x, "figures"))) s <- s + 5
  if (dir.exists(file.path(x, "results"))) s <- s + 2
  s
}

scores <- vapply(candidate_roots, score_root, numeric(1))
project_root <- candidate_roots[which.max(scores)]

out_dir <- file.path(
  project_root,
  "figures",
  "FIGURE 7",
  "Integrative_Model_Nature_v20_donorBlocked"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

W_IN <- 183 / 25.4
H_IN <- 108 / 25.4

pdf_file <- file.path(out_dir, "Fig7_integrativeModel_Nature_v20_donorBlocked.pdf")
png_file <- file.path(out_dir, "Fig7_integrativeModel_Nature_v20_donorBlocked.png")
tif_file <- file.path(out_dir, "Fig7_integrativeModel_Nature_v20_donorBlocked.tif")

# PDF
grDevices::pdf(
  file = pdf_file,
  width = W_IN,
  height = H_IN,
  family = "Helvetica",
  useDingbats = FALSE,
  onefile = FALSE,
  paper = "special"
)
grid.newpage()
grid.draw(fig)
dev.off()

# PNG
grDevices::png(
  filename = png_file,
  width = W_IN,
  height = H_IN,
  units = "in",
  res = 600,
  bg = "white"
)
grid.newpage()
grid.draw(fig)
dev.off()

# TIF
grDevices::tiff(
  filename = tif_file,
  width = W_IN,
  height = H_IN,
  units = "in",
  res = 600,
  compression = "none",
  bg = "white"
)
grid.newpage()
grid.draw(fig)
dev.off()


# Archive README
readme_file <- file.path(out_dir, "README_Fig7_v20_donorBlocked.txt")
writeLines(
  c(
    "Figure 7 v20 — donor-blocked final clean-session lock",
    "",
    "Conceptual structure:",
    "- Native-aligned rescue programme: 81 genes",
    "- Alternative-state rescue programme: 99 genes",
    "- Anatomy-referenced state balance",
    "- 1,373 non-programme genes",
    "- Recurrent model features",
    "- Strong marginal–multivariate concordance",
    "- Nested top-25 compact signature",
    "",
    "Updated quantitative support:",
    "- Full non-programme model: Spearman rho = 0.977",
    "- Full non-programme model: CV R2 = 0.982",
    "- Fixed-prediction within-donor permutation P = 1 x 10^-4",
    "- Nested top-25 reconstruction: Spearman rho = 0.958",
    "- Nested top-25 reconstruction: CV R2 = 0.960",
    "- 15 genes selected in all 200 outer-fold opportunities",
    "",
    "Obsolete elements removed:",
    "- co-expression module branch",
    "- 54 / 300 pairwise-correlation claim",
    "- old rho = 0.39 / R2 = 0.10 metrics",
    "- old nested top-25 rho = 0.24 comparison"
  ),
  readme_file
)

cat(
  "\n============================================\n",
  "FIGURE 7 v20 DONOR-BLOCKED COMPLETE\n",
  "============================================\n\n",
  "PDF:\n", pdf_file, "\n\n",
  "PNG:\n", png_file, "\n\n",
  "TIF:\n", tif_file, "\n"
)
