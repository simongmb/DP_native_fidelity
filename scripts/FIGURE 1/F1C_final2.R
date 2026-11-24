# ============================================================================
# HairRegAI | FIGURE 1 | PANEL C | v3 vertical spacing refinement
# Save as: scripts/FIGURE 1/Fig1_PanelC.R
# In RStudio Console: source("scripts/FIGURE 1/Fig1_PanelC.R")
# Output: figures/FIGURE 1/Panel_C/
# No panel letter; no Figure assembly. Base grid + ragg; vector PDF.
# Existing outputs are NEVER overwritten.
# ============================================================================

# -------------------------- EDITABLE SETTINGS -------------------------------
WIDTH_MM  <- 88.0
HEIGHT_MM <- 64.0
DPI       <- 600L
FONT      <- "Helvetica"
BASENAME  <- "Fig1C_RescueSelection_Nature_final2"

INPUT   <- file.path("results", "GSE44765_rescue",
                     "GSE44765_rescue_analysis_all_genes.csv")
OUT_DIR <- file.path("figures", "FIGURE 1", "Panel_C_Nature_final2")

# True physical font sizes (points), NOT scaled from a large original panel.
PT_HEADING  <- 6.2
PT_NUMBER   <- 7.4
PT_SMALL    <- 5.4
PT_CONTEXT  <- 5.4
PT_OVERLAP  <- 7.3
PT_STRICT   <- 7.8
PT_LEGEND   <- 5.4

COL <- list(
  ink   = "#253945",  muted = "#748A95",  line = "#D8E3E7",
  track = "#EDF2F4",  twoD  = "#819AAA",  threeD = "#B9825E",
  up    = "#B9825E",  down = "#087F79",  other = "#CAD5DB",
  white = "#FFFFFF"
)
# ---------------------------------------------------------------------------

# 1. Source data and exact screening semantics.
if (!file.exists(INPUT)) {
  stop("Missing source file: ", INPUT,
       "\nOpen the HairRegeneration_AI project in RStudio first.")
}
if (dir.exists(OUT_DIR) || file.exists(OUT_DIR)) {
  stop("Previous Panel_C is protected: ", OUT_DIR,
       "\nRename it to Panel_C_previous before generating the new version.")
}
if (!requireNamespace("ragg", quietly = TRUE)) {
  stop("The package 'ragg' is required for 600 dpi PNG and LZW TIFF.\n",
       "Install once with install.packages('ragg').")
}

raw <- utils::read.csv(INPUT, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("GeneSymbol", "Altered_in_2D", "Changed_in_3D",
            "Same_direction", "Closer_to_Intact", "Strict_rescue", "Rescue_type")
if (!all(needed %in% names(raw))) {
  stop("Missing columns: ", paste(setdiff(needed, names(raw)), collapse = ", "))
}
if (anyNA(raw$GeneSymbol) || anyDuplicated(raw$GeneSymbol)) {
  stop("Missing or duplicated gene identifiers. Review original analysis.")
}
flag <- function(z, name) {
  if (anyNA(z)) stop("NA in source flag: ", name)
  if (is.logical(z)) return(z)
  zz <- toupper(trimws(as.character(z)))
  if (!all(zz %in% c("TRUE", "FALSE", "T", "F"))) {
    stop("Invalid TRUE/FALSE values in: ", name)
  }
  zz %in% c("TRUE", "T")
}
A <- flag(raw$Altered_in_2D, "Altered_in_2D")
B <- flag(raw$Changed_in_3D, "Changed_in_3D")
D <- flag(raw$Same_direction, "Same_direction")
C <- flag(raw$Closer_to_Intact, "Closer_to_Intact")
S <- flag(raw$Strict_rescue, "Strict_rescue")

N     <- nrow(raw)
n2    <- sum(A)
n3    <- sum(B)
nboth <- sum(A & B)
ns    <- sum(S)
nu    <- sum(S & raw$Rescue_type == "Restored_upward", na.rm = TRUE)
nd    <- sum(S & raw$Rescue_type == "Restored_downward", na.rm = TRUE)
nother <- nboth - ns

if (any(S & !(A & B & D & C))) {
  stop("Strict rescue flags violate the source screening rules.")
}
if (nu + nd != ns) stop("Strict up/down gene counts do not sum to strict total.")
expected <- c(total = 21355L, twoD = 7804L, threeD = 3143L,
              overlap = 1969L, strict = 1571L, up = 510L,
              down = 1061L, other = 398L)
actual <- c(total = N, twoD = n2, threeD = n3, overlap = nboth,
            strict = ns, up = nu, down = nd, other = nother)
if (!identical(as.integer(actual), as.integer(expected))) {
  stop("Source counts have changed; review input before plotting:\n",
       paste(names(actual), actual, sep = " = ", collapse = "; "))
}
fmt <- function(n) format(n, trim = TRUE, scientific = FALSE, big.mark = ",")
pct <- function(n, denom) sprintf("%.1f%%", 100*n/denom)

# 2. Physical-mm layout. All geometry and labels are drawn at FINAL print size.
mm <- function(v) grid::unit(v, "native")
txt <- function(s, x, y, pt = PT_SMALL, col = COL$ink,
                bold = FALSE, just = "centre") {
  grid::grid.text(s, x = mm(x), y = mm(y), just = just,
                  gp = grid::gpar(fontfamily = FONT, fontsize = pt, col = col,
                                  fontface = if (bold) "bold" else "plain"))
}
rule <- function(x0, y0, x1, y1, col = COL$line, lw = 0.6) {
  grid::grid.lines(x = mm(c(x0, x1)), y = mm(c(y0, y1)),
                   gp = grid::gpar(col = col, lwd = lw, lineend = "round"))
}
box_g <- function(x, y, w, h, fill) {
  grid::grid.rect(x = mm(x), y = mm(y), width = mm(w), height = mm(h),
                  gp = grid::gpar(fill = fill, col = NA))
}

draw_panel <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    xscale = c(0, WIDTH_MM), yscale = c(0, HEIGHT_MM)))
  grid::grid.rect(gp = grid::gpar(fill = COL$white, col = NA))

  # Common denominator
  txt(paste0(fmt(N), " genes analyzed"), 44, 61.0, PT_CONTEXT, COL$muted)

  # Parallel DE sets
  txt("2D-associated alteration", 24.2, 55.8, PT_HEADING, COL$ink, TRUE)
  txt("3D-associated response", 63.8, 55.8, PT_HEADING, COL$ink, TRUE)

  txt(fmt(n2), 24.2, 50.0, PT_NUMBER, COL$twoD, TRUE)
  txt(fmt(n3), 63.8, 50.0, PT_NUMBER, COL$threeD, TRUE)

  txt("Intact DP vs P3 2D", 24.2, 45.6, PT_SMALL, COL$muted)
  txt("P3 3D vs P3 2D", 63.8, 45.6, PT_SMALL, COL$muted)

  # Proportion tracks; denominator = all analyzed genes
  bw <- 29.0
  ty <- 42.1
  bh <- 0.94

  for (cx in c(24.2, 63.8)) {
    box_g(cx, ty, bw, bh, COL$track)
  }

  w2 <- bw * n2 / N
  w3 <- bw * n3 / N

  box_g(24.2 - bw/2 + w2/2, ty, w2, bh, COL$twoD)
  box_g(63.8 - bw/2 + w3/2, ty, w3, bh, COL$threeD)

  txt(pct(n2, N), 24.2, 38.7, PT_SMALL, COL$muted)
  txt(pct(n3, N), 63.8, 38.7, PT_SMALL, COL$muted)

  # Convergence = intersection of the two DE sets
  rule(24.2, 35.6, 44.0, 32.5, lw = 0.7)
  rule(63.8, 35.6, 44.0, 32.5, lw = 0.7)

  txt(fmt(nboth), 44.0, 28.8, PT_OVERLAP, COL$ink, TRUE)
  txt("Significant in both contrasts", 44.0, 24.9, PT_SMALL, COL$muted)

  rule(44.0, 22.9, 44.0, 21.7, lw = 0.7)

  # Strict-rescue subset
  txt(paste0(fmt(ns), " strict rescue genes"), 44.0, 18.9,
      PT_STRICT, COL$ink, TRUE)

  txt(paste0(pct(ns, nboth), " of dual-significant genes"),
      44.0, 15.8, PT_SMALL, COL$muted)

  # Explicit denominator for the three-class bar
  txt(paste0("Classification of ", fmt(nboth), " dual-significant genes"),
      44.0, 12.9, 5.1, COL$muted)

  # Partition of all 1,969 overlap genes
  sx <- 10.5
  sw <- 67.0
  sy <- 10.0
  sh <- 1.65

  wu <- sw * nu / nboth
  wd <- sw * nd / nboth
  wo <- sw * nother / nboth

  box_g(sx + wu/2, sy, wu, sh, COL$up)
  box_g(sx + wu + wd/2, sy, wd, sh, COL$down)
  box_g(sx + wu + wd + wo/2, sy, wo, sh, COL$other)

  # Class labels
  box_g(13.0, 5.5, 1.4, 1.4, COL$up)
  txt(paste0(fmt(nu), " upward"), 15.2, 5.5, PT_LEGEND, COL$ink,
      just = "left")

  box_g(36.0, 5.5, 1.4, 1.4, COL$down)
  txt(paste0(fmt(nd), " downward"), 38.2, 5.5, PT_LEGEND, COL$ink,
      just = "left")

  box_g(64.3, 5.5, 1.4, 1.4, COL$other)
  txt(paste0(fmt(nother), " other"), 66.5, 5.5, PT_LEGEND, COL$ink,
      just = "left")

  grid::popViewport()
}

# 3. Complete export: atomic final folder, no overwrite, validated file headers.
parent_dir <- dirname(OUT_DIR)
dir.create(parent_dir, showWarnings = FALSE, recursive = TRUE)
staging <- tempfile(pattern = ".PanelC_staging_", tmpdir = parent_dir)
if (!dir.create(staging, showWarnings = FALSE)) stop("Cannot create staging folder")

pdf_path <- file.path(staging, paste0(BASENAME, ".pdf"))
png_path <- file.path(staging, paste0(BASENAME, ".png"))
tif_path <- file.path(staging, paste0(BASENAME, ".tif"))
width_in  <- WIDTH_MM/25.4
height_in <- HEIGHT_MM/25.4

export_one <- function(path, open_device) {
  device_open <- FALSE
  tryCatch({
    open_device()
    device_open <- TRUE
    draw_panel()
    grDevices::dev.off()
    device_open <- FALSE
  }, finally = {
    if (device_open) try(grDevices::dev.off(), silent = TRUE)
  })
  if (!file.exists(path) || is.na(file.info(path)$size) ||
      file.info(path)$size < 900) stop("Missing or invalid output: ", path)
}
header <- function(path, bytes) {
  f <- file(path, "rb")
  on.exit(close(f))
  readBin(f, what = "raw", n = bytes)
}

tryCatch({
  # PDF first: native R vector graphics, avoids Cairo/X11 issues.
  export_one(pdf_path, function() grDevices::pdf(
    pdf_path, width = width_in, height = height_in,
    family = FONT, useDingbats = FALSE, onefile = FALSE,
    paper = "special", bg = "white"))
  export_one(png_path, function() ragg::agg_png(
    png_path, width = width_in, height = height_in,
    units = "in", res = DPI, background = "white"))
  export_one(tif_path, function() ragg::agg_tiff(
    tif_path, width = width_in, height = height_in,
    units = "in", res = DPI, compression = "lzw", background = "white"))
  
  if (!identical(header(pdf_path, 5), charToRaw("%PDF-")))
    stop("PDF header invalid")
  if (!identical(as.integer(header(png_path, 8)),
                 c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
    stop("PNG header invalid")
  if (!identical(rawToChar(header(tif_path, 2)), "II") &&
      !identical(rawToChar(header(tif_path, 2)), "MM"))
    stop("TIFF header invalid")
  
  details <- data.frame(
    metric = c("Total genes", "2D-altered", "3D-responsive",
               "Significant in both", "Strict rescue", "Restored upward",
               "Restored downward", "Other dual-significant"),
    genes = as.integer(actual),
    denominator = c(NA, N, N, N, nboth, nboth, nboth, nboth),
    stringsAsFactors = FALSE)
  utils::write.csv(details, file.path(staging, "plot_data.csv"),
                   row.names = FALSE)
  writeLines(c(
    "HairRegAI | Figure 1 Panel C | final-print-size standalone figure",
    paste0("Input: ", normalizePath(INPUT)),
    paste0("Input MD5: ", unname(tools::md5sum(INPUT))),
    "The two DEG sets are PARALLEL, not successive nested screens.",
    "Intersection means Altered_in_2D AND Changed_in_3D.",
    "Strict rescue uses source Strict_rescue flag; rule consistency checked.",
    "Restored upward/downward directions come from source Rescue_type.",
    "510 + 1061 + 398 = 1969. The strip partitions the shared-DE set.",
    "Tall Nature layout with typography normalized below the 8-pt panel-letter scale; no data changed.",
    "Full strict definition is documented in the analysis and figure legend:",
    "FDR < 0.05 and |log2FC| > 0.5 in BOTH contrasts; matched",
    "direction toward intact and closer-to-intact expression.",
    "This is a data-derived screening figure, not independent validation.",
    paste0("Canvas: ", WIDTH_MM, " x ", HEIGHT_MM, " mm"),
    "PDF is vector; PNG/TIFF are 600 dpi; TIFF is LZW compressed.",
    "No panel letter; no full Figure 1 assembly."
  ), con = file.path(staging, "README_provenance.txt"))
  
  if (dir.exists(OUT_DIR) || file.exists(OUT_DIR))
    stop("Final output directory appeared while drawing; cannot overwrite")
  if (!file.rename(staging, OUT_DIR))
    stop("Cannot move completed staging folder to final output folder")
  cat("SUCCESS: Figure 1 Panel C, final size 88 x 64 mm\n")
  cat("Output: ", normalizePath(OUT_DIR), "\n", sep = "")
  cat("Verified: PDF + PNG + LZW TIFF; plots and provenance recorded.\n")
}, error = function(e) {
  message("PANEL C FAILED: ", conditionMessage(e))
  if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE)
  message("Previous Panel C folders were NOT overwritten.")
  stop(e)
})
