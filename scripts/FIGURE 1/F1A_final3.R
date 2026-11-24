# ============================================================================
# HairRegAI | Figure 1, Panel A | Nature-style final-size standalone graphic
# Source: GSE44765, 3 donors x 6 culture conditions (18 expression profiles)
# Save to: scripts/FIGURE 1/Fig1_PanelA.R
# Run in RStudio: source("scripts/FIGURE 1/Fig1_PanelA.R")
# Output to: figures/FIGURE 1/Panel_A/
# PDF (vector) + PNG (600 dpi) + TIFF (600 dpi, LZW)
# No panel letter; no automatic figure assembly; no source data changed.
# ============================================================================

# -------------------------- EDITABLE PARAMETERS -----------------------------
WIDTH_MM  <- 88.0   # Designed for ~one half of a 183-mm-wide two-column figure
HEIGHT_MM <- 44.0
DPI       <- 600L
FONT      <- "Helvetica"

OUT_DIR   <- file.path("figures", "FIGURE 1", "Panel_A_Nature_final3")
META_FILE <- file.path("results", "GSE44765_sample_metadata.csv")
BASENAME  <- "Fig1A_ExperimentalDesign_Nature_final3"

# Font sizes below are ACTUAL points at the exported 88 x 57 mm panel size.
SIZE_HEADING   <- 6.2
SIZE_GROUP     <- 6.4
SIZE_COUNT     <- 5.2
SIZE_CONTRAST  <- 5.8
SIZE_FOOTER    <- 5.0

INK        <- "#253945"
NAVY       <- "#46677A"
MUTED      <- "#748A95"
TEAL       <- "#087F79"
TEAL_LIGHT <- "#E5F3F0"
RULE       <- "#CBD9DE"
PALE       <- "#F0F4F5"
WHITE      <- "#FFFFFF"
# ---------------------------------------------------------------------------

# 1. Check project folder, dependencies and experimental design.
if (!file.exists(META_FILE)) {
  stop("Metadata not found: ", META_FILE,
       "\nOpen HairRegeneration_AI.Rproj in RStudio, then run this script.")
}
if (!requireNamespace("ragg", quietly = TRUE)) {
  stop("Package 'ragg' is required for reliable 600-dpi PNG and LZW TIFF.\n",
       "Install once with install.packages('ragg').")
}

meta <- utils::read.csv(META_FILE, stringsAsFactors = FALSE,
                        check.names = FALSE)
required_cols <- c("sample", "donor", "group")
if (!all(required_cols %in% names(meta))) {
  stop("Metadata needs columns: ", paste(required_cols, collapse = ", "))
}
expected_groups <- c("Intact", "P0", "P1", "P3_2D", "P5_2D", "P3_3D")
expected_donors <- c("D5", "D6", "D7")
observed <- table(factor(meta$group, levels = expected_groups),
                  factor(meta$donor, levels = expected_donors))
if (nrow(meta) != 18L || anyNA(meta[, required_cols]) ||
    anyDuplicated(meta$sample) ||
    !setequal(unique(meta$group), expected_groups) ||
    !setequal(unique(meta$donor), expected_donors) ||
    !all(observed == 1L)) {
  stop("Unexpected metadata. Expected 18 distinct samples: ",
       "each of D5, D6 and D7 represented once in all six conditions.")
}

# Protection: old figures are never overwritten.
if (dir.exists(OUT_DIR) || file.exists(OUT_DIR)) {
  stop("Existing output protected: ", OUT_DIR,
       "\nRename the existing Panel_A directory (e.g. Panel_A_previous) ",
       "before generating this new version.")
}

# 2. Physical-coordinate drawing functions. No ggplot2 or Cairo/X11 required.
# The viewport uses millimetres on BOTH axes, so aspect ratios remain faithful.
mm <- function(x) grid::unit(x, "native")

g_line <- function(x0, y0, x1, y1, color = RULE, width = 0.40) {
  grid::grid.lines(x = mm(c(x0, x1)), y = mm(c(y0, y1)),
                   gp = grid::gpar(col = color, lwd = width,
                                   lineend = "round"))
}

g_text <- function(label, x, y, size = SIZE_GROUP, color = INK,
                   bold = FALSE, just = "centre") {
  grid::grid.text(label, x = mm(x), y = mm(y), just = just,
                  gp = grid::gpar(fontfamily = FONT, fontsize = size,
                                  col = color,
                                  fontface = if (bold) "bold" else "plain"))
}

g_ellipse <- function(cx, cy, width_mm, height_mm, fill, edge = NA,
                      lwd = 0.35, rotate_deg = 0) {
  theta <- seq(0, 2 * pi, length.out = 65)
  ang <- rotate_deg * pi / 180
  ex <- (width_mm / 2) * cos(theta)
  ey <- (height_mm / 2) * sin(theta)
  grid::grid.polygon(
    x = mm(cx + ex * cos(ang) - ey * sin(ang)),
    y = mm(cy + ex * sin(ang) + ey * cos(ang)),
    gp = grid::gpar(fill = fill, col = edge, lwd = lwd))
}

# These motifs are schematics, NOT real microscopy or cell-lineage evidence.
draw_intact <- function(cx, cy) {
  shape <- rbind(
    c(-3.2,  1.0, -25), c( 0.1,  2.2, 18),
    c( 3.1,  1.0, -15), c(-1.6, -2.2, 18),
    c( 2.2, -2.1, -28))
  for (i in seq_len(nrow(shape))) {
    xx <- cx + shape[i, 1]; yy <- cy + shape[i, 2]
    g_ellipse(xx, yy, 5.5, 2.05,
              if (i %% 2) "#DAE7EB" else "#E8EFF1",
              NAVY, 0.36, shape[i, 3])
    g_ellipse(xx, yy, 0.90, 0.66, NAVY)
  }
}

draw_dish <- function(cx, cy, accent = FALSE) {
  col <- if (accent) TEAL else NAVY
  fill <- if (accent) TEAL_LIGHT else PALE
  # Shallow monolayer culture dish, rendered at final printed size.
  g_line(cx-5.1, cy-1.5, cx+5.1, cy-1.5, col, 0.46)
  g_line(cx-4.2, cy-3.1, cx+4.2, cy-3.1, RULE, 0.32)
  g_line(cx-5.1, cy-1.5, cx-4.2, cy-3.1, RULE, 0.32)
  g_line(cx+5.1, cy-1.5, cx+4.2, cy-3.1, RULE, 0.32)
  for (off in c(-3.25, -1.1, 1.1, 3.25)) {
    g_ellipse(cx+off, cy-0.55, 1.65, 1.15, fill, col, 0.32)
    g_ellipse(cx+off, cy-0.55, 0.40, 0.37, col)
  }
}

draw_spheroid <- function(cx, cy) {
  offsets <- rbind(c(-1.8, 1.6), c(0, 2.5), c(2, 1.4),
                   c(-2.5, 0), c(-0.4, 0.2), c(2.5, -0.4),
                   c(-1.8,-1.8), c(0.7,-2.0), c(2,-1.6))
  for (i in seq_len(nrow(offsets))) {
    xx <- cx + offsets[i, 1]; yy <- cy + offsets[i, 2]
    g_ellipse(xx, yy, 3.45, 3.45,
              if (i %% 3 == 0) "#C8E8E2" else TEAL_LIGHT,
              TEAL, 0.32)
    g_ellipse(xx+0.23, yy-0.12, 0.62, 0.62, TEAL)
  }
  g_ellipse(cx, cy, 10.2, 10.1, fill = NA, edge = TEAL, lwd = 0.50)
}

# 3. Plot composition at the REAL target size of 88 x 57 mm.
draw_panel <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    xscale = c(0, WIDTH_MM), yscale = c(0, HEIGHT_MM),
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc")))

  grid::grid.rect(gp = grid::gpar(fill = WHITE, col = NA))

  # Nature-style hierarchy:
  # - all text in dark neutral ink
  # - colour reserved for graphical accents
  # - minimal explanatory text inside the panel

  g_text("Fresh DP",       11.0, 40.5, SIZE_HEADING, INK, TRUE)
  g_text("2D monolayer",   45.0, 40.5, SIZE_HEADING, INK, TRUE)
  g_text("3D spheroid",    78.0, 40.5, SIZE_HEADING, INK, TRUE)

  g_line( 4.0, 38.6, 18.0, 38.6, RULE, 0.36)
  g_line(21.0, 38.6, 68.5, 38.6, RULE, 0.36)
  g_line(71.0, 38.6, 85.0, 38.6, RULE, 0.36)

  # Six experimental states
  centers_2d <- c(27, 39, 51, 63)
  draw_intact(11, 31.7)
  for (i in seq_along(centers_2d)) {
    draw_dish(centers_2d[i], 31.8, accent = (i == 3L))
  }
  draw_spheroid(78, 31.6)

  # Condition labels
  g_text("Intact DP", 11, 23.8, SIZE_GROUP, INK, TRUE)
  labels_2d <- c("P0", "P1", "P3 2D", "P5 2D")
  for (i in seq_along(centers_2d)) {
    g_text(labels_2d[i], centers_2d[i], 23.8,
           SIZE_GROUP, INK, TRUE)
  }
  g_text("P3 3D", 78, 23.8, SIZE_GROUP, INK, TRUE)

  # Replicate count
  for (xx in c(11, centers_2d, 78)) {
    g_text("n = 3", xx, 20.5, SIZE_COUNT, MUTED)
  }

  # Passage-matched comparison; colour retained only as a graphical accent
  g_line(51, 16.6, 51, 15.3, TEAL, 0.50)
  g_line(78, 16.6, 78, 15.3, TEAL, 0.50)
  g_line(51, 15.3, 78, 15.3, TEAL, 0.50)
  g_text("P3 2D vs P3 3D", 64.5, 12.1, SIZE_CONTRAST, INK)

  # Concise donor annotation retained for standalone interpretability.
  g_text(
    "Donors D5, D6 and D7; n = 3 per condition",
    6.5, 6.7,
    size = 5.0,
    color = MUTED,
    bold = FALSE,
    just = "left"
  )

  grid::popViewport()
}

# 4. Export into a temporary staging directory; commit only if all 3 succeed.
parent <- dirname(OUT_DIR)
if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
if (!dir.exists(parent)) stop("Could not create parent output directory: ", parent)
staging <- tempfile(pattern = ".Panel_A_staging_", tmpdir = parent)
if (!dir.create(staging, showWarnings = FALSE)) stop("Cannot create staging directory")

w_in <- WIDTH_MM / 25.4
h_in <- HEIGHT_MM / 25.4
pdf_file <- file.path(staging, paste0(BASENAME, ".pdf"))
png_file <- file.path(staging, paste0(BASENAME, ".png"))
tif_file <- file.path(staging, paste0(BASENAME, ".tif"))

render <- function(open_device) {
  open_device()
  tryCatch(draw_panel(), finally = grDevices::dev.off())
}

tryCatch({
  # Native PDF remains editable vector artwork. No Cairo or X11 required.
  render(function() grDevices::pdf(pdf_file, width = w_in,
                                   height = h_in, family = "Helvetica", useDingbats = FALSE,
                                   onefile = FALSE, paper = "special", bg = WHITE))
  
  # TIFF is truly LZW-compressed when written with ragg::agg_tiff.
  render(function() ragg::agg_png(png_file, width = WIDTH_MM,
                                  height = HEIGHT_MM, units = "mm", res = DPI, background = WHITE))
  render(function() ragg::agg_tiff(tif_file, width = WIDTH_MM,
                                   height = HEIGHT_MM, units = "mm", res = DPI,
                                   compression = "lzw", background = WHITE))
  
  image_files <- c(pdf_file, png_file, tif_file)
  if (!all(file.exists(image_files)) || any(file.info(image_files)$size < 900)) {
    stop("Some image files are missing or too small")
  }
  
  file_header <- function(path, n) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    readBin(con, what = "raw", n = n)
  }
  if (rawToChar(file_header(pdf_file, 5)) != "%PDF-") stop("Invalid PDF output")
  if (!identical(as.integer(file_header(png_file, 8)),
                 c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))) {
    stop("Invalid PNG output")
  }
  tif_header <- rawToChar(file_header(tif_file, 2))
  if (!(tif_header %in% c("II", "MM"))) stop("Invalid TIFF output")
  
  utils::write.csv(meta, file.path(staging, "source_metadata.csv"),
                   row.names = FALSE)
  writeLines(c(
    "Figure 1 Panel A: six experimental culture conditions (schematic)",
    paste0("Input: ", META_FILE),
    "Sample groups: Intact, P0, P1, P3_2D, P5_2D, P3_3D",
    "Donors: D5, D6, D7; each condition has n=3; total 18 profiles",
    "The drawings are schematic, not cell microscopy or lineage evidence.",
    "The bracket identifies P3 2D vs P3 3D for a passage-matched contrast.",
    paste0("Final-size canvas: ", WIDTH_MM, " x ", HEIGHT_MM, " mm"),
    paste0("Text sizes (pt): heading ", SIZE_HEADING,
           "; labels ", SIZE_GROUP, "; counts ", SIZE_COUNT),
    "PDF: vector; PNG/TIFF: 600 dpi; TIFF: LZW; no panel letter.",
    "No source input or existing output folder has been overwritten."
  ), file.path(staging, "README_provenance.txt"))
  
  if (file.exists(OUT_DIR) || dir.exists(OUT_DIR)) {
    stop("Output directory appeared during render; refusing to overwrite")
  }
  if (!file.rename(staging, OUT_DIR)) {
    stop("Could not move validated files from staging directory to: ", OUT_DIR)
  }
  cat("SUCCESS: standalone Panel A exported at final-size 88 x 44 mm.\n")
  cat("Output folder: ", normalizePath(OUT_DIR), "\n", sep = "")
  cat("Verified: vector PDF + 600-dpi PNG + LZW TIFF.\n")
}, error = function(e) {
  if (dir.exists(staging)) {
    message("Temporary troubleshooting files: ", staging)
  }
  stop(conditionMessage(e), call. = FALSE)
})
