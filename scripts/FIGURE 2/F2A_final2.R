# Figure 2 | Panel A | Native human hair-follicle anatomy reference | v5
# Suggested filename: Fig2_PanelA_v2.R
#
# Purpose
# -------
# Introduce GSE272333 as an independent normal human hair-follicle anatomy
# reference used to evaluate native DP fidelity.
#
# This is a schematic/descriptive panel. No statistical result is shown.
#
# Visual language:
# - coordinated with Figure 1
# - Helvetica
# - muted blue-grey context
# - DP highlighted in teal
# - physical-mm layout
# - native grid graphics
# - PDF + 600 dpi PNG/TIFF
#
# Output:
# <project root>/figures/FIGURE 2/Panel_A_v2/

WIDTH_MM  <- 88
HEIGHT_MM <- 56
DPI       <- 600
FONT      <- "Helvetica"
BASE      <- "Fig2A_NativeAnatomy_Nature_final2"

# -----------------------------
# Palette coordinated with Fig. 1
# -----------------------------
INK        <- "#263A45"
MUTED      <- "#728894"
LIGHT_TEXT <- "#8EA1AC"
HAIR_EDGE  <- "#8EA7B4"
HAIR_FILL  <- "#F4F7F8"
INNER_FILL <- "#E4EDF0"
NONDP_FILL <- "#F1F5F6"
NONDP_LINE <- "#B7C6CD"
DP_COLOR   <- "#087D79"
DP_LIGHT   <- "#DCEFED"
GUIDE      <- "#D8E2E7"
WHITE      <- "#FFFFFF"

# -----------------------------
# Locate project root
# -----------------------------
find_root <- function() {
  wd <- normalizePath(getwd(), mustWork = TRUE)
  candidates <- unique(c(
    wd,
    dirname(wd),
    dirname(dirname(wd)),
    dirname(dirname(dirname(wd)))
  ))

  # Preferred: original project name
  hit <- candidates[file.exists(file.path(candidates, "HairRegeneration_AI.Rproj"))]
  if (length(hit)) return(hit[1])

  # Fallback: project-like folder
  hit <- candidates[
    dir.exists(file.path(candidates, "figures")) &
    dir.exists(file.path(candidates, "results"))
  ]
  if (length(hit)) return(hit[1])

  stop(
    "Project root not found.\n",
    "Please run from the HairRegeneration_AI project or one of its subfolders."
  )
}

ROOT <- find_root()
OUT  <- file.path(ROOT, "figures", "FIGURE 2", "Panel_A_Nature_final2")

if (file.exists(OUT) || dir.exists(OUT)) {
  stop(
    "Refusing to overwrite existing output folder:\n", OUT,
    "\nRename/remove it before rerunning."
  )
}

# -----------------------------
# Drawing helpers
# -----------------------------
draw_panel <- function() {
  grid::grid.newpage()

  grid::pushViewport(grid::viewport(
    width  = grid::unit(WIDTH_MM, "mm"),
    height = grid::unit(HEIGHT_MM, "mm")
  ))
  on.exit(grid::popViewport(), add = TRUE)

  text_mm <- function(label, x, y, pt,
                      col = INK,
                      just = "centre",
                      bold = FALSE,
                      rot = 0) {
    grid::grid.text(
      label,
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      just = just,
      rot = rot,
      gp = grid::gpar(
        col = col,
        fontsize = pt,
        fontfamily = FONT,
        fontface = if (bold) "bold" else "plain"
      )
    )
  }

  line_mm <- function(x1, y1, x2, y2,
                      col = GUIDE,
                      lwd = 0.6,
                      lty = 1) {
    grid::grid.lines(
      x = grid::unit(c(x1, x2), "mm"),
      y = grid::unit(c(y1, y2), "mm"),
      gp = grid::gpar(col = col, lwd = lwd, lty = lty)
    )
  }

  roundrect_mm <- function(x, y, w, h,
                           fill = WHITE,
                           border = GUIDE,
                           lwd = 0.6,
                           r = 1.2) {
    grid::grid.roundrect(
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      width  = grid::unit(w, "mm"),
      height = grid::unit(h, "mm"),
      r = grid::unit(r, "mm"),
      gp = grid::gpar(fill = fill, col = border, lwd = lwd)
    )
  }

  circle_mm <- function(x, y, r,
                        fill = WHITE,
                        border = HAIR_EDGE,
                        lwd = 0.55) {
    grid::grid.circle(
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      r = grid::unit(r, "mm"),
      gp = grid::gpar(fill = fill, col = border, lwd = lwd)
    )
  }

  polygon_mm <- function(x, y,
                         fill = WHITE,
                         border = HAIR_EDGE,
                         lwd = 0.65) {
    grid::grid.polygon(
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      gp = grid::gpar(fill = fill, col = border, lwd = lwd)
    )
  }

  xspline_mm <- function(x, y,
                         shape = 1,
                         open = FALSE,
                         fill = NA,
                         border = HAIR_EDGE,
                         lwd = 0.7) {
    grid::grid.xspline(
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      shape = shape,
      open = open,
      gp = grid::gpar(fill = fill, col = border, lwd = lwd)
    )
  }

  # ============================================================
  # Header
  # ============================================================
  text_mm(
    "Normal human hair-follicle anatomy",
    x = 3.0, y = 53.0, pt = 5.6,
    col = MUTED, just = "left", bold = TRUE
  )

  text_mm(
    "GSE272333",
    x = 3.0, y = 49.2, pt = 6.6,
    col = INK, just = "left", bold = TRUE
  )

  text_mm(
    "Independent reference for native DP fidelity",
    x = 3.0, y = 46.3, pt = 5.1,
    col = MUTED, just = "left"
  )

  line_mm(3.0, 44.6, 85.0, 44.6, col = GUIDE, lwd = 0.65)

  # ============================================================
  # Left: refined hair-follicle schematic
  # ============================================================
  cx <- 21.5

  # Outer follicle contour: smoother, tapered, organic silhouette.
  # Drawn with an open xspline so the outline feels less geometric.
  xspline_mm(
    x = c(cx - 4.3, cx - 4.6, cx - 4.8, cx - 5.0, cx - 5.8,
          cx - 5.0, cx - 3.7, cx - 2.1, cx,
          cx + 2.1, cx + 3.7, cx + 5.0, cx + 5.8,
          cx + 5.0, cx + 4.8, cx + 4.6, cx + 4.3),
    y = c(39.8, 36.6, 31.2, 26.0, 21.8,
          18.0, 15.7, 14.7, 14.2,
          14.7, 15.7, 18.0, 21.8,
          26.0, 31.2, 36.6, 39.8),
    shape = 0.6,
    open = TRUE,
    fill = NA,
    border = HAIR_EDGE,
    lwd = 0.85
  )

  # Soft inner epithelial column.
  xspline_mm(
    x = c(cx - 2.15, cx - 2.25, cx - 2.25, cx - 2.45,
          cx - 2.9, cx - 1.9, cx,
          cx + 1.9, cx + 2.9, cx + 2.45,
          cx + 2.25, cx + 2.25, cx + 2.15),
    y = c(39.0, 35.0, 30.0, 25.2,
          21.0, 18.0, 17.0,
          18.0, 21.0, 25.2,
          30.0, 35.0, 39.0),
    shape = 0.65,
    open = TRUE,
    fill = NA,
    border = "#B8C9D1",
    lwd = 0.7
  )

  # Very light central shaft/core fill to add depth without illustration-like styling.
  grid::grid.roundrect(
    x = grid::unit(cx, "mm"),
    y = grid::unit(31.8, "mm"),
    width = grid::unit(3.35, "mm"),
    height = grid::unit(14.8, "mm"),
    r = grid::unit(1.6, "mm"),
    gp = grid::gpar(fill = "#EEF3F5", col = NA)
  )

  # Lower matrix/bulb shading: subtle and layered.
  polygon_mm(
    x = c(cx - 3.0, cx - 2.5, cx - 1.45, cx,
          cx + 1.45, cx + 2.5, cx + 3.0, cx + 2.1,
          cx, cx - 2.1),
    y = c(23.0, 19.8, 17.5, 16.6,
          17.5, 19.8, 23.0, 24.8,
          25.5, 24.8),
    fill = "#E7EFF2", border = "#C5D3D9", lwd = 0.45
  )

  # Native DP: rounded teardrop-like papilla rather than a polygon block.
  xspline_mm(
    x = c(cx - 2.4, cx - 2.0, cx - 1.1, cx,
          cx + 1.1, cx + 2.0, cx + 2.4, cx + 1.6,
          cx, cx - 1.6),
    y = c(17.6, 15.6, 14.6, 14.1,
          14.6, 15.6, 17.6, 18.4,
          18.8, 18.4),
    shape = 0.7,
    open = FALSE,
    fill = DP_LIGHT,
    border = DP_COLOR,
    lwd = 1.0
  )

  # DP cells: small restrained circles.
  circle_mm(cx - 1.05, 16.5, 0.37, fill = DP_COLOR, border = DP_COLOR)
  circle_mm(cx,        15.7, 0.37, fill = DP_COLOR, border = DP_COLOR)
  circle_mm(cx + 1.05, 16.6, 0.37, fill = DP_COLOR, border = DP_COLOR)

  # A few subtle epithelial cells above DP to make the bulb feel biological.
  circle_mm(cx - 1.35, 19.4, 0.27, fill = "#C9D8DE", border = "#C9D8DE")
  circle_mm(cx,        19.0, 0.27, fill = "#C9D8DE", border = "#C9D8DE")
  circle_mm(cx + 1.35, 19.4, 0.27, fill = "#C9D8DE", border = "#C9D8DE")

  # DP callout
  line_mm(cx + 2.15, 16.7, 29.0, 16.7, col = DP_COLOR, lwd = 0.8)
  text_mm("DP", 30.3, 16.7, pt = 5.4,
          col = DP_COLOR, just = "left", bold = TRUE)

  text_mm(
    "hair-follicle compartments",
    cx, 11.0, pt = 5.0, col = LIGHT_TEXT
  )

  # ============================================================
  # Right: native target + compact compartment reference
  # ============================================================
  x_left <- 39.0

  text_mm(
    "Native target",
    x_left, 42.0, pt = 5.2,
    col = MUTED, just = "left", bold = TRUE
  )

  # DP target strip
  roundrect_mm(
    61.0, 36.8, 43.0, 4.8,
    fill = DP_LIGHT, border = DP_COLOR,
    lwd = 0.8, r = 1.25
  )
  text_mm(
    "Dermal papilla (DP)",
    61.0, 36.8, pt = 5.5,
    col = DP_COLOR, bold = TRUE
  )

  text_mm(
    "Non-DP anatomical references",
    x_left, 31.8, pt = 5.2,
    col = MUTED, just = "left", bold = TRUE
  )

  # ---- Group 1: dermal mesenchyme
  text_mm(
    "Dermal mesenchyme",
    x_left, 28.7, pt = 4.85,
    col = LIGHT_TEXT, just = "left", bold = TRUE
  )
  text_mm(
    "DC  Dermal cup",
    x_left, 26.2, pt = 5.25,
    col = INK, just = "left"
  )
  text_mm(
    "CTS  Connective tissue sheath",
    x_left + 21.2, 26.2, pt = 5.05,
    col = INK, just = "left"
  )

  line_mm(x_left, 24.1, 84.4, 24.1, col = GUIDE, lwd = 0.45)

  # ---- Group 2: matrix
  text_mm(
    "Hair matrix",
    x_left, 22.0, pt = 4.85,
    col = LIGHT_TEXT, just = "left", bold = TRUE
  )
  text_mm(
    "HM1  Hair matrix 1",
    x_left, 19.5, pt = 5.25,
    col = INK, just = "left"
  )
  text_mm(
    "HM2  Hair matrix 2",
    x_left + 21.2, 19.5, pt = 5.05,
    col = INK, just = "left"
  )

  line_mm(x_left, 17.4, 84.4, 17.4, col = GUIDE, lwd = 0.45)

  # ---- Group 3: root sheaths
  text_mm(
    "Root sheaths",
    x_left, 15.3, pt = 4.85,
    col = LIGHT_TEXT, just = "left", bold = TRUE
  )

  text_mm(
    "IRS1  Inner root sheath 1",
    x_left, 12.8, pt = 4.95,
    col = INK, just = "left"
  )
  text_mm(
    "ORS1  Outer root sheath 1",
    x_left + 21.2, 12.8, pt = 4.95,
    col = INK, just = "left"
  )
  text_mm(
    "ORS2  Outer root sheath 2",
    x_left, 10.1, pt = 4.95,
    col = INK, just = "left"
  )

  # ============================================================
  # Footer: dataset scale
  # ============================================================
  line_mm(3.0, 7.4, 85.0, 7.4, col = GUIDE, lwd = 0.65)

  text_mm("10 donors",       15.0, 4.45, pt = 5.3, col = INK, bold = TRUE)
  text_mm("110 profiles",    44.0, 4.45, pt = 5.3, col = INK, bold = TRUE)
  text_mm("8 compartments",  73.0, 4.45, pt = 5.3, col = INK, bold = TRUE)

  line_mm(29.5, 3.2, 29.5, 5.7, col = GUIDE, lwd = 0.55)
  line_mm(58.5, 3.2, 58.5, 5.7, col = GUIDE, lwd = 0.55)
}

# -----------------------------
# Export
# -----------------------------
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

STAGE <- paste0(OUT, "_INCOMPLETE_", Sys.getpid())
if (!dir.create(STAGE)) {
  stop("Cannot create staging directory: ", STAGE)
}

ok <- FALSE

tryCatch({

  render <- function(open_device) {
    open_device()
    tryCatch(draw_panel(), finally = grDevices::dev.off())
  }

  # PDF
  render(function() {
    grDevices::pdf(
      file.path(STAGE, paste0(BASE, ".pdf")),
      width  = WIDTH_MM / 25.4,
      height = HEIGHT_MM / 25.4,
      family = FONT,
      useDingbats = FALSE
    )
  })

  # PNG
  render(function() {
    grDevices::png(
      file.path(STAGE, paste0(BASE, ".png")),
      width  = WIDTH_MM / 25.4,
      height = HEIGHT_MM / 25.4,
      units = "in",
      res = DPI,
      bg = "white"
    )
  })

  # TIFF
  if (requireNamespace("ragg", quietly = TRUE)) {
    render(function() {
      ragg::agg_tiff(
        file.path(STAGE, paste0(BASE, ".tif")),
        width  = WIDTH_MM / 25.4,
        height = HEIGHT_MM / 25.4,
        units = "in",
        res = DPI,
        compression = "lzw",
        background = "white"
      )
    })
    tiff_note <- "TIFF: ragg, 600 dpi, LZW compression"
  } else {
    render(function() {
      grDevices::tiff(
        file.path(STAGE, paste0(BASE, ".tif")),
        width  = WIDTH_MM / 25.4,
        height = HEIGHT_MM / 25.4,
        units = "in",
        res = DPI,
        bg = "white"
      )
    })
    tiff_note <- "TIFF: base R, 600 dpi"
  }

  # Provenance
  writeLines(
    c(
      "Figure 2 Panel A v5 — normal human hair-follicle anatomy reference",
      "",
      "Purpose:",
      "Introduce GSE272333 as an independent normal-tissue reference used",
      "to evaluate native dermal papilla fidelity.",
      "",
      "Displayed dataset context:",
      "GSE272333",
      "10 donors | 110 profiles | 8 compartments",
      "DP, DC, CTS, HM1, HM2, IRS1, ORS1, ORS2",
      "",
      "Panel is schematic/descriptive only.",
      "No statistical test is encoded in the artwork.",
      "",
      "Design refinements:",
      "- removed seven oversized comparator boxes",
      "- grouped non-DP references by anatomical class",
      "- simplified DP target emphasis",
      "- retained Figure 1 colour language and physical-mm layout",
      "- typography normalized below the 8-pt panel-letter scale",
      "",
      paste0("Canvas: ", WIDTH_MM, " x ", HEIGHT_MM, " mm"),
      paste0("Raster resolution: ", DPI, " dpi"),
      tiff_note,
      "",
      paste0("Output directory: ", OUT)
    ),
    file.path(STAGE, "README_provenance.txt")
  )

  expected <- c(
    paste0(BASE, ".pdf"),
    paste0(BASE, ".png"),
    paste0(BASE, ".tif"),
    "README_provenance.txt"
  )

  if (!all(file.exists(file.path(STAGE, expected)))) {
    stop("Missing expected output files")
  }

  if (!file.rename(STAGE, OUT)) {
    stop("Cannot finalize output folder")
  }

  ok <- TRUE

}, error = function(e) {
  message("Panel A v2 FAILED: ", conditionMessage(e))
  message("Incomplete staging folder retained: ", STAGE)
  stop(e)
})

if (ok) {
  cat(
    "SUCCESS — Figure 2 Panel A v2 exported to:\n",
    OUT, "\n"
  )
}
