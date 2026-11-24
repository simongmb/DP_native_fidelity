# Figure 2 | Panel B | Direction-aligned native fidelity definition | revised v2
# Suggested filename: F2B_final_revised.R
#
# Purpose
# -------
# Define how native DP fidelity is calculated from an independent normal
# hair-follicle anatomy reference.
#
# Concept:
# 1) For a Restored_upward gene, native support means DP > non-DP.
# 2) For a Restored_downward gene, native support means DP < non-DP.
# 3) Multiply the donor-aware DP-vs-non-DP contrast by the frozen restoration
#    direction sign so that positive values always indicate native-consistent
#    anatomy.
#
# This panel is conceptual/methodological. It contains no gene-specific results.
#
# Output:
# <project root>/figures/FIGURE 2/Panel_B_v1/

WIDTH_MM  <- 88
HEIGHT_MM <- 56
DPI       <- 600
FONT      <- "Helvetica"
BASE      <- "Fig2B_DirectionAlignedFidelity_Nature_final"

# -----------------------------
# Palette coordinated with Figure 1 / Figure 2A
# -----------------------------
INK        <- "#263A45"
MUTED      <- "#728894"
LIGHT_TEXT <- "#8EA1AC"
GUIDE      <- "#D8E2E7"
WHITE      <- "#FFFFFF"

UP_COLOR   <- "#B78361"
UP_LIGHT   <- "#F3E7DF"

DOWN_COLOR <- "#078C82"
DOWN_LIGHT <- "#DDEFEA"

DP_COLOR   <- "#087D79"
DP_LIGHT   <- "#DCEFED"

NONDP_FILL <- "#EEF3F5"
NONDP_EDGE <- "#B8C7CE"

# Shared geometry for the four Native anatomy boxes.
# Keeping identical dimensions and a common baseline prevents visual imbalance
# and leaves a clear gap below the "Native anatomy" labels.
ANATOMY_BOX_W <- 10.8
ANATOMY_BOX_H <- 5.0
ANATOMY_BOX_Y <- 31.4

POS_COLOR  <- "#087D79"
NEG_COLOR  <- "#B78361"

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

  hit <- candidates[file.exists(file.path(candidates, "HairRegeneration_AI.Rproj"))]
  if (length(hit)) return(hit[1])

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
OUT  <- file.path(ROOT, "figures", "FIGURE 2", "Panel_B_Nature_final")

# Safe reruns: an existing output folder is allowed. Only the files generated
# by this script are refreshed; unrelated files in the folder are preserved.
if (file.exists(OUT) && !dir.exists(OUT)) {
  stop(
    "The output path exists but is not a directory:\n", OUT,
    "\nPlease rename that file or change OUT in this script."
  )
}

# -----------------------------
# Drawing
# -----------------------------
draw_panel <- function() {
  grid::grid.newpage()

  grid::pushViewport(grid::viewport(
    width  = grid::unit(WIDTH_MM,  "mm"),
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

  arrow_mm <- function(x1, y1, x2, y2,
                       col = INK,
                       lwd = 0.8,
                       length_mm = 1.5) {
    grid::grid.lines(
      x = grid::unit(c(x1, x2), "mm"),
      y = grid::unit(c(y1, y2), "mm"),
      arrow = grid::arrow(
        length = grid::unit(length_mm, "mm"),
        type = "closed"
      ),
      gp = grid::gpar(col = col, lwd = lwd, fill = col)
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
                        border = INK,
                        lwd = 0.6) {
    grid::grid.circle(
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      r = grid::unit(r, "mm"),
      gp = grid::gpar(fill = fill, col = border, lwd = lwd)
    )
  }

  # ============================================================
  # Header
  # ============================================================
  text_mm(
    "Direction-aligned native fidelity",
    x = 3.0, y = 53.0, pt = 5.6,
    col = MUTED, just = "left", bold = TRUE
  )

  text_mm(
    "The same native-DP direction is made positive for every restored gene",
    x = 3.0, y = 49.4, pt = 5.8,
    col = INK, just = "left", bold = TRUE
  )

  text_mm(
    "DP is compared with matched non-DP follicle compartments within donor",
    x = 3.0, y = 46.5, pt = 5.1,
    col = MUTED, just = "left"
  )

  line_mm(3.0, 44.7, 85.0, 44.7, col = GUIDE, lwd = 0.65)

  # ============================================================
  # Two branches: restored upward / restored downward
  # ============================================================

  # Branch centers
  xL <- 23.0
  xR <- 65.0

  # ----- Restored upward branch
  roundrect_mm(
    xL, 40.2, 31.0, 5.4,
    fill = UP_LIGHT, border = UP_COLOR, lwd = 0.8, r = 1.3
  )
  text_mm(
    "Restored upward",
    xL, 40.2, pt = 5.6,
    col = UP_COLOR, bold = TRUE
  )

  # small expression schematic: DP higher than non-DP
  text_mm("Native anatomy", xL, 35.8, pt = 4.8,
          col = MUTED, bold = TRUE)

  # Native anatomy boxes use identical dimensions and vertical alignment.
  # non-DP (lower expression in the upward-restored branch)
  roundrect_mm(
    xL - 8.5, ANATOMY_BOX_Y, ANATOMY_BOX_W, ANATOMY_BOX_H,
    fill = NONDP_FILL, border = NONDP_EDGE, lwd = 0.55, r = 1.0
  )
  text_mm("non-DP", xL - 8.5, ANATOMY_BOX_Y, pt = 4.8, col = INK)

  # DP (higher expression in the upward-restored branch)
  roundrect_mm(
    xL + 8.5, ANATOMY_BOX_Y, ANATOMY_BOX_W, ANATOMY_BOX_H,
    fill = DP_LIGHT, border = DP_COLOR, lwd = 0.75, r = 1.0
  )
  text_mm("DP", xL + 8.5, ANATOMY_BOX_Y, pt = 5.0,
          col = DP_COLOR, bold = TRUE)

  text_mm("DP > non-DP",
          xL, 27.4, pt = 5.5,
          col = INK, bold = TRUE)

  arrow_mm(xL, 25.4, xL, 22.0, col = UP_COLOR, lwd = 0.8)

  roundrect_mm(
    xL, 19.1, 27.0, 5.2,
    fill = WHITE, border = GUIDE, lwd = 0.6, r = 1.1
  )
  text_mm(
    "direction sign = +1",
    xL, 19.1, pt = 4.9,
    col = UP_COLOR, bold = TRUE
  )

  # ----- Restored downward branch
  roundrect_mm(
    xR, 40.2, 31.0, 5.4,
    fill = DOWN_LIGHT, border = DOWN_COLOR, lwd = 0.8, r = 1.3
  )
  text_mm(
    "Restored downward",
    xR, 40.2, pt = 5.6,
    col = DOWN_COLOR, bold = TRUE
  )

  text_mm("Native anatomy", xR, 35.8, pt = 4.8,
          col = MUTED, bold = TRUE)

  # DP (lower expression in the downward-restored branch)
  roundrect_mm(
    xR - 8.5, ANATOMY_BOX_Y, ANATOMY_BOX_W, ANATOMY_BOX_H,
    fill = DP_LIGHT, border = DP_COLOR, lwd = 0.75, r = 1.0
  )
  text_mm("DP", xR - 8.5, ANATOMY_BOX_Y, pt = 5.0,
          col = DP_COLOR, bold = TRUE)

  # non-DP (higher expression in the downward-restored branch)
  roundrect_mm(
    xR + 8.5, ANATOMY_BOX_Y, ANATOMY_BOX_W, ANATOMY_BOX_H,
    fill = NONDP_FILL, border = NONDP_EDGE, lwd = 0.55, r = 1.0
  )
  text_mm("non-DP", xR + 8.5, ANATOMY_BOX_Y, pt = 4.8, col = INK)

  text_mm("DP < non-DP",
          xR, 27.4, pt = 5.5,
          col = INK, bold = TRUE)

  arrow_mm(xR, 25.4, xR, 22.0, col = DOWN_COLOR, lwd = 0.8)

  roundrect_mm(
    xR, 19.1, 27.0, 5.2,
    fill = WHITE, border = GUIDE, lwd = 0.6, r = 1.1
  )
  text_mm(
    "direction sign = -1",
    xR, 19.1, pt = 4.9,
    col = DOWN_COLOR, bold = TRUE
  )

  # ============================================================
  # Convergence into common native-fidelity score
  # ============================================================

  arrow_mm(xL, 16.5, 38.5, 12.8, col = GUIDE, lwd = 0.8)
  arrow_mm(xR, 16.5, 49.5, 12.8, col = GUIDE, lwd = 0.8)

  # Formula strip
  roundrect_mm(
    44.0, 10.4, 56.0, 6.0,
    fill = "#F8FAFB", border = GUIDE, lwd = 0.65, r = 1.2
  )

  text_mm(
    "Native fidelity = (DP - non-DP) × restoration-direction sign",
    44.0, 10.4, pt = 5.2,
    col = INK, bold = TRUE
  )

  # Interpretation bar
  line_mm(14.0, 5.8, 74.0, 5.8, col = GUIDE, lwd = 0.8)

  # left negative
  line_mm(14.0, 5.8, 44.0, 5.8, col = NEG_COLOR, lwd = 2.1)
  # right positive
  line_mm(44.0, 5.8, 74.0, 5.8, col = POS_COLOR, lwd = 2.1)

  circle_mm(44.0, 5.8, 0.55, fill = WHITE, border = INK, lwd = 0.7)

  text_mm(
    "non-native direction",
    14.0, 3.3, pt = 4.7,
    col = NEG_COLOR, just = "left"
  )

  text_mm(
    "native-consistent direction",
    74.0, 3.3, pt = 4.7,
    col = POS_COLOR, just = "right"
  )

  text_mm(
    "0",
    44.0, 3.3, pt = 4.8,
    col = MUTED
  )
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

  render(function() {
    grDevices::pdf(
      file.path(STAGE, paste0(BASE, ".pdf")),
      width  = WIDTH_MM / 25.4,
      height = HEIGHT_MM / 25.4,
      family = FONT,
      useDingbats = FALSE
    )
  })

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

  writeLines(
    c(
      "Figure 2 Panel B — direction-aligned native fidelity",
      "",
      "Purpose:",
      "Define the donor-aware, direction-aligned native DP fidelity metric.",
      "",
      "Restored_upward:",
      "native support = DP > non-DP; direction sign = +1.",
      "",
      "Restored_downward:",
      "native support = DP < non-DP; direction sign = -1.",
      "",
      "Unified interpretation:",
      "positive native-fidelity values indicate expression aligned with",
      "the expected native DP anatomical direction.",
      "",
      "Formula shown:",
      "Native fidelity = (DP - non-DP) × restoration-direction sign",
      "",
      "This panel is conceptual/methodological only.",
      "Native-anatomy boxes share identical width, height, and vertical alignment.",
      "A clear vertical gap is maintained between the anatomy labels and boxes.",
      "Typography and box geometry were jointly tuned for final-size readability.",
      "No gene-specific result or statistical inference is shown.",
      "",
      paste0("Canvas: ", WIDTH_MM, " x ", HEIGHT_MM, " mm"),
      paste0("Raster resolution: ", DPI, " dpi"),
      tiff_note,
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

  if (dir.exists(OUT)) {
    # Refresh only this script's expected outputs. Do not delete the folder,
    # because it may contain files created by the user or another workflow.
    copied <- file.copy(
      from = file.path(STAGE, expected),
      to = OUT,
      overwrite = TRUE,
      copy.mode = TRUE,
      copy.date = TRUE
    )
    if (!all(copied)) {
      failed <- expected[!copied]
      stop(
        "Could not update the following output file(s): ",
        paste(failed, collapse = ", "),
        "\nCheck whether the files are open or the folder is writable."
      )
    }
    unlink(STAGE, recursive = TRUE, force = TRUE)
  } else {
    if (!file.rename(STAGE, OUT)) {
      stop("Cannot finalize output folder: ", OUT)
    }
  }

  ok <- TRUE

}, error = function(e) {
  message("Panel B FAILED: ", conditionMessage(e))
  message("Incomplete staging folder retained: ", STAGE)
  stop(e)
})

if (ok) {
  cat(
    "SUCCESS — Figure 2 Panel B exported to:\n",
    OUT, "\n"
  )
}
