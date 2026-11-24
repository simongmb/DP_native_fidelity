# ============================================================
# Figure 1 - Panel F (v3): text-only micro-adjustment
# Based on the user's Fig1_PanelF(2).R (Nature-style Panel F v2).
# Original geometry, bars, colors, and frozen statistics preserved.
# Only secondary labels, numerical labels, and explanatory fonts tuned.
#
# Place in: HairRegeneration_AI/scripts/FIGURE 1/
# Run from an open HairRegeneration_AI RStudio project:
# source("scripts/FIGURE 1/Fig1_Panel_F_v3_text_refine.R")
#
# No external input, no file picker, no raw-gene recomputation.
# Produces exactly 5 files. Never overwrites an existing output.
# ============================================================

# --------- 1. FROZEN values from the provided Panel F script ---
N <- 21355L
CLOSER <- 48.8
NOT_CLOSER <- 51.2
MEDIAN_2D <- 0.3602
MEDIAN_3D <- 0.3623
MAX_SCALE <- 0.42
stopifnot(N == 21355L,
          isTRUE(all.equal(CLOSER + NOT_CLOSER, 100)),
          MEDIAN_2D >= 0, MEDIAN_3D >= 0,
          MEDIAN_2D <= MAX_SCALE, MEDIAN_3D <= MAX_SCALE)

# --------- 2. Resolve project, preserve older versions ----------
find_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- unique(c(wd, dirname(wd), dirname(dirname(wd))))
  hits <- candidates[file.exists(file.path(candidates, "HairRegeneration_AI.Rproj"))]
  if (!length(hits)) {
    stop("请先在 RStudio 中打开 HairRegeneration_AI.Rproj，再运行本脚本。")
  }
  hits[1]
}
ROOT <- find_root()
OUT <- file.path(ROOT, "figures", "FIGURE 1", "Panel_F_Nature_final")
BASE <- "Fig1F_GlobalProximity_Nature_final"
W_MM <- 88
H_MM <- 56
DPI <- 600L
FONT <- "Helvetica"

if (file.exists(OUT) || dir.exists(OUT)) {
  stop("为保护旧结果，输出目录已存在，未覆盖：\n", OUT)
}

# --------- 3. Colors and text sizes (pt) -------------------------
COL <- list(
  ink   = "#253946",
  muted = "#738C9B",
  teal  = "#078E86",
  gold  = "#D2A062",
  red   = "#BB6058",
  green = "#2C9B92",
  light = "#E9EFF2",
  rule  = "#D9E2E7"
)

# Relative to v2 (micro-adjustments, not a new layout):
# - percent values remain 24 pt
# - top/body headings unchanged
# - bottom numerical labels 10 -> 12 pt
# - scale note 7.3 -> 9.5 pt
# - bottom conclusion 8.2 -> 10 pt
# - gene count 8.5 -> 10 pt
# - contextual labels 9 -> 9.8 pt; row labels 10 -> 11 pt
PT_PERCENT <- 7.8
PT_MAIN_HEADING <- 6.2
PT_SECTION_HEADING <- 6.0
PT_CONTEXT <- 5.1
PT_GENE_COUNT <- 5.0
PT_ROW_LABEL <- 5.4
PT_VALUE <- 5.4
PT_SCALE_NOTE <- 4.9
PT_CONCLUSION <- 4.9

# --------- 4. Draw one panel using the ORIGINAL v2 geometry -----
draw_panel <- function() {
  grid::grid.newpage()

  grid::pushViewport(
    grid::viewport(
      width = grid::unit(W_MM, "mm"),
      height = grid::unit(H_MM, "mm"),
      x = 0.5,
      y = 0.5
    )
  )

  on.exit(grid::popViewport(), add = TRUE)

  txt <- function(label, x, y, size = 5.2, col = COL$ink,
                  bold = FALSE, just = "left") {
    grid::grid.text(
      label,
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      just = just,
      gp = grid::gpar(
        fontfamily = FONT,
        fontsize = size,
        fontface = if (bold) "bold" else "plain",
        col = col
      )
    )
  }

  rect <- function(x, y, w, h, col) {
    grid::grid.rect(
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      width = grid::unit(w, "mm"),
      height = grid::unit(h, "mm"),
      gp = grid::gpar(fill = col, col = NA)
    )
  }

  rule <- function(x0, x1, y) {
    grid::grid.lines(
      x = grid::unit(c(x0, x1), "mm"),
      y = grid::unit(c(y, y), "mm"),
      gp = grid::gpar(col = COL$rule, lwd = 0.45)
    )
  }

  # -------------------------------------------------------------------------
  # Upper module
  # -------------------------------------------------------------------------

  txt(
    "Fraction of genes closer to intact DP in 3D",
    5.0, 52.2,
    PT_MAIN_HEADING,
    bold = TRUE
  )

  txt(
    "48.8%",
    5.0, 45.2,
    PT_PERCENT,
    COL$teal,
    TRUE
  )

  txt(
    "closer to intact",
    17.9, 45.2,
    PT_CONTEXT,
    COL$teal
  )

  txt(
    "51.2%",
    49.2, 45.2,
    PT_PERCENT,
    COL$gold,
    TRUE
  )

  txt(
    "not closer",
    62.3, 45.2,
    PT_CONTEXT,
    COL$gold
  )

  bar_left <- 5.0
  bar_width <- 78.0
  bar_y <- 39.8
  bar_height <- 4.0
  split <- bar_width * CLOSER / 100

  rect(
    bar_left + split / 2,
    bar_y,
    split,
    bar_height,
    COL$teal
  )

  rect(
    bar_left + split + (bar_width - split) / 2,
    bar_y,
    bar_width - split,
    bar_height,
    COL$gold
  )

  txt(
    "All 21,355 analyzed genes",
    5.0, 35.8,
    PT_GENE_COUNT,
    COL$muted
  )

  rule(
    5.0,
    83.0,
    33.2
  )

  # -------------------------------------------------------------------------
  # Lower module
  # -------------------------------------------------------------------------

  txt(
    "Median absolute expression difference from intact DP",
    5.0, 29.9,
    PT_SECTION_HEADING,
    bold = TRUE
  )

  x0 <- 31.0
  width <- 43.0
  h <- 3.2
  ys <- c(22.0, 15.3)

  vals <- c(
    MEDIAN_2D,
    MEDIAN_3D
  )

  labs <- c(
    "P3 2D",
    "P3 3D"
  )

  fills <- c(
    COL$red,
    COL$green
  )

  for (i in seq_along(ys)) {

    txt(
      labs[i],
      5.0,
      ys[i],
      PT_ROW_LABEL
    )

    rect(
      x0 + width / 2,
      ys[i],
      width,
      h,
      COL$light
    )

    fill_width <- width * vals[i] / MAX_SCALE

    rect(
      x0 + fill_width / 2,
      ys[i],
      fill_width,
      h,
      fills[i]
    )

    txt(
      sprintf("%.4f", vals[i]),
      83.0,
      ys[i],
      PT_VALUE,
      COL$ink,
      TRUE,
      "right"
    )
  }

  txt(
    "Both bars use a common scale (0–0.42)",
    31.0, 10.7,
    PT_SCALE_NOTE,
    COL$muted
  )

  txt(
    "Descriptive gene-level comparison; no global restoration advantage established",
    5.0, 5.2,
    PT_CONCLUSION,
    COL$muted
  )
}

# --------- 5. Export five artifacts; no Cairo dependency --------
parent <- dirname(OUT)
if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(parent)) stop("无法创建输出父目录：", parent)

# Hidden staging folder prevents partially generated final output.
STAGE <- file.path(parent, paste0(".Panel_F_v3_stage_", Sys.getpid()))
if (dir.exists(STAGE)) stop("临时目录已存在：", STAGE)
if (!dir.create(STAGE)) stop("无法创建临时输出目录：", STAGE)

render <- function(path, open_device) {
  open_device()
  closed <- FALSE
  on.exit({
    if (!closed) try(grDevices::dev.off(), silent = TRUE)
  }, add = TRUE)
  draw_panel()
  grDevices::dev.off()
  closed <- TRUE
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 1000) {
    stop("图形输出缺失或文件过小：", path)
  }
}

check_signature <- function(path, expected_raw) {
  actual <- readBin(path, what = "raw", n = length(expected_raw))
  if (!identical(actual, expected_raw)) stop("文件格式签名校验失败：", path)
}

tryCatch({
  pdf_path <- file.path(STAGE, paste0(BASE, ".pdf"))
  png_path <- file.path(STAGE, paste0(BASE, ".png"))
  tif_path <- file.path(STAGE, paste0(BASE, ".tif"))
  csv_path <- file.path(STAGE, "plot_data.csv")
  provenance_path <- file.path(STAGE, "README_provenance.txt")
  
  render(pdf_path, function() grDevices::pdf(
    file = pdf_path,
    width = W_MM / 25.4, height = H_MM / 25.4,
    family = FONT, useDingbats = FALSE, bg = "white"
  ))
  
  # The user's Mac R environment previously failed with Cairo DLL.
  # quartz is native to macOS and was verified in their earlier exports.
  render(png_path, function() grDevices::png(
    filename = png_path, width = W_MM, height = H_MM,
    units = "mm", res = DPI, type = "quartz", bg = "white"
  ))
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    render(tif_path, function() ragg::agg_tiff(
      filename = tif_path,
      width = W_MM, height = H_MM, units = "mm", res = DPI,
      compression = "lzw", background = "white"
    ))
    tiff_note <- "ragg TIFF; LZW compression"
  } else {
    render(tif_path, function() grDevices::tiff(
      filename = tif_path, width = W_MM, height = H_MM,
      units = "mm", res = DPI, type = "quartz", bg = "white"
    ))
    tiff_note <- "macOS Quartz TIFF; no external Cairo dependency"
  }
  
  write.csv(data.frame(
    metric = c("total_genes", "closer_pct", "not_closer_pct",
               "median_abs_diff_P3_2D", "median_abs_diff_P3_3D"),
    value = c(N, CLOSER, NOT_CLOSER, MEDIAN_2D, MEDIAN_3D)
  ), csv_path, row.names = FALSE)
  
  writeLines(c(
    "Figure 1 Panel F - v3 text refinement",
    "Source: user-provided Fig1_PanelF(2).R (frozen descriptive statistics).",
    "Redesigned directly at final 88 x 56 mm panel size; frozen statistics unchanged.",
    "Bar proportions, colors and numerical values preserved; geometry adapted to final panel size.",
    "No raw gene-expression matrix accessed. No inference or efficacy claim.",
    paste("Gene count:", N),
    paste("Closer / not closer:", CLOSER, "/", NOT_CLOSER, "percent"),
    paste("2D / 3D median absolute differences:", MEDIAN_2D, "/", MEDIAN_3D),
    paste("Common scale: 0 to", MAX_SCALE),
    paste("Canvas:", W_MM, "x", H_MM, "mm; raster DPI:", DPI),
    tiff_note,
    "No panel letter; old outputs not modified."
  ), provenance_path, useBytes = TRUE)
  
  check_signature(pdf_path, charToRaw("%PDF-"))
  check_signature(png_path, as.raw(c(137,80,78,71,13,10,26,10)))
  tiff_signature <- readBin(tif_path, "raw", n = 4)
  tiff_ok <- identical(tiff_signature, as.raw(c(73,73,42,0))) ||
    identical(tiff_signature, as.raw(c(77,77,0,42))) ||
    identical(tiff_signature, as.raw(c(73,73,43,0)))
  if (!tiff_ok) stop("TIFF 文件格式校验失败：", tif_path)
  
  expected <- c(paste0(BASE, c(".pdf", ".png", ".tif")),
                "plot_data.csv", "README_provenance.txt")
  if (!setequal(list.files(STAGE), expected) || length(list.files(STAGE)) != 5L) {
    stop("输出文件数量或名称不符合预期（应为5个文件）。")
  }
  if (!file.rename(STAGE, OUT)) {
    stop("无法将临时结果目录设为最终输出目录：", OUT)
  }
  cat("SUCCESS: Panel F Nature final exported (5 files)\n", OUT, "\n", sep = "")
  cat(paste("-", expected, collapse = "\n"), "\n")
}, error = function(e) {
  if (dir.exists(STAGE)) unlink(STAGE, recursive = TRUE, force = TRUE)
  stop("Panel F v3 未完成，旧版结果未受影响：", conditionMessage(e), call. = FALSE)
})
