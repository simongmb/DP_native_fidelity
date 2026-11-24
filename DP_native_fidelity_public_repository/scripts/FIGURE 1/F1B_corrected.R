# ============================================================================
# HairRegAI | FIGURE 1 | PANEL B | Nature-style minimal correction
# Based directly on the original stable Fig1_PanelB.R layout
#
# Key principle:
# - KEEP original 88 x 57 mm geometry
# - KEEP original PCA coordinates and panel positions
# - ONLY make minimal editorial/style corrections
# ============================================================================

# ----------------------------- EDITABLE SETTINGS ----------------------------
WIDTH_MM  <- 88.0
HEIGHT_MM <- 57.0
DPI       <- 600L
FONT      <- "Helvetica"

OUT_DIR   <- file.path("figures", "FIGURE 1", "Panel_B_Nature_final2")
BASENAME  <- "Fig1B_PCA_Nature_final2"
META_FILE <- file.path("results", "GSE44765_sample_metadata.csv")
PCA_FILE  <- file.path("results", "GSE44765_PCA_coordinates.csv")
VAR_FILE  <- file.path("results", "GSE44765_PCA_variance.csv")

# Final-size typography
PT_HEADING <- 6.2
PT_AXIS    <- 6.0
PT_TICK    <- 5.0
PT_LEGEND  <- 5.1
PT_DONOR   <- 5.1
PT_NOTE    <- 5.2

POINT_GLOBAL_MM <- 1.75
POINT_ZOOM_MM   <- 1.85

INK       <- "#253945"
NAVY      <- "#46677A"
MUTED     <- "#748A95"
RULE      <- "#CBD9DE"
TEAL      <- "#087F79"
TEAL_DARK <- "#176176"
AMBER     <- "#B9825E"
PALE_BLUE <- "#A7BAC4"
MID_BLUE  <- "#819AA9"
WHITE     <- "#FFFFFF"

GROUP_ORDER <- c("Intact", "P0", "P1", "P3_2D", "P5_2D", "P3_3D")

GROUP_COL <- c(
  Intact = NAVY,
  P0 = MID_BLUE,
  P1 = PALE_BLUE,
  P3_2D = TEAL,
  P5_2D = AMBER,
  P3_3D = TEAL_DARK
)

GROUP_LABEL <- c(
  Intact = "Intact DP",
  P0 = "P0",
  P1 = "P1",
  P3_2D = "P3 2D",
  P5_2D = "P5 2D",
  P3_3D = "P3 3D"
)

DONOR_ORDER <- c("D5", "D6", "D7")
DONOR_PCH <- c(
  D5 = 21L,
  D6 = 22L,
  D7 = 24L
)

# ORIGINAL STABLE LAYOUT — intentionally unchanged
GLOBAL <- list(
  x = 13.0,
  y = 15.2,
  w = 42.3,
  h = 29.5,
  xlim = c(-160, 65),
  ylim = c(-103, 39)
)

ZOOM <- list(
  x = 64.4,
  y = 15.2,
  w = 19.4,
  h = 29.5,
  xlim = c(15, 54),
  ylim = c(2, 26)
)

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

if (!requireNamespace("ragg", quietly = TRUE)) {
  stop(
    "The 'ragg' package is needed for 600-dpi PNG and LZW TIFF.\n",
    "Install once using install.packages('ragg')."
  )
}

for (p in c(META_FILE, PCA_FILE, VAR_FILE)) {
  if (!file.exists(p)) {
    stop(
      "Missing source file: ",
      p,
      "\nOpen HairRegeneration_AI.Rproj in RStudio first."
    )
  }
}

if (dir.exists(OUT_DIR) || file.exists(OUT_DIR)) {
  stop(
    "Existing output protected. Rename/remove first:\n",
    OUT_DIR
  )
}

meta <- utils::read.csv(
  META_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pca <- utils::read.csv(
  PCA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

vari <- utils::read.csv(
  VAR_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  !all(c("sample", "donor", "group") %in% names(meta)) ||
  !all(c("sample", "PC1", "PC2", "group", "donor") %in% names(pca)) ||
  !all(c("PC", "VariancePercent") %in% names(vari))
) {
  stop("Unexpected input columns.")
}

if (
  nrow(meta) != 18L ||
  nrow(pca) != 18L ||
  anyNA(meta[, c("sample", "donor", "group")]) ||
  anyNA(pca[, c("sample", "donor", "group", "PC1", "PC2")]) ||
  anyDuplicated(meta$sample) ||
  anyDuplicated(pca$sample) ||
  !setequal(meta$sample, pca$sample) ||
  !setequal(unique(meta$group), GROUP_ORDER) ||
  !setequal(unique(meta$donor), DONOR_ORDER) ||
  !setequal(unique(pca$group), GROUP_ORDER) ||
  !setequal(unique(pca$donor), DONOR_ORDER)
) {
  stop(
    "Expected 18 matching profiles, 6 culture conditions and donors D5/D6/D7."
  )
}

meta <- meta[
  match(pca$sample, meta$sample),
  ,
  drop = FALSE
]

if (
  any(meta$group != pca$group) ||
  any(meta$donor != pca$donor)
) {
  stop("PCA sample identities disagree with sample_metadata.csv.")
}

if (
  any(
    table(
      factor(pca$group, levels = GROUP_ORDER),
      factor(pca$donor, levels = DONOR_ORDER)
    ) != 1L
  )
) {
  stop("Expected one sample from each donor in each condition.")
}

if (
  !is.numeric(pca$PC1) ||
  !is.numeric(pca$PC2) ||
  any(!is.finite(pca$PC1)) ||
  any(!is.finite(pca$PC2))
) {
  stop("PCA coordinates must be finite numeric values.")
}

pc1_pct <- vari$VariancePercent[
  as.character(vari$PC) == "PC1"
]

pc2_pct <- vari$VariancePercent[
  as.character(vari$PC) == "PC2"
]

if (
  length(pc1_pct) != 1L ||
  length(pc2_pct) != 1L ||
  any(!is.finite(c(pc1_pct, pc2_pct))) ||
  any(c(pc1_pct, pc2_pct) <= 0) ||
  sum(c(pc1_pct, pc2_pct)) > 100.01
) {
  stop("Unexpected PCA explained variance in source file.")
}

if (
  any(
    pca$PC1 < GLOBAL$xlim[1] |
    pca$PC1 > GLOBAL$xlim[2] |
    pca$PC2 < GLOBAL$ylim[1] |
    pca$PC2 > GLOBAL$ylim[2]
  )
) {
  stop(
    "Some PCA points fall outside the defined global axes."
  )
}

selected_2d <- pca$group %in% c(
  "P0",
  "P1",
  "P3_2D",
  "P5_2D"
)

if (
  sum(selected_2d) != 12L ||
  any(
    pca$PC1[selected_2d] < ZOOM$xlim[1] |
    pca$PC1[selected_2d] > ZOOM$xlim[2] |
    pca$PC2[selected_2d] < ZOOM$ylim[1] |
    pca$PC2[selected_2d] > ZOOM$ylim[2]
  )
) {
  stop(
    "The 12 2D profiles do not fit the predefined zoom window."
  )
}

plot_data <- pca[
  order(
    match(pca$group, GROUP_ORDER),
    match(pca$donor, DONOR_ORDER)
  ),
  ,
  drop = FALSE
]

plot_data$view_inset <-
  plot_data$group %in% c(
    "P0",
    "P1",
    "P3_2D",
    "P5_2D"
  )

# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

mm <- function(v) {
  grid::unit(v, "native")
}

text_mm <- function(
    s,
    x,
    y,
    size = PT_AXIS,
    col = INK,
    just = "centre",
    bold = FALSE,
    rot = 0
) {
  grid::grid.text(
    s,
    x = mm(x),
    y = mm(y),
    just = just,
    rot = rot,
    gp = grid::gpar(
      fontfamily = FONT,
      fontsize = size,
      col = col,
      fontface = if (bold) "bold" else "plain"
    )
  )
}

line_mm <- function(
    x1,
    y1,
    x2,
    y2,
    col = RULE,
    lwd = 0.45,
    lty = "solid"
) {
  grid::grid.lines(
    x = mm(c(x1, x2)),
    y = mm(c(y1, y2)),
    gp = grid::gpar(
      col = col,
      lwd = lwd,
      lty = lty
    )
  )
}

push_data <- function(spec, clip = "on") {
  grid::pushViewport(
    grid::viewport(
      x = mm(spec$x),
      y = mm(spec$y),
      width = mm(spec$w),
      height = mm(spec$h),
      just = c("left", "bottom"),
      xscale = spec$xlim,
      yscale = spec$ylim,
      clip = clip
    )
  )
}

plot_axes <- function(
    spec,
    xticks,
    yticks,
    xlab,
    ylab,
    zoom = FALSE
) {
  push_data(spec, clip = "off")
  
  xl <- spec$xlim
  yl <- spec$ylim
  
  # axes
  grid::grid.lines(
    x = grid::unit(xl, "native"),
    y = grid::unit(rep(yl[1], 2), "native"),
    gp = grid::gpar(
      col = INK,
      lwd = 0.55
    )
  )
  
  grid::grid.lines(
    x = grid::unit(rep(xl[1], 2), "native"),
    y = grid::unit(yl, "native"),
    gp = grid::gpar(
      col = INK,
      lwd = 0.55
    )
  )
  
  # x ticks
  for (v in xticks) {
    grid::grid.lines(
      x = grid::unit(rep(v, 2), "native"),
      y =
        grid::unit(yl[1], "native") +
        grid::unit(c(0, -0.75), "mm"),
      gp = grid::gpar(
        col = INK,
        lwd = 0.38
      )
    )
    
    grid::grid.text(
      as.character(v),
      x = grid::unit(v, "native"),
      y =
        grid::unit(yl[1], "native") -
        grid::unit(1.6, "mm"),
      just = "top",
      gp = grid::gpar(
        fontfamily = FONT,
        fontsize = PT_TICK,
        col = MUTED
      )
    )
  }
  
  # y ticks
  for (v in yticks) {
    grid::grid.lines(
      x =
        grid::unit(xl[1], "native") +
        grid::unit(c(0, -0.75), "mm"),
      y = grid::unit(rep(v, 2), "native"),
      gp = grid::gpar(
        col = INK,
        lwd = 0.38
      )
    )
    
    grid::grid.text(
      as.character(v),
      x =
        grid::unit(xl[1], "native") -
        grid::unit(1.4, "mm"),
      y = grid::unit(v, "native"),
      just = "right",
      gp = grid::gpar(
        fontfamily = FONT,
        fontsize = PT_TICK,
        col = MUTED
      )
    )
  }
  
  grid::popViewport()
  
  text_mm(
    xlab,
    spec$x + spec$w / 2,
    10.3,
    size = PT_AXIS
  )
  
  text_mm(
    ylab,
    if (zoom) 59.0 else 7.5,
    spec$y + spec$h / 2,
    size = PT_AXIS,
    rot = 90
  )
}

draw_points <- function(
    df,
    spec,
    diameter_mm
) {
  push_data(spec, clip = "on")
  
  for (i in seq_len(nrow(df))) {
    
    g <- as.character(df$group[i])
    d <- as.character(df$donor[i])
    
    grid::grid.points(
      x = grid::unit(
        df$PC1[i],
        "native"
      ),
      y = grid::unit(
        df$PC2[i],
        "native"
      ),
      pch = DONOR_PCH[[d]],
      size = grid::unit(
        diameter_mm,
        "mm"
      ),
      gp = grid::gpar(
        fill = GROUP_COL[[g]],
        col = WHITE,
        lwd = 0.55
      )
    )
  }
  
  grid::popViewport()
}

# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------

draw_panel <- function() {
  
  grid::grid.newpage()
  
  grid::pushViewport(
    grid::viewport(
      xscale = c(0, WIDTH_MM),
      yscale = c(0, HEIGHT_MM)
    )
  )
  
  grid::grid.rect(
    gp = grid::gpar(
      fill = WHITE,
      col = NA
    )
  )
  
  # Headers
  text_mm(
    "All conditions",
    13.0,
    52.6,
    size = PT_HEADING,
    bold = TRUE,
    just = "left"
  )
  
  text_mm(
    "2D detail",
    64.4,
    52.6,
    size = PT_HEADING,
    bold = TRUE,
    just = "left"
  )
  
  # Donor shape key
  donor_x <- c(
    35.1,
    43.5,
    51.8
  )
  
  for (i in seq_along(DONOR_ORDER)) {
    
    grid::grid.points(
      x = mm(donor_x[i]),
      y = mm(52.6),
      pch = DONOR_PCH[[DONOR_ORDER[i]]],
      size = grid::unit(
        1.55,
        "mm"
      ),
      gp = grid::gpar(
        fill = MUTED,
        col = WHITE,
        lwd = 0.50
      )
    )
    
    text_mm(
      DONOR_ORDER[i],
      donor_x[i] + 1.50,
      52.6,
      size = PT_DONOR,
      col = MUTED,
      just = "left"
    )
  }
  
  # Global PCA
  plot_axes(
    GLOBAL,
    xticks = c(
      -150,
      -100,
      -50,
      0,
      50
    ),
    yticks = c(
      -100,
      -50,
      0,
      25
    ),
    xlab = sprintf(
      "PC1 (%.1f%%)",
      pc1_pct
    ),
    ylab = sprintf(
      "PC2 (%.1f%%)",
      pc2_pct
    )
  )
  
  # Zoom rectangle
  push_data(
    GLOBAL,
    clip = "on"
  )
  
  grid::grid.rect(
    x = grid::unit(
      mean(ZOOM$xlim),
      "native"
    ),
    y = grid::unit(
      mean(ZOOM$ylim),
      "native"
    ),
    width = grid::unit(
      diff(ZOOM$xlim),
      "native"
    ),
    height = grid::unit(
      diff(ZOOM$ylim),
      "native"
    ),
    gp = grid::gpar(
      fill = NA,
      col = RULE,
      lwd = 0.55,
      lty = "22"
    )
  )
  
  grid::popViewport()
  
  draw_points(
    plot_data,
    GLOBAL,
    POINT_GLOBAL_MM
  )
  
  # Direct labels retained because ORIGINAL geometry has adequate space
  push_data(
    GLOBAL,
    clip = "on"
  )
  
  grid::grid.text(
    "Intact DP",
    x = grid::unit(
      -152,
      "native"
    ),
    y = grid::unit(
      33,
      "native"
    ),
    just = "left",
    gp = grid::gpar(
      fontfamily = FONT,
      fontsize = PT_NOTE,
      fontface = "bold",
      col = INK
    )
  )
  
  grid::grid.text(
    "P3 3D",
    x = grid::unit(
      3,
      "native"
    ),
    y = grid::unit(
      -76,
      "native"
    ),
    just = "left",
    gp = grid::gpar(
      fontfamily = FONT,
      fontsize = PT_NOTE,
      fontface = "bold",
      col = INK
    )
  )
  
  grid::popViewport()
  
  # Same PCA, enlarged 2D region
  plot_axes(
    ZOOM,
    xticks = c(
      20,
      35,
      50
    ),
    yticks = c(
      5,
      15,
      25
    ),
    xlab = "PC1",
    ylab = "PC2",
    zoom = TRUE
  )
  
  draw_points(
    plot_data[
      plot_data$view_inset,
      ,
      drop = FALSE
    ],
    ZOOM,
    POINT_ZOOM_MM
  )
  
  # Condition legend
  legend_text <- unname(
    GROUP_LABEL[GROUP_ORDER]
  )
  
  label_w <- vapply(
    legend_text,
    function(label) {
      
      grid::convertWidth(
        grid::grobWidth(
          grid::textGrob(
            label,
            gp = grid::gpar(
              fontfamily = FONT,
              fontsize = PT_LEGEND
            )
          )
        ),
        "mm",
        valueOnly = TRUE
      )
    },
    numeric(1)
  )
  
  symbol_w <- 1.40
  symbol_to_text <- 0.90
  item_gap <- 2.10
  
  item_w <-
    symbol_w +
    symbol_to_text +
    label_w
  
  total_w <-
    sum(item_w) +
    item_gap *
    (length(item_w) - 1)
  
  if (
    total_w >
    WIDTH_MM - 4
  ) {
    stop(
      sprintf(
        "Legend too wide (%.1f mm) for %.1f mm panel",
        total_w,
        WIDTH_MM
      )
    )
  }
  
  cursor <-
    (WIDTH_MM - total_w) /
    2
  
  for (
    i in seq_along(
      GROUP_ORDER
    )
  ) {
    
    g <- GROUP_ORDER[i]
    
    grid::grid.points(
      x = mm(
        cursor +
          symbol_w / 2
      ),
      y = mm(3.5),
      pch = 21,
      size = grid::unit(
        symbol_w,
        "mm"
      ),
      gp = grid::gpar(
        fill =
          GROUP_COL[[g]],
        col = WHITE,
        lwd = 0.45
      )
    )
    
    text_mm(
      GROUP_LABEL[[g]],
      cursor +
        symbol_w +
        symbol_to_text,
      3.5,
      size = PT_LEGEND,
      col = MUTED,
      just = "left"
    )
    
    cursor <-
      cursor +
      item_w[i] +
      item_gap
  }
  
  grid::popViewport()
}

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

export_panel <- function() {
  
  parent <- dirname(
    OUT_DIR
  )
  
  dir.create(
    parent,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (
    dir.exists(OUT_DIR) ||
    file.exists(OUT_DIR)
  ) {
    stop(
      "Refusing to overwrite existing output: ",
      OUT_DIR
    )
  }
  
  stage <- tempfile(
    pattern =
      ".PanelB_staging_",
    tmpdir = parent
  )
  
  if (
    !dir.create(
      stage,
      showWarnings = FALSE
    )
  ) {
    stop(
      "Cannot create output staging folder"
    )
  }
  
  completed <- FALSE
  
  on.exit({
    
    if (
      !completed &&
      dir.exists(stage)
    ) {
      message(
        "Incomplete staged files retained for diagnosis: ",
        stage
      )
    }
    
  }, add = TRUE)
  
  pdf_file <- file.path(
    stage,
    paste0(
      BASENAME,
      ".pdf"
    )
  )
  
  png_file <- file.path(
    stage,
    paste0(
      BASENAME,
      ".png"
    )
  )
  
  tif_file <- file.path(
    stage,
    paste0(
      BASENAME,
      ".tif"
    )
  )
  
  w_in <- WIDTH_MM /
    25.4
  
  h_in <- HEIGHT_MM /
    25.4
  
  render <- function(
    open_device
  ) {
    
    open_device()
    
    tryCatch(
      draw_panel(),
      finally =
        grDevices::dev.off()
    )
  }
  
  # Vector PDF
  render(
    function() {
      
      grDevices::pdf(
        pdf_file,
        width = w_in,
        height = h_in,
        family = FONT,
        useDingbats = FALSE,
        onefile = FALSE,
        paper = "special",
        bg = WHITE
      )
    }
  )
  
  # PNG 600 dpi
  render(
    function() {
      
      ragg::agg_png(
        png_file,
        width = w_in,
        height = h_in,
        units = "in",
        res = DPI,
        background = WHITE
      )
    }
  )
  
  # TIFF 600 dpi LZW
  render(
    function() {
      
      ragg::agg_tiff(
        tif_file,
        width = w_in,
        height = h_in,
        units = "in",
        res = DPI,
        compression = "lzw",
        background = WHITE
      )
    }
  )
  
  img_files <- c(
    pdf_file,
    png_file,
    tif_file
  )
  
  if (
    !all(
      file.exists(
        img_files
      )
    ) ||
    any(
      file.info(
        img_files
      )$size < 900
    )
  ) {
    stop(
      "An image file is absent or unexpectedly small."
    )
  }
  
  utils::write.csv(
    plot_data,
    file.path(
      stage,
      "plot_data.csv"
    ),
    row.names = FALSE
  )
  
  writeLines(
    c(
      "Figure 1 Panel B | Nature-style minimal correction",
      "Original PCA geometry preserved.",
      "Global view: 18 samples.",
      "Zoom view: 12 2D profiles from SAME PCA.",
      "Donors D5/D6/D7: point shapes.",
      "Six conditions: point fill colours.",
      sprintf(
        "PC1 explains %.8f%%; PC2 explains %.8f%%",
        pc1_pct,
        pc2_pct
      ),
      sprintf(
        "Canvas: %.1f x %.1f mm",
        WIDTH_MM,
        HEIGHT_MM
      ),
      "P3 spheroid label renamed to P3 3D.",
      "Operational note removed.",
      "Original data and PCA coordinates unchanged."
    ),
    file.path(
      stage,
      "README_provenance.txt"
    )
  )
  
  if (
    dir.exists(
      OUT_DIR
    ) ||
    file.exists(
      OUT_DIR
    )
  ) {
    stop(
      "Output folder appeared while rendering."
    )
  }
  
  if (
    !file.rename(
      stage,
      OUT_DIR
    )
  ) {
    stop(
      "Could not move completed output to ",
      OUT_DIR
    )
  }
  
  completed <- TRUE
  
  cat(
    "SUCCESS\n"
  )
  
  cat(
    "Output folder: ",
    normalizePath(
      OUT_DIR
    ),
    "\n",
    sep = ""
  )
}

export_panel()