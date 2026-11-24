# Figure 1 | Panel E | Nature final | 26 representative strict rescue genes (18 samples)
# Save as scripts/FIGURE 1/Fig1_PanelE_Nature_final.R
# Run from HairRegeneration_AI project root:
# source("scripts/FIGURE 1/Fig1_PanelE_Nature_final.R")
# Native grid graphics: no Cairo/X11 or third-party plotting packages.
# This panel is descriptive, NOT independent validation of rescue or hair growth.
# 26 genes are objectively selected by Combined_FDR (13 per rescue direction).

WIDTH_MM <- 88
HEIGHT_MM <- 64
DPI <- 600
N_PER_DIRECTION <- 13L
Z_LIMIT <- 2
FONT <- "Helvetica"
BASE <- "Fig1E_Heatmap_Nature_tall_final4"

# Keep the six conditions in biological progression; donor-matched samples retained.
GROUP_LEVELS <- c("Intact", "P0", "P1", "P3_2D", "P5_2D", "P3_3D")
GROUP_LABELS <- c("Intact", "P0", "P1", "P3 2D", "P5 2D", "P3 3D")
GROUP_COLORS <- c("#405F72", "#829AAA", "#AEC1CA",
                  "#34A49A", "#B78361", "#087D79")
UP_COLOR <- "#B78361"
DOWN_COLOR <- "#078C82"
INK <- "#263A45"
MUTED <- "#728894"

find_root <- function() {
  wd <- normalizePath(getwd(), mustWork=TRUE)
  candidates <- unique(c(wd, dirname(wd), dirname(dirname(wd))))
  hit <- candidates[file.exists(file.path(candidates,"HairRegeneration_AI.Rproj"))]
  if (length(hit)) return(hit[1])
  # Accept project folder even if .Rproj was renamed
  hit <- candidates[dir.exists(file.path(candidates,"results")) &
                      dir.exists(file.path(candidates,"scripts"))]
  if (length(hit)) return(hit[1])
  stop("Open HairRegeneration_AI.Rproj first; project root not found.")
}

ROOT <- find_root()
rpath <- function(...) file.path(ROOT,...)
OUT <- rpath("figures","FIGURE 1","Panel_E_Nature_tall_final4")

if (file.exists(OUT) || dir.exists(OUT)) {
  stop("Refusing to overwrite existing output: ",OUT,
       "\nRename the old Panel_E folder before rerunning.")
}

find_input <- function(paths, name) {
  hits <- paths[file.exists(paths)]
  if (!length(hits)) stop(name," not found. Checked:\n",paste(paths,collapse="\n"))
  hits[1]
}

STRICT_FILE <- find_input(c(
  rpath("results","GSE44765_rescue","GSE44765_STRICT_rescue_genes.csv"),
  rpath("results","GSE44765_rescue","GSE44765_rescue_analysis_all_genes.csv")
), "Strict rescue source")

EXPR_FILE <- find_input(c(
  rpath("results","GSE44765_rescue","GSE44765_gene_expression_RMA.csv.gz"),
  rpath("results","GSE44765_RMA_gene_expression.csv.gz"),
  rpath("results","GSE44765_rescue","GSE44765_RMA_gene_expression.csv.gz")
), "Gene-level RMA expression matrix")

META_FILE <- find_input(c(
  rpath("results","GSE44765_sample_metadata.csv"),
  rpath("results","GSE44765_rescue","GSE44765_sample_metadata.csv")
), "Sample metadata")

strict <- read.csv(STRICT_FILE, check.names=FALSE, stringsAsFactors=FALSE)

for (k in c("GeneSymbol","Rescue_type","Combined_FDR")) {
  if (!(k %in% names(strict))) stop("Missing strict-rescue column: ",k)
}

if ("Strict_rescue" %in% names(strict)) {
  flag <- toupper(as.character(strict$Strict_rescue))
  strict <- strict[flag %in% c("TRUE","T","1"),,drop=FALSE]
}

strict <- strict[!is.na(strict$Combined_FDR) &
                   is.finite(as.numeric(strict$Combined_FDR)),,drop=FALSE]

strict$Combined_FDR <- as.numeric(strict$Combined_FDR)

if (anyDuplicated(strict$GeneSymbol)) stop("Duplicate GeneSymbol in strict list")

up <- strict[strict$Rescue_type=="Restored_upward",,drop=FALSE]
down <- strict[strict$Rescue_type=="Restored_downward",,drop=FALSE]

if (nrow(up)<N_PER_DIRECTION || nrow(down)<N_PER_DIRECTION) {
  stop("Insufficient strict genes by direction. Up=",nrow(up)," Down=",nrow(down))
}

# Deterministic, transparent ranking. No manual cherry-picking.
up <- head(up[order(up$Combined_FDR,up$GeneSymbol),,drop=FALSE],N_PER_DIRECTION)
down <- head(down[order(down$Combined_FDR,down$GeneSymbol),,drop=FALSE],N_PER_DIRECTION)

selected <- c(up$GeneSymbol,down$GeneSymbol)

if (anyDuplicated(selected)) stop("Selected genes are duplicated")

meta <- read.csv(META_FILE,check.names=FALSE,stringsAsFactors=FALSE)

if (!all(c("sample","donor","group") %in% names(meta))) {
  stop("Metadata must contain sample, donor, group")
}

normalize_group <- function(s) {
  x <- gsub("[^A-Z0-9]", "", toupper(trimws(as.character(s))))
  map <- c(INTACT="Intact",INTACTDP="Intact",P0="P0",P1="P1",
           P32D="P3_2D",P52D="P5_2D",P33D="P3_3D",
           P3SPHEROID="P3_3D")
  out <- unname(map[x])
  if (anyNA(out)) stop("Unrecognized sample groups: ",
                       paste(unique(as.character(s[is.na(out)])),collapse=", "))
  out
}

meta$group_norm <- normalize_group(meta$group)
meta$donor <- toupper(trimws(as.character(meta$donor)))

if (!setequal(meta$group_norm,GROUP_LEVELS) || nrow(meta)!=18L ||
    anyDuplicated(meta$sample) || anyDuplicated(paste(meta$group_norm,meta$donor))) {
  stop("Expected 18 unique samples: 3 donor-matched samples x 6 conditions")
}

donors <- sort(unique(meta$donor))

if (length(donors)!=3L || !all(table(meta$group_norm)==3L) ||
    any(vapply(GROUP_LEVELS,function(g) !setequal(meta$donor[meta$group_norm==g],donors),
               logical(1)))) {
  stop("Donor pairing is incomplete")
}

meta <- meta[order(match(meta$group_norm,GROUP_LEVELS),
                   match(meta$donor,donors)),,drop=FALSE]

# Matrix is gene-level RMA log2 expression: first column holds gene IDs.
expr <- read.csv(EXPR_FILE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)

if (anyDuplicated(rownames(expr))) stop("Duplicate expression gene IDs")

if (!all(selected %in% rownames(expr))) {
  stop("Missing selected genes: ",paste(setdiff(selected,rownames(expr)),collapse=", "))
}

if (!all(meta$sample %in% colnames(expr))) {
  stop("Expression columns missing samples: ",
       paste(setdiff(meta$sample,colnames(expr)),collapse=", "))
}

mat <- as.matrix(expr[selected,meta$sample,drop=FALSE])
storage.mode(mat) <- "double"

if (any(!is.finite(mat))) stop("Missing or nonfinite expression values")
if (nrow(mat)!=26L || ncol(mat)!=18L) stop("Unexpected matrix dimensions")

# Per-gene Z-score across all 18 INDIVIDUAL samples (not group averages).
means <- rowMeans(mat)
sds <- apply(mat,1,sd)

if (any(sds<=0)) stop("Zero-variance gene cannot be Z-scored")

z <- sweep(sweep(mat,1,means,"-"),1,sds,"/")

# Preserve the 26 x 18 matrix dimensions during display clipping.
z_display <- z
z_display[z_display < -Z_LIMIT] <- -Z_LIMIT
z_display[z_display >  Z_LIMIT] <-  Z_LIMIT

stopifnot(is.matrix(z_display), identical(dim(z_display), c(26L, 18L)))
stopifnot(identical(rownames(z_display), selected),
          identical(colnames(z_display), as.character(meta$sample)))

# Restrained blue-white-red diverging palette.
palette <- grDevices::colorRampPalette(c("#356B92","#F7F7F5","#C46D5C"))(201)

z_to_col <- function(v) {
  index <- pmax(1L,pmin(201L,as.integer(round((v+Z_LIMIT)/(2*Z_LIMIT)*200))+1L))
  palette[index]
}

# Target placement: 88 x 56 mm in Illustrator at 100% size.
draw_panel <- function() {
  grid::grid.newpage()

  grid::pushViewport(
    grid::viewport(
      width = grid::unit(WIDTH_MM, "mm"),
      height = grid::unit(HEIGHT_MM, "mm")
    )
  )

  on.exit(grid::popViewport(), add = TRUE)

  rect_mm <- function(x, y, w, h, fill, border = NA, lwd = 0.3) {
    grid::grid.rect(
      x = grid::unit(x, "mm"),
      y = grid::unit(y, "mm"),
      width = grid::unit(w, "mm"),
      height = grid::unit(h, "mm"),
      gp = grid::gpar(
        fill = fill,
        col = border,
        lwd = lwd
      )
    )
  }

  text_mm <- function(
    label,
    x,
    y,
    pt,
    col = INK,
    just = "centre",
    bold = FALSE,
    rot = 0
  ) {
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

  # -------------------------------------------------------------------------
  # Heatmap geometry
  #
  # Explicit vertical hierarchy:
  # condition labels  62.8 mm
  # group colour bar  60.4 mm
  # donor labels      58.1 mm
  # heatmap top       56.8 mm
  #
  # This prevents the colour bar from overlapping donor labels.
  # -------------------------------------------------------------------------

  x0 <- 12.8
  x1 <- 75.8

  y0 <- 8.70
  y1 <- 56.80

  cw <- (x1 - x0) / 18
  rh <- (y1 - y0) / 26

  stopifnot(
    cw > 3.4,
    rh > 1.80
  )

  CONDITION_Y <- 62.80
  GROUP_BAR_Y <- 60.40
  DONOR_Y <- 58.10

  # -------------------------------------------------------------------------
  # Six experimental groups in biological order
  # -------------------------------------------------------------------------

  for (g in seq_along(GROUP_LEVELS)) {

    left <- x0 + (g - 1L) * 3 * cw

    rect_mm(
      left + 1.5 * cw,
      GROUP_BAR_Y,
      3 * cw,
      1.15,
      GROUP_COLORS[g]
    )

    text_mm(
      GROUP_LABELS[g],
      left + 1.5 * cw,
      CONDITION_Y,
      5.85,
      INK,
      bold = FALSE
    )

    if (g > 1) {
      grid::grid.lines(
        x = grid::unit(c(left, left), "mm"),
        y = grid::unit(c(y0, y1), "mm"),
        gp = grid::gpar(
          col = "#FFFFFF",
          lwd = 1.5
        )
      )
    }
  }

  # -------------------------------------------------------------------------
  # Donor labels
  # -------------------------------------------------------------------------

  for (j in seq_len(18)) {
    text_mm(
      meta$donor[j],
      x0 + (j - 0.5) * cw,
      DONOR_Y,
      5.50,
      MUTED
    )
  }

  # -------------------------------------------------------------------------
  # Heatmap
  # -------------------------------------------------------------------------

  for (i in seq_len(26)) {

    yy <- y1 - (i - 0.5) * rh

    for (j in seq_len(18)) {

      rect_mm(
        x0 + (j - 0.5) * cw,
        yy,
        cw,
        rh,
        z_to_col(z_display[i, j]),
        border = "#FFFFFF",
        lwd = 0.20
      )
    }

    text_mm(
      selected[i],
      x0 - 0.90,
      yy,
      5.45,
      INK,
      just = "right",
      bold = TRUE
    )

    rect_mm(
      77.55,
      yy,
      0.72,
      rh * 0.84,
      if (i <= 13) UP_COLOR else DOWN_COLOR
    )
  }

  # -------------------------------------------------------------------------
  # Visual divider between restored-upward and restored-downward genes
  # -------------------------------------------------------------------------

  mid <- y1 - 13 * rh

  grid::grid.lines(
    x = grid::unit(c(x0, x1), "mm"),
    y = grid::unit(c(mid, mid), "mm"),
    gp = grid::gpar(
      col = "#FFFFFF",
      lwd = 2.0
    )
  )

  # -------------------------------------------------------------------------
  # Rescue-direction labels
  # -------------------------------------------------------------------------

  text_mm(
    "Restored upward",
    81.7,
    y1 - 6.5 * rh,
    5.80,
    UP_COLOR,
    rot = 270
  )

  text_mm(
    "Restored downward",
    81.7,
    y0 + 6.5 * rh,
    5.80,
    DOWN_COLOR,
    rot = 270
  )

  # -------------------------------------------------------------------------
  # Colour key
  # -------------------------------------------------------------------------

  key_x0 <- 31.8
  key_x1 <- 57.2
  key_y <- 7.55

  for (k in seq_len(201)) {

    xx <- key_x0 + (k - 0.5) * (key_x1 - key_x0) / 201

    rect_mm(
      xx,
      key_y,
      (key_x1 - key_x0) / 201 + 0.035,
      1.05,
      palette[k]
    )
  }

  for (v in c(-2, 0, 2)) {

    xx <- key_x0 + (v + 2) / 4 * (key_x1 - key_x0)

    text_mm(
      as.character(v),
      xx,
      5.90,
      5.65,
      MUTED
    )
  }

  text_mm(
    "Gene-wise Z-score (18 samples; display clipped at ±2)",
    44.5,
    3.30,
    5.55,
    MUTED
  )
}

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

dir.create(dirname(OUT),recursive=TRUE,showWarnings=FALSE)

STAGE <- paste0(OUT,"_INCOMPLETE_",Sys.getpid())

if (!dir.create(STAGE)) {
  stop("Cannot create staging directory: ",STAGE)
}

ok <- FALSE

tryCatch({

  render <- function(open_device) {
    open_device()
    tryCatch(
      draw_panel(),
      finally=grDevices::dev.off()
    )
  }

  # Vector PDF
  render(function() {
    grDevices::pdf(
      file.path(STAGE,paste0(BASE,".pdf")),
      width=WIDTH_MM/25.4,
      height=HEIGHT_MM/25.4,
      family=FONT,
      useDingbats=FALSE
    )
  })

  # 600-dpi PNG
  render(function() {
    grDevices::png(
      file.path(STAGE,paste0(BASE,".png")),
      width=WIDTH_MM/25.4,
      height=HEIGHT_MM/25.4,
      units="in",
      res=DPI,
      type="quartz",
      bg="white"
    )
  })

  # TIFF
  if (requireNamespace("ragg", quietly=TRUE)) {

    render(function() {
      ragg::agg_tiff(
        file.path(STAGE,paste0(BASE,".tif")),
        width=WIDTH_MM/25.4,
        height=HEIGHT_MM/25.4,
        units="in",
        res=DPI,
        compression="lzw",
        background="white"
      )
    })

    tiff_note <- "TIFF: ragg, 600 dpi, LZW compression"

  } else {

    render(function() {
      grDevices::tiff(
        file.path(STAGE,paste0(BASE,".tif")),
        width=WIDTH_MM/25.4,
        height=HEIGHT_MM/25.4,
        units="in",
        res=DPI,
        type="quartz",
        bg="white"
      )
    })

    tiff_note <- "TIFF: Quartz, 600 dpi, compression unspecified (install ragg for LZW)"
  }

  plot_data <- data.frame(
    GeneSymbol=rep(selected,times=18),
    Sample=rep(meta$sample,each=26),
    Donor=rep(meta$donor,each=26),
    Group=rep(meta$group_norm,each=26),
    Expression_log2=as.vector(mat),
    Z_unclipped=as.vector(z),
    Z_display=as.vector(z_display),
    Rescue_type=rep(
      c(
        rep("Restored_upward",13),
        rep("Restored_downward",13)
      ),
      18
    ),
    stringsAsFactors=FALSE
  )

  write.csv(
    plot_data,
    file.path(STAGE,"plot_data.csv"),
    row.names=FALSE
  )

  write.csv(
    data.frame(
      GeneSymbol=selected,
      Rescue_type=c(
        rep("Restored_upward",13),
        rep("Restored_downward",13)
      ),
      Combined_FDR=c(
        up$Combined_FDR,
        down$Combined_FDR
      )
    ),
    file.path(STAGE,"selected_genes.csv"),
    row.names=FALSE
  )

  writeLines(
    c(
      "Figure 1 Panel E — representative strict rescue genes",
      "Descriptive display; no independent validation and no functional proof.",
      paste("Strict input:",STRICT_FILE),
      paste("Expression input:",EXPR_FILE),
      paste("Metadata input:",META_FILE),
      paste("Strict input MD5:",unname(tools::md5sum(STRICT_FILE))),
      "Selection: 13 lowest Combined_FDR genes per strict rescue direction; GeneSymbol tie-break.",
      "Expression: 18 donor-matched samples, not condition averages.",
      "Z score: per gene across 18 samples; sample SD; clipped at [-2,2] for colours only.",
      "Rows ordered by direction then Combined_FDR; columns by condition then donor.",
      tiff_note,
      "Print-size layout: 88 x 64 mm; explicit top hierarchy: condition labels 62.8 mm, group bars 60.4 mm, donor labels 58.1 mm, heatmap top 56.8 mm.",
      "Condition naming unified to P3 3D; all other layout unchanged.",
      "No panel letter; no assembly; original data unchanged.",
      "All graphical coordinates in physical mm; print at 100% to preserve final label sizes."
    ),
    file.path(STAGE,"README_provenance.txt")
  )

  expected <- c(
    paste0(BASE,c(".pdf",".png",".tif")),
    "plot_data.csv",
    "selected_genes.csv",
    "README_provenance.txt"
  )

  if (!all(file.exists(file.path(STAGE,expected)))) {
    stop("Missing expected output")
  }

  sig <- readBin(
    file.path(STAGE,paste0(BASE,".pdf")),
    "raw",
    n=5
  )

  if (!identical(rawToChar(sig),"%PDF-")) {
    stop("PDF header invalid")
  }

  if (any(file.info(file.path(STAGE,expected[1:3]))$size < 100)) {
    stop("Image file unexpectedly small")
  }

  if (!file.rename(STAGE,OUT)) {
    stop("Cannot finalize output folder")
  }

  ok <- TRUE

},error=function(e) {

  message(
    "Panel E FAILED: ",
    conditionMessage(e)
  )

  message(
    "Incomplete staging folder retained for diagnosis: ",
    STAGE
  )

  stop(e)
})

if (ok) {
  cat(
    "SUCCESS — Panel E exported to:\n",
    OUT,
    "\n"
  )
}
