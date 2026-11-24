# Supplementary Figure S1 | GSE44765 | panel B | v4 Nature-style WIDE
# Reproduces the v3 Pearson correlations and grouping, changes presentation only.
# Standalone panel: no panel letter (add 8-pt bold 'b' in Illustrator).
# The RMA matrix and design files are NEVER changed. Existing output is NEVER overwritten.
#
# Place this script in scripts/SUPPLEMENTARY FIGURES/ under HairRegeneration_AI.
# Run from the RStudio project root:
# source("scripts/SUPPLEMENTARY FIGURES/SuppFig_S1_Panel_B_v4_Nature_wide.R")

make_s1_panel_B_v4 <- function() {
  root <- getwd()
  for (pkg in c("ggplot2", "ragg")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("缺少 R 包 ", pkg, "。请运行 install.packages('", pkg, "')")
    }
  }
  suppressPackageStartupMessages(library(ggplot2))

  version <- "S1_Panel_B_v5_Nature_final"
  output_dir <- file.path(root, "figures", "SUPPLEMENTARY FIGURES", version)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("为防止覆盖，输出目录已存在且非空：\n", output_dir,
         "\n如需重运行，请为 version 设定一个新的版本名（不要删除旧结果）。")
  }

  # --- Input discovery: exactly the same two input files as v3 ---
  required <- c("GSE44765_gene_expression_RMA.rds", "GSE44765_design_matrix.csv")
  base_dir <- "GSE44765_rescue"
  candidates <- c(file.path(root, "results", base_dir),
                  file.path(root, base_dir),
                  file.path(root, "scripts", "SUPPLEMENTARY FIGURES", base_dir))
  complete <- function(d) dir.exists(d) && all(file.exists(file.path(d, required)))
  found <- candidates[vapply(candidates, complete, logical(1))]
  if (length(found)) {
    input_dir <- found[1]
  } else {
    zip_candidates <- c(file.path(root, "GSE44765_rescue.zip"),
                        file.path(root, "results", "GSE44765_rescue.zip"),
                        file.path(root, "scripts", "SUPPLEMENTARY FIGURES", "GSE44765_rescue.zip"))
    zip_candidates <- zip_candidates[file.exists(zip_candidates)]
    if (!length(zip_candidates)) {
      stop("找不到 GSE44765_rescue 输入文件夹或 ZIP；请确认已打开项目根目录。")
    }
    tmp <- tempfile("S1B_v4_input_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
    archive_entries <- utils::unzip(zip_candidates[1], list = TRUE)$Name
    # Compatible with ZIP files that do/do not contain a top-level folder.
    prefixes <- c("GSE44765_rescue/", "")
    prefix <- prefixes[vapply(prefixes, function(p) all(paste0(p, required) %in% archive_entries), logical(1))]
    if (!length(prefix)) stop("ZIP 中未找到两个必需输入文件。")
    utils::unzip(zip_candidates[1], files = paste0(prefix[1], required), exdir = tmp)
    dirs <- c(file.path(tmp, base_dir), tmp)
    found <- dirs[vapply(dirs, complete, logical(1))]
    if (!length(found)) stop("ZIP 解压后未找到完整输入文件。")
    input_dir <- found[1]
  }
  message("S1b v4 输入文件夹：", input_dir)

  expr <- as.matrix(readRDS(file.path(input_dir, required[1])))
  storage.mode(expr) <- "numeric"
  design <- read.csv(file.path(input_dir, required[2]), row.names = 1,
                     check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(expr) != 21355L || ncol(expr) != 18L ||
      any(!is.finite(expr)) || anyDuplicated(rownames(expr)) ||
      anyDuplicated(colnames(expr))) {
    stop("表达矩阵不是预期的 21,355 基因 × 18 个唯一样本，或包含非有限值。")
  }
  if (!setequal(rownames(design), colnames(expr))) stop("设计矩阵和表达矩阵的样本名不一致。")
  design <- design[colnames(expr), , drop = FALSE]
  conditions <- c("Intact", "P0", "P1", "P3_2D", "P5_2D", "P3_3D")
  group_cols <- paste0("group", conditions)
  if (!all(c(group_cols, "donorD6", "donorD7") %in% colnames(design))) {
    stop("设计矩阵缺少分组或供体哑变量。")
  }
  group_matrix <- as.matrix(design[, group_cols, drop = FALSE])
  storage.mode(group_matrix) <- "numeric"
  if (any(!is.finite(group_matrix)) || any(!(group_matrix %in% c(0, 1))) ||
      any(rowSums(group_matrix) != 1L)) stop("分组变量应为 0/1，且每个样本恰好属于一组。")
  meta <- data.frame(
    Sample = rownames(design),
    Condition = factor(conditions[max.col(group_matrix)], levels = conditions),
    Donor = factor(ifelse(design$donorD6 == 1, "D6",
                          ifelse(design$donorD7 == 1, "D7", "D5")),
                   levels = c("D5", "D6", "D7")),
    stringsAsFactors = FALSE
  )
  if (!all(table(meta$Condition, meta$Donor) == 1L)) {
    stop("六组 × 三位供体的设计不完整或重复。")
  }
  meta <- meta[order(meta$Condition, meta$Donor), , drop = FALSE]
  expr <- expr[, meta$Sample, drop = FALSE]
  meta$ConditionShort <- factor(as.character(meta$Condition), levels = conditions,
                                labels = c("Intact", "P0", "P1", "P3 2D", "P5 2D", "P3 3D"))
  if (!identical(meta$Sample, paste0("DP", 1:18)) ||
      !identical(as.character(meta$Condition), rep(conditions, each = 3L)) ||
      !identical(as.character(meta$Donor), rep(c("D5", "D6", "D7"), times = 6L))) {
    stop("DP1–DP18 / 条件 / 供体的排序与 v3 不一致，停止以防注释错位。")
  }

  # --- Unchanged scientific computation from v3 ---
  r <- stats::cor(expr, method = "pearson")
  if (any(!is.finite(r)) || !isTRUE(all.equal(r, t(r), tolerance = 1e-12)) ||
      any(abs(diag(r) - 1) > 1e-10)) stop("Pearson 相关矩阵未通过对称性或对角线检查。")
  off_diag <- r[upper.tri(r)]
  message(sprintf("S1b Pearson off-diagonal range: %.5f to %.5f", min(off_diag), max(off_diag)))
  if (min(r) < 0.85 || max(r) > 1 + 1e-10) {
    stop("实际相关系数超出原图固定色标 [0.85, 1.00]，已停止以避免截断颜色。")
  }
  df <- expand.grid(IndexX = 1:18, IndexY = 1:18, KEEP.OUT.ATTRS = FALSE)
  df$Pearson_r <- as.vector(r)
  dp_labels <- paste0("DP", 1:18)

  # Keep all SIX original condition colors and the original Pearson color scale.
  INK <- "#273C4B"
  group_labels <- c("Intact", "P0", "P1", "P3 2D", "P5 2D", "P3 3D")
  group_colors <- c("#6D879B", "#9BA9B6", "#C3A28C", "#51A79E", "#7C9B94", "#177F79")
  group_text_colors <- c("white", INK, INK, "white", INK, "white")
  group_info <- data.frame(Condition = conditions, Display = group_labels,
                           Color = group_colors,
                           FirstArray = seq(1L, 16L, by = 3L),
                           LastArray = seq(3L, 18L, by = 3L))

  # --- Final-size typography: panel letter is added separately in Illustrator. ---
  # Output dimensions are in mm, rather than inches. Place the PDF without rescaling
  # if you want the following actual type sizes at final printed size.
  FIG_WIDTH_MM <- 88
  FIG_HEIGHT_MM <- 62
  BASE_FONT <- "Helvetica"

  p <- ggplot(df, aes(x = IndexX, y = IndexY, fill = Pearson_r)) +
    geom_tile(width = 1, height = 1, linewidth = 0) +
    scale_fill_gradientn(colours = c("#F3F7F8", "#C6E4DF", "#67B7AA", "#087F78"),
                         limits = c(0.85, 1), breaks = c(0.85, 0.90, 0.95, 1),
                         labels = c("0.85", "0.90", "0.95", "1.00"),
                         name = "Pearson r") +
    scale_x_continuous(breaks = 1:18, labels = dp_labels,
                       limits = c(0.5, 18.5), expand = c(0, 0)) +
    scale_y_reverse(breaks = 1:18, labels = dp_labels,
                    limits = c(18.5, -0.75), expand = c(0, 0)) +
    coord_equal() +
    labs(title = "Sample-to-sample expression correlation", x = NULL, y = NULL) +
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barheight = grid::unit(30, "mm"),
      barwidth = grid::unit(2.2, "mm"),
      frame.colour = NA, ticks.colour = INK
    )) +
    theme_minimal(base_family = BASE_FONT, base_size = 6) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(family = BASE_FONT, colour = INK,
                                face = "bold", size = 7,
                                margin = margin(b = 2, unit = "mm")),
      plot.title.position = "plot",
      axis.text.x = element_text(size = 5, colour = INK, angle = 90,
                                 hjust = 1, vjust = 0.5,
                                 margin = margin(t = 0.7, unit = "mm")),
      axis.text.y = element_text(size = 5, colour = INK,
                                 margin = margin(r = 0.7, unit = "mm")),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 5.2, colour = INK,
                                  margin = margin(b = 1, unit = "mm")),
      legend.text = element_text(size = 5.0, colour = INK),
      legend.spacing.x = grid::unit(0.3, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(1.0, 0.8, 0.8, 0.8, unit = "mm")
    )

  # White condition separators run only THROUGH the heatmap, not the top color bar.
  for (boundary in seq(3.5, 15.5, by = 3)) {
    p <- p +
      annotate("segment", x = boundary, xend = boundary, y = 0.5, yend = 18.5,
               colour = "white", linewidth = 0.43) +
      annotate("segment", x = 0.5, xend = 18.5, y = boundary, yend = boundary,
               colour = "white", linewidth = 0.43)
  }

  # A condition strip with equal-width groups, aligned with the data column centers.
  # ggplot 'size' for text is in mm, hence conversion from points.
  point_to_mm <- 25.4 / 72.27
  for (j in seq_along(group_labels)) {
    center <- 3 * j - 1
    p <- p +
      annotate("rect", xmin = center - 1.5, xmax = center + 1.5,
               ymin = -0.50, ymax = 0.24,
               fill = group_colors[j], colour = "white", linewidth = 0.22) +
      annotate("text", x = center, y = -0.13, label = group_labels[j],
               colour = group_text_colors[j], fontface = "bold",
               family = BASE_FONT, size = 5.0 * point_to_mm)
  }

  # No long subtitle or bottom explanatory sentence. Put these in the figure legend:
  # 'Pearson correlations computed using 21,355 genes; 18 arrays ordered by
  #  experimental condition (D5, D6 and D7 within each group).'

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_base <- version
  save_figure <- function(ext) {
    dest <- file.path(output_dir, paste0(output_base, ".", ext))
    if (ext == "pdf") {
      grDevices::pdf(dest, width = FIG_WIDTH_MM / 25.4,
                     height = FIG_HEIGHT_MM / 25.4,
                     family = BASE_FONT, useDingbats = FALSE, bg = "white")
    } else if (ext == "png") {
      ragg::agg_png(dest, width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM,
                    units = "mm", res = 600, background = "white")
    } else if (ext == "tiff") {
      ragg::agg_tiff(dest, width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM,
                     units = "mm", res = 600, compression = "lzw", background = "white")
    } else stop("不支持的输出格式")
    tryCatch(print(p), finally = grDevices::dev.off())
    if (!file.exists(dest) || is.na(file.info(dest)$size) || file.info(dest)$size == 0) {
      stop("输出文件未成功生成：", dest)
    }
    message("已输出：", dest)
  }
  for (ext in c("pdf", "png", "tiff")) save_figure(ext)

  meta$DisplayedLabel <- dp_labels
  write.csv(meta, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
  write.csv(group_info, file.path(output_dir, "condition_annotation_mapping.csv"), row.names = FALSE)
  write.csv(data.frame(Sample = rownames(r), r, check.names = FALSE),
            file.path(output_dir, "Panel_B_Pearson_correlation_matrix.csv"), row.names = FALSE)
  writeLines(c(
    "Supplementary Figure S1b v5 -- Nature-style final, standalone panel",
    paste0("Input folder: ", normalizePath(input_dir, winslash = "/", mustWork = TRUE)),
    "Unchanged v3 analysis: 21,355 gene-level RMA measurements, 18 arrays, Pearson r.",
    "Unchanged order: DP1-DP18 = six conditions x donors D5, D6, D7.",
    sprintf("Observed off-diagonal Pearson r: %.6f to %.6f", min(off_diag), max(off_diag)),
    "Heatmap palette and Pearson scale fixed to 0.85-1.00, as in v3.",
    sprintf("Export size: %s x %s mm (PDF vector, PNG/TIFF 600 dpi)", FIG_WIDTH_MM, FIG_HEIGHT_MM),
    "PDF type sizes refer to this UNRESIZED final size; scaling in Illustrator scales fonts too.",
    "Panel letter b is NOT included. Add it in Illustrator at 8 pt bold.",
    "Subtitle and bottom explanation have been moved to the manuscript legend.",
    "Old outputs and all inputs were not modified.",
    paste0("R: ", R.version.string),
    paste0("ggplot2: ", as.character(utils::packageVersion("ggplot2"))),
    paste0("ragg: ", as.character(utils::packageVersion("ragg")))
  ), file.path(output_dir, "SOURCE_AUDIT.txt"))
  message("Figure S1b v4 完成：", output_dir)
  invisible(output_dir)
}

make_s1_panel_B_v4()
rm(make_s1_panel_B_v4)
