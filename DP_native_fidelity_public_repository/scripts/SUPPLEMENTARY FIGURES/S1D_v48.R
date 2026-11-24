# Supplementary Figure S1 | GSE44765 | standalone panel D v4.2 Nature compact with left-side legend
# This R script draws ONLY panel D (without a panel letter).
# Input files are read-only; output is written to a new dedicated version folder.
# Place this file in scripts/SUPPLEMENTARY FIGURES/ and run from project root:
# source("scripts/SUPPLEMENTARY FIGURES/SuppFig_S1D_v48.R")

make_s1_panel <- function() {
  root <- getwd()
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("缺少 ggplot2；请先运行 install.packages('ggplot2')")
  if (!requireNamespace("ragg", quietly = TRUE))
    stop("缺少 ragg；请先运行 install.packages('ragg')。macOS 的 PNG/TIFF 使用 ragg 导出。")
  output_dir <- file.path(root, "figures", "SUPPLEMENTARY FIGURES", "S1D_v48")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("输出目录已有内容，为防止覆盖请更改本脚本中的版本目录名：\n", output_dir)
  suppressPackageStartupMessages(library(ggplot2))

  # Input search: existing results/ folder preferred; optionally extract ZIP to temp.
  required <- c("GSE44765_gene_expression_RMA.rds",
                "GSE44765_design_matrix.csv",
                "GSE44765_global_rescue_summary.csv")
  basename_dir <- "GSE44765_rescue"
  paths <- c(file.path(root, "results", basename_dir),
             file.path(root, basename_dir),
             file.path(root, "scripts", "SUPPLEMENTARY FIGURES", basename_dir))
  complete <- function(d) dir.exists(d) && all(file.exists(file.path(d, required)))
  paths <- paths[vapply(paths, complete, logical(1))]
  if (length(paths)) {
    input_dir <- paths[1]
  } else {
    zip_paths <- c(file.path(root, "GSE44765_rescue.zip"),
                   file.path(root, "results", "GSE44765_rescue.zip"),
                   file.path(root, "scripts", "SUPPLEMENTARY FIGURES", "GSE44765_rescue.zip"))
    zip_paths <- zip_paths[file.exists(zip_paths)]
    if (!length(zip_paths)) stop("找不到原始 GSE44765_rescue 文件夹或 ZIP。请放在项目根目录。")
    tmp <- tempfile("s1_panel_D_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
    utils::unzip(zip_paths[1], exdir = tmp, files = paste0("GSE44765_rescue/", required))
    folders <- c(file.path(tmp, basename_dir), tmp)
    found <- folders[vapply(folders, complete, logical(1))]
    if (!length(found)) stop("ZIP 内缺少必要的输入文件：\n", paste(required, collapse = "\n"))
    input_dir <- found[1]
  }
  message("Panel D 输入文件夹：", input_dir)
  expr <- as.matrix(readRDS(file.path(input_dir, "GSE44765_gene_expression_RMA.rds")))
  storage.mode(expr) <- "numeric"
  design <- read.csv(file.path(input_dir, "GSE44765_design_matrix.csv"),
                     row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(expr) != 21355L || ncol(expr) != 18L || any(!is.finite(expr)) ||
      anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr)))
    stop("表达矩阵与预期 21,355 基因 × 18 样本不一致。")
  if (!setequal(rownames(design), colnames(expr))) stop("样本名与设计矩阵不一致。")
  design <- design[colnames(expr), , drop = FALSE]
  conditions <- c("Intact", "P0", "P1", "P3_2D", "P5_2D", "P3_3D")
  group_cols <- paste0("group", conditions)
  if (!all(c(group_cols, "donorD6", "donorD7") %in% colnames(design)))
    stop("设计矩阵缺少条件或供体列。")
  group_matrix <- as.matrix(design[, group_cols, drop = FALSE])
  if (any(rowSums(group_matrix) != 1L)) stop("每个样本必须属于一个条件。")
  meta <- data.frame(
    Sample = rownames(design),
    Condition = factor(conditions[max.col(group_matrix)], levels = conditions),
    Donor = factor(ifelse(design$donorD6 == 1, "D6",
                          ifelse(design$donorD7 == 1, "D7", "D5")),
                   levels = c("D5", "D6", "D7")),
    stringsAsFactors = FALSE)
  if (!all(table(meta$Condition, meta$Donor) == 1L))
    stop("六个条件 × 三名供体的样本矩阵不完整。")
  meta <- meta[order(meta$Condition, meta$Donor), , drop = FALSE]
  expr <- expr[, meta$Sample, drop = FALSE]
  meta$ConditionShort <- factor(as.character(meta$Condition),
    levels = conditions, labels = c("Intact", "P0", "P1", "P3 2D", "P5 2D", "P3 3D"))

  INK <- "#273C4B"; DOWN <- "#1C938A"; grey2d <- "#839AAB"
  theme_s1 <- theme_classic(base_family = "Helvetica", base_size = 6) +
    theme(
      plot.title = element_text(color = INK, face = "bold", size = 7,
                                margin = margin(b = 6, unit = "pt")),
      axis.title.x = element_text(color = INK, size = 6,
                                  margin = margin(t = 5, unit = "pt")),
      axis.text = element_text(color = INK, size = 6),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "#899BA5", linewidth = .25),
      axis.ticks.x = element_line(color = "#899BA5", linewidth = .25),
      panel.grid = element_blank(),
      plot.margin = margin(9, 14, 7, 5, unit = "pt")
    )

  get_sample <- function(donor, condition) {
    z <- meta$Sample[meta$Donor == donor & meta$Condition == condition]
    if (length(z) != 1L) stop("供体和条件无法唯一匹配。")
    z
  }
  rows <- lapply(levels(meta$Donor), function(donor) {
    intact <- expr[, get_sample(donor,"Intact")]
    two_d <- expr[, get_sample(donor,"P3_2D")]
    three_d <- expr[, get_sample(donor,"P3_3D")]
    data.frame(Donor=donor,N_Genes=length(intact),
               N_CloserIn3D=sum(abs(three_d-intact)<abs(two_d-intact)),
               Percent_CloserIn3D=100*mean(abs(three_d-intact)<abs(two_d-intact)),
               Median_Abs_Distance_P3_2D=median(abs(two_d-intact)),
               Median_Abs_Distance_P3_3D=median(abs(three_d-intact)))
  })
  donor_df <- do.call(rbind,rows)
  donor_df$Donor <- factor(donor_df$Donor,levels=c("D7","D6","D5"))
  donor_df$PercentLabel <- sprintf("%.1f%%",donor_df$Percent_CloserIn3D)
  long <- rbind(
    data.frame(Donor=donor_df$Donor,Condition="P3 2D",
               MedianDistance=donor_df$Median_Abs_Distance_P3_2D),
    data.frame(Donor=donor_df$Donor,Condition="P3 3D",
               MedianDistance=donor_df$Median_Abs_Distance_P3_3D))
  long$Condition <- factor(long$Condition,levels=c("P3 2D","P3 3D"))

  # Independent reconstruction of the group-mean Figure 1F values.
  means <- sapply(c("Intact","P3_2D","P3_3D"),function(cnd)
    rowMeans(expr[,meta$Sample[meta$Condition==cnd],drop=FALSE]))
  global_recalc <- c(
    Fraction_genes_3D_closer_to_Intact = mean(abs(means[,"P3_3D"]-means[,"Intact"]) <
                                              abs(means[,"P3_2D"]-means[,"Intact"])),
    Median_absolute_distance_P3_2D_to_Intact = median(abs(means[,"P3_2D"]-means[,"Intact"])),
    Median_absolute_distance_P3_3D_to_Intact = median(abs(means[,"P3_3D"]-means[,"Intact"])))
  saved_df <- read.csv(file.path(input_dir,"GSE44765_global_rescue_summary.csv"),
                       check.names=FALSE)
  saved <- setNames(as.numeric(saved_df$Value), saved_df$Metric)
  if (!all(names(global_recalc) %in% names(saved)) ||
      max(abs(global_recalc-saved[names(global_recalc)])) > 1e-7)
    stop("Figure 1F 的汇总指标复算未通过。")

  # Percentage labels denote the fraction of genes closer in 3D,
  # NOT a percentage reduction in the median distance.
  p <- ggplot(long, aes(x = MedianDistance, y = Donor)) +
    geom_segment(data = donor_df, inherit.aes = FALSE,
                 aes(x = Median_Abs_Distance_P3_2D,
                     xend = Median_Abs_Distance_P3_3D, y = Donor, yend = Donor),
                 linewidth = 0.55, colour = "#BDC9D1") +
    geom_point(aes(colour = Condition, shape = Condition),
               size = 2.55, stroke = 0.55) +
    geom_text(data = donor_df, inherit.aes = FALSE,
              aes(x = 0.420, y = Donor, label = PercentLabel),
              hjust = 0, size = 6 / ggplot2::.pt, colour = INK) +
    # Top annotation row: symbols left, percentage-column heading right.
    annotate("point", x = 0.345, y = 3.30, shape = 16,
             size = 1.80, colour = grey2d) +
    annotate("text", x = 0.350, y = 3.30, label = "P3 2D",
             hjust = 0, vjust = 0.5,
             size = 5.5 / ggplot2::.pt, colour = INK) +
    annotate("point", x = 0.378, y = 3.30, shape = 17,
             size = 1.95, colour = DOWN) +
    annotate("text", x = 0.383, y = 3.30, label = "P3 3D",
             hjust = 0, vjust = 0.5,
             size = 5.5 / ggplot2::.pt, colour = INK) +
    annotate("text", x = 0.457, y = 3.30,
             label = "Genes closer in 3D (%)", hjust = 1, vjust = 0.5,
             size = 5.5 / ggplot2::.pt, colour = INK) +
    scale_colour_manual(values = c("P3 2D" = grey2d, "P3 3D" = DOWN),
                        guide = "none") +
    scale_shape_manual(values = c("P3 2D" = 16, "P3 3D" = 17),
                       guide = "none") +
    scale_x_continuous(limits = c(.34, .465),
                       breaks = c(.34, .38, .42, .46),
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    scale_y_discrete(expand = expansion(add = c(.36, .55))) +
    labs(title = "Transcriptome-wide distance from intact DP by donor",
         x = "Median absolute log2-expression distance from intact DP",
         y = NULL) +
    theme_s1

  # Independent outputs. PDF remains editable in Illustrator.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  basename_output <- "S1D_v48"
  emit <- function(ext) {
    file <- file.path(output_dir, paste0(basename_output, ".", ext))
    w <- 88 / 25.4; h <- 70 / 25.4; res <- 600
    if (ext == "pdf") {
      grDevices::pdf(file, width = w, height = h, useDingbats = FALSE, bg = "white")
    } else if (ext == "png") {
      ragg::agg_png(file, width = w, height = h, units = "in", res = res,
                    background = "white")
    } else if (ext == "tiff") {
      ragg::agg_tiff(file, width = w, height = h, units = "in", res = res,
                     compression = "lzw", background = "white")
    }
    tryCatch(print(p, newpage = FALSE), finally = grDevices::dev.off())
    if (!file.exists(file) || file.info(file)$size <= 0) stop("绘图设备未成功生成文件：", file)
    message("成功输出：", file)
  }
  for (ext in c("pdf", "png", "tiff")) emit(ext)
  write.csv(meta, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
  write.csv(donor_df,file.path(output_dir,"Panel_D_global_proximity_by_donor.csv"),row.names=FALSE)
  write.csv(data.frame(Metric=names(global_recalc),Value=unname(global_recalc)),
            file.path(output_dir,"Panel_D_Figure1F_recalculation.csv"),row.names=FALSE)
  writeLines(c(
    "Supplementary Figure S1 panel D | independent output",
    paste0("Input folder: ", normalizePath(input_dir, winslash = "/", mustWork = TRUE)),
    "Gene-level RMA data; 21,355 genes, 18 samples, 3 matched donors, 6 conditions.",
    "Panel D: all genes; per-donor paired intact/2D/3D comparisons are descriptive. Figure 1F metrics recalculated separately from group mean profiles and checked against original CSV.",
    "X axis = median across 21,355 genes of absolute log2-expression difference (condition minus the SAME donor's intact DP).",
    "Circles and triangles show TWO SEPARATE distances to intact DP. Their positions do NOT represent a direct P3 2D versus P3 3D distance.",
    "Lower distance is closer to matched intact DP. These are gene-level medians, NOT PCA/Euclidean transcriptome distances.",
    "Percent closer is the percentage of all genes whose absolute within-donor 3D-intact difference is smaller than their 2D-intact difference; ties are not counted as closer.",
    "Discovery dataset, only 3 matched donors. Descriptive internal assessment: no consistent global restoration advantage has been established.",
    "Figure 1F pooled-group results use donor-averaged profiles; S1D per-donor values use paired samples. These are different summaries.",
     "v4.6 clipping fix: 88 x 70 mm; right-side heading is right-aligned to a safe in-panel anchor so the full label extends leftward and cannot be clipped; data layout unchanged.",
    "Top annotation row is pulled closer to the data; percentage column is shifted left to reduce dead space.",
    "Percentages denote fraction of genes closer in 3D, not median-distance reductions.",
    "The plotted distance is a median of absolute log2-expression differences, not a Euclidean distance.",
    "No panel letter; add 8 pt bold d manually in Illustrator.",
    "All existing source files were read-only; old outputs were not overwritten.",
    paste0("R: ", R.version.string),
    paste0("ggplot2: ", as.character(utils::packageVersion("ggplot2"))),
    paste0("ragg: ", as.character(utils::packageVersion("ragg")))
  ), file.path(output_dir, "SOURCE_AUDIT.txt"))
  message("Panel D 完成，输出目录：", output_dir)
  invisible(output_dir)
}

make_s1_panel()
rm(make_s1_panel)
