# Supplementary Figure S1 | GSE44765 | standalone panel C | v4.1 Nature-style compact-width layout
# This R script draws ONLY panel C (without a panel letter).
# Input files are read-only; output is written to a new dedicated version folder.
# Place this file in scripts/SUPPLEMENTARY FIGURES/ and run from project root:
# source("scripts/SUPPLEMENTARY FIGURES/SuppFig_S1C_v43.R")

make_s1_panel <- function() {
  root <- getwd()
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("缺少 ggplot2；请先运行 install.packages('ggplot2')")
  if (!requireNamespace("ragg", quietly = TRUE))
    stop("缺少 ragg；请先运行 install.packages('ragg')。macOS 的 PNG/TIFF 使用 ragg 导出。")
  output_dir <- file.path(root, "figures", "SUPPLEMENTARY FIGURES", "S1C_v43")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("输出目录已有内容，为防止覆盖请更改本脚本中的版本目录名：\n", output_dir)
  suppressPackageStartupMessages(library(ggplot2))

  # Input search: existing results/ folder preferred; optionally extract ZIP to temp.
  required <- c("GSE44765_gene_expression_RMA.rds",
                "GSE44765_design_matrix.csv",
                "GSE44765_STRICT_rescue_genes.csv")
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
    tmp <- tempfile("s1_panel_C_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
    utils::unzip(zip_paths[1], exdir = tmp, files = paste0("GSE44765_rescue/", required))
    folders <- c(file.path(tmp, basename_dir), tmp)
    found <- folders[vapply(folders, complete, logical(1))]
    if (!length(found)) stop("ZIP 内缺少必要的输入文件：\n", paste(required, collapse = "\n"))
    input_dir <- found[1]
  }
  message("Panel C 输入文件夹：", input_dir)
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

  # Publication typography at native PDF size (88 x 70 mm).
  # Panel letter c is intentionally added separately in Illustrator (8 pt bold).
  INK <- "#273C4B"
  UP <- "#C58D6E"
  DOWN <- "#1C938A"
  theme_s1 <- theme_classic(base_family = "Helvetica", base_size = 6) +
    theme(
      plot.title = element_text(color = INK, face = "bold", size = 7,
                                hjust = 0, margin = margin(b = 6, unit = "pt")),
      axis.title.x = element_text(color = INK, size = 6.5,
                                  margin = margin(t = 5, unit = "pt")),
      axis.text.x = element_text(color = INK, size = 6),
      axis.text.y = element_text(color = INK, size = 6),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "#899BA5", linewidth = 0.25),
      axis.ticks.x = element_line(color = "#899BA5", linewidth = 0.25),
      panel.grid = element_blank(),
      plot.margin = margin(7, 7, 7, 5, unit = "pt")
    )

  strict <- read.csv(file.path(input_dir,"GSE44765_STRICT_rescue_genes.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(nrow(strict) == 1571L, !anyDuplicated(strict$GeneSymbol),
            sum(strict$Rescue_type == "Restored_upward") == 510L,
            sum(strict$Rescue_type == "Restored_downward") == 1061L,
            all(strict$GeneSymbol %in% rownames(expr)))
  get_sample <- function(donor, condition) {
    z <- meta$Sample[meta$Donor == donor & meta$Condition == condition]
    if (length(z)!=1L) stop("供体和条件无法唯一匹配。")
    z
  }
  rows <- list()
  for (donor in levels(meta$Donor)) {
    intact <- expr[,get_sample(donor,"Intact")]
    two_d <- expr[,get_sample(donor,"P3_2D")]
    three_d <- expr[,get_sample(donor,"P3_3D")]
    proximity <- abs(three_d-intact) < abs(two_d-intact)
    direction <- sign(intact-two_d) == sign(three_d-two_d)
    for (arm in c("Restored_upward","Restored_downward")) {
      genes <- strict$GeneSymbol[strict$Rescue_type == arm]
      rows[[length(rows)+1L]] <- data.frame(
        Donor = donor, Arm = arm, N_Predefined = length(genes),
        N_CloserIn3D = sum(proximity[genes]),
        Percent_CloserIn3D = 100 * mean(proximity[genes]),
        N_DirectionMatch = sum(direction[genes]),
        Percent_DirectionMatch = 100 * mean(direction[genes]),
        N_Both = sum(proximity[genes] & direction[genes]),
        Percent_Both = 100 * mean(proximity[genes] & direction[genes]))
    }
  }
  donor_df <- do.call(rbind, rows)
  donor_df$Arm <- factor(donor_df$Arm,
                         levels=c("Restored_upward","Restored_downward"),
                         labels=c("Restored upward","Restored downward"))
  donor_df$Row <- paste(donor_df$Donor,
                         ifelse(donor_df$Arm=="Restored upward","up","down"))
  donor_df$Row <- factor(donor_df$Row,
    levels=rev(c("D5 up","D5 down","D6 up","D6 down","D7 up","D7 down")))
  # Strict rescue gene percentages are descriptive within the discovery cohort.
  # No recomputation or filtering is introduced in v4.
  p <- ggplot(donor_df, aes(x = Percent_CloserIn3D, y = Row, colour = Arm)) +
    geom_segment(aes(x = 0, xend = Percent_CloserIn3D, yend = Row),
                 linewidth = 0.78, alpha = 0.34, show.legend = FALSE) +
    geom_point(size = 2.35, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%", Percent_CloserIn3D)),
              hjust = 0, nudge_x = 1.4, size = 6 / ggplot2::.pt,
              colour = INK, show.legend = FALSE) +
    scale_colour_manual(values = c("Restored upward" = UP,
                                   "Restored downward" = DOWN)) +
    scale_y_discrete(labels = function(x)
      sub(" ", " | ", gsub("up$", "upward", gsub("down$", "downward", x)))) +
    scale_x_continuous(limits = c(0, 106), breaks = c(0, 25, 50, 75, 100),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = "Within-donor proximity for strict rescue genes",
         x = "Strict-rescue genes closer to intact DP in 3D (%)", y = NULL) +
    theme_s1

  # Independent outputs. PDF remains editable in Illustrator.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  basename_output <- "S1C_v43"
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
  write.csv(donor_df, file.path(output_dir, "Panel_C_donor_strict_genes_descriptive.csv"), row.names = FALSE)
  writeLines(c(
    "Supplementary Figure S1 panel C | independent output",
    paste0("Input folder: ", normalizePath(input_dir, winslash = "/", mustWork = TRUE)),
    "Gene-level RMA data; 21,355 genes, 18 samples, 3 matched donors, 6 conditions.",
    "Panel C: donor-wise |3D - intact| < |2D - intact| for strictly selected genes.",
    "The 1,571 strict rescue genes were selected using all three donors (D5/D6/D7); percentages are in-sample descriptive assessments and are NOT independent validation.",
    "v4.2 Nature balanced retains 88 x 70 mm; connecting bars are lighter/thinner, points are slightly stronger, the x-range is tightened to reduce dead space, and percentage labels are pulled closer; gene data, percentages, sample order and palette unchanged.",
    "Native PDF size: 88 x 70 mm. Do not scale without checking final text sizes.",
    "No panel letter; add 8 pt bold c manually in Illustrator.",
    "All existing source files were read-only; old outputs were not overwritten.",
    paste0("R: ", R.version.string),
    paste0("ggplot2: ", as.character(utils::packageVersion("ggplot2"))),
    paste0("ragg: ", as.character(utils::packageVersion("ragg")))
  ), file.path(output_dir, "SOURCE_AUDIT.txt"))
  message("Panel C 完成，输出目录：", output_dir)
  invisible(output_dir)
}

make_s1_panel()
rm(make_s1_panel)
