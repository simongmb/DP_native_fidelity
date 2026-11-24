# Supplementary Figure S1a | GSE44765 | Nature-style refinement (v4.2 TITLE 7PT)
# -----------------------------------------------------------------------------
# Changes from v3 are PRESENTATIONAL only; v4.2 additionally fixes a numeric QC check:
#   1) Shorter title: RMA-normalized expression distributions
#   2) Removes the lengthy in-panel subtitle (details belong in figure legend)
#   3) Removes background gridlines; retains discrete group headers
#   4) Sets text sizes and physical output size for placement in an ~183-mm
#      two-column composite figure (panel A target width = 110 mm)
# IMPORTANT: NO data/quantile/boxplot/color/condition definitions were changed.
# Input files are READ ONLY. v3 and v4 outputs are NEVER overwritten.
# Figure panel letter "a" should be added in Illustrator, not here.
#
# Place this script in: scripts/SUPPLEMENTARY FIGURES/
# Run from the HairRegeneration_AI project root in RStudio:
# source("scripts/SUPPLEMENTARY FIGURES/SuppFig_S1_Panel_A_v4_2_Nature_7pt_title.R")
#
# Expected input at one of:
#   results/GSE44765_rescue/{GSE44765_gene_expression_RMA.rds,
#                           GSE44765_design_matrix.csv}
# OR GSE44765_rescue.zip in the project root / results folder.
# New outputs at:
#   figures/SUPPLEMENTARY FIGURES/S1_Panel_A_v4_2_Nature_7pt_title/
# -----------------------------------------------------------------------------

make_s1_panel_a_v4_2 <- function() {
  root <- getwd()

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("缺少 ggplot2：请先运行 install.packages('ggplot2')")
  }
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("缺少 ragg：请先运行 install.packages('ragg')")
  }
  suppressPackageStartupMessages(library(ggplot2))

  output_dir <- file.path(
    root, "figures", "SUPPLEMENTARY FIGURES", "S1_Panel_A_v4_2_Nature_7pt_title"
  )
  if (dir.exists(output_dir) &&
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("为保护已有结果，v4.2 输出文件夹非空，程序停止：\n", output_dir,
         "\n如果需要再次运行，请先将已有文件夹重命名备份。")
  }

  # --- 1. Locate read-only source files --------------------------------------
  required <- c("GSE44765_gene_expression_RMA.rds",
                "GSE44765_design_matrix.csv")
  input_name <- "GSE44765_rescue"
  candidate_dirs <- c(
    file.path(root, "results", input_name),
    file.path(root, input_name),
    file.path(root, "scripts", "SUPPLEMENTARY FIGURES", input_name)
  )
  has_inputs <- function(d) {
    dir.exists(d) && all(file.exists(file.path(d, required)))
  }
  candidate_dirs <- candidate_dirs[
    vapply(candidate_dirs, has_inputs, logical(1))
  ]

  if (length(candidate_dirs)) {
    input_dir <- candidate_dirs[1]
  } else {
    zip_candidates <- c(
      file.path(root, "GSE44765_rescue.zip"),
      file.path(root, "results", "GSE44765_rescue.zip"),
      file.path(root, "scripts", "SUPPLEMENTARY FIGURES", "GSE44765_rescue.zip")
    )
    zip_candidates <- zip_candidates[file.exists(zip_candidates)]
    if (!length(zip_candidates)) {
      stop("未找到 GSE44765_rescue 文件夹或 ZIP。请确认 RStudio 的工作目录是项目根目录：\n",
           root)
    }
    tmp <- tempfile("s1a_v4_input_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
    utils::unzip(zip_candidates[1], exdir = tmp,
                 files = paste0(input_name, "/", required))
    unzipped_dirs <- c(file.path(tmp, input_name), tmp)
    found <- unzipped_dirs[vapply(unzipped_dirs, has_inputs, logical(1))]
    if (!length(found)) {
      stop("ZIP 中未找到完整输入文件：", paste(required, collapse = ", "))
    }
    input_dir <- found[1]
  }
  message("S1a v4.2 输入位置：", input_dir)

  # --- 2. Read and validate raw expression/design matrices ------------------
  expr <- as.matrix(readRDS(file.path(input_dir, required[1])))
  storage.mode(expr) <- "numeric"
  design <- utils::read.csv(file.path(input_dir, required[2]),
                            row.names = 1, check.names = FALSE,
                            stringsAsFactors = FALSE)

  if (nrow(expr) != 21355L || ncol(expr) != 18L ||
      any(!is.finite(expr)) ||
      anyDuplicated(rownames(expr)) ||
      anyDuplicated(colnames(expr))) {
    stop("表达矩阵检查未通过：应为 21,355 基因 × 18 样本，无缺失或重复 ID。")
  }
  if (!setequal(rownames(design), colnames(expr))) {
    stop("设计矩阵的样本 ID 与表达矩阵的列名不匹配。")
  }
  design <- design[colnames(expr), , drop = FALSE]

  conditions <- c("Intact", "P0", "P1", "P3_2D", "P5_2D", "P3_3D")
  group_cols <- paste0("group", conditions)
  if (!all(c(group_cols, "donorD6", "donorD7") %in% names(design))) {
    stop("设计矩阵缺少 condition 或 donor 字段。")
  }
  group_matrix <- as.matrix(design[, group_cols, drop = FALSE])
  if (anyNA(group_matrix) || any(rowSums(group_matrix) != 1L)) {
    stop("设计矩阵异常：每个样本必须且只能属于一个实验条件。")
  }
  donor_D6 <- design$donorD6 == 1
  donor_D7 <- design$donorD7 == 1
  if (anyNA(donor_D6) || anyNA(donor_D7) || any(donor_D6 & donor_D7)) {
    stop("设计矩阵的 donorD6 / donorD7 标注存在错误。")
  }

  meta <- data.frame(
    Sample = rownames(design),
    Condition = factor(conditions[max.col(group_matrix)], levels = conditions),
    Donor = factor(ifelse(donor_D6, "D6", ifelse(donor_D7, "D7", "D5")),
                   levels = c("D5", "D6", "D7")),
    stringsAsFactors = FALSE
  )
  if (!all(table(meta$Condition, meta$Donor) == 1L)) {
    stop("六个条件 × 三位供体应包含 18 个唯一组合。")
  }
  meta <- meta[order(meta$Condition, meta$Donor), , drop = FALSE]
  expr <- expr[, meta$Sample, drop = FALSE]
  meta$ConditionShort <- factor(
    as.character(meta$Condition), levels = conditions,
    labels = c("Intact", "P0", "P1", "P3 2D", "P5 2D", "P3 3D")
  )

  # --- 3. EXACT v3 quantile definitions, unchanged --------------------------
  quantile_df <- do.call(rbind, lapply(seq_len(nrow(meta)), function(j) {
    q <- stats::quantile(expr[, j], probs = c(.05, .25, .50, .75, .95),
                         names = FALSE)
    data.frame(
      Sample = meta$Sample[j],
      Condition = meta$ConditionShort[j],
      Donor = meta$Donor[j],
      q05 = q[1], lower = q[2], middle = q[3],
      upper = q[4], q95 = q[5]
    )
  }))
  # Data frames are lists internally: is.finite(data.frame) raises
  # "default method not implemented for type 'list'" in R.
  # Check numeric columns explicitly, then convert to a numeric matrix.
  quantile_cols <- c("q05", "lower", "middle", "upper", "q95")
  if (nrow(quantile_df) != 18L ||
      !all(vapply(quantile_df[, quantile_cols, drop = FALSE],
                  is.numeric, logical(1)))) {
    stop("18 个样本的分位数计算失败：数量不符或存在非数值列。")
  }
  quantile_matrix <- as.matrix(quantile_df[, quantile_cols, drop = FALSE])
  if (any(!is.finite(quantile_matrix)) ||
      any(quantile_matrix[, "q05"] > quantile_matrix[, "lower"]) ||
      any(quantile_matrix[, "lower"] > quantile_matrix[, "middle"]) ||
      any(quantile_matrix[, "middle"] > quantile_matrix[, "upper"]) ||
      any(quantile_matrix[, "upper"] > quantile_matrix[, "q95"])) {
    stop("18 个样本的分位数检查失败：存在无效值或分位数顺序异常。")
  }

  # --- 4. Nature-style presentation: v3 palette & data geometry retained ----
  INK <- "#273C4B"
  cond_palette <- c(
    "Intact" = "#6D879B",
    "P0" = "#9BA9B6",
    "P1" = "#C3A28C",
    "P3 2D" = "#51A79E",
    "P5 2D" = "#7C9B94",
    "P3 3D" = "#177F79"
  )

  p <- ggplot(
    quantile_df,
    aes(x = Donor, ymin = q05, lower = lower, middle = middle,
        upper = upper, ymax = q95, fill = Condition)
  ) +
    geom_boxplot(stat = "identity", width = 0.64,
                 colour = INK, linewidth = 0.32) +
    scale_fill_manual(values = cond_palette, guide = "none") +
    facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
    scale_y_continuous(breaks = c(4, 6, 8, 10),
                       expand = expansion(mult = c(0.02, 0.04))) +
    labs(
      title = "RMA-normalized expression distributions",
      x = NULL,
      y = "Normalized log2 expression"
    ) +
    theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(
      plot.title = element_text(color = INK, face = "bold", size = 7.0,
                                hjust = 0, margin = margin(b = 5)),
      plot.subtitle = element_blank(),
      axis.title.y = element_text(color = INK, size = 6.5,
                                  margin = margin(r = 5)),
      axis.title.x = element_blank(),
      axis.text.x = element_text(color = INK, size = 6.0,
                                 margin = margin(t = 3)),
      axis.text.y = element_text(color = INK, size = 6.5,
                                 margin = margin(r = 2)),
      axis.line = element_line(color = INK, linewidth = 0.27),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = INK, linewidth = 0.25),
      axis.ticks.length = grid::unit(1.4, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#F1F5F7", colour = NA),
      strip.text.x = element_text(color = INK, face = "bold", size = 6.5,
                                  margin = margin(t = 3, b = 3)),
      panel.spacing.x = grid::unit(1.5, "mm"),
      plot.margin = margin(t = 6, r = 6, b = 5, l = 4)
    )

  # Title is 7 pt at the exported 110-mm panel width. Avoid scaling after placement.
  # --- 5. Export independent figure files; do not touch v3 -------------------
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  base <- "S1_Panel_A_v4_2_Nature_7pt_title"
  width_mm <- 110  # Intended panel width within a two-column composite figure
  height_mm <- 60
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  emit <- function(ext) {
    file <- file.path(output_dir, paste0(base, ".", ext))
    if (ext == "pdf") {
      grDevices::pdf(file, width = width_in, height = height_in,
                     family = "Helvetica", useDingbats = FALSE, bg = "white")
    } else if (ext == "png") {
      ragg::agg_png(file, width = width_mm, height = height_mm,
                    units = "mm", res = 600, background = "white")
    } else if (ext == "tiff") {
      ragg::agg_tiff(file, width = width_mm, height = height_mm,
                     units = "mm", res = 600, compression = "lzw",
                     background = "white")
    } else {
      stop("不支持的格式：", ext)
    }
    tryCatch(print(p), finally = grDevices::dev.off())
    if (!file.exists(file) || file.info(file)$size <= 0) {
      stop("文件导出失败：", file)
    }
    message("已导出：", file)
  }
  for (ext in c("pdf", "png", "tiff")) emit(ext)

  utils::write.csv(meta, file.path(output_dir, "sample_metadata.csv"),
                   row.names = FALSE)
  utils::write.csv(quantile_df,
                   file.path(output_dir, "Panel_A_expression_quantiles.csv"),
                   row.names = FALSE)
  writeLines(c(
    "Supplementary Figure S1 panel A | v4.2 Nature-style presentation revision",
    paste0("Input directory: ", normalizePath(input_dir, winslash = "/")),
    "21,355 genes; 18 samples; six conditions; three donors D5, D6, D7.",
    "Same data, quantile definitions, group order, and color palette as v3.",
    "Quantiles: 5th, 25th, 50th, 75th, 95th percentiles across genes per sample.",
    "v4.2 adjustments: title 7 pt; strip 6.5 pt; y-axis title 6.5 pt; x tick labels 6 pt; grid removed.",
    paste0("Panel size at export: ", width_mm, " x ", height_mm, " mm."),
    "PDF is vector-based; import/place PDF in Illustrator. Add panel letter a manually.",
    "The RMA distribution check is descriptive, not proof of absence of sample artifacts.",
    "Original data and older figure exports were not modified.",
    paste0("R: ", R.version.string),
    paste0("ggplot2: ", as.character(utils::packageVersion("ggplot2"))),
    paste0("ragg: ", as.character(utils::packageVersion("ragg")))
  ), file.path(output_dir, "SOURCE_AUDIT.txt"))

  message("\nS1a v4.2 完成。将 PDF 置入 Illustrator：\n",
          file.path(output_dir, paste0(base, ".pdf")))
  invisible(output_dir)
}

make_s1_panel_a_v4_2()
rm(make_s1_panel_a_v4_2)
