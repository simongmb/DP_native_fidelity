# ============================================================
# Supplementary Figure S2C | FINAL DATA SOURCE
# Composite native fidelity across rescue-fraction quartiles
#
# DATA SOURCE:
#   results/DP_native_fidelity_reconstruction/
#   02_gene_level_native_fidelity_audit.csv
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

root <- getwd()
OUT_DIR <- file.path(root, "figures", "SUPPLEMENTARY FIGURES", "S2_Panel_C_rescue_quartiles_v5")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

input_file <- file.path(
  root, "results", "DP_native_fidelity_reconstruction",
  "02_gene_level_native_fidelity_audit.csv"
)

if (file.exists(input_file)) {
  dat <- read.csv(input_file, check.names=FALSE, stringsAsFactors=FALSE)
  input_source <- input_file
} else if (exists("audit", envir=.GlobalEnv, inherits=FALSE)) {
  dat <- get("audit", envir=.GlobalEnv)
  input_source <- "GlobalEnv object: audit"
} else {
  stop("找不到 Figure 2 gene-level native-fidelity audit 数据。")
}

required <- c(
  "GeneSymbol","Rescue_fraction","native_fidelity_direction_aligned"
)
miss <- setdiff(required, names(dat))
if (length(miss)) stop("缺少列：", paste(miss, collapse=", "))

d <- data.frame(
  Gene=dat$GeneSymbol,
  RescueFraction=as.numeric(dat$Rescue_fraction),
  CompositeFidelity=as.numeric(dat$native_fidelity_direction_aligned),
  stringsAsFactors=FALSE
)
d <- d[is.finite(d$RescueFraction) & is.finite(d$CompositeFidelity), , drop=FALSE]

# Empirical rescue quartiles; no arbitrary fidelity threshold.
r <- rank(d$RescueFraction, ties.method="average")
q <- pmin(4, pmax(1, ceiling(4*r/max(r))))
d$RescueQuartile <- factor(
  q, levels=1:4,
  labels=c("Q1\nlowest","Q2","Q3","Q4\nhighest")
)

rho <- suppressWarnings(cor(
  d$RescueFraction, d$CompositeFidelity,
  method="spearman", use="complete.obs"
))

sum_out <- do.call(rbind, lapply(split(d,d$RescueQuartile), function(z) {
  data.frame(
    RescueQuartile=as.character(unique(z$RescueQuartile)),
    N=nrow(z),
    MedianRescueFraction=median(z$RescueFraction),
    MedianCompositeFidelity=median(z$CompositeFidelity),
    Q1CompositeFidelity=unname(quantile(z$CompositeFidelity,.25)),
    Q3CompositeFidelity=unname(quantile(z$CompositeFidelity,.75))
  )
}))

INK <- "#273C4B"
FILL <- "#DDEBE8"
EDGE <- "#3D8C82"
POINT <- "#7892A2"

set.seed(20261030)

p <- ggplot(d, aes(RescueQuartile, CompositeFidelity)) +
  geom_hline(yintercept=0, linewidth=0.25, colour="#D2DADF") +
  geom_violin(width=0.72, trim=TRUE, fill=FILL, colour=NA) +
  geom_boxplot(
    width=0.26, outlier.shape=NA,
    fill="white", colour=EDGE, linewidth=0.42
  ) +
  geom_jitter(
    width=0.12, height=0,
    size=0.46, alpha=0.11, colour=POINT
  ) +
  labs(
    title="Stronger rescue does not imply higher native fidelity",
    subtitle=sprintf("Composite native fidelity | Spearman rho = %.3f", rho),
    x="Rescue-fraction quartile",
    y="Composite native-fidelity score"
  ) +
  theme_classic(base_family="Helvetica", base_size=6) +
  theme(
    plot.title=element_text(size=6.2,face="bold",colour=INK,
                            hjust=0.5,margin=margin(b=1.8,"pt")),
    plot.title.position="plot",
    plot.subtitle=element_text(size=4.8,colour="#627982",
                               hjust=0.5,margin=margin(b=4,"pt")),
    axis.title=element_text(size=5.6,colour=INK),
    axis.text=element_text(size=5.0,colour=INK),
    axis.line=element_line(linewidth=0.25,colour="#82949D"),
    axis.ticks=element_line(linewidth=0.25,colour="#82949D"),
    plot.margin=margin(4,5,4,5,"pt")
  )

W_MM <- 110
H_MM <- 58

pdf_file <- file.path(OUT_DIR, "S2_Panel_C_rescue_quartiles_v5.pdf")
grDevices::pdf(pdf_file, width=W_MM/25.4, height=H_MM/25.4,
               family="Helvetica", useDingbats=FALSE, bg="white")
print(p); dev.off()

if (requireNamespace("ragg", quietly=TRUE)) {
  ragg::agg_png(
    file.path(OUT_DIR,"S2_Panel_C_rescue_quartiles_v5.png"),
    width=W_MM,height=H_MM,units="mm",res=600,background="white"
  )
  print(p); dev.off()

  ragg::agg_tiff(
    file.path(OUT_DIR,"S2_Panel_C_rescue_quartiles_v5.tiff"),
    width=W_MM,height=H_MM,units="mm",res=600,
    compression="lzw",background="white"
  )
  print(p); dev.off()
}

write.csv(d, file.path(OUT_DIR,"S2C_plot_data.csv"), row.names=FALSE)
write.csv(sum_out, file.path(OUT_DIR,"S2C_quartile_summary.csv"), row.names=FALSE)

writeLines(c(
  "Supplementary Figure S2C v5",
  paste0("Input source: ", input_source),
  paste0("N genes with finite data: ", nrow(d)),
  "Rescue metric: Rescue_fraction",
  "Composite fidelity: native_fidelity_direction_aligned",
  sprintf("Overall Spearman rho = %.6f", rho)
), file.path(OUT_DIR,"SOURCE_AUDIT.txt"))

message("S2C 完成：", OUT_DIR)
