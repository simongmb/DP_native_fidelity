# ============================================================
# Supplementary Figure S2B | FINAL DATA SOURCE
# Gene-level rescue fraction vs native fidelity across 8 references
#
# DATA SOURCE:
#   results/DP_native_fidelity_reconstruction/
#   02_gene_level_native_fidelity_audit.csv
#
# Columns used:
#   GeneSymbol
#   Rescue_fraction
#   native_fidelity_direction_aligned
#   aligned_DP_vs_DC
#   aligned_DP_vs_CTS
#   aligned_DP_vs_HM1
#   aligned_DP_vs_HM2
#   aligned_DP_vs_IRS1
#   aligned_DP_vs_ORS1
#   aligned_DP_vs_ORS2
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

root <- getwd()
OUT_DIR <- file.path(root, "figures", "SUPPLEMENTARY FIGURES", "S2_Panel_B_gene_scatter_v5")
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
  "GeneSymbol","Rescue_fraction",
  "native_fidelity_direction_aligned",
  "aligned_DP_vs_DC","aligned_DP_vs_CTS",
  "aligned_DP_vs_HM1","aligned_DP_vs_HM2",
  "aligned_DP_vs_IRS1","aligned_DP_vs_ORS1","aligned_DP_vs_ORS2"
)
miss <- setdiff(required, names(dat))
if (length(miss)) stop("缺少列：", paste(miss, collapse=", "))

ref_cols <- c(
  "Composite native fidelity" = "native_fidelity_direction_aligned",
  "DP vs DC"   = "aligned_DP_vs_DC",
  "DP vs CTS"  = "aligned_DP_vs_CTS",
  "DP vs HM1"  = "aligned_DP_vs_HM1",
  "DP vs HM2"  = "aligned_DP_vs_HM2",
  "DP vs IRS1" = "aligned_DP_vs_IRS1",
  "DP vs ORS1" = "aligned_DP_vs_ORS1",
  "DP vs ORS2" = "aligned_DP_vs_ORS2"
)

long <- do.call(rbind, lapply(names(ref_cols), function(lbl) {
  data.frame(
    Gene = dat$GeneSymbol,
    RescueFraction = as.numeric(dat$Rescue_fraction),
    NativeFidelity = as.numeric(dat[[ref_cols[[lbl]]]]),
    Reference = lbl,
    stringsAsFactors = FALSE
  )
}))
long <- long[
  is.finite(long$RescueFraction) & is.finite(long$NativeFidelity),
  , drop=FALSE
]
long$Reference <- factor(long$Reference, levels=names(ref_cols))

stats <- do.call(rbind, lapply(split(long, long$Reference), function(z) {
  data.frame(
    Reference=as.character(unique(z$Reference)),
    n=nrow(z),
    rho=suppressWarnings(cor(z$RescueFraction, z$NativeFidelity,
                             method="spearman", use="complete.obs"))
  )
}))
stats$Reference <- factor(stats$Reference, levels=names(ref_cols))
stats$label <- sprintf("rho = %.3f", stats$rho)

xr <- range(long$RescueFraction, na.rm=TRUE)
stats$x <- xr[1] + 0.04*diff(xr)
stats$y <- vapply(seq_len(nrow(stats)), function(i) {
  z <- long$NativeFidelity[long$Reference == stats$Reference[i]]
  max(z,na.rm=TRUE) - 0.05*diff(range(z,na.rm=TRUE))
}, numeric(1))

INK <- "#273C4B"
POINT <- "#7892A2"
LINE <- "#3D8C82"

p <- ggplot(long, aes(RescueFraction, NativeFidelity)) +
  geom_hline(yintercept=0, linewidth=0.22, colour="#D2DADF") +
  geom_point(size=0.42, alpha=0.17, colour=POINT) +
  geom_smooth(method="lm", formula=y~x, se=FALSE,
              linewidth=0.31, colour=LINE) +
  geom_text(
    data=stats,
    aes(x=x,y=y,label=label),
    inherit.aes=FALSE,
    hjust=0,vjust=1,
    size=4.3/ggplot2::.pt,
    colour=INK
  ) +
  facet_wrap(~Reference, ncol=4, scales="fixed") +
  labs(
    title="Weak rescue-fidelity coupling across references",
    x="Rescue fraction",
    y="Native-fidelity score"
  ) +
  theme_classic(base_family="Helvetica", base_size=6) +
  theme(
    plot.title=element_text(size=6.2, face="bold", colour=INK,
                            hjust=0.5, margin=margin(b=3.0,"pt")),
    plot.title.position="plot",
    axis.title=element_text(size=5.4, colour=INK),
    axis.text=element_text(size=4.3, colour=INK),
    strip.background=element_rect(fill="#F1F4F5", colour=NA),
    strip.text=element_text(size=4.5, face="bold", colour=INK,
                            margin=margin(1.4,1.4,1.4,1.4,"pt")),
    panel.spacing=grid::unit(1.5,"mm"),
    axis.line=element_line(linewidth=0.25, colour="#82949D"),
    axis.ticks=element_line(linewidth=0.25, colour="#82949D"),
    plot.margin=margin(4,4,4,4,"pt")
  )

W_MM <- 103
H_MM <- 58

pdf_file <- file.path(OUT_DIR, "S2_Panel_B_gene_scatter_v5.pdf")
grDevices::pdf(pdf_file, width=W_MM/25.4, height=H_MM/25.4,
               family="Helvetica", useDingbats=FALSE, bg="white")
print(p); dev.off()

if (requireNamespace("ragg", quietly=TRUE)) {
  ragg::agg_png(
    file.path(OUT_DIR,"S2_Panel_B_gene_scatter_v5.png"),
    width=W_MM,height=H_MM,units="mm",res=600,background="white"
  )
  print(p); dev.off()

  ragg::agg_tiff(
    file.path(OUT_DIR,"S2_Panel_B_gene_scatter_v5.tiff"),
    width=W_MM,height=H_MM,units="mm",res=600,
    compression="lzw",background="white"
  )
  print(p); dev.off()
}

write.csv(long, file.path(OUT_DIR,"S2B_plot_data.csv"), row.names=FALSE)
write.csv(stats, file.path(OUT_DIR,"S2B_reference_stats.csv"), row.names=FALSE)

writeLines(c(
  "Supplementary Figure S2B v5",
  paste0("Input source: ", input_source),
  paste0("N source genes: ", nrow(dat)),
  "Rescue metric: Rescue_fraction",
  "Composite fidelity: native_fidelity_direction_aligned",
  paste0(names(ref_cols), " -> ", unname(ref_cols))
), file.path(OUT_DIR,"SOURCE_AUDIT.txt"))

message("S2B 完成：", OUT_DIR)
