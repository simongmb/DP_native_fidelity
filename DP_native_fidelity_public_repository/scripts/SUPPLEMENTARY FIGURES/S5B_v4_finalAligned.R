# ============================================================
# Supplementary Figure S5B v4 — DONOR-BLOCKED RESULTS
# Repeated donor-blocked nested-CV performance across 20 repeats
#
# Input:
#   Fig5C_repeat_metrics.csv from the finalized Figure 5C state-balance reconstruction
#
# This script DOES NOT rerun any model.
# It only visualizes already-computed donor-blocked Fig5C results.
# ============================================================

required_pkgs <- c(
  "ggplot2",
  "readr",
  "dplyr",
  "tidyr"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "\nMissing packages:\n",
      paste(missing_pkgs, collapse = ", "),
      "\n\nInstall with:\ninstall.packages(c(",
      paste0('"', missing_pkgs, '"', collapse = ", "),
      "))\n"
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(
      dirname(
        normalizePath(
          sub("^--file=", "", file_arg[1])
        )
      )
    )
  }

  getwd()
}

start_dir <- get_script_dir()

candidate_roots <- unique(
  c(
    start_dir,
    getwd(),
    dirname(start_dir),
    dirname(dirname(start_dir)),
    dirname(dirname(dirname(start_dir))),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )
)

candidate_roots <- candidate_roots[
  dir.exists(candidate_roots)
]

input_candidates <- character(0)

for (root in candidate_roots) {
  input_candidates <- c(
    input_candidates,
    list.files(
      root,
      pattern = "^Fig5C_repeat_metrics\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
  )
}

input_candidates <- unique(
  input_candidates
)

if (length(input_candidates) == 0) {
  stop(
    "\nCannot find Fig5C_repeat_metrics.csv.\n"
  )
}

score_input <- function(x) {
  score <- 0

  if (
    grepl(
      "stateBalanceReconstruction_v12|Panel_C_stateBalanceReconstruction_v12",
      x,
      ignore.case = TRUE
    )
  ) {
    score <- score + 200
  }

  if (grepl("donor_blocked_v8", x, ignore.case = TRUE)) {
    score <- score + 100
  }

  if (grepl("Panel_C_donor", x, ignore.case = TRUE)) {
    score <- score + 50
  }

  if (grepl("FIGURE[ _]?5", x, ignore.case = TRUE)) {
    score <- score + 10
  }

  score
}

scores <- vapply(
  input_candidates,
  score_input,
  numeric(1)
)

input_file <- input_candidates[
  which.max(scores)
]

cat(
  "\nInput file:\n",
  input_file,
  "\n\n"
)

dat <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

required_cols <- c(
  "Pearson_r",
  "Spearman_rho",
  "CV_R2",
  "MAE",
  "RMSE",
  "Repeat"
)

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "\nQC STOP: required columns missing:\n",
      paste(missing_cols, collapse = ", "),
      "\n"
    )
  )
}

if (nrow(dat) != 20) {
  warning(
    paste0(
      "Expected 20 repeated-CV rows, but found ",
      nrow(dat),
      "."
    )
  )
}

if (anyDuplicated(dat$Repeat) > 0) {
  stop(
    "\nQC STOP: duplicate Repeat values found.\n"
  )
}

plot_dat <- dat %>%
  dplyr::select(
    Repeat,
    Pearson_r,
    Spearman_rho,
    CV_R2,
    MAE
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      Pearson_r,
      Spearman_rho,
      CV_R2,
      MAE
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      Pearson_r = "Pearson r",
      Spearman_rho = "Spearman rho",
      CV_R2 = "CV R2",
      MAE = "MAE"
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Pearson r",
        "Spearman rho",
        "CV R2",
        "MAE"
      )
    )
  )

summary_dat <- plot_dat %>%
  dplyr::group_by(
    Metric
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(Value),
    median = median(Value),
    sd = sd(Value),
    min = min(Value),
    q1 = quantile(Value, 0.25),
    q3 = quantile(Value, 0.75),
    max = max(Value),
    .groups = "drop"
  )

print(
  summary_dat
)

p <- ggplot(
  plot_dat,
  aes(
    x = 1,
    y = Value
  )
) +
  geom_boxplot(
    width = 0.24,
    outlier.shape = NA,
    linewidth = 0.32,
    fill = "white"
  ) +
  geom_jitter(
    width = 0.055,
    height = 0,
    size = 1.45,
    alpha = 0.80
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 2.0
  ) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    title = "Donor-blocked nested-CV performance across 20 repeats",
    subtitle = "Each point represents one repeat comprising 10 leave-one-donor-out outer folds",
    x = NULL,
    y = NULL
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(
      mult = c(0.15, 0.15)
    )
  ) +
  theme_classic(
    base_size = 6
  ) +
  theme(
    plot.title = element_text(
      size = 6.2,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 4.8,
      hjust = 0.5,
      margin = margin(
        b = 4
      )
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 5.2,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 4.7
    ),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = grid::unit(
      0.32,
      "lines"
    ),
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )

project_root <- dirname(
  dirname(
    dirname(
      input_file
    )
  )
)

if (
  !dir.exists(
    file.path(
      project_root,
      "figures"
    )
  )
) {
  project_root <- start_dir
}

out_dir <- file.path(
  project_root,
  "figures",
  "SUPPLEMENTARY FIGURES",
  "S5_Panel_B_repeatedCV_donorBlocked_v4"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(
    out_dir,
    "S5_Panel_B_repeatedCV_donorBlocked_v4.pdf"
  ),
  plot = p,
  width = 100,
  height = 62,
  units = "mm",
  device = "pdf"
)

ggsave(
  file.path(
    out_dir,
    "S5_Panel_B_repeatedCV_donorBlocked_v4.png"
  ),
  plot = p,
  width = 100,
  height = 62,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "S5_Panel_B_repeatedCV_donorBlocked_v4.tiff"
  ),
  plot = p,
  width = 100,
  height = 62,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

readr::write_csv(
  dat,
  file.path(
    out_dir,
    "S5B_repeat_metrics.csv"
  )
)

readr::write_csv(
  summary_dat,
  file.path(
    out_dir,
    "S5B_repeat_metrics_summary.csv"
  )
)

readme <- c(
  "Supplementary Figure S5B v4 — donor-blocked",
  "",
  "Input:",
  input_file,
  "",
  "This script does not rerun models.",
  "It visualizes the 20 repeat-level metrics generated by",
  "the finalized Figure 5C donor-blocked state-balance reconstruction.",
  "",
  "Each plotted point is one complete repeat, not one outer donor fold.",
  "Each repeat comprises 10 leave-one-donor-out outer folds with donor-grouped inner tuning."
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_S5B_donorBlocked_v4.txt"
  )
)

cat(
  "\nS5B donor-blocked plot completed.\n",
  "Output directory:\n",
  out_dir,
  "\n",
  sep = ""
)
