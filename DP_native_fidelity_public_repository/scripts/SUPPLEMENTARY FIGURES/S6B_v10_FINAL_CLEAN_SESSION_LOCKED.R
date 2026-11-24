# ============================================================
# Supplementary Figure S6B v10 — FINAL CLEAN-SESSION LOCKED
# Stability of nested top-25 gene selection across outer CV folds
#
# Input:
#   Frozen clean-session Fig6D v7 selection-frequency output (15 genes selected in all 200 outer folds)
#
# This script DOES NOT rerun any model.
# It visualizes the already-computed nested-CV selection frequencies.
#
# Outputs:
#   figures/SUPPLEMENTARY FIGURES/S6_Panel_B_nestedTop25_v5/
#   PDF / PNG / TIFF
#   plot data CSV
#   summary CSV
#   README
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_pkgs <- c("ggplot2", "readr", "dplyr", "forcats", "scales")

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
      "\n\nInstall them first with:\n",
      "install.packages(c(",
      paste0('"', missing_pkgs, '"', collapse = ", "),
      "))\n"
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(forcats)
  library(scales)
})

# ------------------------------------------------------------
# 1. Locate input automatically by REQUIRED COLUMNS
# ------------------------------------------------------------
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
    getwd(),
    start_dir,
    dirname(getwd()),
    dirname(dirname(getwd())),
    dirname(start_dir),
    dirname(dirname(start_dir)),
    dirname(dirname(dirname(start_dir)))
  )
)

candidate_roots <- candidate_roots[
  dir.exists(candidate_roots)
]

# Prefer the active HairRegeneration_AI project.
root_scores <- vapply(
  candidate_roots,
  function(x) {
    s <- 0

    if (
      grepl(
        "HairRegeneration_AI$",
        basename(x),
        ignore.case = TRUE
      )
    ) s <- s + 50

    if (dir.exists(file.path(x, "figures"))) s <- s + 10
    if (dir.exists(file.path(x, "results"))) s <- s + 10
    if (dir.exists(file.path(x, "data_raw"))) s <- s + 5

    s
  },
  numeric(1)
)

project_root <- candidate_roots[
  which.max(root_scores)
]

cat("\nProject root:\n", project_root, "\n")

required_input_cols <- c(
  "GeneSymbol",
  "Selection_count",
  "Selection_frequency",
  "Mean_abs_train_rho",
  "Median_rank"
)

# FINAL SOURCE LOCK:
# Search all candidate selection-frequency tables, then accept ONLY
# the branch that reproduces the frozen Supplementary Table 6 QC.

all_csv <- list.files(
  project_root,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

schema_hits <- character(0)

if (length(all_csv) > 0) {
  for (f in all_csv) {
    hdr <- tryCatch(
      names(
        readr::read_csv(
          f,
          n_max = 2,
          show_col_types = FALSE,
          progress = FALSE
        )
      ),
      error = function(e) character(0)
    )

    if (
      length(hdr) > 0 &&
      all(required_input_cols %in% hdr)
    ) {
      schema_hits <- c(
        schema_hits,
        f
      )
    }
  }
}

schema_hits <- unique(schema_hits)

EXPECTED_EVER_SELECTED <- 39L
EXPECTED_ALWAYS_SELECTED <- 15L
EXPECTED_TOTAL_OPPORTUNITIES <- 200L

EXPECTED_ALWAYS_GENES <- sort(c(
  "LAMC3", "NFATC4", "JAM2", "SHOX2", "THSD7A", "USP11", "SPON1",
  "CEP126", "VIM", "IFITM1", "CACHD1", "FGF7", "TMEM204", "SCARA5",
  "SERPING1"
))

passes_frozen_qc <- function(tab) {

  if (!is.data.frame(tab)) return(FALSE)

  if (!all(required_input_cols %in% names(tab))) {
    return(FALSE)
  }

  if (nrow(tab) != EXPECTED_EVER_SELECTED) {
    return(FALSE)
  }

  if (anyDuplicated(tab$GeneSymbol) > 0) {
    return(FALSE)
  }

  if (
    any(
      tab$Selection_frequency < 0 |
      tab$Selection_frequency > 1,
      na.rm = TRUE
    )
  ) {
    return(FALSE)
  }

  always_genes <- sort(
    as.character(
      tab$GeneSymbol[
        abs(tab$Selection_frequency - 1) < 1e-12
      ]
    )
  )

  if (
    length(always_genes) != EXPECTED_ALWAYS_SELECTED ||
    !identical(always_genes, EXPECTED_ALWAYS_GENES)
  ) {
    return(FALSE)
  }

  denom_est <- tab %>%
    dplyr::filter(
      Selection_frequency > 0
    ) %>%
    dplyr::mutate(
      denom =
        Selection_count /
        Selection_frequency
    ) %>%
    dplyr::pull(
      denom
    )

  if (
    length(denom_est) == 0 ||
    !all(is.finite(denom_est))
  ) {
    return(FALSE)
  }

  total_opportunities <- round(
    median(
      denom_est
    )
  )

  total_opportunities ==
    EXPECTED_TOTAL_OPPORTUNITIES
}

validated_files <- character(0)

if (length(schema_hits) > 0) {
  for (f in schema_hits) {

    tmp <- tryCatch(
      readr::read_csv(
        f,
        show_col_types = FALSE,
        progress = FALSE
      ),
      error = function(e) NULL
    )

    if (
      !is.null(tmp) &&
      passes_frozen_qc(tmp)
    ) {
      validated_files <- c(
        validated_files,
        f
      )
    }
  }
}

validated_files <- unique(
  validated_files
)

input_file <- NA_character_
dat_from_env <- NULL

if (length(validated_files) > 0) {

  hit_scores <- vapply(
    validated_files,
    function(x) {
      s <- 0

      if (
        grepl(
          "Fig6D|Figure.?6|FIGURE[ _]?6",
          x,
          ignore.case = TRUE
        )
      ) s <- s + 50

      if (
        grepl(
          "top25|top_25|nestedTop25",
          x,
          ignore.case = TRUE
        )
      ) s <- s + 30

      if (
        grepl(
          "selection",
          x,
          ignore.case = TRUE
        )
      ) s <- s + 10

      if (
        grepl(
          "frequency",
          x,
          ignore.case = TRUE
        )
      ) s <- s + 10

      s
    },
    numeric(1)
  )

  input_file <- validated_files[
    which.max(
      hit_scores
    )
  ]

  cat(
    "\nFrozen Table 6 branch selected:\n",
    input_file,
    "\n"
  )

} else {

  # Fallback to GlobalEnv, but only if an object passes the same
  # frozen Supplementary Table 6 QC.
  env_names <- ls(
    envir = .GlobalEnv
  )

  env_hits <- list()

  for (nm in env_names) {

    obj <- tryCatch(
      get(
        nm,
        envir = .GlobalEnv
      ),
      error = function(e) NULL
    )

    if (
      is.data.frame(obj) &&
      passes_frozen_qc(obj)
    ) {
      env_hits[[nm]] <- obj
    }
  }

  if (length(env_hits) > 0) {

    picked_name <- names(env_hits)[1]

    dat_from_env <-
      env_hits[[picked_name]]

    cat(
      "\nFrozen Table 6 GlobalEnv object selected:\n",
      picked_name,
      "\n"
    )
  }
}

if (
  is.na(input_file) &&
  is.null(dat_from_env)
) {
  stop(
    paste0(
      "\nQC STOP\n\n",
      "No selection-frequency table reproduces the frozen Supplementary Table 6 branch.\n",
      "Required QC:\n",
      "- 39 genes selected at least once\n",
      "- 200 outer-fold opportunities\n",
      "- 15 genes selected in 100% of opportunities\n",
      "- exact frozen 14-gene set\n\n",
      "The script will not export a figure from a historical or mismatched branch.\n"
    )
  )
}

# ------------------------------------------------------------
# 2. Read + QC
# ------------------------------------------------------------
dat <- if (
  !is.null(dat_from_env)
) {
  tibble::as_tibble(
    dat_from_env
  )
} else {
  readr::read_csv(
    input_file,
    show_col_types = FALSE
  )
}

required_cols <- c(
  "GeneSymbol",
  "Selection_count",
  "Selection_frequency",
  "Mean_abs_train_rho",
  "Median_rank"
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
      "\n\nAvailable columns:\n",
      paste(names(dat), collapse = ", "),
      "\n"
    )
  )
}

if (anyDuplicated(dat$GeneSymbol) > 0) {
  stop(
    "\nQC STOP: duplicate GeneSymbol entries detected.\n"
  )
}

if (any(dat$Selection_frequency < 0 | dat$Selection_frequency > 1, na.rm = TRUE)) {
  stop(
    "\nQC STOP: Selection_frequency must be between 0 and 1.\n"
  )
}

# Infer total number of outer-fold opportunities from count/frequency
denom_est <- dat %>%
  filter(
    Selection_frequency > 0
  ) %>%
  mutate(
    denom = Selection_count / Selection_frequency
  ) %>%
  pull(denom)

total_opportunities <- round(
  median(denom_est)
)

cat(
  "Estimated total outer-fold selection opportunities:",
  total_opportunities,
  "\n"
)

# ------------------------------------------------------------
# 3. Choose top recurrent genes for display
# ------------------------------------------------------------
top_n_display <- 25

plot_dat <- dat %>%
  arrange(
    desc(Selection_frequency),
    Median_rank,
    desc(Mean_abs_train_rho)
  ) %>%
  slice_head(
    n = top_n_display
  ) %>%
  mutate(
    GeneSymbol = forcats::fct_reorder(
      GeneSymbol,
      Selection_frequency
    )
  )

# ------------------------------------------------------------
# 4. Summary
# ------------------------------------------------------------
summary_dat <- tibble::tibble(
  Metric = c(
    "Genes ever selected",
    "Total outer-fold selection opportunities",
    "Genes selected in >= 50% of opportunities",
    "Genes selected in >= 80% of opportunities",
    "Genes selected in >= 90% of opportunities",
    "Maximum selection frequency"
  ),
  Value = c(
    nrow(dat),
    total_opportunities,
    sum(dat$Selection_frequency >= 0.50),
    sum(dat$Selection_frequency >= 0.80),
    sum(dat$Selection_frequency >= 0.90),
    max(dat$Selection_frequency)
  )
)

print(summary_dat)

# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------
p <- ggplot(
  plot_dat,
  aes(
    x = Selection_frequency,
    y = GeneSymbol
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Selection_frequency,
      yend = GeneSymbol
    ),
    linewidth = 0.32
  ) +
  geom_point(
    size = 1.65
  ) +
  geom_vline(
    xintercept = c(0.50, 0.80),
    linetype = "dashed",
    linewidth = 0.24
  ) +
  scale_x_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    title = "Nested top-25 selection identifies recurrent genes",
    subtitle = paste0(
      "Top ",
      top_n_display,
      " recurrent genes across ",
      total_opportunities,
      " outer-fold selection opportunities"
    ),
    x = "Selection frequency across outer CV folds",
    y = NULL
  ) +
  theme_classic(
    base_size = 6
  ) +
  theme(
    plot.title = element_text(
      size = 6.0,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    plot.subtitle = element_text(
      size = 4.6,
      hjust = 0.5,
      margin = margin(
        b = 4
      )
    ),
    axis.title.x = element_text(
      size = 5.2,
      margin = margin(
        t = 4
      )
    ),
    axis.text.x = element_text(
      size = 4.6
    ),
    axis.text.y = element_text(
      size = 4.5
    ),
    plot.margin = margin(
      5,
      10,
      5,
      5
    )
  )

# ------------------------------------------------------------
# 6. Locate project root
# ------------------------------------------------------------
project_root <- start_dir

root_candidates <- unique(
  c(
    start_dir,
    dirname(start_dir),
    dirname(dirname(start_dir)),
    getwd(),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )
)

root_candidates <- root_candidates[
  dir.exists(root_candidates)
]

for (r in root_candidates) {
  if (
    dir.exists(file.path(r, "figures")) ||
    dir.exists(file.path(r, "results")) ||
    grepl(
      "HairRegeneration_AI$",
      basename(r),
      ignore.case = TRUE
    )
  ) {
    project_root <- r
    break
  }
}

out_dir <- file.path(
  project_root,
  "figures",
  "SUPPLEMENTARY FIGURES",
  "S6_Panel_B_nestedTop25_FINAL_CLEAN_SESSION_v10"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 7. Save figure
# ------------------------------------------------------------
pdf_file <- file.path(
  out_dir,
  "S6_Panel_B_nestedTop25_FINAL_CLEAN_SESSION_v10.pdf"
)

png_file <- file.path(
  out_dir,
  "S6_Panel_B_nestedTop25_FINAL_CLEAN_SESSION_v10.png"
)

tiff_file <- file.path(
  out_dir,
  "S6_Panel_B_nestedTop25_FINAL_CLEAN_SESSION_v10.tiff"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 58,
  height = 66,
  units = "mm",
  device = "pdf"
)

ggsave(
  filename = png_file,
  plot = p,
  width = 58,
  height = 66,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = p,
  width = 58,
  height = 66,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ------------------------------------------------------------
# 8. Save tables
# ------------------------------------------------------------
readr::write_csv(
  plot_dat %>%
    mutate(
      GeneSymbol = as.character(GeneSymbol)
    ),
  file.path(
    out_dir,
    "S6B_plot_data.csv"
  )
)

readr::write_csv(
  summary_dat,
  file.path(
    out_dir,
    "S6B_summary.csv"
  )
)

readr::write_csv(
  dat,
  file.path(
    out_dir,
    "S6B_all_selection_frequencies.csv"
  )
)

# ------------------------------------------------------------
# 9. README
# ------------------------------------------------------------
readme_lines <- c(
  "Supplementary Figure S6B v10 — FINAL CLEAN-SESSION LOCKED",
  "",
  "Purpose:",
  "Visualize the stability of genes selected by the nested top-25 procedure with feature ranking confined to each outer training fold.",
  "",
  paste0(
    "Input: ",
    if (
      !is.null(dat_from_env)
    ) {
      "GlobalEnv data.frame"
    } else {
      normalizePath(input_file)
    }
  ),
  "",
  paste0(
    "Genes ever selected: ",
    nrow(dat)
  ),
  paste0(
    "Estimated total outer-fold selection opportunities: ",
    total_opportunities
  ),
  "",
  "Display:",
  paste0(
    "- Top ",
    top_n_display,
    " recurrent genes ranked by Selection_frequency"
  ),
  "- Dashed vertical reference lines at 50% and 80%",
  "",
  "Important:",
  "This script does NOT rerun any prediction model.",
  "It visualizes selection-frequency results already generated by the nested top-25 analysis.",
  "",
  "Output directory:",
  normalizePath(
    out_dir,
    mustWork = FALSE
  )
)

writeLines(
  readme_lines,
  file.path(
    out_dir,
    "README_S6B_nestedTop25_FINAL_CLEAN_SESSION_v10.txt"
  )
)

cat(
  "\nDone.\n\n",
  "Output directory:\n",
  out_dir,
  "\n\nMain file:\n",
  pdf_file,
  "\n"
)
