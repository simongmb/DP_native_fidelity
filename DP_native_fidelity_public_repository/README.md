# Anatomy-referenced transcriptomic fidelity in human dermal papilla cells

Reproducibility code archive for the manuscript:

**Anatomy-referenced transcriptomics separates molecular reversal from native identity restoration in human dermal papilla cells**

This repository reproduces the transcriptomic rescue, anatomy-referenced native-fidelity, programme-definition, donor-aware sensitivity, model-reconstruction, and figure-generation analyses reported in the manuscript.

## Study overview

The analysis integrates two independent public human transcriptomic resources:

1. **GSE44765** — donor-matched human dermal papilla (DP) culture series used to quantify culture-associated transcriptional loss and 3D-associated rescue.
2. **GSE272333** — laser-capture microdissection atlas of normal human hair-follicle compartments used independently to quantify direction-aligned native anatomical fidelity.

The central analysis distinguishes **how much expression reverses** (rescue magnitude) from **whether the recovered direction is concordant with native DP anatomy** (native fidelity).

## Repository layout

```text
.
├── README.md
├── RUN_ORDER.md
├── DATA_INPUTS.md
├── REPOSITORY_NOTES.md
├── BEFORE_PUBLIC_RELEASE.md
├── CITATION.cff
├── MANIFEST_SHA256.tsv
├── environment/
│   ├── requirements.txt
│   ├── R_CRAN_packages.txt
│   ├── R_Bioconductor_packages.txt
│   └── ENVIRONMENT_NOTES.md
├── data/
│   └── README.md
└── scripts/
    ├── 01_GSE44765_preprocess.R
    ├── 02_GSE44765_rescue_signature.R
    ├── 03_GSE272333_native_fidelity.py
    ├── 04_freeze_restoration_programmes_final.R
    ├── 05_build_GSE272333_subexpr.R
    ├── 06_sensitivity_81_vs_105.R
    ├── 07_donor_aware_association_sensitivity_FINAL.R
    ├── load_subexpr_clean.R
    ├── FIGURE 1/ ... FIGURE 7/
    ├── SUPPLEMENTARY FIGURES/
    └── _archive_supporting/
```

## Primary analysis order

Run the core scripts from the project root in numerical order:

```bash
Rscript scripts/01_GSE44765_preprocess.R
Rscript scripts/02_GSE44765_rescue_signature.R
python scripts/03_GSE272333_native_fidelity.py --project-root .
Rscript scripts/04_freeze_restoration_programmes_final.R
Rscript scripts/05_build_GSE272333_subexpr.R
Rscript scripts/06_sensitivity_81_vs_105.R
Rscript scripts/07_donor_aware_association_sensitivity_FINAL.R
```

See `RUN_ORDER.md` for figure scripts and additional details.

## Expected manuscript-final checkpoints

The core scripts contain hard or explicit QC checks for the principal manuscript-final quantities, including:

- 1,571 strict-rescue genes from GSE44765.
- 1,553 strict-rescue genes represented in the GSE272333 anatomical-reference matrix.
- Rescue fraction versus native fidelity: Spearman rho approximately -0.108.
- 81-gene Native-aligned rescue programme.
- 105-gene threshold-defined low-fidelity branch; a secondary proliferation-associated filter yields the 99-gene Alternative-state set used in main displays.
- 1,373 non-programme genes after exclusion of the displayed 81 + 99 programme genes.
- Donor-blocked reconstruction uses held-out **donors**, not held-out genes.

## Important terminology note

Several historical internal variable names in the frozen analysis scripts use `heldout` or `heldout_genes` for the 1,373-gene set remaining after programme-gene exclusion. These names are retained to avoid altering the finalized analysis code. In the manuscript and public interpretation, these are **non-programme predictor genes**. In donor-blocked cross-validation, the held-out unit is the **donor**.

## Public data and redistribution

The original GEO source files are not redistributed in this archive. Download the public source data directly from NCBI GEO and place them as described in `DATA_INPUTS.md` / `data/README.md`.

## Archived supporting analyses

`scripts/_archive_supporting/` contains exploratory or cross-context analyses retained for auditability. They are **not independent external validation of the complete anatomy-referenced discovery pipeline** and should not be interpreted as such.

## Reproducibility status

- The Python native-fidelity script has been syntax-checked in a clean Python environment.
- R scripts were preserved from the submitted analysis archive. Exact R package versions were not encoded in the source archive supplied for this public-package build; required package names are listed under `environment/`.
- Local RStudio/Mac session artifacts (`.RData`, `.Rhistory`, `.DS_Store`, `__MACOSX`) were removed from this public package.

## Data and code citation

Add the final repository DOI and manuscript DOI to `CITATION.cff` and `BEFORE_PUBLIC_RELEASE.md` before freezing the public release.


## Reuse and collaboration

This repository is made publicly accessible primarily for transparency, peer review,
and reproducibility assessment. Reuse of the code or analytical framework in a new
research project, derivative work, redistribution, publication, or commercial activity
requires prior written permission from the corresponding author.

Any authorized academic use must cite the associated publication and this repository.
Researchers interested in adapting the framework to new datasets are encouraged to
contact the authors regarding permission and potential collaboration. Authorship in any
subsequent publication should be based on substantive contribution and the authorship
criteria of the target journal.

See `LICENSE` for the full terms and `CITATION.cff` for citation metadata.
