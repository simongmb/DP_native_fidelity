# Run order

Run commands from the repository/project root so that relative `data_raw/`, `results/`, and `figures/` paths resolve consistently.

## A. Core analysis

1. `scripts/01_GSE44765_preprocess.R`
   - Preprocesses GSE44765 raw Affymetrix CEL files by RMA.
   - Generates the expression object used by the rescue analysis.

2. `scripts/02_GSE44765_rescue_signature.R`
   - Defines culture-associated changes, 3D-associated responses, strict rescue, and rescue fraction.

3. `scripts/03_GSE272333_native_fidelity.py`
   - Reconstructs native anatomical fidelity using the independent GSE272333 atlas.
   - Produces the 1,553-gene native-fidelity audit table and QC output.

4. `scripts/04_freeze_restoration_programmes_final.R`
   - Freezes the Native-aligned and Alternative-state programme definitions from the finalized audit table.

5. `scripts/05_build_GSE272333_subexpr.R`
   - Reconstructs the 1,553 x 110 expression matrix used by downstream Figures 4–6.

6. `scripts/06_sensitivity_81_vs_105.R`
   - Sensitivity analysis retaining the full 105-gene threshold-defined low-fidelity branch.

7. `scripts/07_donor_aware_association_sensitivity_FINAL.R`
   - Donor-aware within-donor permutation sensitivity analyses for programme-score and non-programme gene associations.

## B. Main figures

After the required upstream result files have been produced, run scripts within:

- `scripts/FIGURE 1/`
- `scripts/FIGURE 2/`
- `scripts/FIGURE 3/`
- `scripts/FIGURE 4/`
- `scripts/FIGURE 5/`
- `scripts/FIGURE 6/`
- `scripts/FIGURE 7/`

The filenames contain version labels retained from the finalized analysis archive. They are preserved rather than renamed to avoid breaking internal source/output references.

## C. Supplementary figures

Run the corresponding scripts under:

- `scripts/SUPPLEMENTARY FIGURES/`

## D. Archived analyses

`scripts/_archive_supporting/` contains exploratory/cross-context analyses retained for provenance and auditability. These analyses are not part of the independent validation claim for the main anatomy-referenced state axis.
