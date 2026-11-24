# Clean-session batch correction for final figure scripts

## Purpose
Remove dependence on a pre-existing RStudio Global Environment / legacy `.RData`
for figure scripts that use the GSE272333 `subexpr` matrix.

## Standard applied
Each corrected script now:
1. locates the project root using its existing logic;
2. sources `scripts/load_subexpr_clean.R`;
3. obtains `subexpr` from the clean loader;
4. hard-stops unless `subexpr` is exactly **1553 genes x 110 profiles**;
5. leaves the downstream plotting/statistical logic unchanged.

## Included
- FIGURE 4/Fig4A_wide_v5_clean_donorAware.R
- FIGURE 4/Fig4C_v3_clean.R
- FIGURE 4/Fig4E_swapped_v3_clean.R
- FIGURE 4/Fig4F_tall_v4_clean.R
- FIGURE 5/Fig5A_wide_v4_clean.R
- FIGURE 5/Fig5B_v4_clean.R
- SUPPLEMENTARY FIGURES/S4B_v4_clean.R

Fig4A is the already corrected donor-aware version; it was not reverted.

## Important
Two scripts identified during the manuscript audit were not present as accessible
standalone files in the extracted script set available in this session:
- FIGURE 5/Fig5C_v8_donor_blocked.R
- FIGURE 6/Fig6A_v5_donor_blocked_complete.R

They should receive the same clean-session loader/QC pattern before the final public
code archive is frozen. Their current scientific outputs do not need to be changed
merely because of this reproducibility cleanup.

## Expected console line
`Clean-session subexpr dimensions: 1553 x 110`
