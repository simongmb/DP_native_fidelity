#!/usr/bin/env python3
"""
03_GSE272333_native_fidelity.py

Reconstruct direction-aligned native anatomical fidelity for the 1,571
GSE44765 strict-rescue genes using the independent GSE272333 human hair-
follicle anatomical atlas.

Final frozen method used in the manuscript:
  1. Author-provided GSE272333 normalized counts.
  2. log2(x + 1).
  3. Gene-wise z-standardization across all 110 atlas profiles using the
     population standard deviation (ddof=0).
  4. Within each donor, raw fidelity = mean(DP) - mean(all available non-DP).
  5. Average donor-specific differences across donors.
  6. Direction align by GSE44765 restoration direction:
       Restored_upward   -> +1
       Restored_downward -> -1
  7. Repeat donor-aware DP-minus-comparator calculation for each of
     DC, CTS, HM1, HM2, IRS1, ORS1 and ORS2.
  8. Anatomical consistency = number of the seven aligned contrasts > 0.

This script does NOT define the 81/99 programmes. It produces the finalized
1,553-gene audit table consumed by 04_freeze_restoration_programmes.R.

Required inputs (searched under project root):
  - results/GSE44765_rescue/GSE44765_STRICT_rescue_genes.csv
  - GSE272333 normalized-count matrix (*.tab.gz / *.tsv.gz / *.txt.gz)
  - a two-column Ensembl-to-HGNC mapping table containing Ensembl_ID and
    Approved_Symbol (recommended location: results/GSE272333/
    GSE272333_gene_ID_mapping.csv)

Outputs:
  results/DP_native_fidelity_reconstruction/
    02_gene_level_native_fidelity_audit.csv
    00_native_fidelity_QC.csv

The script contains hard QC stops for the manuscript-final counts and the
final rescue-vs-fidelity correlation.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import numpy as np
import pandas as pd

COMPARATORS = ["DC", "CTS", "HM1", "HM2", "IRS1", "ORS1", "ORS2"]
DIRECTION_SIGN = {"Restored_upward": 1.0, "Restored_downward": -1.0}
EXPECTED_STRICT = 1571
EXPECTED_MAPPED = 1557
EXPECTED_ANALYZABLE = 1553
EXPECTED_RHO = -0.107691
RHO_TOL = 5e-4


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Reconstruct GSE272333 native fidelity.")
    p.add_argument("--project-root", type=Path, default=None)
    p.add_argument("--strict", type=Path, default=None,
                   help="GSE44765_STRICT_rescue_genes.csv")
    p.add_argument("--normalized-counts", type=Path, default=None,
                   help="GSE272333 author-provided normalized-count matrix")
    p.add_argument("--mapping", type=Path, default=None,
                   help="Ensembl-to-HGNC mapping CSV/TSV")
    p.add_argument("--output-dir", type=Path, default=None)
    return p.parse_args()


def find_project_root(start: Path) -> Path:
    candidates = [start.resolve(), *start.resolve().parents]
    for root in candidates:
        if (root / "results").exists() or (root / "scripts").exists():
            return root
    return start.resolve()


def first_existing(candidates: list[Path]) -> Path | None:
    for p in candidates:
        if p.exists() and p.is_file():
            return p
    return None


def find_one(root: Path, patterns: list[str], label: str) -> Path:
    hits: list[Path] = []
    for pat in patterns:
        hits.extend(root.rglob(pat))
    hits = sorted({p.resolve() for p in hits if p.is_file() and "__MACOSX" not in p.parts})
    if len(hits) == 1:
        return hits[0]
    if len(hits) == 0:
        raise FileNotFoundError(f"Cannot find {label}. Searched: {patterns}")
    raise RuntimeError(
        f"Found multiple candidate files for {label}; pass an explicit command-line path:\n" +
        "\n".join(str(x) for x in hits)
    )


def read_table(path: Path) -> pd.DataFrame:
    lower = path.name.lower()
    if lower.endswith((".tsv", ".tsv.gz", ".tab", ".tab.gz", ".txt", ".txt.gz")):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def strip_ensembl_version(x: pd.Series) -> pd.Series:
    return x.astype(str).str.replace(r"\.\d+$", "", regex=True)


def donor_diff_for_comparator(long: pd.DataFrame, comp: str) -> pd.Series:
    sub = long[long["compartment"].isin(["DP", comp])]
    piv = (
        sub.groupby(["GeneSymbol", "donor", "compartment"], observed=True)["z"]
        .mean()
        .unstack("compartment")
    )
    piv = piv.dropna(subset=["DP", comp])
    delta = piv["DP"] - piv[comp]
    return delta.groupby(level=0).mean()


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    root = args.project_root.resolve() if args.project_root else find_project_root(script_dir)

    strict_path = args.strict.resolve() if args.strict else first_existing([
        root / "results" / "GSE44765_rescue" / "GSE44765_STRICT_rescue_genes.csv",
    ])
    if strict_path is None:
        strict_path = find_one(root, ["GSE44765_STRICT_rescue_genes.csv"], "strict-rescue table")

    norm_path = args.normalized_counts.resolve() if args.normalized_counts else first_existing([
        root / "data" / "GSE272333_normalized_counts.tab.gz",
        root / "results" / "GSE272333" / "GSE272333_normalized_counts.tab.gz",
        root / "GSE272333_normalized_counts.tab.gz",
    ])
    if norm_path is None:
        norm_path = find_one(
            root,
            ["*GSE272333*normalized*count*.tab.gz", "*GSE272333*normalized*count*.tsv.gz",
             "*GSE272333*normalized*count*.txt.gz", "*GSE272333*normalized*count*.csv.gz"],
            "GSE272333 normalized-count matrix",
        )

    mapping_path = args.mapping.resolve() if args.mapping else first_existing([
        root / "results" / "GSE272333" / "GSE272333_gene_ID_mapping.csv",
        root / "data" / "GSE272333_gene_ID_mapping.csv",
    ])
    if mapping_path is None:
        mapping_path = find_one(
            root,
            ["*GSE272333*gene*ID*mapping*.csv", "*gene_ID_mapping_audit.csv"],
            "GSE272333 Ensembl-to-HGNC mapping table",
        )

    out_dir = args.output_dir.resolve() if args.output_dir else root / "results" / "DP_native_fidelity_reconstruction"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Project root:       {root}")
    print(f"Strict rescue:      {strict_path}")
    print(f"Normalized counts:  {norm_path}")
    print(f"Gene mapping:       {mapping_path}")
    print(f"Output directory:   {out_dir}\n")

    strict = read_table(strict_path)
    require_columns(strict, ["GeneSymbol", "Rescue_fraction", "Rescue_type"], "Strict-rescue table")
    strict = strict.drop_duplicates("GeneSymbol").copy()
    if len(strict) != EXPECTED_STRICT:
        raise RuntimeError(f"QC FAIL: expected {EXPECTED_STRICT} strict-rescue genes, found {len(strict)}")

    mapping = read_table(mapping_path)
    # tolerate common names but normalize to the two final names
    ren = {}
    for c in mapping.columns:
        cl = c.lower()
        if cl in {"ensembl_id", "ensembl", "gene_id", "ensembl_gene_id"}:
            ren[c] = "Ensembl_ID"
        elif cl in {"approved_symbol", "genesymbol", "gene_symbol", "symbol", "hgnc_symbol"}:
            ren[c] = "Approved_Symbol"
    mapping = mapping.rename(columns=ren)
    require_columns(mapping, ["Ensembl_ID", "Approved_Symbol"], "Gene mapping table")
    mapping = mapping[["Ensembl_ID", "Approved_Symbol"]].dropna().copy()
    mapping["Ensembl_ID"] = strip_ensembl_version(mapping["Ensembl_ID"])
    mapping["Approved_Symbol"] = mapping["Approved_Symbol"].astype(str).str.strip()
    mapping = mapping.drop_duplicates("Ensembl_ID", keep="first")

    mapped_symbols = set(mapping.loc[mapping["Approved_Symbol"].isin(strict["GeneSymbol"]), "Approved_Symbol"])
    if len(mapped_symbols) != EXPECTED_MAPPED:
        raise RuntimeError(
            f"QC FAIL: expected {EXPECTED_MAPPED} strict-rescue symbols represented in mapping, "
            f"found {len(mapped_symbols)}"
        )

    expr = read_table(norm_path)
    if expr.shape[1] < 2:
        raise RuntimeError("Normalized-count matrix has fewer than two columns.")
    id_col = expr.columns[0]
    expr = expr.copy()
    expr[id_col] = strip_ensembl_version(expr[id_col])
    sample_cols = list(expr.columns[1:])
    if len(sample_cols) != 110:
        raise RuntimeError(f"QC FAIL: expected 110 atlas profiles, found {len(sample_cols)}")

    meta_rows = []
    for s in sample_cols:
        parts = str(s).split("_")
        if len(parts) < 2:
            raise RuntimeError(f"Cannot parse compartment/donor from sample name: {s}")
        meta_rows.append((s, parts[0], parts[1]))
    meta = pd.DataFrame(meta_rows, columns=["sample", "compartment", "donor"])
    missing_comps = [c for c in ["DP", *COMPARATORS] if c not in set(meta["compartment"])]
    if missing_comps:
        raise RuntimeError(f"QC FAIL: missing expected compartments: {missing_comps}")

    x = expr.merge(mapping, left_on=id_col, right_on="Ensembl_ID", how="left")
    x = x[x["Approved_Symbol"].isin(set(strict["GeneSymbol"]))].copy()
    # If multiple Ensembl rows map to one HGNC symbol, arithmetic-average normalized counts.
    x[sample_cols] = x[sample_cols].apply(pd.to_numeric, errors="coerce")
    x = x.groupby("Approved_Symbol", as_index=False)[sample_cols].mean()

    if x["Approved_Symbol"].nunique() != EXPECTED_ANALYZABLE:
        raise RuntimeError(
            f"QC FAIL: expected {EXPECTED_ANALYZABLE} analyzable strict-rescue genes in normalized matrix, "
            f"found {x['Approved_Symbol'].nunique()}"
        )

    g = np.log2(x.set_index("Approved_Symbol")[sample_cols] + 1.0)
    mu = g.mean(axis=1)
    sd = g.std(axis=1, ddof=0).replace(0, np.nan)
    z = g.sub(mu, axis=0).div(sd, axis=0)
    if z.isna().all(axis=1).any():
        bad = z.index[z.isna().all(axis=1)].tolist()
        raise RuntimeError(f"QC FAIL: zero-variance genes after transformation: {bad[:20]}")

    long = (
        z.reset_index()
        .melt(id_vars="Approved_Symbol", var_name="sample", value_name="z")
        .rename(columns={"Approved_Symbol": "GeneSymbol"})
        .merge(meta, on="sample", how="left", validate="many_to_one")
    )

    # Composite donor-aware DP-minus-all-non-DP fidelity.
    donor_stats = (
        long.assign(is_dp=long["compartment"].eq("DP"))
        .groupby(["GeneSymbol", "donor", "is_dp"], observed=True)["z"]
        .mean()
        .unstack("is_dp")
    )
    require_dp = {False, True}
    if not require_dp.issubset(set(donor_stats.columns)):
        raise RuntimeError("QC FAIL: donor-level table lacks DP or non-DP values.")
    raw_fidelity = (donor_stats[True] - donor_stats[False]).groupby(level=0).mean()

    strict_idx = strict.set_index("GeneSymbol")
    direction = strict_idx["Rescue_type"].map(DIRECTION_SIGN)
    if direction.isna().any():
        bad = direction.index[direction.isna()].tolist()
        raise RuntimeError(f"QC FAIL: unrecognized Rescue_type values for genes: {bad[:20]}")
    aligned_fidelity = raw_fidelity * direction.reindex(raw_fidelity.index)

    comp_df = pd.DataFrame(index=aligned_fidelity.index)
    for comp in COMPARATORS:
        raw_delta = donor_diff_for_comparator(long, comp)
        comp_df[f"aligned_DP_vs_{comp}"] = raw_delta * direction.reindex(raw_delta.index)

    consistency_n = (comp_df > 0).sum(axis=1)

    audit = strict_idx.loc[aligned_fidelity.index].copy()
    audit["native_fidelity_raw_DP_minus_nonDP"] = raw_fidelity.reindex(audit.index)
    audit["native_fidelity_direction_aligned"] = aligned_fidelity.reindex(audit.index)
    audit["anatomical_consistency_n_of_7"] = consistency_n.reindex(audit.index).astype(int)
    audit = audit.join(comp_df)
    audit = audit.reset_index()

    rescue_q3 = float(audit["Rescue_fraction"].quantile(0.75))
    fidelity_q1 = float(audit["native_fidelity_direction_aligned"].quantile(0.25))
    fidelity_q3 = float(audit["native_fidelity_direction_aligned"].quantile(0.75))

    # These flags are annotations only. Programme freezing is performed by script 04.
    audit["high_rescue"] = audit["Rescue_fraction"] > rescue_q3
    audit["high_native_fidelity"] = audit["native_fidelity_direction_aligned"] > fidelity_q3
    audit["low_native_fidelity"] = audit["native_fidelity_direction_aligned"] < fidelity_q1
    audit["anatomy_consistent_candidate"] = audit["high_rescue"] & audit["high_native_fidelity"] & (audit["anatomical_consistency_n_of_7"] >= 5)
    audit["non_native_candidate_pre_exclusion"] = audit["high_rescue"] & audit["low_native_fidelity"] & (audit["anatomical_consistency_n_of_7"] <= 2)

    rho = float(audit[["Rescue_fraction", "native_fidelity_direction_aligned"]].corr(method="spearman").iloc[0, 1])
    if not math.isfinite(rho) or abs(rho - EXPECTED_RHO) > RHO_TOL:
        raise RuntimeError(f"QC FAIL: expected Spearman rho ~ {EXPECTED_RHO:.6f}, reconstructed {rho:.6f}")

    audit_file = out_dir / "02_gene_level_native_fidelity_audit.csv"
    audit.to_csv(audit_file, index=False)

    qc = pd.DataFrame({
        "metric": [
            "strict_rescue_n", "mapped_strict_rescue_n", "analyzable_n", "atlas_profiles_n",
            "rescue_fraction_Q3", "native_fidelity_Q1", "native_fidelity_Q3",
            "rescue_vs_native_fidelity_spearman_rho"
        ],
        "value": [
            len(strict), len(mapped_symbols), len(audit), len(sample_cols), rescue_q3,
            fidelity_q1, fidelity_q3, rho
        ]
    })
    qc_file = out_dir / "00_native_fidelity_QC.csv"
    qc.to_csv(qc_file, index=False)

    print("\nQC PASS")
    print(f"Strict rescue genes:      {len(strict)}")
    print(f"Mapped to HGNC:           {len(mapped_symbols)}")
    print(f"Analyzable atlas genes:   {len(audit)}")
    print(f"Rescue Q3:                {rescue_q3:.9f}")
    print(f"Native fidelity Q1/Q3:    {fidelity_q1:.9f} / {fidelity_q3:.9f}")
    print(f"Spearman rho:             {rho:.6f}")
    print(f"\nWrote:\n  {audit_file}\n  {qc_file}")


if __name__ == "__main__":
    main()
