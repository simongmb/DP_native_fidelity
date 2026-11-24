# Required public data inputs

## GSE44765 — dermal papilla culture series

Primary use: culture-associated transcriptional loss, 3D-associated response, strict-rescue definition, and rescue fraction.

Expected input for `01_GSE44765_preprocess.R`:

```text
data_raw/GSE44765/*.CEL.gz
```

The finalized script expects **18 CEL files** corresponding to the donor-matched culture series.

Source: NCBI Gene Expression Omnibus accession **GSE44765**.

## GSE272333 — human hair-follicle anatomical reference

Primary use: direction-aligned native anatomical fidelity and downstream anatomy-referenced state analyses.

Required inputs include:

1. The author-provided normalized-count matrix for GSE272333.
2. An Ensembl-to-HGNC mapping table with columns equivalent to:
   - `Ensembl_ID`
   - `Approved_Symbol`

Recommended mapping-table location:

```text
results/GSE272333/GSE272333_gene_ID_mapping.csv
```

The Python script can also locate compatible normalized-count/mapping files by filename pattern or accept explicit paths:

```bash
python scripts/03_GSE272333_native_fidelity.py \
  --project-root . \
  --normalized-counts /path/to/GSE272333_normalized_counts.tab.gz \
  --mapping /path/to/GSE272333_gene_ID_mapping.csv
```

The finalized analysis expects **110 atlas profiles** across DP and seven non-DP comparator compartments.

Source: NCBI Gene Expression Omnibus accession **GSE272333**.

## Other GEO accessions

Additional accessions referenced by scripts under `_archive_supporting/` are exploratory/cross-context analyses and are not required to reconstruct the main discovery pipeline.
