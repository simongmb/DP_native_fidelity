# Environment notes

## Python

`03_GSE272333_native_fidelity.py` uses Python syntax requiring Python 3.10+ and imports:

- numpy
- pandas

Install with:

```bash
python -m pip install -r environment/requirements.txt
```

## R

The R scripts use a mixture of CRAN and Bioconductor packages listed in:

- `R_CRAN_packages.txt`
- `R_Bioconductor_packages.txt`

Many core scripts contain their own package-install checks. For a publication-grade frozen environment, authors may optionally initialize `renv` and record a lockfile on the machine used for the final rerun.

Exact historical package versions could not be reconstructed reliably from the supplied script archive alone and therefore are not fabricated here.
