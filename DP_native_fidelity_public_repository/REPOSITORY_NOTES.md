# Repository notes

## What was changed when preparing this public package

The scientific scripts were preserved from the supplied analysis archive. Public-release cleanup was limited to repository hygiene and documentation:

- Removed `__MACOSX` metadata directories.
- Removed `.DS_Store` files.
- Removed `.RData` and `.Rhistory` session-state files.
- Replaced the old internal clean-session README with `scripts/CLEAN_SESSION_AUDIT_NOTE.md`.
- Added repository-level documentation, dependency lists, citation metadata, and SHA-256 manifest.

No numerical result, statistical threshold, programme membership, or modelling procedure was intentionally changed during repository packaging.

## Legacy `heldout` naming

Some finalized R scripts use internal object names such as `heldout_genes`. In those scripts this denotes the genes remaining after removing programme-defining genes. It does **not** mean that genes were the cross-validation hold-out unit. The donor-blocked modelling holds out donors.

The legacy names are retained in the frozen code for traceability and to minimize risk of introducing a code change during public-release packaging.
