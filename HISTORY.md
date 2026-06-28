# History

## 0.5.0 - 2026-06-28

Initial PyPI-oriented release for motif-based protein function inference.

Highlights:

- Added `prosig scan` for sequence and FASTA motif scanning with GO-set
  inference from motif-cluster weights.
- Added complete runtime library resolution with `--library-dir`, current
  working-directory discovery, and packaged default fallback.
- Added `prosig build-library --package` to prepare fixed-name core runtime
  library artifacts for distribution.
- Added motif-cluster score board generation with Jeffreys-prior weighting,
  positive-weight storage, and internal calibration metadata.
- Added synthetic cluster GO signatures in `clusters_meta.tsv` under
  `composed_go`.
- Added `prosig inspect function` cluster ID resolution and
  `prosig inspect cluster` reports for cluster metadata and identifying motifs.
- Added multiprocessing motif scanning with sparse motif-hit output.
- Added freshness checks for GO clustering and motif feature extraction.
