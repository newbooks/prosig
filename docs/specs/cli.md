# CLI Spec

## Goal

Provide a command and subcommand interface for ProSig workflows.

The CLI should be a thin user-facing layer over package modules. Command modules should parse arguments, validate user intent, and call implementation functions; they should not contain motif discovery, signature scanning, or prediction logic directly.

## CLI Framework

Use Typer for the command system.

Rationale:

- Typer has been used successfully in related prior work.
- It supports command and subcommand apps with `typer.Typer()`.
- Type hints can define CLI argument and option parsing.
- It keeps command definitions readable as the command tree grows.

## Package Layout

Recommended initial layout:

```text
src/prosig/
  cli/
    __init__.py
    app.py
    setup_data.py
    discover.py
    annotate.py
  discovery/
  prediction/
  signatures/
  io/
```

## Entry Point

Expose a `prosig` executable through `pyproject.toml`:

```toml
[project.scripts]
prosig = "prosig.cli.app:main"
```

Typer should be tracked as a runtime dependency:

```toml
dependencies = [
  "typer>=0.12",
]
```

## Initial Commands

Start with these top-level subcommands:

- `prosig setup-data`: download and cache external data for offline use.
- `prosig build-library`: build the minimized GO graph, adjustable function
  clusters, and motif library.
- `prosig scan`: scan one sequence or a FASTA file and infer motif-supported
  GO sets.
- `prosig inspect`: inspect GO terms, accessions, motifs, similarity scores,
  and other diagnostic artifacts.
- `prosig discover`: discover sequence signatures from positive and background or negative sequence sets.
- `prosig annotate`: scan sequences against a signature library and predict
  protein function from signature hits.
- `prosig signatures`: inspect, validate, convert, or summarize signature libraries.

## Design Rules

- Keep CLI commands small and focused.
- Put reusable behavior in package modules outside `prosig.cli`.
- Use explicit names for input and output options.
- Prefer file paths and formats that can be used in reproducible scripts.
- Return non-zero exit codes for invalid inputs or failed workflows.

## `scan`

`scan` performs motif-based function inference from existing build-library
artifacts. Exactly one query source is required:

```text
prosig scan --seq MSEQUENCE
prosig scan --fasta queries.fasta
```

For each query, the command loads the complete runtime library, scans
`prosig_motifs.tsv`, looks up matching motifs in
`motif_cluster_scoreboard.pkl`, and reports inferred GO sets with
motif-cluster weight at or above `--min-weight` (`2.0` by default). The default
screen and JSON report includes the top 5 inferred GO sets; use `--top-n N` to
change the limit, or `--top-n 0` to report all inferred GO sets. Cluster
metadata from `clusters_meta.tsv` provides composed GO terms and descriptions.
Each prediction includes the strongest contributing motif signature.

Runtime library selection is all-or-nothing. `--library-dir DIR` must point to a
directory containing every core runtime file. When `--library-dir` is omitted,
`scan` uses the current working directory if any core library file is present;
otherwise it falls back to packaged defaults. The core files are
`prosig_motifs.tsv`, `motif_cluster_scoreboard.pkl`,
`motif_cluster_scoreboard_meta.json`, `clusters_meta.tsv`, `go_graph.pkl`, and
`accession_mf_go.tsv`.

When `motif_cluster_scoreboard_meta.json` contains calibration records, `scan`
also reports a calibrated confidence reference: the observed `set_accuracy` at
the highest calibration threshold less than or equal to the prediction weight.
This is not a per-cluster probability; it is an empirical calibration summary
for predictions at that weight scale.

Use `--json-out PATH` to write the same report as JSON instead of printing the
human-readable report.
- Keep output formats stable once documented.

## `inspect` Command Group

`prosig inspect` is for diagnostics rather than production analysis. It should
make intermediate artifacts easy to verify before they feed clustering, motif
discovery, or prediction.

Initial implemented commands:

- `prosig inspect go-summary`: report GO artifact metadata and IC coverage
  from the resolved runtime library. `--go-graph PATH` overrides the library
  GO graph for explicit artifact diagnostics.
- `prosig inspect go-term GO:0005524`: report one MF GO term from the resolved
  runtime library. `--ancestors` includes the term itself plus all ancestors.
  `--go-graph PATH` overrides the library GO graph.
- `prosig inspect go-sim GO:0005524 GO:0004672`: report the Lin similarity
  score from the resolved runtime library. `--verbose` or `-v` adds term
  descriptions,
  common ancestors, selected MICA, IC values, status, reason, formula, and a
  compact GO path. The path uses Unicode tree connectors by default; use
  `--tree-style ascii` for ASCII output. `--go-graph PATH` overrides the
  library GO graph.
- `prosig inspect go-set-sim`, `prosig inspect function`, and
  `prosig inspect cluster` use the same complete runtime library resolution as
  `scan`; pass `--library-dir DIR` to override the selected library.

Planned diagnostic commands include accession lookup, motif lookup, motif
summary, GO term-set similarity, and cluster/member inspection.

## `build-library` GO Clustering

`prosig build-library` includes GO accession clustering as part of library
construction rather than exposing a separate top-level clustering command.

The clustering workflow consumes `go_graph.pkl` and
`accession_mf_go.tsv`, builds a sparse GO-set similarity kNN graph, and runs
Leiden community detection to produce freshness-managed intermediates:

```text
leiden_clusters.tsv
leiden_clusters_meta.tsv
leiden_clusters_stats.json
```

It then freshness-manages complete-linkage refinement and writes:

```text
clusters.tsv
clusters_meta.tsv
clusters_stats.json
```

`--min-cluster-similarity FLOAT` controls the final all-pairs similarity floor
and defaults to `0.25`. `--leiden-cluster-out` and `--cluster-out` control the
intermediate and final membership paths respectively.

After final refinement, `build-library` synthesizes cluster-level GO signatures
into the final metadata file:

```text
clusters_meta.tsv
```

The metadata column is `composed_go`. It contains up to 10 semicolon-separated
GO IDs per cluster. The signature step consumes final cluster membership,
`accession_mf_go.tsv`, and `go_graph.pkl`; votes equally by accession;
propagates direct MF terms to ancestors; scores candidates as `support × IC`;
and writes only the selected GO IDs.

`build-library` scans the motif library against final cluster member sequences
to produce sparse motif hit features, then builds motif-based
function-prediction weights:

```text
motif_features.tsv
motif_cluster_scoreboard.pkl
motif_cluster_scoreboard_meta.json
```

Relevant options:

```text
--motif-hits PATH
--motif-scoreboard-out PATH
--motif-scoreboard-meta-out PATH
--motif-scoreboard-min-cluster-size INT
--motif-scoreboard-min-support INT
--motif-scan-processes INT
--package
--package-dir PATH
```

`--motif-hits` is the sparse feature output path and scoreboard input path. The
motif pattern source is the `--motif-out` library, whose `prosig_pattern`
column is compiled and scanned. Motif scanning uses 8 worker processes by
default, controlled by `--motif-scan-processes`, and parent-process progress
logs are emitted as accession chunks complete. The score board step skips
clusters smaller than 10 by default, skips motif-cluster pairs with `TP < 5`,
stores only positive weights in the pickle artifact, and logs/internal-metadata
reports calibration top-1, top-3, set accuracy, average prediction count, and
coverage at weight thresholds 2.0 through 8.0 in 0.5 increments.

`--package` copies the core runtime library artifacts from the working
directory into the package default library directory so installed users can run
`scan` and runtime `inspect` commands without first running `setup-data` and
`build-library`. `--package-dir PATH` overrides the packaging target, mainly for
maintainer workflows and tests. Packaging uses fixed filenames and fails if any
core runtime file is missing.

Clustering graph, Leiden, matrix, cache, and candidate-filter parameters live
in `cluster_config.yaml`, created from a packaged starter template when
missing. The same config also controls the stats output path and clustering
progress log interval. Because the config is a clustering input, editing it
invalidates the cluster outputs.

## Open Questions

- Which command should be implemented first after `setup-data`?
- What should be the first canonical signature library file format?
- Should command output default to human-readable text, structured JSON/TSV, or both?
