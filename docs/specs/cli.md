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
- Keep output formats stable once documented.

## `inspect` Command Group

`prosig inspect` is for diagnostics rather than production analysis. It should
make intermediate artifacts easy to verify before they feed clustering, motif
discovery, or prediction.

Initial implemented commands:

- `prosig inspect go-summary --go-graph go_graph.pkl`: report GO artifact
  metadata and IC coverage.
- `prosig inspect go-term GO:0005524 --go-graph go_graph.pkl`: report one MF
  GO term. `--ancestors` includes the term itself plus all ancestors.
- `prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl`:
  report the Lin similarity score. `--verbose` or `-v` adds term descriptions,
  common ancestors, selected MICA, IC values, status, reason, formula, and a
  compact GO path. The path uses Unicode tree connectors by default; use
  `--tree-style ascii` for ASCII output.

Planned diagnostic commands include accession lookup, motif lookup, motif
summary, GO term-set similarity, and cluster/member inspection.

## `build-library` GO Clustering

`prosig build-library` includes GO accession clustering as part of library
construction rather than exposing a separate top-level clustering command.

The clustering workflow consumes `go_graph.pkl` and
`accession_mf_go.tsv`, builds a sparse GO-set similarity kNN graph, and runs
Leiden community detection to produce:

```text
go_clusters.tsv
go_clusters_stats.json
```

The command rebuilds these cluster artifacts only when they are missing,
older than their dependencies, or explicitly forced. See
`go_accession_clustering.md` for the algorithm, output formats, CLI options,
dependency policy, and tests.

Clustering graph, Leiden, matrix, cache, and candidate-filter parameters live
in `cluster_config.yaml`, created from a packaged starter template when
missing. The same config also controls the stats output path and clustering
progress log interval. Because the config is a clustering input, editing it
invalidates the cluster outputs.

## Open Questions

- Which command should be implemented first after `setup-data`?
- What should be the first canonical signature library file format?
- Should command output default to human-readable text, structured JSON/TSV, or both?
