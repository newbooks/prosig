# Implementation Spec: `prosig inspect`

## Goal

Add a diagnostic command group for inspecting ProSig artifacts before later
GO similarity, GO clustering, motif discovery, and prediction workflows depend
on them.

`inspect` is not a replacement for production commands. It is a stable place to
answer questions such as:

- Does this GO term exist in the MF graph?
- What IC and ancestors does this GO term have?
- What Lin similarity and MICA are produced for two GO terms?
- Which accession or motif records will later contribute to a workflow?

## Command Group

```text
prosig inspect
```

The command group should use Typer and keep CLI code thin. Diagnostics should
call reusable package modules such as `prosig.go.similarity`.

## Implemented Commands

### `go-summary`

```text
prosig inspect go-summary --go-graph go_graph.pkl
```

Reports:

- artifact path
- namespace
- schema version
- number of terms
- number of terms with IC
- creation date when available

### `go-term`

```text
prosig inspect go-term GO:0005524 --go-graph go_graph.pkl
```

Reports one MF term:

- GO ID
- name
- frequency
- IC
- depth
- parents
- children

Options:

- `--ancestors`: include ancestors plus the term itself.
- `--json`: emit JSON diagnostic output.

### `go-sim`

```text
prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl
```

By default, reports only the Lin similarity score with up to four decimal
places. If the score is unavailable, reports `NA`.

```text
0.7314
```

Verbose mode:

```text
prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl --verbose
prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl -v
prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl -v --tree-style ascii
```

Reports Lin similarity diagnostics:

- queried GO IDs
- similarity or `None`
- input GO term descriptions
- compact GO path from the Molecular Function root through the selected
  shared path and out to each query term
- MICA
- IC values
- common ancestors with descriptions
- formula and numeric substitution used to derive the score
- status
- unavailable reason, when applicable
- common ancestors with valid IC

Options:

- `--json`: emit JSON diagnostic output.
- `--tree-style unicode|ascii`: choose Unicode or ASCII connectors for the
  compact GO path in verbose output. Unicode is the default.

## Planned Commands

Future diagnostic commands should be added under this command group:

- accession lookup against an accession-to-MF-GO artifact
- motif summary and motif lookup
- term-set similarity
- cluster summary and cluster member inspection
- explanation views connecting motif hits to predicted functions

Output should default to compact, tab-separated human-readable diagnostics, with
JSON available where nested data is useful.
