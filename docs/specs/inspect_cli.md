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

### `go-set-sim`

```text
prosig inspect go-set-sim "(GO:0005524;GO:0004672)" "(GO:0005515)"
prosig inspect go-set-sim "GO:0005524;GO:0004672" "GO:0005515"
prosig inspect go-set-sim GO:0005524,GO:0004672 GO:0005515
prosig inspect go-set-sim P00533 Q9SVY5 --accession-go accession_mf_go.tsv
prosig inspect go-set-sim "(GO:0005524;GO:0004672)" Q9SVY5 --accession-go accession_mf_go.tsv
```

This command computes AMB similarity between two MF GO term sets. Direct sets
use semicolon- or comma-separated GO IDs, with optional parentheses. Accession
inputs resolve through a separate accession-to-MF-GO TSV artifact. Mixed
direct-set/accession inputs are allowed. See `go_set_similarity.md`.

Shell note: direct GO sets must be quoted or escaped in common shells because
parentheses and semicolons are shell syntax. Comma-separated sets without
parentheses usually do not require quotes.

Verbose output includes query labels, accession-expanded GO sets inline, one GO
term description section, A-to-B and B-to-A best-match rows using GO IDs only,
directional means, and the AMB formula. Best-match edge scores are shown with
four decimal places:

```text
B query: A0A024B7W1 (GO:0003724;GO:0003725;GO:0003968)
```

```text
GO:0005524 --0.2211--> GO:0004672
```

### `function`

```text
prosig inspect function P00533 --accession-go accession_mf_go.tsv
prosig inspect function "GO:0004672;GO:0005524"
prosig inspect function cluster_0008 --cluster-meta clusters_meta.tsv
```

This command composes a concise Molecular Function description from one of
three query forms:

- accession: resolved through `accession_mf_go.tsv`
- direct GO set: semicolon- or comma-separated GO IDs
- cluster ID: `cluster_` followed by digits, resolved through the
  `composed_go` column in `clusters_meta.tsv`

Cluster metadata defaults to `clusters_meta.tsv` in the working directory.

### `cluster`

```text
prosig inspect cluster cluster_0008 \
  --cluster-meta clusters_meta.tsv \
  --motif-scoreboard motif_cluster_scoreboard.pkl \
  --motif-library prosig_motifs.tsv
```

This command reports one functional cluster and the positive motif-cluster
weights that identify it. The default output is a readable report; use `--json`
for structured output.

Reports:

- cluster ID
- cluster size
- synthetic GO terms from `composed_go`
- composed description, using `composed_description` when available or deriving
  one from `composed_go` and `go_graph.pkl`
- motifs with positive stored weights for the cluster, including motif ID,
  motif description, TP/FP/FN/TN, odds ratio, and weight

The motif truth table is rendered as in-cluster versus outside-cluster counts:

```text
              In cluster_0008      Outside cluster_0008
-------------------------------------------------------
Motif present              TP                         FP
Motif absent               FN                         TN
```

`TP + FN` should equal the reported cluster size. If not, the command prints a
note because the scoreboard was likely built from different cluster artifacts
than the current `clusters_meta.tsv`.

Options:

- `--json`: emit JSON diagnostic output.

## Planned Commands

Future diagnostic commands should be added under this command group:

- accession lookup against an accession-to-MF-GO artifact
- motif summary and motif lookup
- cluster member inspection
- explanation views connecting motif hits to predicted functions

Output should default to compact, tab-separated human-readable diagnostics, with
JSON available where nested data is useful.
