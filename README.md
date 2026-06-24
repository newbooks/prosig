# ProSig

ProSig: Protein Signature Discovery and Function Inference

## Command Plan

ProSig will expose a small command set that separates setup workflows from
routine analysis workflows:

- `prosig setup-data`: download and cache external data for offline use.
- `prosig build-library`: build the minimized GO graph, adjustable Leiden
  function clusters, and customizable motif library.
- `prosig inspect`: inspect ProSig artifacts and diagnostic calculations,
  including GO terms and Lin similarity scores.
- `prosig discover`: discover discriminative motifs from grouped function
  clusters and background sequences.
- `prosig annotate`: scan sequence(s), report motif hits, and infer sequence
  function from those motif hits as prediction evidence.

`setup-data` and `build-library` are expected to be run less often than
`discover` and `annotate`. Clustering is treated as part of library construction
rather than a separate top-level command, because function clusters are a
prerequisite for the motif library. Function prediction is treated as part of
annotation, because predictions should be reported together with the motif scan
hits that justify them.

GO evidence-code filtering in `build-library` is intended to be used with the
reviewed Swiss-Prot accession file `uniprot_sprot.dat.gz`. The excluded GO
assignment evidence codes are maintained in `src/prosig/go/build.py`; applying
the same exclusion rule to unreviewed annotation sources requires a separate
review.

`build-library` also writes `accession_mf_go.tsv`, a headerless two-column TSV
mapping primary Swiss-Prot accessions to semicolon-separated high-quality direct
MF GO terms for later diagnostics and GO set similarity.

`build-library` skips derived artifacts that are newer than their dependencies.
Use `--force` or `-f` to rebuild them anyway.

Planned GO accession clustering for `build-library` is specified in
`docs/specs/go_accession_clustering.md`. It builds a sparse GO-set similarity
k-nearest-neighbor graph from `go_graph.pkl` and `accession_mf_go.tsv`, runs
Leiden community detection, and writes `go_clusters.tsv` only when the cluster
artifact is missing or older than its dependencies.

## Diagnostic Inspection

`prosig inspect` is a diagnostic command group for checking intermediate
artifacts before clustering, motif discovery, or prediction work depends on
them. The initial GO commands use only the Molecular Function namespace and only
the Lin semantic similarity score:

```text
prosig inspect go-summary --go-graph go_graph.pkl
prosig inspect go-term GO:0005524 --go-graph go_graph.pkl --ancestors
prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl
prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl --verbose
prosig inspect go-sim GO:0005524 GO:0004672 --go-graph go_graph.pkl -v --tree-style ascii
prosig inspect go-set-sim "(GO:0005524;GO:0004672)" Q9SVY5 --go-graph go_graph.pkl --accession-go accession_mf_go.tsv
prosig inspect go-set-sim "GO:0005524;GO:0004672" Q9SVY5 --go-graph go_graph.pkl --accession-go accession_mf_go.tsv
```

The inspect surface is intended to grow with artifact diagnostics for
accessions, motifs, clustering inputs, and standalone similarity calculations.

## Project Structure

- `AGENT.md`: working instructions for coding agents.
- `docs/todos/`: project TODOs and backlog.
- `docs/prosig_motifs.md`: user guide for ProSig motif syntax and motif libraries.
- `docs/specs/`: implementation specs for motif discovery, ProSig motif handling, function prediction, and related workflows.
- `docs/decisions/`: durable technical and scientific decision records.
- `src/prosig/`: package source code.
- `tests/`: test notes and future test suite.
- `data/`: local data workspace; large datasets should remain untracked.
- `notebooks/`: exploratory analysis.
- `scripts/`: operational helper scripts.

## Initial TODOs

- [ ] Fetch dependencies for offline use.
- [ ] Implement a STREME-like discriminative protein motif discovery module for ProSig using k-mer enumeration, Fisher exact enrichment, motif generalization, and optional PWM refinement.
