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

GO accession clustering builds a sparse GO-set similarity k-nearest-neighbor
graph, writes freshness-managed Leiden artifacts (`leiden_clusters.tsv` and
`leiden_clusters_meta.tsv`), then freshness-manages complete-linkage refinement
of those communities. Final outputs are `clusters.tsv` and
`clusters_meta.tsv`. Use `--min-cluster-similarity` to set the required
all-pairs similarity floor; the default is `0.25`.

After final clustering, `build-library` also synthesizes cluster-level GO MF
signatures from `clusters.tsv`, `accession_mf_go.tsv`, and `go_graph.pkl`.
The `clusters_meta.tsv` output includes a `composed_go` column with up to 10
semicolon-separated representative GO terms per cluster. Terms are selected by
`support × IC`, but only GO IDs are written. Each accession contributes one
equal vote per propagated GO term; no additional evidence-code weighting is
applied because the retained Swiss-Prot annotations have already been filtered
during library construction.

`build-library` scans `prosig_motifs.tsv` against `accession.fasta` for final
cluster members and writes the sparse motif hit table `motif_features.tsv`.
The scan uses 8 worker processes by default; override with
`--motif-scan-processes`. It then builds a pickled motif-cluster prediction
score board. The score board ignores clusters with fewer than 10 members,
ignores motif-cluster pairs with support below 5, and stores only positive
motif-cluster weights. Weights are log2 enrichment scores computed with a
Jeffreys prior pseudocount of 0.5, so zero-background hits remain finite and
support-sensitive. Metadata is written alongside the pickle with counts for
ignored combinations, stored weights, and internal calibration at motif-weight
thresholds 2.0 through 8.0, including top-1, top-3, set accuracy, average
prediction count, and coverage.

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
prosig inspect function cluster_0008 --go-graph go_graph.pkl --cluster-meta clusters_meta.tsv
```

`build-library` also derives `accession.fasta` and
`accession.fasta.idx` from the configured `uniprot_sprot.dat.gz`. FASTA labels
are primary Swiss-Prot accessions, and the index stores byte offsets for fast
internal sequence lookup. The pair is rebuilt when either artifact is missing
or older than the Swiss-Prot source. The index records FASTA size and
modification time so edited FASTA files are rejected until the pair is rebuilt.

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
