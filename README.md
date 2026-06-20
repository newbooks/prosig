# ProSig

ProSig: Protein Signature Discovery and Function Inference

## Command Plan

ProSig will expose a small command set that separates setup workflows from
routine analysis workflows:

- `prosig setup-data`: download and cache external data for offline use.
- `prosig build-library`: build the minimized GO graph, adjustable Leiden
  function clusters, and customizable motif library.
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

## Project Structure

- `AGENT.md`: working instructions for coding agents.
- `docs/todos/`: project TODOs and backlog.
- `docs/specs/`: implementation specs for motif discovery, function prediction, and related workflows.
- `docs/decisions/`: durable technical and scientific decision records.
- `src/prosig/`: package source code.
- `tests/`: test notes and future test suite.
- `data/`: local data workspace; large datasets should remain untracked.
- `notebooks/`: exploratory analysis.
- `scripts/`: operational helper scripts.

## Initial TODOs

- [ ] Fetch dependencies for offline use.
- [ ] Implement a STREME-like discriminative protein motif discovery module for ProSig using k-mer enumeration, Fisher exact enrichment, motif generalization, and optional PWM refinement.
