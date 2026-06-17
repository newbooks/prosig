# ProSig

ProSig: Protein Signature Discovery and Function Inference

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
