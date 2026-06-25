# Specs

Implementation specs live here. Keep each spec focused on one subsystem or workflow.

Current specs:

- `cli.md`: Typer-based command and subcommand interface.
- `setup_data.md`: `setup-data` external resource retrieval command.
- `go_pkl.md`: minimal GO Molecular Function graph and IC artifact.
- `go_similarity.md`: Molecular Function Lin similarity over `go_graph.pkl`.
- `go_set_similarity.md`: AMB similarity between Molecular Function GO term
  sets and accession-resolved GO profiles.
- `go_accession_clustering.md`: Leiden clustering over a sparse GO-set
  similarity kNN graph, integrated into `build-library` with freshness checks.
- `go_cluster_complete_linkage_refinement.md`: complete-linkage second-stage
  refinement that enforces a final all-pairs similarity floor within each
  Leiden community.
- `go_mf_natural_language_composer.md`: rule-based natural-language summaries
  from Molecular Function GO term sets.
- `inspect_cli.md`: diagnostic `prosig inspect` command group.
- `motif_discovery.md`: protein signature discovery workflow.
- `prosig_motif_implementation.md`: ProSig motif format, PROSITE translation,
  and scanning behavior to preserve or port.
- `function_prediction.md`: function prediction from signatures.

Recommended structure:

- Problem and goals
- Inputs and outputs
- Data model
- Algorithm or workflow
- Edge cases
- Tests and validation
- Open questions
