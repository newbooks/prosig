# Function Prediction Spec

## Goal

Predict protein function from discovered or curated sequence signatures.

## Inputs

- Protein sequences.
- Signature library.
- Signature-to-function annotations.
- Optional taxonomic, family, or domain metadata.

## Outputs

- Predicted functions.
- Supporting signature hits.
- Confidence or ranking scores.
- Explanation metadata suitable for inspection.

## Initial Workflow Notes

- Scan sequences against a signature library.
- Aggregate signature hits into function evidence.
- Rank candidate functions by evidence strength.
- Report predictions with traceable signature support.

## Motif-Cluster Score Board

The baseline motif-based prediction model uses a pickled motif-cluster score
board built from sparse motif hits and final functional clusters.

Inputs:

- `clusters.tsv`: final `member_id` / `cluster_id` table.
- `prosig_motifs.tsv`: motif library; the `prosig_pattern` column is compiled
  and scanned.
- `accession.fasta` and `accession.fasta.idx`: sequence artifacts.
- `motif_features.tsv`: generated sparse binary motif hit table with
  `accession` and `motif_id` columns. Row presence means motif presence.

Outputs:

- `motif_features.tsv`: sparse positive accession-motif hit features.
- `motif_cluster_scoreboard.pkl`: machine-readable positive motif-cluster
  weights.
- `motif_cluster_scoreboard_meta.json`: build statistics and ignored-combo
  counts.

Filters:

- ignore clusters with fewer than 10 members;
- ignore motif-cluster pairs with `TP < 5`;
- skip zero or negative motif-cluster weights in the pickle artifact.

Weight:

```text
cluster_frequency = TP / (TP + FN)
background_frequency = FP / (FP + TN)
weight = log2(cluster_frequency / background_frequency)
```

The pickle stores only positive weights, nested by motif ID then cluster ID, so
prediction can retrieve all cluster evidence for each matched motif.

## Open Questions

- Should prediction use rule-based scoring first, supervised learning first, or both?
- How should conflicting signatures be resolved?
- Which metrics define success for the first benchmark?
