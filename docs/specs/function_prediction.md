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
cluster_frequency = (TP + 0.5) / (TP + FN + 1)
background_frequency = (FP + 0.5) / (FP + TN + 1)
weight = log2(cluster_frequency / background_frequency)
```

The `0.5` pseudocount is a Jeffreys prior for binary motif presence. It prevents
`FP = 0` from producing infinite weights while still allowing cluster-specific
motifs to receive large positive scores when support is strong.

The pickle stores only positive weights, nested by motif ID then cluster ID, so
prediction can retrieve all cluster evidence for each matched motif.

## Internal Calibration

At the end of `build-library`, the motif score board build performs an internal
calibration from the same motif-accession scan and final cluster assignments.
Calibration evaluates motif inference as cluster prediction at these minimum
weight thresholds:

```text
2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
```

Accessions in clusters smaller than `--motif-scoreboard-min-cluster-size` are
excluded from calibration. Accessions without cluster IDs, such as singletons
not present in `clusters.tsv`, are also excluded.

For each eligible accession and threshold:

1. collect all motif hits for the accession;
2. collect all clusters linked to those motifs with `weight >= threshold`;
3. rank predicted clusters by the maximum motif-cluster weight contributed by
   any hit motif, descending, with cluster ID as a deterministic tie-breaker;
4. mark the accession as covered if at least one cluster is predicted;
5. evaluate whether the true cluster is the top-1 prediction, in the top-3
   predictions, or anywhere in the predicted set.

Reported values:

```text
coverage = covered_accessions / eligible_accessions
top1_accuracy = top1_correct_accessions / covered_accessions
top3_accuracy = top3_correct_accessions / covered_accessions
set_accuracy = set_correct_accessions / covered_accessions
avg_predictions = total_predicted_clusters / covered_accessions
```

If no accessions are covered at a threshold, the accuracy fields are reported as
`null` in `motif_cluster_scoreboard_meta.json` and `NA` in logs. `set_accuracy`
treats the multi-cluster motif inference as successful when the true cluster is
present anywhere in the predicted cluster set. The same calibration records are
logged by `build-library` and written under `stats.calibration` in
`motif_cluster_scoreboard_meta.json`.

## Open Questions

- Should prediction use rule-based scoring first, supervised learning first, or both?
- How should conflicting signatures be resolved?
- Which metrics define success for the first benchmark?
