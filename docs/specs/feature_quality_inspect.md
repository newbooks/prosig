# Implementation Spec: `prosig inspect feature`

## Goal

Add a diagnostic command that evaluates generic per-member feature values
against GO-derived ProSig functional clusters, writes feature quality scores,
and renders score card plots.

The command should use GO-based functional similarity only. Do not include the
old EC-prefix similarity fallback from PClass.

## Command

```text
prosig inspect feature
```

Default invocation:

```bash
prosig inspect feature
```

is equivalent to:

```bash
prosig inspect feature \
  --feature-file feature_values.tsv \
  --output-file feature_scores.tsv
```

## CLI Options

| Option | Default | Description |
|---|---:|---|
| `--feature-file PATH` | `feature_values.tsv` | Input TSV containing member IDs, cluster IDs, and numeric feature columns. |
| `--output-file PATH` | `feature_scores.tsv` | Output TSV containing feature quality metric scores. |
| `--go-graph PATH` | resolved runtime library | Optional compact GO graph pickle. |
| `--accession-go PATH` | resolved runtime library | Optional accession-to-MF-GO TSV. |
| `--library-dir PATH` | existing inspect behavior | Runtime library directory fallback for GO artifacts. |
| `--min-cluster-size INT` | `10` | Minimum number of members required for a cluster to be evaluated. This has a hard lower limit of `10`; smaller requested values are rejected. Clusters smaller than the accepted threshold are skipped and logged. |
| `--format png\|svg` | `png` | Score card image format. |
| `--feature NAME` | none | Optional repeated feature selector for score card rendering. |

The command should print the paths it writes:

```text
Wrote feature scores: feature_scores.tsv
Wrote score cards to: score_cards
```

## Input Format

`feature_values.tsv` is a tab-separated file with required columns:

```text
member_id	cluster_id	feature_a	feature_b	...
P1	cluster_0001	0.12	5.0
P2	cluster_0001	0.18	4.7
Q1	cluster_0002	1.30	2.1
```

Rules:

- `member_id` is an accession or member identifier.
- `cluster_id` is a GO-derived ProSig cluster identifier, such as `cluster_0001`.
- Every column other than `member_id` and `cluster_id` is treated as one numeric feature.
- At least one feature column is required.
- Feature columns must contain numeric values in every cell. Blank cells are not
  allowed.
- Raise a user-facing error if any feature cell is blank, non-numeric, or
  non-finite.
- Cluster size is the number of unique `member_id` values for each `cluster_id`.
- Clusters with fewer than `--min-cluster-size` unique members are skipped before evaluation.
- Skipped clusters must be logged at info level with the threshold and skipped
  cluster count. Debug-level logs must include skipped cluster IDs and their
  unique-member counts.
- At least two clusters must remain after filtering. If fewer than two clusters
  remain, raise a user-facing error.
- Output feature order follows input column order.

## Score Output

Write `feature_scores.tsv` with exactly these columns:

```text
feature	compactness	separation	gradient	specificity
```

One row is written per feature column from the input. Float values should be
written with four decimals.

Example:

```text
feature	compactness	separation	gradient	specificity
feature_a	1.0000	0.8000	0.6000	0.4000
feature_b	0.7142	0.5000	0.3000	0.1250
```

## GO Functional Similarity

Use GO-based similarity for all clusters.

By default, both `go_graph.pkl` and `accession_mf_go.tsv` are resolved through
the existing runtime library protocol: use a complete runtime library in the
working directory when present, otherwise use `--library-dir` when provided,
otherwise use the package-shipped runtime library. The runtime library is
all-or-nothing; a partial local library should fail rather than mixing local and
packaged artifacts. `--go-graph` and `--accession-go` may be used as explicit
per-file overrides.

The accession GO artifact must be a tab-separated file that can be read by
ProSig's `load_accession_mf_go_terms` helper and must provide MF GO terms for
retained `member_id` values.

For each retained cluster:

1. Resolve each `member_id` to its MF GO terms from `accession_mf_go.tsv`.
2. Fail clearly if any retained member has no GO terms.
3. Choose one centroid accession per cluster:
   - For each accession in the cluster, compute mean AMB GO-set similarity to
     every other accession in the same cluster.
   - Select the accession with the highest mean similarity.
   - Tie-break by choosing the lexicographically smallest accession string.
4. Define cluster-cluster functional similarity as AMB GO-set similarity between
   the centroid accessions for the two clusters.
5. If a GO-set similarity is unavailable, use `0.0`.
6. Clamp usable similarity values to `[0, 1]`.

The implementation should use ProSig's current GO APIs, including
`GoSimilarity`, `load_accession_mf_go_terms`, and existing AMB GO-set
similarity behavior.

## Metric Algorithms

Use the same four scalar feature quality algorithms previously used in PClass.

### Compactness

Compactness rewards low within-cluster variance relative to global variance.

For one feature:

1. Use the dense numeric feature values validated during input loading.
2. Keep clusters with enough values for this metric.
3. Compute global sample variance across all retained values.
4. For each cluster:

```text
cluster_score = 1 - cluster_sample_variance / global_sample_variance
```

5. Clip each cluster score to `[0, 1]`.
6. Return the median cluster score.
7. If the global feature variance is effectively zero, return `0.0`.

### Separation

Separation measures the fraction of cluster pairs with separated means.

For each cluster pair:

```text
abs(mean_i - mean_j) > std_i + std_j
```

The score is:

```text
separable_pairs / total_pairs
```

The threshold is strict: equality is not separable.

### Gradient

Gradient measures whether feature distance increases as GO functional distance
increases.

For each cluster pair:

```text
functional_distance = 1 - go_similarity(cluster_i, cluster_j)
feature_distance = abs(mean_i - mean_j)
```

Compute Spearman correlation between functional distances and feature
distances. Negative correlations are clipped to `0.0`.

If all functional distances are identical, return `NaN`. If all feature
distances are identical, return `0.0`.

The evaluator should normally avoid `NaN` scores by requiring at least two
qualified clusters and by handling constant feature-distance cases explicitly.
If a metric still returns `NaN`, preserve it in `feature_scores.tsv`.

### Specificity

Specificity rewards features with a strong local gap between at least one
cluster and its nearest feature-space neighbors.

For each cluster pair:

```text
gap = max(0, abs(mean_i - mean_j) - (std_i + std_j))
```

Neighborhood size defaults to:

```text
round(0.05 * cluster_count)
```

clamped to at least `1` and at most `cluster_count - 1`.

The score is:

```text
specificity = specificity_raw / (specificity_raw + global_std)
```

If global standard deviation is effectively zero, return `0.0`.

## Score Cards

After writing `feature_scores.tsv`, render score cards into `score_cards/`.

Metric order is fixed:

```python
["compactness", "separation", "gradient", "specificity"]
```

Behavior:

- If no `--feature` option is passed, write one score card per feature:

```text
score_cards/<feature>_score_card.png
```

- If one or more `--feature` values are passed, render selected features
  together:

```text
score_cards/selected_features_score_card.png
```

- The score card output directory is always `score_cards/`.
- Existing score card files may be overwritten, including when different
  feature names slugify to the same filename stem.
- Plot values should be clipped to `[0, 1]`.
- If a metric value is `NaN`, plot it as `0.0`.
- Display raw metric values formatted to four decimals.
- Use the existing PClass score card style unless ProSig already has a newer
  plotting convention.

Feature names used in generated filenames should be slugified by replacing
characters outside `[A-Za-z0-9_.-]` with underscores.

## Recommended Code Structure

Add reusable modules instead of putting metric logic in the CLI layer:

```text
src/prosig/features/__init__.py
src/prosig/features/quality.py
src/prosig/features/evaluation.py
src/prosig/features/score_card.py
```

Wire the command in:

```text
src/prosig/cli/inspect.py
```

The CLI should stay thin:

1. Resolve input, output, GO graph, and accession GO paths.
2. Call `evaluate_feature_file(...)`.
3. Call `write_feature_quality_score_cards(...)`.
4. Print written output paths.

The evaluation step must apply `--min-cluster-size` before any metric
calculation. The default threshold is `10`, and `10` is a hard lower limit.
Values below `10` must be rejected. Only clusters with at least the accepted
number of unique `member_id` values contribute to compactness, separation,
gradient, specificity, GO centroid selection, and score card outputs. Clusters
below the threshold are not errors; they are skipped and logged. Evaluation
requires at least two qualified clusters after filtering.

## Error Handling

Report user-facing errors for:

- Missing `feature_values.tsv`.
- Missing required columns: `member_id`, `cluster_id`.
- No feature columns.
- Blank, non-numeric, or non-finite feature value.
- `--min-cluster-size < 10`.
- Fewer than two clusters remain after cluster-size filtering.
- Missing GO graph.
- Missing accession GO file.
- Retained members missing MF GO terms.
- GO similarity outside `[0, 1]`.
- Selected score-card feature not present in `feature_scores.tsv`.
- Unsupported score card format.

## Tests

Add test coverage for:

- Default command reads `feature_values.tsv`, writes `feature_scores.tsv`, and
  creates `score_cards/`.
- Input validation for required columns.
- Input validation for missing feature columns.
- Input validation for blank, non-numeric, and non-finite feature values.
- Error when `--min-cluster-size` is smaller than `10`.
- Cluster-size filtering uses unique `member_id` counts, not row counts.
- Error when fewer than two clusters remain after filtering.
- Centroid tie-breaking chooses the lexicographically smallest accession.
- GO-only centroid similarity path.
- Default GO graph and accession GO inputs are resolved from the same runtime
  library protocol.
- `NaN` metric values plot as `0.0`.
- Output score columns and deterministic feature order.
- Score card files are created as PNG by default.
- Selected features render into one combined score card.
- No EC-specific fallback is used.
