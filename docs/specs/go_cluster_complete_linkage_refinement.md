# Complete-Linkage GO Cluster Refinement

## Status

Implemented as the freshness-managed second clustering stage in
`prosig build-library`.

## Problem

The current workflow constructs a sparse GO-profile similarity graph and uses
Leiden community detection. An edge is admitted only when its AMB Lin
similarity is at least `min_similarity`, but that constraint applies only to
direct graph neighbors.

Community membership is transitive. A Leiden community may therefore contain a
chain such as:

```text
A -- B -- C -- D
```

where every retained edge meets `min_similarity`, while the unconnected pair
`A, D` has low or zero similarity. Large observed communities with
`sim_min = 0.0` demonstrate this behavior.

The downstream goal is stronger than graph community detection: every pair of
accessions in a final cluster should meet a configured similarity floor.

## Current Singleton Behavior

Current behavior distinguishes isolates from singleton communities:

- An accession with no graph edge meeting `min_similarity` is removed before
  Leiden and omitted from `leiden_clusters.tsv`.
- Leiden receives only accessions incident to at least one retained edge.
- If Leiden assigns one of those active accessions to a one-member community,
  that singleton is retained in `leiden_clusters.tsv`.
- `leiden_clusters_meta.tsv` reports retained singleton communities with
  `size = 1` and `NA` for `sim_ave`, `sim_min`, and `sim_max`.

The refinement revision should preserve this distinction. A singleton produced
by complete-linkage splitting remains a final cluster; an initial no-edge
isolate remains excluded. This avoids silently losing active accessions merely
because a strict compactness constraint separates them from their Leiden
community.

## Goals

- Preserve sparse kNN graph construction and Leiden as the scalable first
  stage.
- Split each Leiden community until every pair in every final non-singleton
  cluster has similarity greater than or equal to
  `min_cluster_similarity`.
- Prevent single-linkage-style chaining in final output.
- Keep deterministic cluster membership and stable cluster IDs.
- Retain current output formats, freshness behavior, and isolate policy.
- Report enough pre- and post-refinement statistics to tune thresholds.

## Non-Goals

- Do not replace Leiden with global hierarchical clustering.
- Do not require every accession to have `min_similarity` edges to every other
  member of its Leiden community.
- Do not merge accessions from different Leiden communities during refinement.
- Do not silently approximate complete linkage when a community is large.

## Defaults

```yaml
min_similarity: 0.50
min_cluster_similarity: 0.25
```

The two thresholds have different scopes:

- `min_similarity` is the minimum AMB Lin similarity for a sparse kNN graph
  edge. It controls graph quality and first-stage connectivity.
- `min_cluster_similarity` is the minimum AMB Lin similarity for every unique
  accession pair in a final cluster. It controls final cluster diameter.

`min_similarity = 0.50` is a reasonable initial edge threshold because it
removes weak graph links before community detection while retaining enough
connectivity for Leiden to identify broad candidate communities.

`min_cluster_similarity = 0.25` is a reasonable conservative initial refinement
threshold. It guarantees removal of zero- and very-low-similarity pairings
without immediately requiring every final pair to meet the stricter edge
threshold. It is intentionally permissive: a value of `0.25` guarantees
compactness only at that level and should not be interpreted as strong
functional equivalence.

It is valid for `min_cluster_similarity` to be lower than `min_similarity`.
The edge threshold controls which local relationships build the sparse graph;
the cluster threshold controls the maximum pairwise distance allowed after
Leiden. The values need not be ordered by validation.

The initial defaults should be evaluated empirically after implementation.
Compare:

- pre- and post-refinement cluster counts and size distributions
- number and fraction of final singletons
- final `sim_min` and `sim_ave` distributions
- number of Leiden communities split
- largest final cluster

If final clusters remain functionally broad, raise
`min_cluster_similarity`, with `0.35`, `0.40`, and `0.50` as useful calibration
points. If singleton production is excessive, lower it.

## Command-Line Configuration

The refinement threshold is intentionally independent of
`cluster_config.yaml`, which controls the persisted Leiden stage:

```bash
prosig build-library --min-cluster-similarity 0.25
```

Validation:

```text
0 < min_cluster_similarity <= 1
```

The threshold is inclusive. A pair with similarity exactly equal to
`min_cluster_similarity` may remain in the same cluster.

Complete-linkage refinement is freshness-managed by `build-library`. Changing
this option does not rebuild current Leiden artifacts, but it does invalidate
the refinement outputs because `clusters_stats.json` records
`min_cluster_similarity`.

Refinement writes `clusters.tsv`, `clusters_meta.tsv`, and `clusters_stats.json`
when any output is missing, when any dependency is newer, when `--force` is
used, or when existing stats were built with a different
`min_cluster_similarity`.

Refinement dependencies are:

- `go_graph.pkl`
- `accession_mf_go.tsv`
- `leiden_clusters.tsv`

## Revised Algorithm

### Stage 1: Sparse Graph and Leiden

Keep the current implementation:

1. Clean accession GO profiles.
2. Build candidate postings.
3. Build a sparse kNN graph using only edges with
   `similarity >= min_similarity`.
4. Exclude accessions incident to no retained edge.
5. Run weighted Leiden over active accessions.

The result is a set of coarse candidate communities.

### Stage 2: Complete-Linkage Refinement

Refine each Leiden community independently.

For a community with members `m[0:n]`:

1. If `n == 1`, retain it unchanged.
2. Compute AMB Lin similarity for every unique member pair.
3. Convert similarity to dissimilarity:

   ```text
   distance(a, b) = 1.0 - similarity(a, b)
   ```

   Clamp numerical similarity values to `[0.0, 1.0]`. Treat an unavailable
   pair score consistently with metadata generation as similarity `0.0`.

4. Store distances in SciPy condensed-distance order.
5. Run:

   ```python
   from scipy.cluster.hierarchy import fcluster, linkage

   hierarchy = linkage(distances, method="complete")
   labels = fcluster(
       hierarchy,
       t=1.0 - min_cluster_similarity,
       criterion="distance",
   )
   ```

6. Convert each flat label to a refined community.

Complete linkage defines distance between two candidate clusters as the maximum
distance between their members. Cutting at
`1 - min_cluster_similarity` therefore limits final cluster diameter and gives
the required postcondition:

```text
for every final cluster C:
    for every distinct a, b in C:
        similarity(a, b) >= min_cluster_similarity
```

Validate this postcondition while computing final metadata. Allow only a small
floating-point comparison tolerance, such as `1e-12`. A violation indicates an
implementation error and should fail rather than silently write an invalid
cluster.

### Stage 3: Stable Final IDs

Do not preserve preliminary Leiden IDs in the public output.

After all communities are refined:

1. Sort members within each final cluster.
2. Sort final clusters by first accession, then size, then full member tuple.
3. Assign `cluster_0001`, `cluster_0002`, and so on.

This preserves deterministic final IDs even if SciPy flat labels differ.

## Why Complete Linkage

Complete linkage directly represents the required all-pairs constraint because
its inter-cluster distance is the farthest member-pair distance.

Rejected alternatives:

- Single linkage permits chaining and reproduces the current failure mode.
- Average linkage can retain a low-similarity pair when the overall average is
  high.
- A medoid-radius constraint guarantees closeness to a representative but not
  between every member pair.
- Increasing Leiden resolution may reduce community sizes but does not enforce
  a pairwise similarity floor.

Leiden remains valuable as a scalable pre-partitioning stage. Running complete
linkage independently inside each Leiden community avoids a global all-pairs
comparison across Swiss-Prot and prevents refinement from joining unrelated
coarse communities.

## Complexity and Memory

For a Leiden community of size `n`, complete linkage requires:

```text
pair scores: n * (n - 1) / 2
time:        O(n^2)
memory:      O(n^2)
```

SciPy accepts a condensed vector of `n * (n - 1) / 2` distances and implements
complete linkage in `O(n^2)` time and memory.

Approximate condensed-vector sizes using `float64`:

```text
n = 3,020: 4,558,690 values, about 36.5 MB
n = 2,367: 2,800,161 values, about 22.4 MB
```

SciPy requires additional working memory, so these numbers are lower bounds.
Process Leiden communities sequentially and release each condensed matrix before
processing the next one. Do not retain matrices for all communities.

The current metadata writer already evaluates every pair in each final cluster.
The implementation should avoid unnecessary duplicate scoring where practical:

- use the existing profile-pair cache during condensed-matrix construction
- compute refined-cluster metadata from the same parent-community pair scores,
  or ensure the cache is reused if metadata is computed afterward

SciPy is a clustering dependency:

```toml
dependencies = [
  "scipy>=1.11",
]
```

Do not implement a silent approximate fallback for memory pressure because that
would weaken the all-pairs guarantee. Raise a clear error identifying the
community size and required pair count.

## Outputs

`clusters.tsv` contains final refined cluster IDs, not preliminary Leiden IDs.
The independently reusable first-stage membership remains in
`leiden_clusters.tsv`.

`clusters_meta.tsv` uses:

```text
cluster_id	sim_ave	sim_min	sim_max	size	composed_go
```

For every non-singleton row:

```text
sim_min >= min_cluster_similarity
```

Singleton rows retain `NA` similarity fields.

`clusters_stats.json` contains:

```json
{
  "algorithm": "go_set_similarity_knn_leiden_complete",
  "min_similarity": 0.5,
  "min_cluster_similarity": 0.25,
  "leiden_clusters": 12,
  "refined_clusters": 47,
  "leiden_clusters_split": 8,
  "refined_singletons": 9,
  "refinement_pairs_scored": 7358851,
  "pre_refinement_cluster_size_max": 3020,
  "post_refinement_cluster_size_max": 410
}
```

Keep the existing `clusters` and cluster-size summary fields, but define them as
final post-refinement values.

## Logging

Log at INFO level:

- start and completion of complete-linkage refinement
- current Leiden community number and size for long-running communities
- pair count and approximate condensed-vector size before allocation
- number of refined clusters produced from each split community
- total Leiden clusters, refined clusters, split communities, and final
  singletons
- pre- and post-refinement maximum cluster sizes
- configured `min_similarity` and `min_cluster_similarity`

Use the existing progress interval for pairwise matrix construction.

## API Changes

Extend `GoClusteringConfig` and `cluster_accessions_by_go`:

```python
min_cluster_similarity: float = 0.25
```

Recommended internal helpers:

```python
def refine_go_clusters_complete_linkage(
    *,
    active_accessions: list[str],
    cluster_by_accession: dict[str, str],
    accession_terms: dict[str, tuple[str, ...]],
    go_index: FastGoSimilarityIndex,
    min_cluster_similarity: float,
    profile_pair_cache: ProfilePairCache | None,
    lin_similarity_matrix: np.ndarray | None,
    progress_interval_seconds: float,
) -> dict[str, str]:
    ...


def complete_linkage_labels(
    condensed_distances: np.ndarray,
    *,
    member_count: int,
    min_cluster_similarity: float,
) -> np.ndarray:
    ...
```

The first helper groups by Leiden community, refines communities sequentially,
and returns deterministic final membership.

## Edge Cases

### Empty Active Graph

Write header-only cluster outputs and zero refinement counts.

### Leiden Singleton

Retain it unchanged. No pairwise matrix is needed.

### Refinement Singleton

Retain it as a final cluster with `NA` similarity metadata.

### Identical Profiles

All pair distances are zero and complete linkage keeps them together for every
valid threshold.

### Threshold of 1.0

Only members with pairwise similarity `1.0` may share a final cluster.

### Unavailable Similarity

Treat it as similarity `0.0`, matching current metadata behavior. Such a pair
cannot coexist in a final cluster for any valid positive threshold.

### Numerical Boundaries

Use double-precision condensed distances. Similarity exactly equal to the
threshold qualifies. Verify final `sim_min` with a `1e-12` tolerance.

## Tests

Add unit tests for:

- command-line default and range validation
- conversion from similarity to condensed distance
- a chain `A-B-C-D` that Leiden places together but complete linkage splits
- a cluster whose every pair equals the threshold and remains intact
- a pair immediately below the threshold that is split
- zero-similarity pairs never sharing a final cluster
- unavailable similarities treated as zero
- Leiden singletons retained
- refinement-produced singletons retained
- deterministic final IDs independent of SciPy label values
- every final non-singleton satisfying
  `sim_min >= min_cluster_similarity`
- stats containing pre- and post-refinement counts

Add integration tests for:

- `build-library` exposing the command-line default
- changing the threshold rerunning refinement without rebuilding Leiden
- metadata rows satisfying the configured minimum
- large synthetic communities being processed one at a time

## Acceptance Criteria

- Leiden configuration uses `min_similarity: 0.50`; the CLI default for
  `--min-cluster-similarity` is `0.25`.
- No final non-singleton cluster has `sim_min` below
  `min_cluster_similarity`, apart from the documented floating-point tolerance.
- Initial no-edge isolates remain excluded.
- Leiden and refinement-produced singletons remain represented.
- Final cluster IDs are deterministic.
- Existing cluster and metadata schemas remain consumable.
- Stats distinguish Leiden communities from final refined clusters.

## References

- SciPy complete-linkage definition and complexity:
  <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html>
- SciPy flat-cluster distance threshold:
  <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html>
- Leiden API and resolution behavior:
  <https://leidenalg.readthedocs.io/en/stable/reference.html>
