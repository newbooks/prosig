# Standalone Swiss-Prot GO Accession Clustering

## Status

Implementation specification only. No implementation is included with this
document.

## Current Behavior

`clusters.tsv` is currently produced by the integrated `prosig build-library`
workflow in `src/prosig/cli/build_library.py`:

1. `write_accession_mf_go_tsv` extracts primary Swiss-Prot accessions and their
   high-quality Molecular Function (MF) GO annotations into
   `accession_mf_go.tsv`.
2. `cluster_accessions_by_go` builds a sparse AMB Lin GO-similarity kNN graph,
   runs Leiden, and writes the intermediate `leiden_clusters.tsv`.
3. `refine_go_clusters_complete_linkage` refines each Leiden community and
   writes `clusters.tsv`, `clusters_meta.tsv`, and `clusters_stats.json`.

The current refinement guarantees the configured all-pairs similarity floor
inside each output cluster. Its default is `0.25`, not `0.8`. It does not apply
a minimum final cluster size, write an orphan table, or merge compatible
clusters across Leiden community boundaries.

## Goal

Add an operational script that clusters Swiss-Prot primary accessions by MF GO
set similarity without running unrelated motif-library build stages. Retain
Leiden for scalable coarse community detection and complete linkage for strict
within-cluster refinement.

Proposed entry point:

```text
scripts/cluster_swissprot_by_go.py
```

This remains a script rather than a new top-level `prosig cluster` command,
consistent with the current product command plan. The core orchestration and
validation should live in `prosig.go.clustering` so the script is a thin CLI
wrapper and behavior is directly testable.

## Required Guarantees and Interpretation

Let `S(a, b)` be the existing AMB Lin similarity between the valid MF GO
profiles of accessions `a` and `b`, and let `t` be `--similarity-threshold`.

### 1. Within-cluster floor

For every retained cluster `C`:

```text
for every distinct a, b in C: S(a, b) >= t
```

The default is `t = 0.8`. The comparison is inclusive, with a numerical
tolerance of `1e-12`. Missing/unavailable pair scores count as `0.0`.

### 2. Across-cluster separation

For every pair of distinct retained clusters `C1`, `C2`, every cross-cluster
accession pair must have similarity strictly below the threshold:

```text
for every distinct retained C1, C2:
    for every a in C1 and b in C2: S(a, b) < t
```

Thus, an accession pair scoring exactly `t` may occur within one retained
cluster but may not occur across two retained clusters. Missing/unavailable
pair scores count as `0.0` and therefore satisfy cross-cluster separation.

Guarantees 1 and 2 together are not always feasible for every input accession
because thresholded GO similarity need not be transitive. For example, if
`S(A,B) >= t`, `S(B,C) >= t`, and `S(A,C) < t`, placing all three together
violates guarantee 1, while separating either above-threshold pair across
clusters violates guarantee 2. At least one conflicting accession must be
reported as an orphan. The algorithm must never weaken either guarantee merely
to retain more accessions.

### 3. Minimum retained size

Only clusters containing at least `--min-cluster-size` accessions are written
to `clusters.tsv`. Default: `10`. Accessions in smaller candidate clusters are
orphans; small clusters must not be silently discarded.

### 4. Complete orphan accounting

Every unique input accession must occur in exactly one of these states:

- one retained cluster in `clusters.tsv`; or
- `orphan_accessions.tsv`.

This includes accessions rejected during GO-profile validation, accessions
with no qualifying graph edge, complete-linkage singletons, and members of
clusters below the minimum size.

## Inputs

Required:

```text
go_graph.pkl
accession_mf_go.tsv
```

`go_graph.pkl` is the compact MF graph and information-content artifact made by
`prosig build-library`. `accession_mf_go.tsv` is the existing headerless input:

```text
P00533\tGO:0004672;GO:0005524
Q9SVY5\tGO:0000002
```

The script clusters the accessions present in this file; it does not download
Swiss-Prot or rebuild GO artifacts. A user starting from `uniprot_sprot.dat.gz`
must first run the relevant `prosig build-library` preparation or
`prosig setup-data` workflow.

Input validation:

- accession IDs must be non-empty;
- duplicate accession rows are combined into one sorted, deduplicated GO set;
- GO tokens must use `GO:NNNNNNN` syntax;
- alternate/obsolete IDs are handled exactly as in the existing GO artifact
  and similarity loader;
- non-MF or unknown GO terms are excluded from the valid scoring profile but
  retained for orphan reporting as input annotations;
- malformed rows fail with a path and line-number diagnostic.

## Command-Line Interface

Proposed usage:

```bash
python scripts/cluster_swissprot_by_go.py \
  --accession-go accession_mf_go.tsv \
  --go-graph go_graph.pkl \
  --clusters-out clusters.tsv \
  --orphans-out orphan_accessions.tsv \
  --similarity-threshold 0.8 \
  --min-cluster-size 10
```

Options:

```text
--accession-go PATH
    Input accession/MF-GO TSV. Default: accession_mf_go.tsv.

--go-graph PATH
    ProSig GO graph and IC pickle. Default: go_graph.pkl.

--clusters-out PATH
    Final retained cluster membership. Default: clusters.tsv.

--orphans-out PATH
    Accessions not in retained clusters. Default: orphan_accessions.tsv.

--similarity-threshold FLOAT
    Inclusive all-pairs AMB Lin similarity floor. Default: 0.8.
    Validation: 0 < value <= 1.

--min-cluster-size INTEGER
    Minimum number of accessions in a retained cluster. Default: 10.
    Validation: value >= 1.

--leiden-clusters-out PATH
    Optional persisted coarse Leiden memberships. Default:
    leiden_clusters.tsv beside --clusters-out.

--meta-out PATH
    Final cluster metrics. Default: <clusters stem>_meta.tsv.

--stats-out PATH
    Run statistics and provenance. Default: <clusters stem>_stats.json.

--cluster-config PATH
    Existing Leiden/kNN tuning config. Default: cluster_config.yaml.

--force
    Replace output artifacts even when they are current.
```

The script must refuse output-path collisions, create parent directories, and
write each artifact atomically through a temporary file followed by replace.

## Outputs

### `clusters.tsv`

Retain the existing downstream-compatible schema:

```text
member_id\tcluster_id
P00533\tcluster_0001
Q9SVY5\tcluster_0001
```

Rules:

- one row per retained accession;
- rows sorted by accession ID;
- members sorted within cluster before ID assignment;
- final clusters sorted by first accession, then size, then full member tuple;
- IDs assigned as `cluster_0001`, `cluster_0002`, ...;
- every cluster has at least `min_cluster_size` members;
- every non-singleton cluster has `sim_min >= similarity_threshold`.

### `orphan_accessions.tsv`

Required schema:

```text
accession_id\tgo_term
P12345\tGO:0004672
P12345\tGO:0005524
Q99999\t
```

Use long form: one row per orphan accession/GO-term pair. This preserves the
requested two-column schema for accessions with multiple GO annotations. Write
one row with an empty `go_term` when an accession has no input GO annotation.
Sort by `accession_id`, then `go_term`, and deduplicate rows.

The file includes all original input GO terms, not only the terms retained in
the valid MF scoring profile. This makes filtering auditable.

### Metadata and statistics

`clusters_meta.tsv` retains the current columns:

```text
cluster_id\tsim_ave\tsim_min\tsim_max\tsize\tcomposed_go
```

`clusters_stats.json` adds:

- `similarity_threshold`;
- `min_cluster_size`;
- input, valid-profile, active-graph, pre-size-filter, retained, and orphan
  accession counts;
- Leiden, refined, reconciled, retained, and undersized cluster counts;
- counts of cross-Leiden merges, threshold-conflict removals, and
  reconciliation pair scores;
- orphan reason counts;
- paths and SHA-256 digests for inputs;
- algorithm/library versions and Leiden seed;
- parameters used from `cluster_config.yaml`.

Orphan reasons belong in statistics, not the required orphan TSV schema.
Internally assign exactly one precedence-ordered reason:

1. `no_valid_mf_go_profile`;
2. `no_qualifying_similarity_edge`;
3. `cross_cluster_threshold_conflict`;
4. `below_min_cluster_size`.

## Algorithm

### Stage 1: Load and validate profiles

Load both raw annotations for reporting and valid MF profiles for scoring.
Build the existing fast GO similarity index and precomputed Lin term matrix.

### Stage 2: Leiden coarse partition

Reuse `cluster_accessions_by_go`:

1. construct the filtered candidate index;
2. construct the sparse weighted kNN graph;
3. keep edges meeting the config's graph `min_similarity`;
4. run deterministic Leiden with seed `0`.

The graph edge threshold remains independent of the final
`similarity_threshold`. To avoid excluding pairs needed by the final `0.8`
clustering, validate that graph `min_similarity <= similarity_threshold`.
Candidate-search recall remains a scalability heuristic and must be reported;
Leiden alone cannot prove global cross-cluster separation.

### Stage 3: Complete-linkage refinement

Reuse the existing complete-linkage implementation inside each Leiden
community, but pass `similarity_threshold` (default `0.8`) rather than the
integrated build's current default `0.25`.

Validate every refined cluster's all-pairs postcondition before proceeding.

### Stage 4: Global threshold reconciliation

The current implementation stops at Leiden boundaries, so it cannot enforce
guarantee 2. Reconcile the refined clusters globally.

Conceptually, form a threshold graph whose vertices are active accessions and
whose edges are pairs with `S(a, b) >= similarity_threshold`. Guarantees 1 and
2 mean that every connected component retained from this graph must be a
clique. Turning an arbitrary graph into disjoint cliques by removing the fewest
vertices is the cluster-vertex-deletion problem; exact maximum-retention
optimization is not required by this spec. The implementation must use the
following deterministic heuristic and report its removals as orphans:

1. Score cross-cluster accession pairs. The candidate index may determine the
   order in which pairs are examined, but it may not be used to assume that an
   unexamined pair is below the threshold.
2. If two clusters have every within- and cross-cluster pair at or above the
   threshold, merge them. Choose the eligible union with the highest minimum
   cross-pair similarity; break ties using the sorted full member tuples.
3. Otherwise, identify all cross-cluster violation edges, meaning accession
   pairs in different clusters with similarity at or above the threshold.
4. Orphan the accession incident to the largest number of violation edges.
   Break ties by lower average similarity to its current cluster, then by
   lexicographically larger accession ID. Removing an accession cannot violate
   the existing within-cluster floor.
5. Remove empty clusters and repeat merging and conflict removal until no
   cross-cluster violation edge remains.
6. Exhaustively audit every pair of accessions in different surviving clusters
   and fail before writing outputs if any similarity is at or above the
   threshold.

The exhaustive scoring and audit can be expensive, up to quadratic in the
number of active accessions. Log progress at the configured interval and record
the number of accession-pair scores. Do not silently skip or sample the audit.
Future approximate modes require a separate option and must not claim
guarantee 2.

### Stage 5: Minimum-size filter and orphan assignment

Apply `min_cluster_size` only after reconciliation. Retain qualifying clusters;
send every member of an undersized cluster to the orphan set. Do not attempt to
attach an orphan to a retained cluster unless the entire union satisfies the
threshold; such compatible unions should already have been found during
reconciliation.

After the size filter, re-audit both the within-cluster floor and strict
cross-cluster separation among retained clusters. Removing undersized clusters
cannot introduce a violation, but the audit protects the public postconditions
against implementation regressions.

### Stage 6: Stable IDs and atomic writes

Assign final IDs only after all merging and filtering. Validate complete input
partitioning and all guarantees before atomically replacing any public output.

## Freshness and Relationship to `build-library`

The standalone script reuses integrated clustering code but has independent
defaults and outputs. It must not invoke the whole `build-library` command.

Freshness dependencies:

- `go_graph.pkl`;
- `accession_mf_go.tsv`;
- `cluster_config.yaml`;
- all CLI parameters affecting results.

An output set is current only when all expected outputs exist and the stats
file records matching input digests and parameters. A partial output set is
rebuilt as a unit. `--force` always rebuilds.

Do not change `prosig build-library` defaults as part of this work. A later,
explicit decision may align its `--min-cluster-similarity` default with `0.8`
or add minimum-size/orphan outputs to the integrated workflow.

## Implementation Changes

Expected files:

- add `scripts/cluster_swissprot_by_go.py` as a thin argument parser and logger
  setup;
- extend `src/prosig/go/clustering.py` with a reusable orchestration result,
  global threshold reconciliation, deterministic conflict removal, orphan
  classification, audits, and atomic writers;
- add focused tests in `tests/test_go_clustering.py`;
- add script CLI tests in
  `tests/test_cluster_swissprot_by_go_script.py`;
- document usage in `scripts/README.md` after implementation.

Prefer extracting shared helpers over changing the semantics of
`cluster_accessions_by_go` or `refine_go_clusters_complete_linkage`, because
the existing `build-library` tests and artifacts depend on their current
stage-specific behavior.

## Tests and Acceptance Criteria

Unit and integration coverage must include:

1. defaults are threshold `0.8` and minimum size `10`;
2. values outside `0 < threshold <= 1` and sizes below `1` fail;
3. a Leiden community containing one below-threshold pair is split;
4. every retained cluster passes an independently computed all-pairs audit;
5. two fully compatible clusters from different Leiden communities are merged;
6. every accession pair drawn from different retained clusters has similarity
   strictly below the threshold;
7. the non-transitive A/B/C example deterministically orphans at least one
   conflicting accession and satisfies both similarity guarantees;
8. clusters of sizes `9` and `10` are orphaned and retained respectively at the
   default minimum;
9. invalid profiles, no-edge accessions, refined singletons, and undersized
   clusters all appear in `orphan_accessions.tsv`;
10. multi-term orphans produce one sorted row per GO term; no-term orphans
    produce one blank-term row;
11. no accession occurs in both outputs, and the output union equals the unique
    input accession set;
12. duplicate input rows and GO terms are deterministically deduplicated;
13. reruns produce byte-identical TSVs and stable IDs;
14. partial or stale output sets rebuild together, while current outputs skip;
15. a failure before commit leaves the previous complete output set intact;
16. existing `build-library` clustering tests remain unchanged and pass.

## Confirmed Threshold Semantics

- Within a retained cluster, every distinct pair has similarity greater than
  or equal to the threshold.
- Across two retained clusters, every pair has similarity strictly lower than
  the threshold.
- Accessions that prevent both conditions from holding are reported as
  orphans according to the deterministic reconciliation policy.
