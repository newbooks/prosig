# Implementation Spec: GO Accession Clustering

## Goal

Cluster Swiss-Prot primary accessions by Molecular Function GO-set similarity
and write a stable `member_id` / `cluster_id` table for downstream motif
discovery.

This adapts the relevant design from `../pclass`:

- `../pclass/src/pclass/datasets/go_clusters.py`
- `../pclass/src/pclass/commands/research.py`, command `cluster_by_go`
- `../pclass/Feature_Evaluation_Guide.md`, section `cluster_by_go`

ProSig integration differs from pclass in one important way: clustering should
not be exposed as a separate top-level command initially. It should be part of
`prosig build-library` and should run only when the cluster artifact is missing
or stale relative to its dependencies.

## External Algorithm and API References

- Leiden algorithm:
  <https://arxiv.org/abs/1810.08473>
  - Traag, Waltman, and van Eck introduced Leiden as an improvement over
    Louvain with connected-community guarantees and a refinement phase.
- `leidenalg` Python API:
  <https://leidenalg.readthedocs.io/en/stable/reference.html>
  - `find_partition(graph, partition_type, weights, n_iterations,
    max_comm_size, seed, **kwargs)` is the relevant entry point.
  - `RBConfigurationVertexPartition` accepts `resolution_parameter`, matching
    the pclass implementation.
- GO ontology downloads:
  <https://geneontology.org/docs/download-ontology/>
  - `go-basic` is acyclic and recommended for most GO annotation tools.
- GO annotations:
  <https://geneontology.org/docs/go-annotations/>
  - Positive annotations propagate upward; ProSig already materializes direct
    MF annotations in `accession_mf_go.tsv` and propagated counts/IC in
    `go_graph.pkl`.
- GO evidence codes:
  <https://geneontology.org/docs/guide-go-evidence-codes/>
  - ProSig should keep using the existing high-quality evidence filtering in
    `src/prosig/go/build.py`; clustering should consume the filtered artifact,
    not re-parse evidence policy.

## Inputs

Required clustering inputs:

```text
go_graph.pkl
accession_mf_go.tsv
```

`go_graph.pkl` is the compact Molecular Function graph and IC artifact produced
by `build_go_pkl`.

`accession_mf_go.tsv` is a headerless two-column TSV:

```text
P00533	GO:0004672;GO:0005524
Q9SVY5	GO:0000002;GO:0000003
```

Rules:

- column 1 is the Swiss-Prot primary accession
- column 2 is semicolon-separated direct high-quality MF GO terms
- terms absent from `go_graph.pkl` are ignored before clustering
- profiles with no valid GO terms are excluded from graph construction

## Outputs

Primary output:

```text
go_clusters.tsv
```

Format:

```text
member_id	cluster_id
P00533	cluster_0001
Q9SVY5	cluster_0001
```

Rules:

- include a header
- `member_id` is a primary accession
- `cluster_id` is deterministic: `cluster_0001`, `cluster_0002`, ...
- accessions with no positive GO-similarity edge are omitted
- rows are written in active accession order, sorted by accession

Recommended secondary output:

```text
go_clusters_stats.json
go_clusters_meta.tsv
```

Fields:

```json
{
  "algorithm": "go_set_similarity_knn_leiden",
  "similarity": "lin_amb",
  "partition": "RBConfigurationVertexPartition",
  "resolution": 1.0,
  "neighbors": 10,
  "seed": 0,
  "min_informative_ic": 0.5,
  "max_posting_fraction": 0.05,
  "max_posting_size": 0,
  "input_accessions": 123,
  "clustered_accessions": 120,
  "excluded_accessions": 3,
  "informative_terms_before_filtering": 3456,
  "informative_terms_after_filtering": 1234,
  "posting_cap": 7,
  "fallback_accessions_after_filtering": 2,
  "edges": 456,
  "clusters": 12,
  "cluster_size_min": 2,
  "cluster_size_mean": 10.0,
  "cluster_size_median": 8.0,
  "cluster_size_max": 31,
  "lin_matrix": {
    "dtype": "float32",
    "storage": "memory"
  },
  "profile_cache": {
    "budget_mb": 128,
    "max_entries": 262144,
    "entries": 2345,
    "hits": 500,
    "misses": 1000,
    "evictions": 0
  },
  "dependencies": {
    "go_graph": "go_graph.pkl",
    "accession_go": "accession_mf_go.tsv"
  }
}
```

The stats file is not required for downstream behavior, but it gives users a
stable diagnostic artifact and supports future stale-build checks that include
parameter provenance.

Cluster metadata output:

```text
cluster_id	incluster_sim	composed_description
cluster_0001	0.8475	
cluster_0002	NA	
```

Rules:

- `cluster_id`: cluster identifier from `go_clusters.tsv`
- `incluster_sim`: average AMB Lin similarity over all unique accession pairs
  in the cluster, excluding self-pairs
- singleton clusters use `NA`
- unavailable accession-pair similarities contribute `0.0`
- `composed_description` is reserved for the future GO description composer and
  is currently written as an empty value

## Build-Library Integration

`prosig build-library` includes clustering options:

```bash
prosig build-library \
  --go-obo go-basic.obo \
  --swissprot uniprot_sprot.dat.gz \
  --go-out go_graph.pkl \
  --cluster-out go_clusters.tsv \
  --cluster-config cluster_config.yaml
```

Recommended options:

```text
--cluster-out PATH
    Path to write GO accession clusters. Default: go_clusters.tsv.

--cluster-config PATH
    Path to GO clustering config. Created from the starter template when
    missing. Default: cluster_config.yaml.

--force, -f
    Rebuild derived build-library artifacts even when outputs are newer than
    dependencies.
```

The default should build clusters when needed, because the README and command
plan treat clustering as part of library construction.

## Cluster Config

`cluster_config.yaml` is a user-editable flat YAML file created from
`src/prosig/data/cluster_config.yaml.template` when missing. `--force` must not
overwrite an existing config file.

Default template:

```yaml
stats_file: go_clusters_stats.json
meta_file: go_clusters_meta.tsv
neighbors: 10
resolution: 1.0
progress_interval_seconds: 60.0
term_cache_size_mb: 256
profile_cache_size_mb: 128
min_informative_ic: 0.5
max_posting_fraction: 0.05
max_posting_size: 0
```

Fields:

```text
stats_file
    Path to write clustering stats/provenance.

meta_file
    Path to write per-cluster metadata.

neighbors
    Number of nearest GO-similarity neighbors used to build the sparse graph.
    Default: 10. Must be at least 1.

resolution
    Leiden resolution parameter. Default: 1.0. Must be greater than 0.

progress_interval_seconds
    Seconds between long-running kNN progress logs. Must be > 0.

term_cache_size_mb
    Deprecated. Clustering now uses the precomputed Lin matrix instead.

profile_cache_size_mb
    Approximate GO profile-pair AMB cache size in MB. Use 0 to disable.

min_informative_ic
    Minimum ancestor IC for candidate-index terms. Must be >= 0.

max_posting_fraction
    Drop candidate-index terms whose accession posting list covers more than
    this fraction of cleaned accessions. Must be in (0, 1].

max_posting_size
    Optional absolute posting-list cap. Use 0 to derive the cap from
    max_posting_fraction.
```

## Freshness Rules

Clustering should run when any of these conditions is true:

- `--force` is set
- `cluster_out` does not exist
- `cluster_out` exists but is older than `go_graph.pkl`
- `cluster_out` exists but is older than `accession_mf_go.tsv`
- `cluster_out` exists but is older than `cluster_config.yaml`
- the configured `stats_file` does not exist
- the configured `stats_file` exists but is older than `go_graph.pkl`,
  `accession_mf_go.tsv`, or `cluster_config.yaml`
- the configured `meta_file` does not exist
- the configured `meta_file` exists but is older than `go_graph.pkl`,
  `accession_mf_go.tsv`, or `cluster_config.yaml`

Clustering should be skipped when all enabled cluster outputs exist and are at
least as new as every dependency.

Dependency timestamps should be checked after `build_go_pkl` completes, because
that step writes both `go_graph.pkl` and `accession_mf_go.tsv`.

Implementation helper:

```python
def artifact_is_stale(output: Path, dependencies: Iterable[Path]) -> bool:
    if not output.exists():
        return True
    output_mtime = output.stat().st_mtime
    return any(dep.stat().st_mtime > output_mtime for dep in dependencies)
```

If a dependency is missing, fail with a clear diagnostic rather than silently
skipping clustering.

## Algorithm

### 1. Load and Clean Profiles

Load `accession_mf_go.tsv` using the existing parser from
`prosig.go.similarity`:

```python
load_accession_go_terms(path) -> dict[str, tuple[str, ...]]
```

Build a fast scalar similarity index from `go_graph.pkl`. Do not use the
diagnostic `GoSimilarity.set_lin_amb()` path inside clustering hot loops.

The fast index should port the pclass design:

```text
GO ID -> integer term ID
term ID -> IC
term ID -> ancestor bit mask including self
term ID -> ancestors sorted by descending IC
```

Scalar Lin should use integer IDs, bit-mask intersection, and descending-IC
ancestor scans to find the MICA without constructing dataclass diagnostics.
Diagnostic methods in `prosig.go.similarity` should remain available for
interactive inspection; clustering should use the scalar-only API.

Recommended placement:

```text
src/prosig/go/similarity.py
    FastGoSimilarityIndex
    BoundedTermPairCache
    BoundedProfilePairCache
    build_fast_go_similarity_index(...)
    lin_fast(...)
    set_lin_amb_fast_for_valid_profiles(...)

src/prosig/go/clustering.py
    candidate indexing, kNN edge construction, Leiden orchestration
```

This keeps one source of truth for GO similarity while separating clustering
workflow logic from diagnostic formatting.

For each accession:

- deduplicate GO terms
- keep terms present in `go_graph.pkl`
- keep terms with usable IC
- sort terms for deterministic profiles
- drop accessions with an empty cleaned profile

### 2. Build Candidate Index

Avoid dense all-vs-all scoring. For each accession profile, collect informative
candidate terms:

```text
informative_terms(accession) =
  every ancestor of every profile term whose IC is >= min_informative_ic
```

Default `min_informative_ic` is `0.5`. With natural-log IC this removes very
broad ancestors with frequency greater than about `exp(-0.5) = 0.6065`, while
keeping moderately common terms that may still connect related functions. This
is intentionally conservative; users can raise it for faster, more specific
graphs or lower it when too many accessions become isolated.

Build an inverted index:

```text
informative GO term -> accession indices containing that informative term
```

After building the raw inverted index, apply broad-ancestor filtering:

```text
posting_cap = cluster_max_posting_size
if posting_cap <= 0:
    posting_cap = ceil(n_cleaned_accessions * cluster_max_posting_fraction)

drop informative terms whose posting-list size > posting_cap
```

Default `cluster_max_posting_fraction` is `0.05`, so an informative ancestor can
nominate at most 5% of cleaned accessions. This prevents high-level GO terms
from dominating candidate generation. The cap must never drop all informative
terms for an accession if the accession had candidates before filtering; in
that fallback case, keep the accession's rarest available informative term.

This mirrors pclass's inverted-index approach, but adds a default broad-term
cap to reduce unnecessary candidate scoring.

### 3. Build Sparse kNN Edges

For each accession:

1. Add exact same-profile candidates first with similarity `1.0` when the
   profile has at least one positive-IC term.
2. Use the informative-term inverted index to find additional candidates.
3. Score candidates with the fast scalar profile similarity API:

```text
similarity = go_set_similarity_for_valid_profiles_fast(profile_a, profile_b)
```

4. Treat unavailable similarity as no edge.
5. Keep the top `k = neighbors` positive candidates, sorting by:

```text
similarity descending, accession ascending
```

6. Store edges as undirected pairs with the maximum observed weight.

The graph edge weight is the GO-set similarity in `[0.0, 1.0]`.

### 3a. Lin Matrix and Bounded Similarity Cache

Clustering precomputes an in-memory dense `float32` GO term-pair Lin matrix by
filling the upper triangle and mirroring each score. Unavailable scores are
stored as `NaN` so real zero-valued similarities remain distinguishable.

Lin matrix:

```text
shape: n_valid_go_terms x n_valid_go_terms
dtype: float32
storage: memory
```

Profile-pair cache:

```text
key: canonical pair of sorted valid GO profiles
value: AMB Lin similarity or None
default budget: 128 MB
```

Use approximate entry sizes to convert MB budgets to max entries. A profile
cache size of `0` disables that cache. The stats output should include the Lin
matrix path and profile-cache budgets, max entries, final entries, hits, misses,
and evictions.

The profile cache is useful because Swiss-Prot contains repeated or near-common
GO profiles, and broad candidate terms can cause the same profile pairs to be
requested many times.

### 4. Remove Isolates

Only accessions incident to at least one positive edge enter the Leiden graph.

Log the count of excluded accessions. Do not write singleton clusters for
isolated accessions in the initial implementation, because pclass excluded
accessions with no determinable positive GO-similarity edges and downstream
motif discovery needs clusters with comparative support.

### 5. Run Leiden

Use `igraph` plus `leidenalg`:

```python
partition = leidenalg.find_partition(
    graph,
    leidenalg.RBConfigurationVertexPartition,
    weights="weight",
    resolution_parameter=resolution,
    seed=0,
)
```

Rationale:

- Leiden is designed for community detection on weighted graphs.
- The refinement phase avoids the disconnected-community failure mode of
  Louvain.
- `RBConfigurationVertexPartition` preserves the pclass behavior and exposes a
  tunable `resolution_parameter`.
- `seed=0` gives deterministic output for stable build artifacts.

### 6. Assign Stable Cluster IDs

Convert raw Leiden membership to deterministic cluster IDs:

1. group accessions by membership ID
2. sort members within each group
3. sort communities by first accession, then cluster size
4. assign `cluster_0001`, `cluster_0002`, ...

This removes dependence on implementation-specific Leiden community labels.

## Package Layout

Recommended implementation module:

```text
src/prosig/go/clustering.py
```

Recommended public entry point:

```python
def cluster_accessions_by_go(
    accession_go_file: str | Path,
    *,
    go_graph_file: str | Path = "go_graph.pkl",
    output_file: str | Path = "go_clusters.tsv",
    stats_file: str | Path | None = "go_clusters_stats.json",
    meta_file: str | Path | None = "go_clusters_meta.tsv",
    resolution: float = 1.0,
    neighbors: int = 10,
    term_cache_size_mb: int = 256,
    profile_cache_size_mb: int = 128,
    min_informative_ic: float = 0.5,
    max_posting_fraction: float = 0.05,
    max_posting_size: int = 0,
    progress_interval_seconds: float = 60.0,
) -> GoClusteringResult:
    ...
```

Recommended dataclass:

```python
@dataclass(frozen=True)
class GoClusteringResult:
    output_file: Path
    stats_file: Path | None
    meta_file: Path | None
    input_accessions: int
    clustered_accessions: int
    excluded_accessions: int
    edges: int
    clusters: int
```

Keep CLI-specific timestamp decisions in `src/prosig/cli/build_library.py` or a
small build helper. Keep clustering algorithm code outside `prosig.cli`.

## Dependencies

Runtime dependencies needed for clustering:

```toml
dependencies = [
  "igraph>=0.11",
  "leidenalg>=0.10",
  "numpy>=1.26",
]
```

If these are later considered too heavy for base installs, move them to an
optional extra:

```toml
[project.optional-dependencies]
cluster = [
  "igraph>=0.11",
  "leidenalg>=0.10",
]
```

If optional dependencies are used, `build-library` should raise a clear error
when clustering is needed but the extra is missing:

```text
GO clustering requires igraph and leidenalg. Install with: pip install -e .[cluster]
```

## Logging

Log at INFO level:

- whether cluster artifacts are current or will be rebuilt
- GO term Lin matrix build/load progress
- profile-pair cache budget when enabled
- candidate index progress for long runs
- informative-term count before and after broad-ancestor filtering
- posting-list cap and fallback accession count
- retained kNN edge count
- excluded no-edge accession count
- final clustered accessions, cluster count, `k`, and resolution
- cluster size min, mean, median, and max
- output paths

Progress logging should default to every 60 seconds to keep large Swiss-Prot
runs observable without flooding logs.

## Edge Cases

### Empty Input

Write a header-only `go_clusters.tsv` and stats with zero counts.

### One Valid Accession

Write a header-only `go_clusters.tsv`. A single accession has no positive edge
and should be excluded under the initial no-isolate rule.

### No Positive Edges

Write a header-only `go_clusters.tsv`, log zero clusters, and write stats.

### Invalid Parameters

Raise `ValueError` before doing work:

- `resolution <= 0`
- `neighbors < 1`
- `term_cache_size_mb < 0`
- `profile_cache_size_mb < 0`
- `min_informative_ic < 0`
- `max_posting_fraction <= 0` or `max_posting_fraction > 1`
- `max_posting_size < 0`
- `progress_interval_seconds <= 0`

The CLI should convert these into `typer.BadParameter`.

### Missing GO Terms

Silently drop individual terms absent from `go_graph.pkl`, matching existing
GO-set similarity behavior. Report aggregate counts in stats if easy.

### Missing IC

Drop terms without usable IC from clustering profiles. Accessions that become
empty after this filtering are excluded.

### Identical GO Profiles

Connect identical profiles deterministically with similarity `1.0` before
searching broader candidates. This avoids unnecessary repeated Lin scoring.

### Broad-Ancestor Filtering Removes All Candidates

If broad-ancestor filtering removes every informative term that could nominate
candidates for an accession, keep the rarest available informative term for that
accession as a fallback. This preserves graph coverage while still removing
high-fanout ancestors for ordinary cases.

### Cache Limits

The bounded profile cache may evict scores during long runs. Eviction must not
change results, only runtime. Cache key canonicalization must make profile-pair
lookups symmetric.

### Ties

Break candidate ties by accession ID and cluster ID assignment by sorted member
lists. Do not rely on dict insertion from non-sorted sources.

## Tests

Add tests for `prosig.go.clustering`:

- parameter validation
- fast scalar Lin against diagnostic `GoSimilarity.lin`
- fast scalar AMB against diagnostic `GoSimilarity.set_lin_amb`
- dense float32 Lin matrix values against scalar Lin
- bounded profile-pair cache hit, miss, and eviction behavior
- accession GO parsing through the existing two-column artifact
- valid profile filtering against a tiny `go_graph.pkl`
- broad-ancestor filtering by IC and posting-list cap
- broad-ancestor rarest-term fallback
- kNN edges for identical profiles
- kNN edges for partially similar profiles using AMB Lin
- accessions with no positive edges are excluded
- deterministic `cluster_0001` assignment independent of Leiden raw IDs
- header-only output for no edges
- stats JSON contents

Add CLI tests for `build-library`:

- first run builds `go_clusters.tsv`
- second run skips clustering when cluster outputs are newer than
  `go_graph.pkl` and `accession_mf_go.tsv`
- missing `cluster_out` triggers clustering
- older `cluster_out` triggers clustering
- `--force` triggers clustering
- invalid cluster parameters produce Typer errors

Use tiny OBO/Swiss-Prot fixtures similar to `tests/test_build_library.py`.
Mock or monkeypatch the Leiden call in CLI freshness tests where the exact
partition is not under test.

## Open Questions

- Should cluster stats include hashes of dependencies, not just mtimes?
- Should isolated accessions eventually be written as singleton clusters for
  recall, or remain omitted for motif-discovery precision?
- Should `igraph`/`leidenalg` be base dependencies or a `[cluster]` extra?
- Should resolution profiles be supported later to help choose the configured
  `resolution`?
- Should kNN edge construction be parallelized with multiprocessing after the
  single-process implementation is validated?
