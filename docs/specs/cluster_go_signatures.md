# Cluster GO Signature Synthesis

Status: implemented

## Purpose

Generate a compact Molecular Function GO signature for each functional cluster
produced by `prosig build-library`.

Each functional cluster contains one or more accessions. Each accession can have
one or more direct MF GO annotations in `accession_mf_go.tsv`. The cluster-level
signature should summarize the dominant shared molecular functions while
preserving informative, specific GO terms where possible.

This design is ported from pclass `synthetic_go`, with one intentional policy
change: ProSig treats all retained Swiss-Prot accessions as equally weighted.
Evidence codes are already filtered during `build-library`, and no additional
evidence-code weighting is applied during cluster signature synthesis.

Reference implementation and design notes:

- `../pclass/src/pclass/workflows/synthetic_go.py`
- `../pclass/src/pclass/commands/research.py`, command `synthetic_go`
- `../pclass/docs/science/synthetic_go.md`

## Inputs

`clusters.tsv`

: Final complete-linkage-refined cluster membership table from
  `build-library`.

  ```tsv
  member_id	cluster_id
  P00533	cluster_0001
  Q9SVY5	cluster_0001
  ```

  Requirements:

  - `member_id` is a primary accession.
  - `cluster_id` is the final refined cluster ID.
  - Each accession appears in at most one final cluster.

`accession_mf_go.tsv`

: Headerless two-column TSV written by `build-library`.

  ```tsv
  P00533	GO:0004714;GO:0005524
  Q9SVY5	GO:0004672;GO:0005524
  ```

  Requirements:

  - Column 1 is accession.
  - Column 2 is semicolon-separated direct MF GO IDs.
  - Evidence codes are not present and are not used for weighting.

`go_graph.pkl`

: Compact ProSig GO artifact containing MF term records, ancestors, IC values,
  depth, names, and semantic roles.

## Outputs

Default output from `build-library` is embedded in the cluster metadata report:

```text
clusters_meta.tsv
```

Format:

```tsv
cluster_id	sim_ave	sim_min	sim_max	size	composed_go
cluster_0001	0.84750	0.72000	1.00000	8	GO:0004672;GO:0005524
cluster_0002	NA	NA	NA	1	GO:0016491
```

Terms are selected using raw `support × IC` scores, but only GO IDs are written
to `composed_go`. Terms are ordered from highest selected score to lowest
selected score and separated by `;`.

## Build-library integration

`prosig build-library` should synthesize cluster GO signatures after final
complete-linkage refinement has written `clusters.tsv`.

Default artifacts:

```text
clusters.tsv
clusters_meta.tsv
clusters_stats.json
```

The synthesis step depends on:

- `clusters.tsv`
- `accession_mf_go.tsv`
- `go_graph.pkl`

The current implementation writes signatures as part of metadata generation.
Complete-linkage refinement refreshes `composed_go` whenever final cluster
metadata is rebuilt.

## Algorithm

### 1. Load cluster members

Read `clusters.tsv` with `csv.DictReader`.

Required columns:

- `member_id`
- `cluster_id`

Reject duplicate `member_id` values, because downstream motif discovery and
cluster-level annotation assume one final cluster per accession.

Sort cluster IDs before writing output for deterministic results.

### 2. Load accession GO terms

Load `accession_mf_go.tsv` into:

```python
dict[str, tuple[str, ...]]
```

Parse semicolon-separated GO IDs and remove duplicate GO IDs within each
accession while preserving deterministic ordering.

If a cluster member is absent from `accession_mf_go.tsv`, count it as
`missing_members` and exclude it from the annotated-member denominator.

If a member is present and has at least one parsed GO term, include it in the
annotated-member denominator even if every term is later skipped because it is
missing from `go_graph.pkl` or lacks usable IC.

### 3. Propagate each accession's terms to ancestors

For each accession, create a set of supported terms:

```text
direct MF terms + all ancestors of those terms
```

Each accession contributes at most one vote to a propagated GO term. This is
the key equal-weighting rule: accessions vote, repeated terms do not, and
evidence codes do not change vote weight.

Skip and count terms that are not present in `go_graph.pkl`.

Skip and count terms that do not have numeric IC.

Skip terms with `ic <= 0`, including the MF root, because they do not produce an
informative signature score.

### 4. Compute support and score

For each candidate GO term:

```text
support(term) = accessions supporting term / annotated cluster members
score(term) = support(term) × IC(term)
```

Discard candidates with:

```text
support < min_support
```

Default:

```text
min_support = 0.1
```

Rationale: the threshold removes accidental or very rare terms while allowing
large clusters to retain real minority subfunctions. The default should be
treated as a reporting threshold, not a clustering constraint.

### 5. Rank candidates

Sort candidate terms by:

1. score descending
2. support descending
3. GO ID ascending

This makes output stable when scores tie.

### 6. Suppress redundant parent/child terms

GO terms are hierarchical. Reporting both a broad parent and its more specific
child can be redundant.

Use the pclass conditional parent-suppression rule:

```text
coverage = support(descendant) / support(parent)

if coverage >= parent_coverage_cutoff:
    keep descendant and suppress parent
else:
    keep parent and suppress descendant
```

Default:

```text
parent_coverage_cutoff = 0.8
```

This keeps a specific term when it explains most of the parent support, but
keeps the parent when child support is too narrow.

### 7. Apply relative-drop selection

After redundancy suppression, walk ranked terms in order. Always keep the first
term if one exists. For later terms, stop before the first candidate where:

```text
candidate.score / previous_selected.score < relative_drop_cutoff
```

Default:

```text
relative_drop_cutoff = 0.5
```

Then truncate to:

```text
max_terms = 10
```

This keeps the signature concise and avoids reporting a long tail of weak
terms.

### 8. Write output

For each cluster, write:

```tsv
cluster_id	mf_terms
```

Format `mf_terms` as:

```text
GO:0004672(2.4151);GO:0005524(1.0523)
```

Recommended numeric format:

```python
f"{score:.4f}"
```

If a cluster has no usable synthesized terms, write the cluster row with an
empty `mf_terms` value and increment `empty_clusters`.

## Public API

Implementation lives in `src/prosig/go/clustering.py` because signatures are
currently written only as part of cluster metadata generation. The helper uses
the already-built `FastGoSimilarityIndex`, so `build-library` does not need a
second GO pickle read.

## Validation

Parameter validation:

- `max_terms >= 1`
- `0.0 <= min_support <= 1.0`
- `relative_drop_cutoff >= 0.0`
- `0.0 <= parent_coverage_cutoff <= 1.0`

Input validation:

- `clusters.tsv` must contain `member_id` and `cluster_id`.
- Duplicate cluster members are rejected.
- Missing GO artifact fields should produce clear `ValueError` messages.

## Edge cases

No clusters

: Write header-only metadata.

Singleton clusters

: Process normally. A singleton can produce a signature from that accession's
  propagated GO terms.

Members missing from `accession_mf_go.tsv`

: Exclude from denominator and count in stats.

Terms absent from `go_graph.pkl`

: Skip and count unique missing terms.

Terms lacking usable IC

: Skip and count unique missing-IC terms.

Only root or zero-IC terms

: Write an empty `composed_go` value for that cluster.

## Tests

Add unit tests for:

- equal accession voting when one accession has repeated or multiple direct
  terms;
- ancestor propagation;
- support denominator handling for missing members;
- missing GO terms and missing IC counts;
- `support × IC` ranking;
- parent suppression when descendant coverage is above and below `0.8`;
- relative-drop cutoff;
- deterministic tie ordering by GO ID;
- singleton cluster signature output;
- duplicate `member_id` rejection;
- `build-library` metadata output with `composed_go`.

## Justification

The cluster signature is a reporting artifact, not a clustering input. It should
therefore summarize what the final clusters contain without changing membership.

Using equal accession votes is appropriate for ProSig because retained
annotations come from reviewed Swiss-Prot records after evidence-code filtering.
Applying additional evidence weights at this stage would double-count evidence
quality and make cluster summaries harder to interpret.

The `support × IC` score gives a defensible balance:

- support favors functions shared by many cluster members;
- IC favors specific, informative GO terms;
- ancestor propagation prevents overly specific direct annotations from hiding
  shared broader functions;
- parent suppression avoids redundant signatures.

The defaults are intentionally conservative and match the pclass behavior unless
ProSig has a specific reason to diverge.
