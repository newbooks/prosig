# Implementation Spec: GO Lin Similarity

## Goal

Compute GO-based similarity in ProSig using only Molecular Function terms from
`go_graph.pkl` and only the Lin semantic similarity score.

The runtime implementation should be reusable by diagnostics, standalone
similarity calculation, clustering, and later prediction modules. A class-based
index is preferred because it centralizes artifact validation, term lookup,
ancestor handling, and future caches while keeping CLI code thin.

## Runtime Artifact

Input:

```text
go_graph.pkl
```

The artifact is the ProSig build output documented in `go_pkl.md`:

```python
{
    "meta": {"namespace": "molecular_function", ...},
    "terms": {
        "GO:0005524": {
            "name": "ATP binding",
            "parents": [...],
            "children": [...],
            "ancestors": set(...),
            "depth": 2,
            "freq": 0.044,
            "ic": 3.12,
        },
    },
}
```

No accession-to-GO mapping is required for term-to-term similarity.

## Lin Score

For two MF terms `t1` and `t2`:

```text
Lin(t1, t2) = 2 * IC(MICA) / (IC(t1) + IC(t2))
```

`MICA` is the most informative common ancestor, meaning the common ancestor
with the highest IC. Each term must be included in its own ancestor set so that
identical terms with positive IC return `1.0`.

## Required Behavior

The implementation must:

- Validate that the artifact namespace is `molecular_function` when metadata is
  present.
- Use the existing `ic` values from the artifact; do not recompute IC.
- Treat `ancestors` as ancestors excluding self and add the queried term during
  similarity calculation.
- Return unavailable results with explicit reasons for missing terms, missing
  IC, no common ancestor, or zero IC denominator.
- Never convert unavailable similarity to `0.0`.

## Public API

Initial module:

```text
prosig.go.similarity
```

Initial class:

```python
GoSimilarity.from_pickle(path)
GoSimilarity.term(go_id)
GoSimilarity.ancestors_including_self(go_id)
GoSimilarity.find_mica(go1, go2)
GoSimilarity.lin(go1, go2)
GoSimilarity.lin_with_details(go1, go2)
```

The detailed result should include:

- queried GO IDs
- similarity or `None`
- MICA GO ID
- IC values for both terms and the MICA
- status: `ok` or `unavailable`
- reason when unavailable
- common ancestors used for diagnostics

## Future Scaling

The class can later add:

- pairwise score caches
- GO term-set similarity using Average Maximum Best-match; see
  `go_set_similarity.md`
- accession-to-term resolution
- vectorized or matrix output for clustering

Those extensions should reuse the same MF-only term lookup and Lin pairwise
calculation instead of reimplementing GO graph traversal.
