# Implementation Spec: GO Set Similarity

## Goal

Compute similarity between two sets of Molecular Function GO terms using the
Average Maximum Best-match (AMB) method over pairwise Lin scores.

This spec adapts the relevant design from `../pclass`:

- `../pclass/src/pclass/annotations/go_similarity.py`
- `../pclass/tests/test_annotations_go_similarity.py`
- `../pclass/docs/science/go_sim_lin.md`

ProSig keeps the scope narrower than pclass for now:

- Molecular Function terms only.
- Lin score only.
- `go_graph.pkl` is the source of graph topology and IC values.
- Unavailable pairwise scores are excluded from AMB aggregation, never treated
  as `0.0`.

## User Inputs

Each side of a set-similarity query can be supplied in one of two forms:

1. A direct GO term set.
2. An accession whose MF GO terms are resolved from an accession-to-GO artifact.

### Direct GO Term Sets

Use a semicolon- or comma-separated list, with optional parentheses:

```text
(GO:0005524;GO:0004672)
(GO:0005515)
GO:0005524;GO:0004672
GO:0005524,GO:0004672
```

Rules:

- Strip surrounding whitespace.
- Parentheses are optional, but if one parenthesis is present both are required.
- Split on semicolons. Commas are also tolerated for compatibility.
- Strip whitespace around individual terms.
- Deduplicate while preserving input order.
- Reject malformed GO IDs. A valid token matches:

```text
GO:\d{7}
```

Examples:

```text
GO:0005524; GO:0004672
(GO:0005515)
```

Invalid examples:

```text
(GO:0005524;GO:0004672
(GO:00000012)
()
```

### Accession Inputs

An input without direct-set parentheses is treated as an accession:

```text
P00533
Q9SVY5
```

The two query sides may use different input modes:

- set vs set is valid
- accession vs accession is valid
- set vs accession is valid

## Accession-to-GO Artifact

GO graph artifacts intentionally do not store accession annotations. Set
similarity by accession requires a separate accession-to-MF-GO artifact.

Initial TSV format produced by `prosig build-library`:

```text
P00533	GO:0004672;GO:0005524
Q9SVY5	GO:0000002;GO:0000003
```

The file is headerless two-column TSV:

- column 1: primary Swiss-Prot accession
- column 2: semicolon-separated high-quality direct MF GO terms

Parsing rules:

- Read as tab-delimited rows.
- Require at least two columns per non-empty row.
- Extract GO IDs from column 2 with `GO:\d{7}`.
- Deduplicate extracted terms while preserving artifact order.
- If an accession is missing or resolves to zero GO terms, fail with a clear
  diagnostic such as `No GO terms found for accession(s): P00533`.

This artifact should be built separately from `go_graph.pkl` because term
similarity only requires graph topology and IC values.

## Algorithm

Given two cleaned term sets:

```text
A = {a1, a2, ..., am}
B = {b1, b2, ..., bn}
```

First remove terms not present in the MF GO graph:

```text
valid_A = terms in A present in go_graph.pkl
valid_B = terms in B present in go_graph.pkl
```

Track missing terms for diagnostics. Do not treat missing terms as zero-score
matches.

For every pair `(ai, bj)` in `valid_A x valid_B`, compute the pairwise Lin
score with the efficient scalar method:

```text
score(ai, bj) = GoSimilarity.lin(ai, bj)
```

Ignore pairwise results where the score is `None`.

For each term, keep its best valid score in the opposite set:

```text
best_A_to_B(ai) = max score(ai, bj) over bj in valid_B
best_B_to_A(bj) = max score(ai, bj) over ai in valid_A
```

Then compute directional means:

```text
mean_A_to_B = mean(best_A_to_B values)
mean_B_to_A = mean(best_B_to_A values)
```

The final AMB score is:

```text
AMB(A, B) = (mean_A_to_B + mean_B_to_A) / 2
```

Example from pclass-style tests:

```text
A = (GO:0000002, GO:0000003)
B = (GO:0000002)

Lin(GO:0000002, GO:0000002) = 1.0
Lin(GO:0000003, GO:0000002) = 0.4

A->B best scores = {GO:0000002: 1.0, GO:0000003: 0.4}
B->A best scores = {GO:0000002: 1.0}

mean_A_to_B = 0.7
mean_B_to_A = 1.0
AMB = 0.85
```

## Edge Cases

### Empty Cleaned Set

If either side has no terms remaining after graph filtering:

```text
similarity = None
status = unavailable
reason = empty_cleaned_set
```

### No Valid Pairwise Similarity

If both sides have graph terms but every pairwise Lin score is unavailable:

```text
similarity = None
status = unavailable
reason = no_valid_pairwise_similarity
```

### Missing GO Terms

Missing GO terms are reported in diagnostics:

```text
missing_terms1 = (...)
missing_terms2 = (...)
```

They do not participate in directional means and are not counted as zeroes.

### Missing IC

Terms with missing IC may be present in the graph, but pairwise scores involving
them may be unavailable. They are ignored unless they have at least one valid
pairwise match.

### Duplicate Input Terms

Deduplicate terms before scoring. Preserve first-seen order for diagnostics.

### Symmetry

The final AMB score is symmetric by definition. Detailed directional best-match
tables are directional and should preserve the input-side labels.

## Public API

Extend `prosig.go.similarity`.

Recommended dataclasses:

```python
@dataclass(frozen=True)
class GoBestMatch:
    source: str
    target: str
    score: float


@dataclass(frozen=True)
class GoSetSimilarityResult:
    query1: str
    query2: str
    terms1: tuple[str, ...]
    terms2: tuple[str, ...]
    valid_terms1: tuple[str, ...]
    valid_terms2: tuple[str, ...]
    similarity: float | None
    status: str
    reason: str
    best_matches_1_to_2: tuple[GoBestMatch, ...]
    best_matches_2_to_1: tuple[GoBestMatch, ...]
    missing_terms1: tuple[str, ...]
    missing_terms2: tuple[str, ...]
```

Recommended methods on `GoSimilarity`:

```python
GoSimilarity.set_lin_amb(terms1, terms2) -> float | None
GoSimilarity.set_lin_amb_with_details(terms1, terms2, query1="", query2="") -> GoSetSimilarityResult
```

The scalar method should use `GoSimilarity.lin()` so it benefits from the
efficient term-pair calculation. Detailed mode can additionally collect
best-match rows and term metadata for explanation output.

Recommended standalone helpers:

```python
parse_go_term_set("(GO:0005524,GO:0004672)") -> tuple[str, ...]
load_accession_go_terms(path) -> dict[str, tuple[str, ...]]
lookup_accession_go_terms(path, accessions) -> dict[str, tuple[str, ...]]
resolve_go_set_inputs(query1, query2, accession_file) -> tuple[tuple[str, ...], tuple[str, ...]]
```

## CLI

Implemented diagnostic inspect command:

```text
prosig inspect go-set-sim QUERY1 QUERY2 \
  --go-graph go_graph.pkl \
  --accession-go accession_mf_go.tsv
```

Examples:

```text
prosig inspect go-set-sim "(GO:0005524;GO:0004672)" "(GO:0005515)"
prosig inspect go-set-sim P00533 Q9SVY5 --accession-go accession_mf_go.tsv
prosig inspect go-set-sim "(GO:0005524;GO:0004672)" Q9SVY5 --accession-go accession_mf_go.tsv
```

Default output:

```text
0.8500
```

Use the same score formatting policy as `go-sim`: up to four decimal places,
or `NA` when unavailable.

Verbose output:

```text
prosig inspect go-set-sim P00533 Q9SVY5 -v
```

Verbose output should include:

- resolved query labels and term sets
- ignored missing terms
- A-to-B best-match table
- B-to-A best-match table
- directional means
- final AMB formula
- unavailable reason when applicable

JSON output should serialize the detailed result dataclass.

## Performance Expectations

For set sizes `m` and `n`, AMB requires `m * n` pairwise Lin scores. The simple
score path should:

- use `GoSimilarity.lin()` rather than `lin_with_details()`
- avoid term description formatting
- avoid graph path rendering
- optionally use a term-pair score cache for repeated clustering workloads

Future clustering code should reuse the same AMB/Lin semantics, but not the
diagnostic implementation path. Clustering should use a separate fast scalar
index with bounded term-pair and bounded profile-pair caches; the diagnostic
`GoSimilarity` API remains optimized for clarity and detailed output.

## Tests

Required focused tests:

- set-vs-set parser accepts parenthesized GO terms
- parser rejects malformed or empty sets
- accession resolver loads required TSV columns
- accession resolver rejects missing accessions
- AMB score uses directional best matches
- missing graph terms are ignored and reported
- empty cleaned set returns `empty_cleaned_set`
- all unavailable pairwise scores return `no_valid_pairwise_similarity`
- simple CLI output uses scalar set method
- verbose CLI output reports directional best matches and final AMB score

## Open Questions

- Should the diagnostic set-similarity API expose the fast scalar index, or keep
  it private to clustering workloads?
