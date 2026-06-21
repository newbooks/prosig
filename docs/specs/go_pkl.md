# Implementation Spec: Build Minimal GO MF Graph with IC (`go_graph.pkl`)

## Goal

Create a compact runtime artifact, `go_graph.pkl`, containing only the Gene Ontology Molecular Function graph structure and information content values needed for function similarity calculations.

This file should **not** contain accession-to-GO mappings. Accession annotations are stored separately and are not required for function similarity once term IC values have been computed.

## Inputs

### 1. GO OBO file

Recommended:

```text
go-basic.obo (assume it exists in working directory)
```

Used to obtain:

```text
GO term ID
GO term name
namespace
parent relationships
obsolete status
```

Only terms in the `molecular_function` namespace should be kept.

### 2. Swiss-Prot GO annotation source

Used only during artifact construction to compute MF term frequencies and IC values.

Acceptable sources include:

```text
uniprot_sprot.dat.gz
Swiss-Prot accession-to-GO mapping TSV
GOA filtered to Swiss-Prot reviewed proteins
```

Only Molecular Function GO terms should be used.

## Output

### File

```text
go_graph.pkl
```

### Intended Runtime Use

The file should support:

```text
GO term lookup
ancestor lookup
lowest common ancestor / most informative common ancestor
Lin similarity
GO term depth lookup
IC lookup
```

It should not support accession lookup.

## Output Structure

Recommended top-level structure:

```python
{
    "meta": {...},
    "terms": {...}
}
```

### Diagnostic Files

The build also writes a pretty JSON diagnostic view of the pickle artifact:

```text
go_graph.json
```

This file must include a top-level `_comment` field stating that the content is
for diagnostic use only. The runtime artifact remains `go_graph.pkl`.

The build also writes excluded MF annotation diagnostics:

```text
excluded_mf_annotations.tsv
```

This file contains only Molecular Function GO annotations excluded from IC
calculation because the evidence code is in `EXCLUDED_EVIDENCE`. Biological
Process and Cellular Component annotations must not be written to this file.
The format is:

```text
accession	go_term	evidence
...
```

The build also writes primary accession-to-MF-GO terms for later accession
diagnostics and GO set similarity:

```text
accession_mf_go.tsv
```

This file contains only primary Swiss-Prot accessions with at least one
high-quality Molecular Function GO annotation. GO terms are direct annotations
from `uniprot_sprot.dat.gz`, not propagated ancestors. The file is headerless
two-column TSV:

```text
accession	mf_go_terms
A0A023FBW4	GO:0019958
A0A024B7W1	GO:0008289;GO:0034062;GO:0060090;GO:0140272
```

Excluded evidence codes are omitted. Terms are sorted within each accession.

## `meta` Fields

```python
"meta": {
    "schema_version": "1.0",
    "namespace": "molecular_function",
    "source_obo": "go-basic.obo",
    "annotation_source": "Swiss-Prot",
    "ic_formula": "-log(freq)",
    "frequency_denominator": "accessions with at least one valid MF graph term",
    "propagated_counts": True,
    "obsolete_terms_removed": True,
    "n_terms": 0,
    "n_accessions_provided": 0,
    "n_accessions_with_hq_mf_go": 0,
    "n_accessions_with_any_mf_go": 0,
    "n_hq_mf_go_assignments_not_in_graph": 0,
    "n_hq_mf_go_assignments_obsolete": 0,
    "created_at": "YYYY-MM-DD"
}
```

## `terms` Structure

Each key is a GO term ID.

```python
"terms": {
    "GO:0003674": {
        "name": "molecular_function",
        "parents": [],
        "children": ["GO:0003824", "GO:0005488", ...],
        "ancestors": set(),
        "depth": 0,
        "count": 123456,
        "freq": 1.0,
        "ic": 0.0
    },

    "GO:0005524": {
        "name": "ATP binding",
        "parents": ["GO:0005488"],
        "children": [...],
        "ancestors": {"GO:0003674", "GO:0005488"},
        "depth": 2,
        "count": 5432,
        "freq": 0.044,
        "ic": 3.12
    }
}
```

## Required Term Fields

| Field       | Type        | Description                                                       |
| ----------- | ----------- | ----------------------------------------------------------------- |
| `name`      | `str`       | Human-readable GO term name                                       |
| `parents`   | `list[str]` | Direct parent GO terms within MF                                  |
| `children`  | `list[str]` | Direct child GO terms within MF                                   |
| `ancestors` | `set[str]`  | All ancestor terms within MF                                      |
| `depth`     | `int`       | Shortest distance from MF root                                    |
| `count`     | `int`       | Number of annotated Swiss-Prot accessions propagated to this term |
| `freq`      | `float`     | `count / accessions with at least one valid MF graph term`         |
| `ic`        | `float | None` | `-log(freq)` for counted terms; `None` for zero-count terms     |

## Construction Algorithm

### Step 1: Parse OBO

Load the OBO file and keep only terms where:

```text
namespace == molecular_function
is_obsolete != true
```

Store:

```text
GO ID
name
parents from is_a relationships
```

Use `is_a` relationships only. Do not include `part_of` relationships in the
minimal graph.

### Step 2: Build MF Graph

For each retained MF term:

```text
parents = direct MF parents
children = reverse parent links
```

Remove parent links pointing to terms outside MF or obsolete terms.

### Step 3: Compute Ancestors and Depth

For each term:

```text
ancestors = all recursive parents
depth = shortest path length from GO:0003674
```

Recommended root:

```text
GO:0003674  molecular_function
```

Terms not connected to the MF root should be excluded or flagged.

### Step 4: Load Swiss-Prot MF Annotations

Read accession-to-GO annotations from the selected Swiss-Prot source.

When the source is `uniprot_sprot.dat.gz`, follow `docs/specs/uniprot_sprot_extraction.md`:

```text
primary accession -> set of high-quality Molecular Function GO terms
```

Use only the first accession from the combined `AC` lines for each Swiss-Prot
entry. For the reviewed Swiss-Prot source `uniprot_sprot.dat.gz`, keep only
`DR   GO;` records where the namespace is `F` and the evidence code is not in:

```python
EXCLUDED_EVIDENCE = {"ND", "NAS"}
```

Secondary accessions, non-MF GO terms, and excluded evidence codes must not
contribute to IC counts.

For each accession:

```text
direct_terms = directly annotated MF GO terms
```

Ignore:

```text
non-MF terms
obsolete terms
terms missing from the graph
```

### Step 5: Propagate Counts

For each accession, create:

```text
propagated_terms = direct_terms + all ancestors of direct_terms
```

Count each accession at most once per GO term.

Example:

```text
P12345 -> GO:0005524 ATP binding

Count incremented for:
GO:0005524
GO:0005488
GO:0003674
```

### Step 6: Compute Frequencies and IC

Let:

```python
n_accessions_provided = number of primary Swiss-Prot accessions parsed from the annotation source
n_accessions_with_hq_mf_go = number of accessions with at least one high-quality MF GO annotation before graph filtering
ic_denominator = number of accessions with at least one valid MF term present in the graph
n_hq_mf_go_assignments_not_in_graph = number of high-quality MF GO annotation assignments whose GO IDs are not present in the retained MF graph
n_hq_mf_go_assignments_obsolete = number of high-quality MF GO annotation assignments whose GO IDs are obsolete
```

For each term:

```python
freq = count / ic_denominator
ic = -log(freq)
```

For terms with zero count:

```python
count = 0
freq = 0.0
ic = None
```

Do not assign infinite IC by default.

## Similarity Support

The resulting file should allow downstream functions such as:

```python
def resnik(term_a, term_b, go):
    common = go[term_a]["ancestors"] | {term_a}
    common &= go[term_b]["ancestors"] | {term_b}
    return max(go[t]["ic"] for t in common if go[t]["ic"] is not None)


def lin(term_a, term_b, go):
    mica_ic = resnik(term_a, term_b, go)
    ic_a = go[term_a]["ic"]
    ic_b = go[term_b]["ic"]
    if mica_ic is None or ic_a is None or ic_b is None or ic_a + ic_b == 0:
        return None
    return 2 * mica_ic / (ic_a + ic_b)
```

## Validation Checks

After building `go_graph.pkl`, report:

```text
number of MF terms
number of root terms
number of provided accessions
number of accessions with any parsed MF annotation
number of accessions used for IC
number of accessions skipped because no valid MF term remained after graph filtering
number of HQ MF GO assignments skipped because they were not in the MF graph
number of terms with count > 0
number of terms with IC value
maximum depth
top 20 most frequent MF terms
top 20 highest-IC terms with nonzero count
```

Sanity checks:

```text
GO:0003674 should have freq = 1.0 and ic = 0.0
No obsolete terms should remain
All parents and children should exist in terms
Every child should list the parent reciprocally
No accession_to_terms field should exist in go_graph.pkl
```

## Recommended CLI

```bash
prosig build-library \
    --go-obo go-basic.obo \
    --swissprot uniprot_sprot.dat.gz \
    --go-out go_graph.pkl
```

Optional arguments:

```bash
--write-report go_pkl_report.txt
```

The command always builds the Molecular Function namespace, uses `is_a`
relationships only, computes IC with the natural logarithm, and keeps zero-count
terms with `freq = 0.0` and `ic = None`.

## Design Notes

`go_graph.pkl` is a runtime graph and IC artifact. It should remain compact.

Do not store:

```text
accession_to_terms
raw OBO records
raw Swiss-Prot records
evidence codes
full GOA rows
sequence data
```

Those belong in separate preprocessing artifacts.

The accession annotation file should be loaded only when evaluating motif quality, not when computing term-to-term function similarity.
