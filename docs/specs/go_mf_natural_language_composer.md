# GO Molecular Function Natural-Language Composer

## Goal

Implement a rule-based composer that turns a set of Molecular Function GO terms for a UniProt accession into a concise natural-language function description.

Example input:

```text
P0DO87  GO:0004497;GO:0005506;GO:0016705;GO:0020037
```

Example output:

```text
P0DO87 is annotated as a heme- and iron-binding monooxygenase with oxidoreductase activity.
```

The composer should not simply concatenate GO term names. It should classify GO terms by semantic role, remove redundant ancestor terms, choose a primary functional head, convert binding terms into modifiers, and compose a readable sentence.

---

## Inputs

1. accession: str
2. go_terms: list[str]
3. go_graph: mapping from GO ID to metadata

Each GO term metadata entry should provide:

```python
{
    "id": "GO:0004497",
    "name": "monooxygenase activity",
    "namespace": "molecular_function",
    "parents": set[str],
    "ancestors": set[str],
    "depth": int | None,
    "ic": float | None,
}
```

---

## Outputs

Return a structured object:

```python
{
    "accession": "P0DO87",
    "sentence": "...",
    "head_term": {...},
    "modifiers": [...],
    "supporting_terms": [...],
    "dropped_terms": [...]
}
```

---

## Core Algorithm

### Step 1: Resolve GO Terms

- Look up GO IDs in go_graph.
- Skip missing terms with warning.
- Skip non-MF terms.

### Pre-existing GO Role Knowledge

Semantic role inference is a build-time preprocessing step in
`prosig build-library`. The final composer should consume role metadata already
compiled into `go_graph.pkl`; it should not infer every GO term role from
scratch at description time.

The composer should still treat the GO Molecular Function graph itself as the
first source of semantic role knowledge. In `work/go_graph.pkl` built on
2026-06-21, the root term `GO:0003674 molecular_function` has 33 depth-1 child
branches. These depth-1 terms are useful high-level activity categories:

```text
GO:0140657 ATP-dependent activity
GO:0140691 RNA folding chaperone
GO:0016209 antioxidant activity
GO:0005488 binding
GO:0038024 cargo receptor activity
GO:0003824 catalytic activity
GO:0003774 cytoskeletal motor activity
GO:0102910 dirigent protein activity
GO:0009055 electron transfer activity
GO:0140522 fusogenic activity
GO:0140223 general transcription initiation factor activity
GO:0180020 membrane bending activity
GO:0140912 membrane destabilizing activity
GO:0180024 membrane grommet activity
GO:0060090 molecular adaptor activity
GO:0140104 molecular carrier activity
GO:0098772 molecular function regulator activity
GO:0140313 molecular sequestering activity
GO:0141047 molecular tag activity
GO:0140489 molecular template activity
GO:0060089 molecular transducer activity
GO:0045735 nutrient reservoir activity
GO:0140911 pore-forming activity
GO:0044183 protein folding chaperone
GO:0140776 protein-containing complex destabilizing activity
GO:0140777 protein-containing complex stabilizing activity
GO:0005048 signal sequence receptor activity
GO:0005198 structural molecule activity
GO:0090729 toxin activity
GO:0140110 transcription regulator activity
GO:0180051 translation factor activity
GO:0045182 translation regulator activity
GO:0005215 transporter activity
```

Depth-1 terms do not fully solve role classification. Some biological concepts
span branches. For example, terms containing `receptor` appear under molecular
transducer, binding, catalytic, regulator, cargo receptor, and transporter
branches. Therefore role inference must combine GO-anchor ancestry with
conservative keyword rules and an explicit unknown-term audit.

The current `go_graph.pkl` artifact stores GO IDs, names, parent links,
children, ancestors, depth, IC, frequency, and counts. It does not store GO
definitions, synonyms, comments, xrefs, or chemical reaction participants. If
the composer needs definitions or synonyms for better wording, extend the GO
build artifact to retain those fields from `go-basic.obo` or use an external GO
lookup source.

### Build-Time Semantic Role Inference

`prosig build-library` enriches every non-root GO term in `go_graph.pkl` with
semantic role metadata.

Input:

```text
go_graph.pkl
role_map.yaml
```

CLI option:

```text
prosig build-library --role-map role_map.yaml
```

Default behavior:

- Load `role_map.yaml` from the working directory.
- If the file does not exist, create it from
  `docs/specs/templates/role_map.yaml.template`.
- Compile semantic roles into the `go_graph.pkl` term records.
- Write unknown-role audit rows to `go_terms_unknown_role.txt` next to
  `go_graph.pkl`.

For each GO term except `GO:0003674 molecular_function`:

1. Collect the term itself plus all ancestors.
2. Sort those GO IDs from high IC to low IC. Terms with missing IC sort last.
3. Apply Layer 1 anchor mapping from `role_map.yaml`.
4. If Layer 1 fails, apply Layer 2 keyword or regex matching.
5. If both layers fail, assign role `unknown`.
6. Write unknown terms to:

```text
go_terms_unknown_role.txt
```

Format:

```text
GO:0000000: term name
```

The enriched `go_graph.pkl` term record should add a property such as:

```python
"semantic_role": {
    "role": "catalytic",
    "priority": 100,
    "source": "anchor",
    "matched": "GO:0003824",
}
```

Allowed `source` values:

```text
anchor
keyword
unknown
```

For unknown terms:

```python
"semantic_role": {
    "role": "unknown",
    "priority": 0,
    "source": "unknown",
    "matched": None,
}
```

Build logs should report at minimum:

```text
Loading GO semantic role map: role_map.yaml
Assigning GO semantic roles to <n> non-root GO terms
Applying Layer 1 GO anchor/ancestor role matching
Applying Layer 2 keyword role matching to remaining terms
Processed <n> GO terms for semantic role assignment
GO semantic role layer summary:
  total non-root terms = <n>
  anchor assigned      = <n>
  keyword assigned     = <n>
  unknown              = <n>
GO semantic role stats: binding_cofactor=<n>; catalytic=<n>; ...; unknown=<n>
Wrote GO terms with unknown semantic role: go_terms_unknown_role.txt
```

Role stats should list each assigned role category, with `unknown` last. The
unknown count can be zero when every non-root GO term was assigned by an anchor
or keyword rule. A zero unknown count does not prove the role map is
biologically perfect; it only means the current map covered every term through
the implemented rules.

Legacy shorthand:

```text
GO semantic roles: anchor=<n>; keyword=<n>; unknown=<n>
```

If useful, also log the path to `go_terms_unknown_role.txt`.

Layer 1 anchor matching:

- Build an anchor-to-role index from `roles.*.anchors`.
- Walk the IC-ranked `self + ancestors` list.
- The first GO ID found in the anchor index assigns the role.
- If the same GO ID is configured under multiple roles, prefer the role with
  higher configured priority and log or warn about the ambiguous anchor.
- Record the matched anchor GO ID in `semantic_role.matched`.

Layer 2 keyword or regex matching:

- Match against the GO term name.
- Use case-insensitive matching.
- Exact substring keyword matching is sufficient for the starter map.
- The YAML can later grow an explicit `regex:` list per role if substring
  matching is not expressive enough.
- If multiple keyword roles match, prefer the highest configured priority.
- Record the matched keyword or regex in `semantic_role.matched`.

### Role Map YAML Template

The starter role map is stored as:

```text
docs/specs/templates/role_map.yaml.template
```

This file is a template, not the active runtime configuration. A future
`build-library` creation option should copy or materialize this template to a
user-editable `role_map.yaml`, which can then be amended through curation
cycles. The final reviewed version can later move to a permanent package data
location.

Starter content:

```yaml
roles:

  catalytic:
    priority: 100
    anchors:
      - GO:0003824
    keywords:
      - "catalytic activity"
      - "oxidoreductase activity"
      - "transferase activity"
      - "hydrolase activity"
      - "lyase activity"
      - "isomerase activity"
      - "ligase activity"
      - "kinase activity"
      - "phosphatase activity"
      - "monooxygenase activity"
      - "helicase activity"
      - "ATPase activity"

  binding:
    priority: 20
    anchors:
      - GO:0005488

  binding_cofactor:
    priority: 40
    keywords:
      - "heme binding"
      - "iron ion binding"
      - "zinc ion binding"
      - "magnesium ion binding"
      - "calcium ion binding"
      - "ATP binding"
      - "NAD binding"
      - "FAD binding"
      - "FMN binding"

  binding_nucleic_acid:
    priority: 45
    keywords:
      - "DNA binding"
      - "RNA binding"
      - "nucleic acid binding"

  binding_generic:
    priority: 20
    keywords:
      - "protein binding"
      - "binding"

  transporter:
    priority: 90
    anchors:
      - GO:0005215
    keywords:
      - "transporter activity"
      - "transmembrane transporter activity"
      - "channel activity"

  receptor:
    priority: 80
    anchors:
      - GO:0004872
    keywords:
      - "receptor activity"
      - "signaling receptor activity"

  transcription_factor:
    priority: 85
    anchors:
      - GO:0003700
    keywords:
      - "DNA-binding transcription factor activity"
      - "transcription factor activity"

  regulator:
    priority: 70
    anchors:
      - GO:0098772

  structural:
    priority: 60
    anchors:
      - GO:0005198
    keywords:
      - "structural molecule activity"
      - "structural constituent"

  motor:
    priority: 65
    anchors:
      - GO:0003774
```

Layer summary:

```text
Layer 1: GO anchor/ancestor role map
Layer 2: keyword/regex matching
Layer 3: unknown audit report
```

If the role YAML changes during curation, rebuild `go_graph.pkl` so the compiled
semantic role property stays synchronized with the artifact.

### Step 2: Remove Uninformative Terms

Always drop:

```text
GO:0003674 molecular_function
```

Drop generic ancestors when a more specific descendant is present.

### Step 3: Ancestor Pruning

If term A is an ancestor of term B and both are present:

- Keep B as candidate.
- Remove A from head selection.
- Optionally retain A as supporting context.

### Step 4: Semantic Role Classification

At description time, read semantic role metadata from each GO term record. Do
not repeat the full anchor/keyword inference unless the role property is absent.
If the property is absent, the composer may fall back to the same three-layer
inference as a diagnostic compatibility path.

Initial roles:

```text
catalytic
binding
binding_cofactor
binding_nucleic_acid
binding_generic
transporter
receptor
regulator
structural
transcription_factor
motor
unknown
```

Examples:

```text
heme binding                 -> binding
ATP binding                  -> binding
protein kinase activity      -> catalytic
DNA helicase activity        -> catalytic
transporter activity         -> transporter
receptor activity            -> receptor
```

### Step 5: Head-Term Selection

Priority should come from the compiled semantic role metadata. The starting
priority values are defined in `role_map.yaml`, initially created from
`docs/specs/templates/role_map.yaml.template`.

Initial head-role preference:

```text
catalytic
transporter
receptor
transcription_factor
regulator
structural
motor
binding_cofactor
binding_nucleic_acid
binding_generic
unknown
```

Within a role:

1. Highest IC
2. Greatest depth
3. Longest name

### Step 6: Convert Head to Noun Phrase

Examples:

```text
monooxygenase activity -> monooxygenase
protein kinase activity -> protein kinase
DNA helicase activity -> DNA helicase
ATPase activity -> ATPase
```

### Step 7: Convert Binding Terms to Modifiers

Examples:

```text
heme binding -> heme-binding
iron ion binding -> iron-binding
ATP binding -> ATP-binding
DNA binding -> DNA-binding
RNA binding -> RNA-binding
```

### Step 8: Filter Weak Binding Modifiers

Low-priority modifiers:

```text
binding
protein binding
ion binding
metal ion binding
small molecule binding
```

Prefer specific cofactors and ligands.

### Step 9: Merge Modifiers

Examples:

```text
heme-binding
iron-binding
```

becomes:

```text
heme- and iron-binding
```

Examples:

```text
ATP-binding
```

stays:

```text
ATP-binding
```

### Step 10: Compose Final Sentence

Template:

```text
{accession} is annotated as a {modifier_phrase} {head_phrase}.
```

With supporting context:

```text
{accession} is annotated as a {modifier_phrase} {head_phrase} with {supporting_phrase} activity.
```

---

## Examples

### Example 1

Input:

```text
GO:0004497 monooxygenase activity
GO:0005506 iron ion binding
GO:0016705 oxidoreductase activity
GO:0020037 heme binding
```

Output:

```text
P0DO87 is annotated as a heme- and iron-binding monooxygenase with oxidoreductase activity.
```

### Example 2

Input:

```text
GO:0004672 protein kinase activity
GO:0005524 ATP binding
GO:0000287 magnesium ion binding
GO:0016740 transferase activity
```

Output:

```text
PXXXXX is annotated as an ATP- and magnesium-binding protein kinase.
```

### Example 3

Input:

```text
GO:0003677 DNA binding
GO:0003700 DNA-binding transcription factor activity
GO:0043565 sequence-specific DNA binding
```

Output:

```text
PXXXXX is annotated as a sequence-specific DNA-binding transcription factor.
```

---

## CLI

```bash
prosig inspect accession P0DO87
```

Options:

```bash
--format text|json
--style conservative|direct|detailed
--max-modifiers 3
--show-dropped
--show-roles
--accession-go accession_mf_go.tsv
--go-graph go_graph.pkl
```

`prosig inspect accession` should initially be diagnostic. It can show the
resolved MF GO terms, role assignments, dropped ancestor terms, selected head
term, modifiers, supporting terms, and the composed sentence. If this becomes a
production function-summary command later, add a separate product command or
document the decision.

---

## Design Constraints

1. Do not use an LLM in production composition.
2. Do not concatenate all GO names.
3. Prefer specific descendants over generic ancestors.
4. Prefer catalytic/transporter/receptor terms as heads.
5. Treat binding terms as modifiers.
6. Produce conservative biological language.
7. Return structured debug information.

---

## Suggested Files

```text
src/prosig/go_describe.py
src/prosig/go_roles.py
src/prosig/go_phrase_rules.py
tests/test_go_describe.py
```

---

## Future Extensions

- Evidence-aware wording
- IC-aware ranking improvements
- GO-slim fallback
- Domain-specific phrase overrides
- Multi-sentence summaries

---

## Reference Sources for Role Knowledge

- Gene Ontology overview:
  <https://geneontology.org/docs/ontology-documentation/>
  - Defines the three GO aspects.
  - Explains that Molecular Function terms represent molecular-level
    activities and are commonly named with `activity`.
  - Documents GO term elements such as ID, name, aspect, definition,
    relationships, synonyms, comments, obsolete status, and xrefs.
- GO relations documentation:
  <https://geneontology.org/docs/ontology-relations/>
  - Explains the graph model, parent/child terminology, `is_a`, `part_of`,
    `has_part`, and regulation relations.
  - Notes which relations are safe for grouping annotations.
- GO ontology downloads:
  <https://geneontology.org/docs/download-ontology/>
  - Describes `go-basic`, the acyclic ontology version recommended for most
    GO-based annotation tools.
  - Provides OBO, JSON, and OWL download locations.
- AmiGO term pages:
  <https://amigo.geneontology.org/amigo/term/GO:0004672>
  - Useful for reviewing individual GO term definitions, synonyms, parents,
    children, graph neighborhoods, mappings, and annotations.
- GO annotations:
  <https://geneontology.org/docs/go-annotations/>
  - Explains how gene products are linked to GO terms.
- GO evidence codes:
  <https://geneontology.org/docs/guide-go-evidence-codes/>
  - Relevant for future evidence-aware wording.
- Rhea:
  <https://www.rhea-db.org/>
  - Curated reaction knowledgebase used by UniProtKB and useful for enzymatic
    molecular function wording when GO terms carry reaction xrefs.
