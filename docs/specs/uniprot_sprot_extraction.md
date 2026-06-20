# Parsing `uniprot_sprot.dat.gz` for Primary Accessions and High-Quality MF GO Terms

## Goal

Parse `uniprot_sprot.dat.gz` and extract a mapping from each Swiss-Prot **primary accession** to its associated high-quality **Molecular Function** GO terms.

The output should use only:

```text
primary accession -> set of GO molecular function terms
```

Secondary accessions should be ignored.

## Input File

```text
uniprot_sprot.dat.gz
```

This is the UniProtKB/Swiss-Prot flat-file format. It is a gzipped text file containing many protein entries.

Each entry is separated by:

```text
//
```

Example entry fragment:

```text
ID   1433B_HUMAN             Reviewed;         246 AA.
AC   P31946; Q53XZ2; Q96QU6;
DT   21-JUL-1986, integrated into UniProtKB/Swiss-Prot.
DE   RecName: Full=14-3-3 protein beta/alpha;
GN   Name=YWHAB;
OS   Homo sapiens (Human).
DR   GO; GO:0005634; C:nucleus; IDA:UniProtKB.
DR   GO; GO:0005829; C:cytosol; IDA:UniProtKB.
DR   GO; GO:0005515; F:protein binding; IPI:UniProtKB.
DR   GO; GO:0050815; F:phosphoserine residue binding; IDA:UniProtKB.
//
```

## Relevant Line Types

### `ID` line

Provides entry name and review status.

Example:

```text
ID   1433B_HUMAN             Reviewed;         246 AA.
```

This parser does not need the `ID` line unless validation/debug output is desired.

### `AC` line

Provides accession numbers.

Example:

```text
AC   P31946; Q53XZ2; Q96QU6;
```

The first accession listed is the **primary accession**:

```text
P31946
```

The remaining accessions are secondary accessions:

```text
Q53XZ2
Q96QU6
```

Use only the primary accession as the key.

### `DR   GO` lines

Provide GO cross-references.

Example:

```text
DR   GO; GO:0005515; F:protein binding; IPI:UniProtKB.
```

Format:

```text
DR   GO; <GO_ID>; <Namespace>:<GO term name>; <Evidence>:<Source>.
```

Fields:

```text
GO_ID      = GO:0005515
Namespace  = F
GO name    = protein binding
Evidence   = IPI
Source     = UniProtKB
```

For this project, keep only namespace:

```text
F
```

where `F` means Molecular Function.

Ignore:

```text
P  Biological Process
C  Cellular Component
```

## Multi-Line `AC` Handling

Some entries contain multiple `AC` lines.

Example:

```text
AC   Q9XYZ1; Q9XYZ2; Q9XYZ3;
AC   Q9XYZ4; Q9XYZ5;
```

Treat all `AC` lines in one entry as one combined accession list:

```text
Q9XYZ1; Q9XYZ2; Q9XYZ3; Q9XYZ4; Q9XYZ5;
```

Then split by semicolon:

```text
Q9XYZ1
Q9XYZ2
Q9XYZ3
Q9XYZ4
Q9XYZ5
```

Use only the first accession:

```text
Q9XYZ1
```

Do not create records for secondary accessions.

## Multi-Line Entry Handling

Read the file entry-by-entry. An entry ends when the line is exactly:

```text
//
```

Within each entry:

1. collect all `AC` lines
2. collect all `DR   GO` lines
3. extract the first accession from the combined `AC` list
4. extract qualifying GO terms
5. emit mapping if the accession has at least one qualifying GO term

## Evidence Code Exclusion Rule

When parsing the reviewed Swiss-Prot accession file `uniprot_sprot.dat.gz`,
exclude only these evidence codes:

```python
EXCLUDED_EVIDENCE = {"ND", "NAS"}
```

All other evidence codes are retained. This source-specific policy is tied to
reviewed Swiss-Prot records and should not be applied blindly to unreviewed
annotation sources.

## GO Filtering Rules

For each `DR   GO` line, keep the GO term only if all conditions are true:

```text
line starts with "DR   GO;"
namespace is "F"
evidence code is not in EXCLUDED_EVIDENCE
GO ID is present
```

Example kept line:

```text
DR   GO; GO:0005515; F:protein binding; IPI:UniProtKB.
```

Example excluded because Cellular Component:

```text
DR   GO; GO:0005634; C:nucleus; IDA:UniProtKB.
```

Example excluded because Biological Process:

```text
DR   GO; GO:0006468; P:protein phosphorylation; IDA:UniProtKB.
```

Example excluded because low-quality/automatic evidence:

```text
DR   GO; GO:0005524; F:ATP binding; IEA:InterPro.
```

## Recommended Parser Behavior

Return:

```python
dict[str, set[str]]
```

Example:

```python
{
    "P31946": {"GO:0005515", "GO:0050815"},
    "P00533": {"GO:0004672"},
}
```

Do not store:

```text
secondary accessions
GO term names
GO evidence codes
GO sources
P namespace terms
C namespace terms
entries without qualifying MF GO terms
```

## Pseudocode

```python
import gzip

EXCLUDED_EVIDENCE = {"ND", "NAS"}

def parse_swissprot_mf_go(path: str) -> dict[str, set[str]]:
    accession_to_go = {}

    with gzip.open(path, "rt") as f:
        entry_lines = []

        for line in f:
            line = line.rstrip("\n")

            if line == "//":
                accession, mf_terms = parse_entry(entry_lines)

                if accession and mf_terms:
                    accession_to_go[accession] = mf_terms

                entry_lines = []
            else:
                entry_lines.append(line)

    return accession_to_go


def parse_entry(lines: list[str]) -> tuple[str | None, set[str]]:
    accessions = []
    mf_terms = set()

    for line in lines:
        if line.startswith("AC"):
            # Example:
            # AC   P31946; Q53XZ2; Q96QU6;
            ac_text = line[5:].strip()
            accessions.extend(
                x.strip()
                for x in ac_text.split(";")
                if x.strip()
            )

        elif line.startswith("DR   GO;"):
            # Example:
            # DR   GO; GO:0005515; F:protein binding; IPI:UniProtKB.
            parts = [x.strip() for x in line.split(";")]

            if len(parts) < 4:
                continue

            go_id = parts[1]
            namespace_and_name = parts[2]
            evidence_and_source = parts[3]

            if not go_id.startswith("GO:"):
                continue

            if not namespace_and_name.startswith("F:"):
                continue

            evidence = evidence_and_source.split(":", 1)[0].strip()

            if evidence in EXCLUDED_EVIDENCE:
                continue

            mf_terms.add(go_id)

    primary_accession = accessions[0] if accessions else None

    return primary_accession, mf_terms
```

## Validation Checks

After parsing, report:

```text
number of Swiss-Prot entries read
number of entries with AC lines
number of entries with any GO annotation
number of entries with qualifying MF GO terms
number of unique primary accessions emitted
number of unique MF GO terms emitted
top 20 most frequent MF GO terms
```

Sanity checks:

```text
No accession key should contain semicolon
No secondary accession should be emitted as a separate key
All GO terms should start with "GO:"
All retained GO lines should have namespace F
No retained evidence code should be in EXCLUDED_EVIDENCE
```

## Notes

This parser intentionally keeps the output minimal. It is suitable for building a separate accession-to-MF-GO artifact or for computing propagated MF term counts during `go_graph.pkl` construction.

The `go_graph.pkl` runtime artifact should not include this accession mapping. It should only include GO graph topology and IC values.
