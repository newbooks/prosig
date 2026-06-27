# ProSig Motifs

ProSig motifs are compact amino-acid sequence patterns used to describe
interpretable protein traits. A motif can describe an exact residue sequence, a
small set of allowed residues, a residue that must be absent, or a fixed or
variable gap. ProSig uses these motifs to scan protein sequences and convert
matches into numeric features such as motif count and match position.

This page describes the ProSig motif format for motif libraries built from
external sources, discovered signatures, or user-defined rules.

## Motif Pattern Syntax

Patterns are written as adjacent tokens with no separator between them.

```text
C??CH
[AG]CCH
N!P[ST]!P
C??C?(10,15)H??H
[KR]EL>
```

The basic tokens are:

| Token | Meaning |
| --- | --- |
| `A` | Exact amino-acid residue. |
| `?` | Any single residue. |
| `[AG]` | Either `A` or `G`. |
| `!A` | Any residue except `A`. |
| `![AG]` | Any residue except `A` or `G`. |
| `?(3)` | Any run of exactly 3 residues. |
| `?(2,5)` | Any run of 2 through 5 residues, inclusive. |
| `<` | N-terminal anchor. The motif must start at the first residue. |
| `>` | C-terminal anchor. The motif must end at the last residue. |

Standard one-letter amino-acid codes are accepted as exact residues:

```text
A C D E F G H I K L M N P Q R S T V W Y
```

The reference implementation also accepts any uppercase alphabetic residue
symbol when compiling a ProSig motif. This allows `U` and `O` when they appear
in source motifs, but ProSig should validate allowed sequence alphabets
explicitly when the production sequence policy is finalized.

## Examples

An exact motif:

```text
GG
```

matches every `GG` occurrence. Overlapping matches count separately, so `GGGG`
contains three `GG` starts.

A wildcard motif:

```text
C??CH
```

matches `C`, then any two residues, then `C`, then `H`.

A common glycosylation-style pattern:

```text
N!P[ST]!P
```

matches `N`, then any residue except `P`, then `S` or `T`, then any residue
except `P`.

A variable-gap motif:

```text
C??C?(10,15)H??H
```

matches two cysteines separated by two arbitrary residues, then a gap of 10 to
15 residues, then `H`, two arbitrary residues, and `H`.

A C-terminal motif:

```text
[KR]EL>
```

matches `KEL` or `REL` only when the match reaches the C terminus.

## Residue Macros

ProSig defines a small set of residue-set macros. Macros are enclosed in braces
and can be negated with `!`.

| Macro | Expansion | Meaning |
| --- | --- | --- |
| `{+}` | `[KR]` | Positively charged residues. |
| `{-}` | `[DE]` | Negatively charged residues. |
| `{p}` | `[STNQ]` | Uncharged polar residues. |
| `{h}` | `[AILMVFWY]` | Hydrophobic residues. |

Examples:

```text
{-}      D or E
!{-}     any residue except D or E
```

## Terminal Alternatives

The current implementation allows `>` inside a bracket set to mean "one of
these residues, or the C terminus." For example:

```text
F[IVFY]G[LM]M[G>]
```

matches either a final `G` after `M`, or the motif ending immediately after
`M`.

## Motif Library File

A ProSig motif library is a UTF-8 tab-separated file. Comment lines begin with
`#`. The first non-comment line is the header.

```text
# ProSig motif library
name	prosite_ac	description	prosite_pattern	prosig_pattern	status
N_GLYCOSYLATION	PS00001	N-glycosylation site	N-{P}-[ST]-{P}	N!P[ST]!P	prosite
```

Required columns:

| Column | Meaning |
| --- | --- |
| `name` | Unique motif name. This becomes the output `motif_id`. |
| `prosite_ac` | Source PROSITE accession when available. Leave empty for non-PROSITE motifs. |
| `description` | Human-readable motif description. |
| `prosite_pattern` | Original source pattern when available, without the terminal PROSITE period. Leave empty for new native ProSig motifs. |
| `prosig_pattern` | Pattern in ProSig motif syntax. |
| `status` | Motif source/status label. Use `prosite` for motifs translated from PROSITE and `prosig` for native ProSig motifs. Other values are preserved as provided. |

When users define new motifs directly, they should provide a unique `name`, a
clear `description`, a valid `prosig_pattern`, and `status` set to `prosig`.
Source-specific columns can be empty if the motif did not come from PROSITE.

## Match Features

When a motif is scanned against a sequence for function prediction, ProSig uses
binary motif detection. Scanning stops for that accession-motif pair as soon as
the first positive hit is found. Motifs ending with the C-terminal anchor `>`
are scanned from the possible C-terminal suffix positions instead of from every
N-terminal start position.

The current feature output is sparse: it writes only positive accession-motif
pairs. Row presence means the motif is present for that accession.

```text
accession	motif_id
```

If no match is found, no sparse output row is written.

## Expanding The Motif Language

The ProSig syntax is intentionally small. Users and implementers can extend it
in several directions, but each extension should define both user syntax and
matching behavior before it is added to a production library.

Possible extensions include:

- New residue macros, such as domain-specific biochemical groups.
- Dense output mode that reports zero-hit accession-motif pairs.
- Repeated motif groups, for example `(RGG){3}`.
- Motif density rules, for example "at least 5 basic residues in 8 positions."
- Structural tracks such as predicted helix, strand, disorder, or solvent
  exposure.
- Scored motifs such as position weight matrices, profiles, or HMMs.

For ProSig, extensions should be added as separate documented feature versions
instead of silently changing the meaning of existing motif libraries.
