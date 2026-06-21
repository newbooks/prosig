# ProSig Motif Format And Translation Implementation Spec

## Purpose

This spec defines the ProSig motif format and PROSITE translation behavior.
The behavior was audited against the reference repository at `../pclass`, but
the public format and implementation target in this repository should be named
ProSig.

Reference material inspected:

- `../pclass/docs/science/motif_language.md`
- `../pclass/docs/motif_feature_extraction_spec.md`
- `../pclass/src/pclass/data_sources/prosite.py`
- `../pclass/src/pclass/features/motifs.py`
- `../pclass/tests/test_features_motifs.py`
- `../pclass/tests/test_commands.py`

## Scope

The implementation surface has three related parts:

1. ProSig motif pattern syntax.
2. PROSITE `PATTERN` entry translation into a motif-library TSV.
3. Motif library loading, scanning, and sparse feature extraction.

ProSig should keep this behavior stable for ProSig motif libraries while
leaving room for future ProSig-specific motif discovery output.

## Pattern Language

Implement the following tokens:

| Syntax | Behavior |
| --- | --- |
| `A` | Match exact uppercase residue `A`. |
| `?` | Match one uppercase residue. |
| `[AG]` | Match one residue from the bracket set. |
| `!A` | Match one uppercase residue that is not `A`. |
| `![AG]` | Match one uppercase residue that is not in the bracket set. |
| `{+}` | Macro for `[KR]`. |
| `{-}` | Macro for `[DE]`. |
| `{p}` | Macro for `[STNQ]`. |
| `{h}` | Macro for `[AILMVFWY]`. |
| `!{-}` | Negated macro; match one residue not in the macro expansion. |
| `?(n)` | Match exactly `n` uppercase residues. |
| `?(m,n)` | Match `m` through `n` uppercase residues, inclusive, using the shortest valid span from a fixed start. |
| `<` | N-terminal anchor. |
| `>` | C-terminal anchor. |
| `[G>]` | Match `G` or the C terminus. |
| `[?>]` | Match any one residue or the C terminus. |

Validation rules:

- Wildcard run bounds must be non-negative integers.
- Wildcard range minimum must not exceed maximum.
- Bracket residue sets must not be empty.
- Macros must be known.
- Exclusions may target exact residues, bracket sets, or macros.
- Exclusions must not contain `>` or `?`.
- Unsupported tokens should fail fast during library loading.

Audited reference implementation detail:

- The reference implementation compiles motifs to Python regular expressions.
- It uses `[A-Z]` as the single-residue class.
- It uses non-greedy regex quantifiers for variable wildcard runs.
- It scans each sequence start position with `compiled_pattern.match(sequence, start)`.

## Motif Library TSV

The canonical library file is a UTF-8 TSV with optional comment lines beginning
with `#`. The first non-comment line must contain at least:

```text
name	description	prosig_pattern	status
```

The PROSITE-generated form uses this full header:

```text
name	prosite_ac	description	prosite_pattern	prosig_pattern	status
```

Loading rules:

- Ignore lines beginning with `#`.
- Require `name`, `description`, and `prosig_pattern`.
- For import compatibility only, a loader may accept legacy `pclass_pattern`
  when `prosig_pattern` is absent, then normalize it to `prosig_pattern`
  internally.
- Require `status`.
- Treat `status` as a source/status label, not as an active/inactive filter.
- Use `status = "prosite"` for motifs translated from PROSITE.
- Use `status = "prosig"` for motifs provided natively by ProSig.
- Preserve any other non-empty status value as provided.
- Require non-empty `name` and `prosig_pattern` for every row.
- Preserve motif order from the library file.
- Compile all motifs during loading and fail fast on invalid syntax.

## PROSITE Entry Selection

Translate only PROSITE flat-file entries whose `ID` line has `PATTERN` as the
second semicolon-separated field.

Parsing rules:

- Entries end at `//`.
- Count every non-empty entry as a total entry.
- Count every `ID ...; PATTERN.` entry as a pattern entry.
- Skip selected `PATTERN` entries that have no `PA` line.
- Parse `name` from the first `ID` field.
- Parse `prosite_ac` from `AC`, removing a trailing semicolon.
- Parse `description` from `DE`, stripping one terminal `.` if present.
- Build `prosite_pattern` by stripping each `PA` line and concatenating all
  `PA` fragments in entry order.
- Strip one terminal `.` from the concatenated `prosite_pattern` before writing
  the motif library row.

Output rules:

- Write `# ProSig motif library` as the first line.
- Write the full PROSITE header.
- Write one row only for selected `PATTERN` entries that have at least one `PA`
  line and translate successfully.
- Set `status = "prosite"` for every motif written by PROSITE translation.
- Do not write failed or unsupported translations to the motif file.
- Log failed or unsupported translations with motif name, PROSITE accession when
  available, source pattern, and error reason.

The command-level summary should report:

- total entries read
- PATTERN entries found
- entries successfully translated
- entries converted to ProSig macros
- ambiguous residue codes translated
- entries skipped because they have no `PA` line
- entries not written because translation failed or the pattern is unsupported
- output file path

## PROSITE Translation Rules

Before translation, trim whitespace from the full source pattern and remove a
single trailing `.` if present.

Translate source tokens as follows:

| PROSITE source | ProSig target |
| --- | --- |
| `-` | Omit separator. |
| `<` | Preserve N-terminal anchor. |
| `>` | Preserve C-terminal anchor. |
| `x` | `?` |
| `X` | `?` and count one ambiguous-code translation. |
| `x(n)` | `?(n)`, except `x(1)` becomes `?`. |
| `x(m,n)` | `?(m,n)`. |
| `B` | `[DN]` |
| `Z` | `[EQ]` |
| `J` | `[LI]` |
| `U` | `U` |
| `O` | `O` |
| `[ST]` | `[ST]`, unless it exactly matches a macro expansion. |
| `[DE]` or `[ED]` | `{-}` |
| `[KR]` or `[RK]` | `{+}` |
| `[STNQ]` in any order | `{p}` |
| `[AILMVFWY]` in any order | `{h}` |
| `{P}` | `!P` |
| `{EDRK}` | `![EDRK]` after ambiguous-code expansion. |

Ambiguous codes inside bracket sets or exclusion groups must be expanded before
macro detection:

| Source code | Set expansion |
| --- | --- |
| `B` | `DN` |
| `Z` | `EQ` |
| `J` | `LI` |
| `X` | `ACDEFGHIKLMNPQRSTVWYUO` |
| `U` | `U` |
| `O` | `O` |

Set expansion rules:

- Preserve first-seen residue order after expansion.
- Remove duplicates.
- If `X` appears in an allowed residue set, translate the whole set to `?`.
- If `X` appears in a terminal allowed set containing `>`, translate to `[?>]`.
- If `X` appears in an exclusion group, fail translation because "any residue
  except any residue" is not representable in the current language.

PROSITE token repetitions:

- Bracket sets, exclusions, and exact residues may be followed by `(n)`.
- Such repetitions must use a single positive integer.
- Range repetitions on non-wildcard tokens are invalid.
- Translate repetition by repeating the translated token text `n` times.

Examples:

```text
PROSITE: C-x(2)-C-x(10,15)-H-x(2)-H.
ProSig:  C?(2)C?(10,15)H?(2)H

PROSITE: N-{P}-[ST]-{P}.
ProSig:  N!P[ST]!P

PROSITE: B-[DE]-N-{P}-[ST]-{P}.
ProSig:  [DN]{-}N!P[ST]!P
```

## Motif Matching

Use sparse, deterministic feature extraction.

For each sequence and motif:

1. If the sequence is empty, return `count = 0`, `first_position = 0`,
   `last_position = 0`, `match_fraction = "NA"`.
2. Test the motif at every zero-based start index from `0` through
   `len(sequence) - 1`.
3. Record one match for each start index that matches.
4. For variable-length motifs, record the shortest valid span for that start.
5. Allow overlapping matches.

Match coordinates:

- `start_n` is 1-based from the N terminus.
- `end_n` is 1-based inclusive from the N terminus.
- For a regex backend, Python `match.end()` is a zero-based exclusive offset,
  but for non-empty residue-consuming matches its numeric value equals the
  1-based inclusive coordinate of the last consumed residue. Store that value as
  `end_n`; do not expose zero-based exclusive coordinates in the ProSig data
  model.

Feature definitions:

```text
count = number of matches
first_position = min(start_n), or 0 if no matches
last_start_n = max(start_n)
last_position = sequence_length - last_start_n + 1, or 0 if no matches
match_fraction = count / sequence_length, or "NA" for empty sequences
```

Output rows:

- Write only sequence-motif pairs where `count > 0`.
- Omit zero-hit sequence-motif pairs.
- Omit accessions with no motif hits.
- Order rows by cluster/member input order, then motif library order.

Output header:

```text
accession	motif_id	count	first_position	last_position	match_fraction
```

## Command Behavior To Port

Reference commands in `pclass_tools`:

```bash
pclass_tools prosite_to_pclass \
  --prosite-file prosite.dat \
  --output-file prosig_motifs.tsv

pclass_tools motif_extraction clusters.tsv \
  --sequence-file accession_seq.tsv \
  --motif-file prosig_motifs.tsv \
  --output-file motif_features.tsv \
  --progress-interval 60

pclass_tools accession_motif_stat \
  --motif-features-file motif_features.tsv
```

In ProSig, these behaviors likely belong under:

- `prosig build-library` for motif library construction and PROSITE translation.
- `prosig annotate` for sequence scanning and motif-hit reporting.

Do not add separate top-level `scan` or `predict` commands without a decision
record, because `AGENT.md` currently folds scanning and prediction into
`annotate`.

## Tests Codex Should Preserve Or Add

Port focused tests for:

- PROSITE `PATTERN` entry selection and skip behavior for no-`PA` entries.
- PROSITE translation writes `status = "prosite"` for every translated motif.
- Failed PROSITE translations are logged and omitted from the output motif file.
- Native ProSig motifs can use `status = "prosig"`.
- Other non-empty status values are preserved and do not suppress loading.
- Translation of ambiguous exact residues: `B`, `Z`, `J`, `X`, `U`, `O`.
- Translation of ambiguous residues inside bracket sets.
- Macro conversion for `[DE]`, `[KR]`, `[STNQ]`, and `[AILMVFWY]` in any order.
- Failure for malformed groups and exclusions containing `X`.
- Exact, wildcard, alternative, exclusion, macro, and anchor motif compilation.
- Terminal alternatives such as `[G>]` and `[?>]`.
- `MotifMatch.end_n` is stored as a 1-based inclusive coordinate.
- Overlapping match counting with `GG` against `GGGG`.
- Variable-length shortest-span policy with a pattern like `A?(0,4)B`.
- Sparse output order and omission of zero-hit rows.
- Missing sequence warning behavior, if the same cluster/sequence workflow is
  retained.
- Motif support ranking by unique accession count.

## Observed Inconsistencies

1. `../pclass/docs/science/motif_language.md` supports `<`, `>`, and terminal
   alternatives in the translation table and examples, but its token summary
   omits `<`, `>`, `[G>]`, and `[?>]`. The implementation supports them.

2. `../pclass/docs/motif_feature_extraction_spec.md` says progress logs should
   include the missing sequence count so far. The implementation logs accessions
   scanned, total accessions, and percent; missing accessions are logged only as
   an end warning.

3. The user-facing motif language lists a standard amino-acid alphabet, while
   the compiler accepts any uppercase `A-Z` as an exact residue and uses `[A-Z]`
   for wildcard matching. ProSig needs an explicit production sequence alphabet
   policy.

## Completeness Assessment

The inspected reference material is complete enough to implement a compatible
first ProSig motif parser, PROSITE translator, motif-library loader, scanner,
and sparse feature writer under ProSig names.

Information that remains incomplete for a production ProSig design:

- Exact command placement and user CLI names in ProSig.
- Whether ProSig should keep the full PROSITE TSV schema for discovered motifs
  or introduce a source-neutral schema.
- Formal sequence alphabet validation, including whether `U`, `O`, `B`, `Z`,
  `J`, and unknown residues are valid in query sequences.
- Dense output mode for zero-hit motifs.
- Versioning of motif-library syntax and future extensions.
- Semantics for future scored, structural, repeated-group, or density motifs.
- Whether `end_n` should become a public span feature.
