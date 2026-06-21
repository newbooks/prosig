"""Translate PROSITE sequence patterns into ProSig motif libraries."""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

PROSIG_MOTIF_HEADER = [
    "name",
    "prosite_ac",
    "description",
    "prosite_pattern",
    "prosig_pattern",
    "status",
]

RESIDUE_SET_MACROS = {
    frozenset("KR"): "{+}",
    frozenset("DE"): "{-}",
    frozenset("STNQ"): "{p}",
    frozenset("AILMVFWY"): "{h}",
}

AMBIGUOUS_RESIDUE_TRANSLATIONS = {
    "B": "[DN]",
    "Z": "[EQ]",
    "J": "[LI]",
    "X": "?",
    "U": "U",
    "O": "O",
}

AMBIGUOUS_RESIDUE_SET_EXPANSIONS = {
    "B": "DN",
    "Z": "EQ",
    "J": "LI",
    "X": "ACDEFGHIKLMNPQRSTVWYUO",
    "U": "U",
    "O": "O",
}


@dataclass(frozen=True)
class PrositePatternEntry:
    """A PROSITE PATTERN entry selected for ProSig motif translation."""

    name: str
    prosite_ac: str
    description: str
    prosite_pattern: str


@dataclass(frozen=True)
class ProSigMotifRow:
    """A row in a ProSig motif library TSV."""

    name: str
    prosite_ac: str
    description: str
    prosite_pattern: str
    prosig_pattern: str
    status: str


@dataclass(frozen=True)
class PrositeReadResult:
    """Selected PROSITE PATTERN entries and parser counts."""

    entries: list[PrositePatternEntry]
    total_entries: int
    pattern_entries: int
    skipped_pattern_entries_without_pa: int


@dataclass(frozen=True)
class PatternTranslation:
    """Translated ProSig pattern with conversion counters."""

    pattern: str
    macro_conversions: int
    ambiguous_code_translations: int


@dataclass(frozen=True)
class PrositeMotifTranslationStats:
    """Summary statistics for PROSITE-to-ProSig motif translation."""

    total_entries: int
    pattern_entries: int
    skipped_pattern_entries_without_pa: int
    translated_entries: int
    failed_entries: int
    macro_converted_entries: int
    ambiguous_codes_translated: int


@dataclass(frozen=True)
class PrositeMotifTranslationResult:
    """Written output path and summary statistics for motif translation."""

    output_file: Path
    stats: PrositeMotifTranslationStats


def write_prosig_motif_library(
    prosite_file: str | Path = "prosite.dat",
    output_file: str | Path = "prosig_motifs.tsv",
    *,
    logger: logging.Logger | None = None,
) -> PrositeMotifTranslationResult:
    """Translate PROSITE PATTERN entries into a ProSig motif library TSV.

    Failed or unsupported translations are logged and omitted from the motif
    file. Every row written by this translator has ``status`` set to
    ``prosite``.
    """

    log = logger or logging.getLogger(__name__)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    read_result = read_prosite_patterns(prosite_file)
    rows, translated_entries, failed_entries, macro_entries, ambiguous_count = (
        translate_prosite_entries(read_result.entries, logger=log)
    )

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# ProSig motif library\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(PROSIG_MOTIF_HEADER)
        for row in rows:
            writer.writerow(
                [
                    row.name,
                    row.prosite_ac,
                    row.description,
                    row.prosite_pattern,
                    row.prosig_pattern,
                    row.status,
                ]
            )

    return PrositeMotifTranslationResult(
        output_file=output_path,
        stats=PrositeMotifTranslationStats(
            total_entries=read_result.total_entries,
            pattern_entries=read_result.pattern_entries,
            skipped_pattern_entries_without_pa=read_result.skipped_pattern_entries_without_pa,
            translated_entries=translated_entries,
            failed_entries=failed_entries,
            macro_converted_entries=macro_entries,
            ambiguous_codes_translated=ambiguous_count,
        ),
    )


def read_prosite_patterns(path: str | Path) -> PrositeReadResult:
    """Read PROSITE flat-file entries and select PATTERN entries with PA lines."""

    entries: list[PrositePatternEntry] = []
    current: dict[str, list[str] | str] = {"pa": []}
    total_entries = 0
    pattern_entries = 0
    skipped_without_pa = 0

    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "//":
                total_entries, pattern_entries, skipped_without_pa = _finish_entry(
                    current,
                    entries,
                    total_entries,
                    pattern_entries,
                    skipped_without_pa,
                )
                current = {"pa": []}
                continue

            if len(line) < 5:
                continue

            key = line[:2]
            value = line[5:].strip()
            if key == "ID":
                current["id"] = value
            elif key == "AC":
                current["ac"] = value.rstrip(";").strip()
            elif key == "DE":
                current["de"] = value
            elif key == "PA":
                pa_lines = current["pa"]
                if isinstance(pa_lines, list):
                    pa_lines.append(value.strip())

    total_entries, pattern_entries, skipped_without_pa = _finish_entry(
        current,
        entries,
        total_entries,
        pattern_entries,
        skipped_without_pa,
    )

    return PrositeReadResult(
        entries=entries,
        total_entries=total_entries,
        pattern_entries=pattern_entries,
        skipped_pattern_entries_without_pa=skipped_without_pa,
    )


def translate_prosite_entries(
    entries: list[PrositePatternEntry],
    *,
    logger: logging.Logger | None = None,
) -> tuple[list[ProSigMotifRow], int, int, int, int]:
    """Translate selected PROSITE entries into ProSig motif rows."""

    log = logger or logging.getLogger(__name__)
    rows: list[ProSigMotifRow] = []
    translated_entries = 0
    failed_entries = 0
    macro_converted_entries = 0
    ambiguous_codes_translated = 0

    for entry in entries:
        try:
            translation = translate_prosite_pattern_with_counts(entry.prosite_pattern)
        except ValueError as exc:
            failed_entries += 1
            log.warning(
                "Skipping unsupported PROSITE motif %s (%s): %s; pattern=%s",
                entry.name,
                entry.prosite_ac or "no accession",
                exc,
                entry.prosite_pattern,
            )
            continue

        translated_entries += 1
        ambiguous_codes_translated += translation.ambiguous_code_translations
        if translation.macro_conversions > 0:
            macro_converted_entries += 1
        rows.append(
            ProSigMotifRow(
                name=entry.name,
                prosite_ac=entry.prosite_ac,
                description=entry.description,
                prosite_pattern=entry.prosite_pattern,
                prosig_pattern=translation.pattern,
                status="prosite",
            )
        )

    return (
        rows,
        translated_entries,
        failed_entries,
        macro_converted_entries,
        ambiguous_codes_translated,
    )


def translate_prosite_pattern(pattern: str) -> str:
    """Translate one PROSITE pattern string into ProSig motif syntax."""

    return translate_prosite_pattern_with_counts(pattern).pattern


def translate_prosite_pattern_with_counts(pattern: str) -> PatternTranslation:
    """Translate one PROSITE pattern and return conversion counters."""

    source = pattern.strip()
    if source.endswith("."):
        source = source[:-1]

    output: list[str] = []
    macro_conversions = 0
    ambiguous_code_translations = 0
    index = 0
    while index < len(source):
        char = source[index]
        if char == "-":
            index += 1
            continue
        if char in "<>":
            output.append(char)
            index += 1
            continue
        if char == "[":
            token, index = _read_group(source, index, "[", "]")
            repeat, index = _read_repeat(source, index)
            translated_token, macro_count, ambiguous_count = _translate_residue_set(
                token
            )
            output.append(translated_token * repeat)
            macro_conversions += macro_count * repeat
            ambiguous_code_translations += ambiguous_count * repeat
            continue
        if char == "{":
            token, index = _read_group(source, index, "{", "}")
            repeat, index = _read_repeat(source, index)
            translated_token, ambiguous_count = _translate_exclusion(token)
            output.append(translated_token * repeat)
            ambiguous_code_translations += ambiguous_count * repeat
            continue
        if char in {"x", "X"}:
            if char == "X":
                ambiguous_code_translations += 1
            repeat_text, next_index = _read_optional_parenthesized_text(
                source, index + 1
            )
            if repeat_text is None:
                output.append("?")
            elif "," in repeat_text:
                minimum, maximum = _read_range(repeat_text)
                output.append(f"?({minimum},{maximum})")
            else:
                repeat = _read_positive_int(repeat_text)
                output.append("?" if repeat == 1 else f"?({repeat})")
            index = next_index
            continue
        if char.isalpha():
            repeat, index = _read_repeat(source, index + 1)
            translated_token, ambiguous_count = _translate_residue_code(char)
            output.append(translated_token * repeat)
            ambiguous_code_translations += ambiguous_count * repeat
            continue

        raise ValueError(f"unsupported PROSITE pattern character: {char!r}")

    return PatternTranslation(
        pattern="".join(output),
        macro_conversions=macro_conversions,
        ambiguous_code_translations=ambiguous_code_translations,
    )


def _finish_entry(
    fields: dict[str, list[str] | str],
    entries: list[PrositePatternEntry],
    total_entries: int,
    pattern_entries: int,
    skipped_without_pa: int,
) -> tuple[int, int, int]:
    if not _has_entry_content(fields):
        return total_entries, pattern_entries, skipped_without_pa

    total_entries += 1
    if _is_pattern_entry(fields):
        pattern_entries += 1
        entry = _entry_from_fields(fields)
        if entry is None:
            skipped_without_pa += 1
        else:
            entries.append(entry)

    return total_entries, pattern_entries, skipped_without_pa


def _entry_from_fields(
    fields: dict[str, list[str] | str],
) -> PrositePatternEntry | None:
    id_value = fields.get("id")
    pa_lines = fields.get("pa", [])
    if not isinstance(id_value, str) or not isinstance(pa_lines, list) or not pa_lines:
        return None

    id_parts = [part.strip().rstrip(".") for part in id_value.split(";")]
    if len(id_parts) < 2 or id_parts[1] != "PATTERN":
        return None

    return PrositePatternEntry(
        name=id_parts[0],
        prosite_ac=str(fields.get("ac", "")),
        description=_strip_terminal_period(str(fields.get("de", ""))),
        prosite_pattern=_strip_terminal_period("".join(pa_lines)),
    )


def _has_entry_content(fields: dict[str, list[str] | str]) -> bool:
    return any(value for value in fields.values())


def _is_pattern_entry(fields: dict[str, list[str] | str]) -> bool:
    id_value = fields.get("id")
    if not isinstance(id_value, str):
        return False
    id_parts = [part.strip().rstrip(".") for part in id_value.split(";")]
    return len(id_parts) >= 2 and id_parts[1] == "PATTERN"


def _strip_terminal_period(value: str) -> str:
    value = value.strip()
    if value.endswith("."):
        return value[:-1]
    return value


def _read_group(source: str, start: int, opener: str, closer: str) -> tuple[str, int]:
    end = source.find(closer, start + 1)
    if end == -1:
        raise ValueError(f"unclosed group starting with {opener!r}")
    return source[start : end + 1], end + 1


def _translate_residue_set(token: str) -> tuple[str, int, int]:
    residues = token[1:-1]
    if not residues:
        raise ValueError("empty residue set")

    expanded_residues, ambiguous_count, contains_any_residue = _expand_residue_codes(
        residues
    )
    if contains_any_residue:
        if ">" in residues:
            return "[?>]", 0, ambiguous_count
        return "?", 0, ambiguous_count

    macro = RESIDUE_SET_MACROS.get(frozenset(expanded_residues))
    if macro is not None:
        return macro, 1, ambiguous_count
    if expanded_residues == list(residues):
        return token, 0, ambiguous_count
    return f"[{''.join(expanded_residues)}]", 0, ambiguous_count


def _translate_residue_code(residue: str) -> tuple[str, int]:
    translation = AMBIGUOUS_RESIDUE_TRANSLATIONS.get(residue)
    if translation is None:
        return residue, 0
    return translation, 1


def _translate_exclusion(token: str) -> tuple[str, int]:
    residues = token[1:-1]
    if not residues:
        raise ValueError("empty exclusion group")

    expanded_residues, ambiguous_count, contains_any_residue = _expand_residue_codes(
        residues
    )
    if contains_any_residue:
        raise ValueError("exclusion containing X cannot be represented")
    if len(expanded_residues) == 1:
        return f"!{expanded_residues[0]}", ambiguous_count
    return f"![{''.join(expanded_residues)}]", ambiguous_count


def _expand_residue_codes(residues: str) -> tuple[list[str], int, bool]:
    expanded_residues: list[str] = []
    ambiguous_count = 0
    contains_any_residue = False

    for residue in residues:
        expansion = AMBIGUOUS_RESIDUE_SET_EXPANSIONS.get(residue)
        if expansion is None:
            if residue not in expanded_residues:
                expanded_residues.append(residue)
            continue

        ambiguous_count += 1
        if residue == "X":
            contains_any_residue = True
        for expanded_residue in expansion:
            if expanded_residue not in expanded_residues:
                expanded_residues.append(expanded_residue)

    return expanded_residues, ambiguous_count, contains_any_residue


def _read_repeat(source: str, index: int) -> tuple[int, int]:
    repeat_text, next_index = _read_optional_parenthesized_text(source, index)
    if repeat_text is None:
        return 1, next_index
    if "," in repeat_text:
        raise ValueError("token repetitions must use a single integer")
    return _read_positive_int(repeat_text), next_index


def _read_optional_parenthesized_text(
    source: str,
    index: int,
) -> tuple[str | None, int]:
    if index >= len(source) or source[index] != "(":
        return None, index
    end = source.find(")", index + 1)
    if end == -1:
        raise ValueError("unclosed parenthesized expression")
    return source[index + 1 : end].strip(), end + 1


def _read_range(text: str) -> tuple[int, int]:
    parts = [part.strip() for part in text.split(",")]
    if len(parts) != 2:
        raise ValueError(f"invalid range expression: {text!r}")
    minimum = _read_non_negative_int(parts[0])
    maximum = _read_non_negative_int(parts[1])
    if minimum > maximum:
        raise ValueError(f"invalid range with minimum greater than maximum: {text!r}")
    return minimum, maximum


def _read_positive_int(text: str) -> int:
    value = _read_non_negative_int(text)
    if value < 1:
        raise ValueError(f"expected a positive integer: {text!r}")
    return value


def _read_non_negative_int(text: str) -> int:
    if not text.isdecimal():
        raise ValueError(f"expected a non-negative integer: {text!r}")
    return int(text)
