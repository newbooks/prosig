"""ProSig motif library loading, matching, and sparse feature writing."""

from __future__ import annotations

import csv
import logging
import multiprocessing as mp
import os
import re
import stat
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from re import Pattern

from prosig.sequences import indexed_fasta_sequence

PROSIG_MOTIF_COLUMNS = {"name", "description", "prosig_pattern", "status"}
MOTIF_FEATURE_HEADER = [
    "accession",
    "motif_id",
]
MOTIF_FEATURE_COMPLETION_MARKER = "# completed\ttrue"
RESIDUE_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
RESIDUE_CLASS = "[A-Z]"
MACRO_RESIDUES = {
    "+": "KR",
    "-": "DE",
    "p": "STNQ",
    "h": "AILMVFWY",
}
MOTIF_SCAN_PROGRESS_INTERVAL_SECONDS = 60.0
DEFAULT_MOTIF_SCAN_PROCESSES = 8
DEFAULT_MOTIF_SCAN_CHUNK_SIZE = 5_000
MIN_LITERAL_PREFILTER_LENGTH = 2

_WORKER_MOTIFS: tuple[Motif, ...] = ()
_WORKER_FASTA_FILE = ""
_WORKER_FASTA_INDEX_FILE = ""
_WORKER_SHARD_DIR = ""


@dataclass(frozen=True)
class Motif:
    """One compiled ProSig motif."""

    name: str
    description: str
    pattern: str
    status: str
    compiled_pattern: Pattern[str]
    c_terminal_anchor: bool = False
    min_width: int = 0
    max_width: int | None = None
    literal_prefilter: str | None = None


@dataclass(frozen=True)
class MotifMatch:
    """A motif match with 1-based inclusive N-terminal coordinates."""

    start_n: int
    end_n: int


@dataclass(frozen=True)
class MotifFeatureResult:
    """Sparse features extracted from one sequence for one motif."""

    count: int
    first_position: int
    last_position: int
    match_fraction: float | str


@dataclass(frozen=True)
class MotifFeatureWriteResult:
    """Summary counts for motif feature extraction."""

    output_file: Path
    cluster_rows: int
    motifs: int
    accessions_scanned: int
    missing_sequences: int
    feature_rows: int
    motifs_with_hits: int


@dataclass(frozen=True)
class _AccessionScanTask:
    chunk_index: int
    member_ids: tuple[str, ...]


@dataclass(frozen=True)
class _AccessionScanResult:
    chunk_index: int
    shard_file: str
    feature_rows: int
    accessions_scanned: int
    missing_sequences: int
    missing_preview: tuple[str, ...]
    motifs_with_hits: tuple[str, ...]


def write_motif_features(
    *,
    cluster_file: str | Path,
    motif_file: str | Path = "prosig_motifs.tsv",
    fasta_file: str | Path = "accession.fasta",
    fasta_index_file: str | Path = "accession.fasta.idx",
    output_file: str | Path = "motif_features.tsv",
    processes: int = DEFAULT_MOTIF_SCAN_PROCESSES,
    progress_interval_seconds: float = MOTIF_SCAN_PROGRESS_INTERVAL_SECONDS,
    logger: logging.Logger | None = None,
) -> MotifFeatureWriteResult:
    """Scan cluster member sequences and write sparse motif features."""
    if processes < 1:
        raise ValueError("motif scan processes must be at least 1")
    if progress_interval_seconds <= 0.0:
        raise ValueError("progress log interval must be greater than 0")

    log = logger or logging.getLogger(__name__)
    cluster_path = Path(cluster_file)
    motif_path = Path(motif_file)
    fasta_path = Path(fasta_file)
    fasta_index_path = Path(fasta_index_file)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = _temporary_path(output_path)

    member_ids = tuple(_iter_cluster_member_ids(cluster_path))
    cluster_rows = len(member_ids)
    motifs = read_prosig_motif_library(motif_path)
    chunks = tuple(_chunk_member_ids(member_ids, DEFAULT_MOTIF_SCAN_CHUNK_SIZE))
    worker_count = min(processes, len(chunks)) if chunks else 1
    missing_sequences = 0
    missing_preview: list[str] = []
    accessions_scanned = 0
    motifs_with_hits: set[str] = set()
    feature_rows = 0
    accessions_completed = 0
    chunks_completed = 0
    pending_results: dict[int, _AccessionScanResult] = {}
    next_write_index = 1
    shard_paths: set[Path] = set()

    log.info(
        "Scanning %s cluster member rows against %s ProSig motifs using %s process(es)",
        f"{cluster_rows:,}",
        f"{len(motifs):,}",
        f"{worker_count:,}",
    )
    try:
        with temp_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(MOTIF_FEATURE_HEADER)
            last_log_time = time.monotonic()
            tasks = [
                _AccessionScanTask(
                    chunk_index=chunk_index,
                    member_ids=chunk_member_ids,
                )
                for chunk_index, chunk_member_ids in enumerate(chunks, start=1)
            ]
            for result in _iter_accession_scan_results(
                tasks,
                worker_count,
                motifs=motifs,
                fasta_file=fasta_path,
                fasta_index_file=fasta_index_path,
                shard_dir=output_path.parent,
            ):
                result_path = Path(result.shard_file)
                shard_paths.add(result_path)
                pending_results[result.chunk_index] = result
                chunks_completed += 1
                accessions_completed += len(tasks[result.chunk_index - 1].member_ids)
                feature_rows += result.feature_rows
                motifs_with_hits.update(result.motifs_with_hits)
                if result.missing_sequences:
                    missing_sequences += result.missing_sequences
                    missing_preview.extend(
                        result.missing_preview[: 20 - len(missing_preview)]
                    )
                if result.accessions_scanned:
                    accessions_scanned += result.accessions_scanned

                now = time.monotonic()
                if chunks_completed == len(chunks) or _should_log_progress(
                        last_log_time,
                        interval_seconds=progress_interval_seconds,
                ):
                    last_log_time = now
                    _log_motif_scan_progress(
                        log,
                        accessions_completed=accessions_completed,
                        total_accessions=cluster_rows,
                    )

                while next_write_index in pending_results:
                    ordered_result = pending_results.pop(next_write_index)
                    ordered_path = Path(ordered_result.shard_file)
                    with ordered_path.open("r", encoding="utf-8", newline="") as shard:
                        for line in shard:
                            handle.write(line)
                    ordered_path.unlink(missing_ok=True)
                    shard_paths.discard(ordered_path)
                    next_write_index += 1
            handle.write(f"{MOTIF_FEATURE_COMPLETION_MARKER}\n")
        temp_path.chmod(_output_mode(output_path))
        os.replace(temp_path, output_path)
    finally:
        temp_path.unlink(missing_ok=True)
        for shard_path in shard_paths:
            shard_path.unlink(missing_ok=True)

    if missing_sequences:
        log.warning(
            "%s accessions from %s had no sequence in %s: %s",
            f"{missing_sequences:,}",
            cluster_path,
            fasta_path,
            ", ".join(missing_preview),
        )
    log.info(
        "Wrote motif features to %s with %s positive motif-accession rows; "
        "motifs with hits=%s",
        output_path,
        f"{feature_rows:,}",
        f"{len(motifs_with_hits):,}",
    )
    return MotifFeatureWriteResult(
        output_file=output_path,
        cluster_rows=cluster_rows,
        motifs=len(motifs),
        accessions_scanned=accessions_scanned,
        missing_sequences=missing_sequences,
        feature_rows=feature_rows,
        motifs_with_hits=len(motifs_with_hits),
    )


def read_prosig_motif_library(path: str | Path) -> list[Motif]:
    """Read and compile motifs from a ProSig motif library TSV."""
    rows = _read_motif_rows(Path(path))
    motifs: list[Motif] = []
    seen_names: set[str] = set()
    for row_number, row in rows:
        name = str(row.get("name", "")).strip()
        pattern = str(row.get("prosig_pattern", "")).strip()
        description = str(row.get("description", "")).strip()
        status = str(row.get("status", "")).strip()
        if not name:
            raise ValueError(f"{path} row {row_number} is missing motif name")
        if name in seen_names:
            raise ValueError(
                f"{path} row {row_number} has duplicate motif name: {name}"
            )
        if not pattern:
            raise ValueError(f"{path} row {row_number} is missing prosig_pattern")
        if not status:
            raise ValueError(f"{path} row {row_number} is missing status")
        seen_names.add(name)
        min_width, max_width = motif_width_bounds(pattern)
        motifs.append(
            Motif(
                name=name,
                description=description,
                pattern=pattern,
                status=status,
                compiled_pattern=compile_prosig_motif(pattern),
                c_terminal_anchor=pattern.endswith(">"),
                min_width=min_width,
                max_width=max_width,
                literal_prefilter=motif_literal_prefilter(pattern),
            )
        )
    return motifs


def _iter_accession_scan_results(
    tasks: list[_AccessionScanTask],
    worker_count: int,
    *,
    motifs: list[Motif],
    fasta_file: Path,
    fasta_index_file: Path,
    shard_dir: Path,
):
    if worker_count <= 1:
        _initialize_accession_scan_worker(
            tuple(motifs),
            str(fasta_file),
            str(fasta_index_file),
            str(shard_dir),
        )
        for task in tasks:
            yield _scan_accession_chunk_to_shard(task)
        return

    with mp.Pool(
        processes=worker_count,
        initializer=_initialize_accession_scan_worker,
        initargs=(
            tuple(motifs),
            str(fasta_file),
            str(fasta_index_file),
            str(shard_dir),
        ),
    ) as pool:
        yield from pool.imap_unordered(
            _scan_accession_chunk_to_shard,
            tasks,
            chunksize=1,
        )


def _initialize_accession_scan_worker(
    motifs: tuple[Motif, ...],
    fasta_file: str,
    fasta_index_file: str,
    shard_dir: str,
) -> None:
    global _WORKER_MOTIFS
    global _WORKER_FASTA_FILE
    global _WORKER_FASTA_INDEX_FILE
    global _WORKER_SHARD_DIR
    _WORKER_MOTIFS = motifs
    _WORKER_FASTA_FILE = fasta_file
    _WORKER_FASTA_INDEX_FILE = fasta_index_file
    _WORKER_SHARD_DIR = shard_dir


def _scan_accession_chunk_to_shard(
    task: _AccessionScanTask,
) -> _AccessionScanResult:
    descriptor, shard_name = tempfile.mkstemp(
        prefix=f".motif_features.chunk.{task.chunk_index:06d}.",
        suffix=".tsv",
        dir=_WORKER_SHARD_DIR,
    )
    os.close(descriptor)
    shard_path = Path(shard_name)
    feature_rows = 0
    accessions_scanned = 0
    missing_sequences = 0
    missing_preview: list[str] = []
    motifs_with_hits: set[str] = set()
    try:
        with shard_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            for member_id in task.member_ids:
                sequence = indexed_fasta_sequence(
                    member_id,
                    _WORKER_FASTA_FILE,
                    _WORKER_FASTA_INDEX_FILE,
                )
                if sequence is None:
                    missing_sequences += 1
                    if len(missing_preview) < 20:
                        missing_preview.append(member_id)
                    continue
                accessions_scanned += 1
                for motif in _WORKER_MOTIFS:
                    if motif_present(sequence, motif):
                        motifs_with_hits.add(motif.name)
                        writer.writerow([member_id, motif.name])
                        feature_rows += 1
    except Exception:
        shard_path.unlink(missing_ok=True)
        raise

    return _AccessionScanResult(
        chunk_index=task.chunk_index,
        shard_file=str(shard_path),
        feature_rows=feature_rows,
        accessions_scanned=accessions_scanned,
        missing_sequences=missing_sequences,
        missing_preview=tuple(missing_preview),
        motifs_with_hits=tuple(sorted(motifs_with_hits)),
    )


def compile_prosig_motif(pattern: str) -> Pattern[str]:
    """Compile a ProSig motif pattern into a Python regular expression."""
    regex_parts: list[str] = []
    index = 0
    while index < len(pattern):
        char = pattern[index]
        if char == "<":
            regex_parts.append(r"\A")
            index += 1
        elif char == ">":
            regex_parts.append(r"\Z")
            index += 1
        elif char == "?":
            part, index = _parse_wildcard(pattern, index)
            regex_parts.append(part)
        elif char == "[":
            residues, index = _parse_residue_set(pattern, index)
            regex_parts.append(_residue_set_regex(residues))
        elif char == "!":
            part, index = _parse_exclusion(pattern, index)
            regex_parts.append(part)
        elif char == "{":
            residues, index = _parse_macro(pattern, index)
            regex_parts.append(_residue_set_regex(residues))
        elif char.isalpha() and char.isupper():
            regex_parts.append(re.escape(char))
            index += 1
        else:
            raise ValueError(
                f"Unsupported ProSig motif token {char!r} at position {index + 1}"
            )
    return re.compile("".join(regex_parts))


def find_motif_matches(sequence: str, motif: Motif) -> list[MotifMatch]:
    """Find motif matches, allowing overlap and one shortest match per start."""
    matches: list[MotifMatch] = []
    for start_0 in range(len(sequence)):
        match = motif.compiled_pattern.match(sequence, start_0)
        if match is None:
            continue
        matches.append(MotifMatch(start_n=start_0 + 1, end_n=match.end()))
    return matches


def motif_present(sequence: str, motif: Motif) -> bool:
    """Return whether a motif is present, stopping at the first match."""
    if motif.literal_prefilter is not None and motif.literal_prefilter not in sequence:
        return False
    if motif.c_terminal_anchor:
        return _c_terminal_motif_present(sequence, motif)
    for start_0 in range(len(sequence)):
        if motif.compiled_pattern.match(sequence, start_0) is not None:
            return True
    return False


def motif_literal_prefilter(pattern: str) -> str | None:
    """Return a required literal substring suitable for safe prefiltering."""
    longest_literal = ""
    current_literal = ""
    index = 0
    while index < len(pattern):
        char = pattern[index]
        if char.isalpha() and char.isupper():
            current_literal += char
            index += 1
            continue
        if len(current_literal) > len(longest_literal):
            longest_literal = current_literal
        current_literal = ""
        index = _skip_variable_width_token(pattern, index)

    if len(current_literal) > len(longest_literal):
        longest_literal = current_literal
    if len(longest_literal) < MIN_LITERAL_PREFILTER_LENGTH:
        return None
    return longest_literal


def motif_width_bounds(pattern: str) -> tuple[int, int | None]:
    """Return minimum and maximum residue width consumed by a ProSig motif."""
    minimum = 0
    maximum = 0
    index = 0
    while index < len(pattern):
        char = pattern[index]
        if char in "<>":
            index += 1
        elif char == "?":
            min_part, max_part, index = _wildcard_width(pattern, index)
            minimum += min_part
            maximum += max_part
        elif char == "[":
            residues, index = _parse_residue_set(pattern, index)
            if ">" in residues:
                maximum += 1
            else:
                minimum += 1
                maximum += 1
        elif char == "!":
            _part, index = _parse_exclusion(pattern, index)
            minimum += 1
            maximum += 1
        elif char == "{":
            _residues, index = _parse_macro(pattern, index)
            minimum += 1
            maximum += 1
        elif char.isalpha() and char.isupper():
            minimum += 1
            maximum += 1
            index += 1
        else:
            raise ValueError(
                f"Unsupported ProSig motif token {char!r} at position {index + 1}"
            )
    return minimum, maximum


def _c_terminal_motif_present(sequence: str, motif: Motif) -> bool:
    if not sequence:
        return False

    sequence_length = len(sequence)
    latest_start = sequence_length - motif.min_width
    if latest_start < 0:
        return False
    if motif.max_width is None:
        earliest_start = 0
    else:
        earliest_start = max(0, sequence_length - motif.max_width)

    for start_0 in range(latest_start, earliest_start - 1, -1):
        if motif.compiled_pattern.match(sequence, start_0) is not None:
            return True
    return False


def extract_motif_features(sequence: str, motif: Motif) -> MotifFeatureResult:
    """Extract count, first position, last position, and match fraction."""
    sequence_length = len(sequence)
    if sequence_length == 0:
        return MotifFeatureResult(0, 0, 0, "NA")

    count = 0
    first_position = 0
    last_start_n = 0
    for start_0 in range(sequence_length):
        match = motif.compiled_pattern.match(sequence, start_0)
        if match is None:
            continue
        start_n = start_0 + 1
        count += 1
        if first_position == 0:
            first_position = start_n
        last_start_n = start_n

    if count == 0:
        return MotifFeatureResult(0, 0, 0, 0.0)

    return MotifFeatureResult(
        count=count,
        first_position=first_position,
        last_position=sequence_length - last_start_n + 1,
        match_fraction=count / sequence_length,
    )


def _chunk_member_ids(
    member_ids: tuple[str, ...],
    chunk_size: int,
):
    for start in range(0, len(member_ids), chunk_size):
        yield member_ids[start : start + chunk_size]


def _iter_cluster_member_ids(cluster_file: Path):
    with cluster_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        if "member_id" not in fieldnames:
            raise ValueError(f"{cluster_file} is missing required column: member_id")
        for row in reader:
            member_id = str(row.get("member_id", "")).strip()
            if member_id:
                yield member_id


def _count_cluster_member_rows(cluster_file: Path) -> int:
    return sum(1 for _member_id in _iter_cluster_member_ids(cluster_file))


def _read_motif_rows(path: Path) -> list[tuple[int, dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        data_lines = [
            (line_number, line)
            for line_number, line in enumerate(handle, start=1)
            if not line.startswith("#")
        ]
    if not data_lines:
        raise ValueError(f"{path} does not contain a motif TSV header")
    header = data_lines[0][1].rstrip("\n").split("\t")
    missing_columns = PROSIG_MOTIF_COLUMNS - set(header)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"{path} is missing required column(s): {missing}")
    reader = csv.DictReader([line for _, line in data_lines], delimiter="\t")
    return [
        (line_number, {key: value or "" for key, value in row.items()})
        for (line_number, _), row in zip(data_lines[1:], reader, strict=False)
    ]


def motif_features_complete(path: str | Path) -> bool:
    """Return whether a motif feature TSV has the end-of-file completion marker."""
    feature_path = Path(path)
    try:
        with feature_path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - 4096))
            tail = handle.read().decode("utf-8", errors="replace")
    except FileNotFoundError:
        return False
    return any(
        line.strip() == MOTIF_FEATURE_COMPLETION_MARKER
        for line in tail.splitlines()
    )


def _parse_wildcard(pattern: str, index: int) -> tuple[str, int]:
    next_index = index + 1
    if next_index >= len(pattern) or pattern[next_index] != "(":
        return RESIDUE_CLASS, next_index
    close_index = pattern.find(")", next_index)
    if close_index == -1:
        raise ValueError(f"Unclosed wildcard run at position {index + 1}")
    quantifier = pattern[next_index + 1 : close_index]
    numbers = quantifier.split(",")
    token = pattern[index : close_index + 1]
    if len(numbers) == 1:
        length = _parse_nonnegative_int(numbers[0], token=token)
        return f"{RESIDUE_CLASS}{{{length}}}", close_index + 1
    if len(numbers) == 2:
        minimum = _parse_nonnegative_int(numbers[0], token=token)
        maximum = _parse_nonnegative_int(numbers[1], token=token)
        if minimum > maximum:
            raise ValueError(f"Invalid wildcard range at position {index + 1}")
        return f"{RESIDUE_CLASS}{{{minimum},{maximum}}}?", close_index + 1
    raise ValueError(f"Invalid wildcard run at position {index + 1}")


def _wildcard_width(pattern: str, index: int) -> tuple[int, int, int]:
    next_index = index + 1
    if next_index >= len(pattern) or pattern[next_index] != "(":
        return 1, 1, next_index
    close_index = pattern.find(")", next_index)
    if close_index == -1:
        raise ValueError(f"Unclosed wildcard run at position {index + 1}")
    quantifier = pattern[next_index + 1 : close_index]
    numbers = quantifier.split(",")
    token = pattern[index : close_index + 1]
    if len(numbers) == 1:
        length = _parse_nonnegative_int(numbers[0], token=token)
        return length, length, close_index + 1
    if len(numbers) == 2:
        minimum = _parse_nonnegative_int(numbers[0], token=token)
        maximum = _parse_nonnegative_int(numbers[1], token=token)
        if minimum > maximum:
            raise ValueError(f"Invalid wildcard range at position {index + 1}")
        return minimum, maximum, close_index + 1
    raise ValueError(f"Invalid wildcard run at position {index + 1}")


def _skip_variable_width_token(pattern: str, index: int) -> int:
    char = pattern[index]
    if char in "<>":
        return index + 1
    if char == "?":
        _min_part, _max_part, next_index = _wildcard_width(pattern, index)
        return next_index
    if char == "[":
        _residues, next_index = _parse_residue_set(pattern, index)
        return next_index
    if char == "!":
        _part, next_index = _parse_exclusion(pattern, index)
        return next_index
    if char == "{":
        _residues, next_index = _parse_macro(pattern, index)
        return next_index
    raise ValueError(
        f"Unsupported ProSig motif token {char!r} at position {index + 1}"
    )


def _parse_exclusion(pattern: str, index: int) -> tuple[str, int]:
    next_index = index + 1
    if next_index >= len(pattern):
        raise ValueError(f"Dangling exclusion token at position {index + 1}")
    char = pattern[next_index]
    if char == "[":
        residues, next_index = _parse_residue_set(pattern, next_index)
    elif char == "{":
        residues, next_index = _parse_macro(pattern, next_index)
    elif char.isalpha() and char.isupper():
        residues = char
        next_index += 1
    else:
        raise ValueError(f"Unsupported exclusion token at position {index + 1}")
    if ">" in residues:
        raise ValueError("Terminal anchors are not supported inside exclusion tokens")
    if "?" in residues:
        raise ValueError("Wildcard exclusions cannot match any residue")
    return f"(?![{re.escape(residues)}]){RESIDUE_CLASS}", next_index


def _parse_residue_set(pattern: str, index: int) -> tuple[str, int]:
    close_index = pattern.find("]", index + 1)
    if close_index == -1:
        raise ValueError(f"Unclosed residue set at position {index + 1}")
    residues = pattern[index + 1 : close_index]
    if not residues:
        raise ValueError(f"Empty residue set at position {index + 1}")
    unsupported = set(residues) - set(RESIDUE_ALPHABET) - {">", "?"}
    if unsupported:
        bad = "".join(sorted(unsupported))
        raise ValueError(f"Unsupported residue set token(s): {bad}")
    return residues, close_index + 1


def _parse_macro(pattern: str, index: int) -> tuple[str, int]:
    close_index = pattern.find("}", index + 1)
    if close_index == -1:
        raise ValueError(f"Unclosed macro at position {index + 1}")
    macro_name = pattern[index + 1 : close_index]
    try:
        residues = MACRO_RESIDUES[macro_name]
    except KeyError as exc:
        raise ValueError(f"Unknown macro {{{macro_name}}}") from exc
    return residues, close_index + 1


def _residue_set_regex(residues: str) -> str:
    if "?" in residues and ">" in residues:
        return f"(?:{RESIDUE_CLASS}|\\Z)"
    if "?" in residues:
        return RESIDUE_CLASS
    if ">" in residues:
        residue_options = residues.replace(">", "")
        if not residue_options:
            return r"\Z"
        return f"(?:[{re.escape(residue_options)}]|\\Z)"
    return f"[{re.escape(residues)}]"


def _parse_nonnegative_int(value: str, *, token: str) -> int:
    if not value.isdecimal():
        raise ValueError(f"Invalid non-negative integer in token {token!r}")
    return int(value)


def _format_motif_result(result: MotifFeatureResult) -> list[str]:
    return [
        str(result.count),
        str(result.first_position),
        str(result.last_position),
        _format_match_fraction(result.match_fraction),
    ]


def _format_match_fraction(value: float | str) -> str:
    if isinstance(value, str):
        return value
    return f"{value:.12g}"


def _should_log_progress(last_log_time: float, *, interval_seconds: float) -> bool:
    return time.monotonic() - last_log_time >= interval_seconds


def _log_motif_extraction_progress(
    logger: logging.Logger,
    *,
    accessions_completed: int,
    total_accessions: int,
) -> None:
    logger.info(
        "Motif scan progress: %s/%s accessions completed",
        f"{accessions_completed:,}",
        f"{total_accessions:,}",
    )


_log_motif_scan_progress = _log_motif_extraction_progress


def _temporary_path(destination: Path) -> Path:
    descriptor, name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
    )
    os.close(descriptor)
    return Path(name)


def _output_mode(destination: Path) -> int:
    if destination.exists():
        return stat.S_IMODE(destination.stat().st_mode)
    current_umask = os.umask(0)
    os.umask(current_umask)
    return 0o666 & ~current_umask
