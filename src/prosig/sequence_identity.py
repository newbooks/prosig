"""Protein sequence identity calculated by EMBOSS Needle global alignment."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from collections.abc import Callable
from pathlib import Path


IDENTITY_PATTERN = re.compile(r"^# Identity:\s*(\d+)\s*/\s*(\d+)\b", re.MULTILINE)
NeedleRunner = Callable[[Path, Path], str]


def write_fasta(path: Path, accession: str, sequence: str) -> None:
    """Write one protein sequence as wrapped FASTA."""
    wrapped_sequence = "\n".join(
        sequence[offset : offset + 60] for offset in range(0, len(sequence), 60)
    )
    path.write_text(f">{accession}\n{wrapped_sequence}\n", encoding="ascii")


def run_needle(first_fasta: Path, second_fasta: Path) -> str:
    """Run EMBOSS Needle in explicit protein mode and return its output."""
    if shutil.which("needle") is None:
        raise RuntimeError("EMBOSS needle was not found on PATH")

    command = [
        "needle",
        "-asequence",
        str(first_fasta),
        "-bsequence",
        str(second_fasta),
        "-sprotein1",
        "-sprotein2",
        "-stdout",
        "-auto",
    ]
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        detail = (
            exc.stderr.strip()
            or exc.stdout.strip()
            or f"exit status {exc.returncode}"
        )
        raise RuntimeError(f"needle failed: {detail}") from exc
    return result.stdout


def sequence_identity(
    first_sequence: str,
    second_sequence: str,
    *,
    needle_runner: NeedleRunner | None = None,
) -> float:
    """Globally align two protein sequences and return exact identity fraction."""
    if not first_sequence or not second_sequence:
        raise ValueError("sequences must not be empty")
    if not first_sequence.isalpha() or not second_sequence.isalpha():
        raise ValueError("sequences must contain letters only")

    runner = needle_runner or run_needle
    with tempfile.TemporaryDirectory(prefix="prosig-needle-") as temp_dir:
        temp_path = Path(temp_dir)
        first_fasta = temp_path / "first.fasta"
        second_fasta = temp_path / "second.fasta"
        write_fasta(first_fasta, "first", first_sequence.upper())
        write_fasta(second_fasta, "second", second_sequence.upper())
        output = runner(first_fasta, second_fasta)

    match = IDENTITY_PATTERN.search(output)
    if match is None:
        raise RuntimeError("could not find identity in needle output")
    identical, alignment_length = map(int, match.groups())
    if alignment_length == 0:
        raise RuntimeError("needle returned a zero-length alignment")
    return identical / alignment_length
