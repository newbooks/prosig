from __future__ import annotations

from pathlib import Path

import pytest

from prosig.sequence_identity import run_needle, sequence_identity


def test_sequence_identity_parses_needle_identity() -> None:
    def fake_needle(first: Path, second: Path) -> str:
        assert first.read_text(encoding="ascii") == ">first\nAAAA\n"
        assert second.read_text(encoding="ascii") == ">second\nAAAT\n"
        return "# Identity:       3/4 (75.0%)\n"

    assert sequence_identity("aaaa", "aaat", needle_runner=fake_needle) == 0.75


@pytest.mark.parametrize("first,second", [("", "AAAA"), ("AAAA", "")])
def test_sequence_identity_rejects_empty_sequences(first: str, second: str) -> None:
    with pytest.raises(ValueError, match="must not be empty"):
        sequence_identity(first, second, needle_runner=lambda _a, _b: "")


def test_sequence_identity_rejects_unparseable_output() -> None:
    with pytest.raises(RuntimeError, match="could not find identity"):
        sequence_identity("AAAA", "AAAT", needle_runner=lambda _a, _b: "bad output")


def test_run_needle_reports_missing_program(monkeypatch) -> None:
    monkeypatch.setattr("prosig.sequence_identity.shutil.which", lambda _name: None)

    with pytest.raises(RuntimeError, match="EMBOSS needle was not found on PATH"):
        run_needle(Path("first.fasta"), Path("second.fasta"))
