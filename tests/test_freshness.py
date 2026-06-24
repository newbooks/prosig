import os
from pathlib import Path

import pytest

from prosig.io.freshness import artifact_is_stale


def test_artifact_is_stale_when_output_missing(tmp_path: Path) -> None:
    dependency = tmp_path / "source.txt"
    dependency.write_text("source", encoding="utf-8")

    assert artifact_is_stale(tmp_path / "missing.txt", [dependency])


def test_artifact_is_stale_when_dependency_is_newer(tmp_path: Path) -> None:
    output = tmp_path / "output.txt"
    dependency = tmp_path / "source.txt"
    output.write_text("output", encoding="utf-8")
    dependency.write_text("source", encoding="utf-8")
    os.utime(output, (100, 100))
    os.utime(dependency, (200, 200))

    assert artifact_is_stale(output, [dependency])


def test_artifact_is_current_when_output_is_newer(tmp_path: Path) -> None:
    output = tmp_path / "output.txt"
    dependency = tmp_path / "source.txt"
    output.write_text("output", encoding="utf-8")
    dependency.write_text("source", encoding="utf-8")
    os.utime(dependency, (100, 100))
    os.utime(output, (200, 200))

    assert not artifact_is_stale(output, [dependency])


def test_artifact_is_stale_when_force_is_true(tmp_path: Path) -> None:
    output = tmp_path / "output.txt"
    dependency = tmp_path / "source.txt"
    output.write_text("output", encoding="utf-8")
    dependency.write_text("source", encoding="utf-8")
    os.utime(dependency, (100, 100))
    os.utime(output, (200, 200))

    assert artifact_is_stale(output, [dependency], force=True)


def test_artifact_is_stale_requires_existing_dependencies(tmp_path: Path) -> None:
    output = tmp_path / "output.txt"
    output.write_text("output", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="Dependency file not found"):
        artifact_is_stale(output, [tmp_path / "missing.txt"])
