"""Artifact freshness checks for build workflows."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path


def artifact_is_stale(
    output: Path,
    dependencies: Iterable[Path],
    *,
    force: bool = False,
) -> bool:
    """Return whether output should be rebuilt from dependencies."""
    if force or not output.exists():
        return True

    output_mtime = output.stat().st_mtime
    for dependency in dependencies:
        if not dependency.exists():
            raise FileNotFoundError(f"Dependency file not found: {dependency}")
        if dependency.stat().st_mtime > output_mtime:
            return True
    return False
