"""Runtime library artifact resolution for ProSig commands."""

from __future__ import annotations

import shutil
from contextlib import ExitStack
from dataclasses import dataclass, field
from importlib.resources import as_file, files
from pathlib import Path
from types import TracebackType

CORE_LIBRARY_FILES = (
    "prosig_motifs.tsv",
    "motif_cluster_scoreboard.pkl",
    "motif_cluster_scoreboard_meta.json",
    "clusters_meta.tsv",
    "go_graph.pkl",
    "accession_mf_go.tsv",
)


@dataclass(frozen=True)
class ResolvedLibrary:
    """Resolved all-or-nothing runtime library artifact paths."""

    source: str
    directory: Path
    files: dict[str, Path]
    _resource_stack: ExitStack | None = field(
        default=None,
        repr=False,
        compare=False,
    )

    def path(self, filename: str) -> Path:
        return self.files[filename]

    def close(self) -> None:
        """Release any temporary importlib resource extraction."""
        if self._resource_stack is not None:
            self._resource_stack.close()

    def __enter__(self) -> ResolvedLibrary:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()


def package_default_library_dir() -> Path:
    """Return the source/package default runtime library directory."""
    return Path(__file__).resolve().parent / "default"


def resolve_core_library(library_dir: Path | None = None) -> ResolvedLibrary:
    """Resolve core runtime library files with all-or-nothing semantics."""
    if library_dir is not None:
        return _resolve_directory(library_dir, source="user-specified")

    cwd = Path.cwd()
    cwd_files = [cwd / filename for filename in CORE_LIBRARY_FILES]
    if any(path.exists() for path in cwd_files):
        return _resolve_directory(cwd, source="current-directory")

    return _resolve_packaged_default()


def package_core_library(
    *,
    source_dir: Path,
    target_dir: Path | None = None,
) -> ResolvedLibrary:
    """Copy core runtime artifacts from a build directory into package data."""
    target = target_dir or package_default_library_dir()
    source_library = _resolve_directory(source_dir, source="package-source")
    target.mkdir(parents=True, exist_ok=True)
    for filename in CORE_LIBRARY_FILES:
        shutil.copy2(source_library.path(filename), target / filename)
    return _resolve_directory(target, source="packaged-target")


def _resolve_directory(directory: Path, *, source: str) -> ResolvedLibrary:
    resolved_dir = directory.resolve()
    missing = [
        filename
        for filename in CORE_LIBRARY_FILES
        if not (resolved_dir / filename).is_file()
    ]
    if missing:
        missing_text = ", ".join(missing)
        raise FileNotFoundError(
            f"library directory {resolved_dir} is missing required file(s): "
            f"{missing_text}"
        )
    return ResolvedLibrary(
        source=source,
        directory=resolved_dir,
        files={
            filename: resolved_dir / filename
            for filename in CORE_LIBRARY_FILES
        },
    )


def _resolve_packaged_default() -> ResolvedLibrary:
    resource = files("prosig.library.default")
    stack = ExitStack()
    try:
        package_dir = stack.enter_context(as_file(resource))
        resolved = _resolve_directory(package_dir, source="packaged-default")
        return ResolvedLibrary(
            source=resolved.source,
            directory=resolved.directory,
            files=resolved.files,
            _resource_stack=stack,
        )
    except Exception:
        stack.close()
        raise
