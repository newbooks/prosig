from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import prosig.library as library


class _TemporaryPackagedResource:
    def __init__(self) -> None:
        self._temporary_directory: TemporaryDirectory | None = None

    def __enter__(self) -> Path:
        self._temporary_directory = TemporaryDirectory()
        package_dir = Path(self._temporary_directory.name)
        for filename in library.CORE_LIBRARY_FILES:
            (package_dir / filename).write_text("", encoding="utf-8")
        return package_dir

    def __exit__(self, *args) -> None:
        assert self._temporary_directory is not None
        self._temporary_directory.cleanup()


def test_packaged_default_resource_paths_remain_alive_until_library_close(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(library, "files", lambda package: object())
    monkeypatch.setattr(
        library,
        "as_file",
        lambda resource: _TemporaryPackagedResource(),
    )

    resolved = library.resolve_core_library()

    assert resolved.source == "packaged-default"
    assert resolved.directory.exists()
    assert all(
        resolved.path(filename).is_file()
        for filename in library.CORE_LIBRARY_FILES
    )

    package_dir = resolved.directory
    resolved.close()

    assert not package_dir.exists()
