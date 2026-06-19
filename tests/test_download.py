from pathlib import Path

import pytest

from prosig.io.download import (
    MIN_THREADED_BYTES,
    DownloadError,
    DownloadMetadata,
    download_file,
    should_use_threaded_download,
)


class FakeHeaders:
    def __init__(self, values: dict[str, str] | None = None) -> None:
        self.values = values or {}

    def get(self, key: str) -> str | None:
        return self.values.get(key)


class FakeResponse:
    def __init__(
        self,
        chunks: list[bytes],
        *,
        status: int = 200,
        headers: dict[str, str] | None = None,
        fail_after_reads: int | None = None,
    ) -> None:
        self.chunks = chunks
        self.status = status
        self.headers = FakeHeaders(headers)
        self.fail_after_reads = fail_after_reads
        self.read_count = 0

    def __enter__(self) -> "FakeResponse":
        return self

    def __exit__(self, *args: object) -> None:
        return None

    def read(self, chunk_size: int) -> bytes:
        self.read_count += 1
        should_fail = (
            self.fail_after_reads is not None
            and self.read_count > self.fail_after_reads
        )
        if should_fail:
            raise OSError("connection dropped")
        if not self.chunks:
            return b""
        return self.chunks.pop(0)


def test_download_file_overwrites_destination_atomically(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "nested" / "artifact.tsv"
    destination.parent.mkdir()
    destination.write_text("old", encoding="utf-8")
    progress_calls: list[tuple[int, int | None]] = []

    def fake_urlopen(request: object, timeout: float) -> FakeResponse:
        return FakeResponse(
            [b"new", b" content"],
            headers={"Content-Length": "11"},
        )

    monkeypatch.setattr("prosig.io.download.urlopen", fake_urlopen)

    result = download_file(
        "https://example.test/artifact.tsv",
        destination,
        chunk_size=3,
        progress=lambda written, total: progress_calls.append((written, total)),
    )

    assert destination.read_text(encoding="utf-8") == "new content"
    assert not destination.with_name("artifact.tsv.part").exists()
    assert result.destination == destination
    assert result.bytes_written == 11
    assert result.content_length == 11
    assert progress_calls == [(11, 11)]
    assert result.threaded is False


def test_download_file_keeps_existing_destination_when_download_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "artifact.json"
    destination.write_text("old", encoding="utf-8")

    def fake_urlopen(request: object, timeout: float) -> FakeResponse:
        return FakeResponse([b"partial"], fail_after_reads=1)

    monkeypatch.setattr("prosig.io.download.urlopen", fake_urlopen)

    with pytest.raises(DownloadError):
        download_file("https://example.test/artifact.json", destination)

    assert destination.read_text(encoding="utf-8") == "old"
    assert not destination.with_name("artifact.json.part").exists()


def test_download_file_raises_on_http_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "missing.pkl"

    def fake_urlopen(request: object, timeout: float) -> FakeResponse:
        return FakeResponse([], status=404)

    monkeypatch.setattr("prosig.io.download.urlopen", fake_urlopen)

    with pytest.raises(DownloadError, match="HTTP status 404"):
        download_file("https://example.test/missing.pkl", destination)

    assert not destination.exists()


def test_download_file_rejects_incomplete_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "artifact.tsv"
    destination.write_text("old", encoding="utf-8")

    def fake_urlopen(request: object, timeout: float) -> FakeResponse:
        return FakeResponse(
            [b"short"],
            headers={"Content-Length": "10"},
        )

    monkeypatch.setattr("prosig.io.download.urlopen", fake_urlopen)

    with pytest.raises(DownloadError, match="expected 10"):
        download_file("https://example.test/artifact.tsv", destination)

    assert destination.read_text(encoding="utf-8") == "old"
    assert not destination.with_name("artifact.tsv.part").exists()


def test_should_use_threaded_download_for_large_ranged_file() -> None:
    metadata = DownloadMetadata(
        content_length=MIN_THREADED_BYTES + 1,
        accepts_ranges=True,
    )

    assert should_use_threaded_download(metadata, threads=16)


def test_should_not_use_threaded_download_without_range_support() -> None:
    metadata = DownloadMetadata(
        content_length=MIN_THREADED_BYTES + 1,
        accepts_ranges=False,
    )

    assert not should_use_threaded_download(metadata, threads=16)


def test_download_file_throttles_progress_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "artifact.tsv"
    progress_calls: list[tuple[int, int | None]] = []
    times = iter([0.0, 10.0, 31.0, 32.0])

    def fake_urlopen(request: object, timeout: float) -> FakeResponse:
        return FakeResponse(
            [b"abc", b"def", b"ghi"],
            headers={"Content-Length": "9"},
        )

    monkeypatch.setattr("prosig.io.download.urlopen", fake_urlopen)

    download_file(
        "https://example.test/artifact.tsv",
        destination,
        chunk_size=3,
        progress=lambda written, total: progress_calls.append((written, total)),
        clock=lambda: next(times),
    )

    assert progress_calls == [(6, 9), (9, 9)]
