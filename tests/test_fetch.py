import logging
from pathlib import Path

from typer.testing import CliRunner

from prosig.cli.app import app
from prosig.cli.fetch import FETCH_SOURCES, FetchSource, fetch_sources, format_progress
from prosig.io.download import DownloadResult


def test_fetch_help_includes_threads_option() -> None:
    result = CliRunner().invoke(app, ["fetch", "-h"])

    assert result.exit_code == 0
    assert "--threads" in result.stdout
    assert "--force" in result.stdout
    assert "--dry-run" in result.stdout


def test_default_sources_include_prosite_dat() -> None:
    assert FetchSource(
        description="PROSITE",
        url="https://ftp.expasy.org/databases/prosite/prosite.dat",
        destination="prosite.dat",
    ) in FETCH_SOURCES


def test_fetch_skips_existing_destination_without_force(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger="prosig")
    (tmp_path / "go.obo").write_text("existing", encoding="utf-8")
    downloads = []

    fetch_sources(
        force=False,
        dry_run=False,
        threads=16,
        sources=[
            FetchSource(
                description="GO Graph",
                url="https://example.test/go.obo",
                destination="go.obo",
            )
        ],
        downloader=lambda *args, **kwargs: downloads.append((args, kwargs)),
    )

    assert downloads == []
    assert caplog.messages == [
        "Skipped go.obo: already exists. "
        "Use --force to overwrite existing files."
    ]


def test_fetch_force_downloads_existing_destination(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger="prosig")
    destination = tmp_path / "go.obo"
    destination.write_text("existing", encoding="utf-8")
    calls = []

    def fake_downloader(url: str, path: Path, **kwargs: object) -> DownloadResult:
        calls.append((url, path, kwargs))
        path.write_text("new", encoding="utf-8")
        return DownloadResult(
            url=url,
            destination=path,
            bytes_written=3,
            content_length=3,
            threaded=False,
        )

    fetch_sources(
        force=True,
        dry_run=False,
        threads=4,
        sources=[
            FetchSource(
                description="GO Graph",
                url="https://example.test/go.obo",
                destination="go.obo",
            )
        ],
        downloader=fake_downloader,
    )

    assert destination.read_text(encoding="utf-8") == "new"
    assert calls[0][0] == "https://example.test/go.obo"
    assert calls[0][1] == destination
    assert calls[0][2]["threads"] == 4
    assert "Download threads requested: 4" in caplog.messages
    assert (
        "Downloaded go.obo: 3 bytes using single-threaded download"
        in caplog.messages
    )


def test_fetch_dry_run_does_not_download(tmp_path: Path, monkeypatch, caplog) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger="prosig")
    downloads = []

    fetch_sources(
        force=False,
        dry_run=True,
        threads=16,
        sources=[
            FetchSource(
                description="Swiss-Prot fasta",
                url="https://example.test/uniprot_sprot.fasta.gz",
                destination="uniprot_sprot.fasta.gz",
            )
        ],
        downloader=lambda *args, **kwargs: downloads.append((args, kwargs)),
    )

    assert downloads == []
    assert caplog.messages == [
        "Would download Swiss-Prot fasta: "
        "https://example.test/uniprot_sprot.fasta.gz -> uniprot_sprot.fasta.gz"
    ]


def test_format_progress_with_known_total() -> None:
    assert format_progress(3_299_328, 126_917_632) == (
        "Downloaded 3,222 KB / 123,943 KB (2.6%)"
    )


def test_format_progress_without_known_total() -> None:
    assert format_progress(3_299_328, None) == "Downloaded 3,222 KB"
