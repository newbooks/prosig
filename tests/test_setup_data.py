import logging
from pathlib import Path

from typer.testing import CliRunner

from prosig.cli.app import app
from prosig.cli.setup_data import (
    SETUP_DATA_SOURCES,
    SetupDataSource,
    format_progress,
    setup_data_sources,
)
from prosig.io.download import DownloadResult


def test_setup_data_help_includes_threads_option() -> None:
    result = CliRunner().invoke(app, ["setup-data", "-h"])

    assert result.exit_code == 0
    assert "--threads" in result.stdout
    assert "--force" in result.stdout
    assert "--dry-run" in result.stdout


def test_default_sources_include_prosite_dat() -> None:
    assert SetupDataSource(
        description="PROSITE",
        url="https://ftp.expasy.org/databases/prosite/prosite.dat",
        destination="prosite.dat",
    ) in SETUP_DATA_SOURCES


def test_default_sources_do_not_include_swissprot_fasta() -> None:
    assert all(
        source.destination != "uniprot_sprot.fasta.gz"
        for source in SETUP_DATA_SOURCES
    )


def test_setup_data_skips_existing_destination_without_force(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger="prosig")
    (tmp_path / "go-basic.obo").write_text("existing", encoding="utf-8")
    downloads = []

    setup_data_sources(
        force=False,
        dry_run=False,
        threads=16,
        sources=[
            SetupDataSource(
                description="GO Graph",
                url="https://example.test/go-basic.obo",
                destination="go-basic.obo",
            )
        ],
        downloader=lambda *args, **kwargs: downloads.append((args, kwargs)),
    )

    assert downloads == []
    assert caplog.messages == [
        "Skipped go-basic.obo: already exists. "
        "Use --force to overwrite existing files."
    ]


def test_setup_data_force_downloads_existing_destination(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger="prosig")
    destination = tmp_path / "go-basic.obo"
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

    setup_data_sources(
        force=True,
        dry_run=False,
        threads=4,
        sources=[
            SetupDataSource(
                description="GO Graph",
                url="https://example.test/go-basic.obo",
                destination="go-basic.obo",
            )
        ],
        downloader=fake_downloader,
    )

    assert destination.read_text(encoding="utf-8") == "new"
    assert calls[0][0] == "https://example.test/go-basic.obo"
    assert calls[0][1] == destination
    assert calls[0][2]["threads"] == 4
    assert "Download threads requested: 4" in caplog.messages
    assert (
        "Downloaded go-basic.obo: 3 bytes using single-threaded download"
        in caplog.messages
    )


def test_setup_data_dry_run_does_not_download(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO, logger="prosig")
    downloads = []

    setup_data_sources(
        force=False,
        dry_run=True,
        threads=16,
        sources=[
            SetupDataSource(
                description="Swiss-Prot GO",
                url="https://example.test/uniprot_sprot.dat.gz",
                destination="uniprot_sprot.dat.gz",
            )
        ],
        downloader=lambda *args, **kwargs: downloads.append((args, kwargs)),
    )

    assert downloads == []
    assert caplog.messages == [
        "Would download Swiss-Prot GO: "
        "https://example.test/uniprot_sprot.dat.gz -> uniprot_sprot.dat.gz"
    ]


def test_format_progress_with_known_total() -> None:
    assert format_progress(3_299_328, 126_917_632) == (
        "Downloaded 3,222 KB / 123,943 KB (2.6%)"
    )


def test_format_progress_without_known_total() -> None:
    assert format_progress(3_299_328, None) == "Downloaded 3,222 KB"
