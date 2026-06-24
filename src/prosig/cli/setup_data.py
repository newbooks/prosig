from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import typer

from prosig.cli.logging import get_logger
from prosig.io.download import DEFAULT_THREADS, DownloadResult, download_file


@dataclass(frozen=True)
class SetupDataSource:
    description: str
    url: str
    destination: str


SETUP_DATA_SOURCES = [
    SetupDataSource(
        description="GO Graph",
        url="https://current.geneontology.org/ontology/go-basic.obo",
        destination="go-basic.obo",
    ),
    SetupDataSource(
        description="Swiss-Prot GO",
        url=(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/complete/uniprot_sprot.dat.gz"
        ),
        destination="uniprot_sprot.dat.gz",
    ),
    SetupDataSource(
        description="PROSITE",
        url="https://ftp.expasy.org/databases/prosite/prosite.dat",
        destination="prosite.dat",
    ),
]

Downloader = Callable[..., DownloadResult]


def setup_data(
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Overwrite existing files.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show what would be fetched without writing files.",
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS,
        "--threads",
        min=1,
        help="Number of download threads to use when supported.",
    ),
) -> None:
    """Download external files required by ProSig workflows."""
    setup_data_sources(force=force, dry_run=dry_run, threads=threads)


def setup_data_sources(
    *,
    force: bool,
    dry_run: bool,
    threads: int,
    sources: list[SetupDataSource] | None = None,
    downloader: Downloader = download_file,
) -> list[DownloadResult]:
    if threads <= 0:
        raise ValueError("threads must be greater than zero")

    logger = get_logger()
    results = []
    for source in sources or SETUP_DATA_SOURCES:
        destination = Path.cwd() / source.destination
        if destination.exists() and not force:
            logger.info(
                f"Skipped {destination.name}: already exists. "
                "Use --force to overwrite existing files."
            )
            continue

        if dry_run:
            action = "Would overwrite" if destination.exists() else "Would download"
            logger.info(
                f"{action} {source.description}: {source.url} -> {destination.name}"
            )
            continue

        logger.info(
            f"Downloading {source.description}: {source.url} -> {destination.name}"
        )
        logger.info(f"Download threads requested: {threads}")
        result = downloader(
            source.url,
            destination,
            threads=threads,
            progress=_progress_logger(logger.info),
        )
        mode = "multi-threaded" if result.threaded else "single-threaded"
        logger.info(
            f"Downloaded {destination.name}: {result.bytes_written:,} bytes "
            f"using {mode} download"
        )
        results.append(result)

    return results


def _progress_logger(
    log_message: Callable[[str], None],
) -> Callable[[int, int | None], None]:
    def log(bytes_written: int, content_length: int | None) -> None:
        log_message(format_progress(bytes_written, content_length))

    return log


def format_progress(bytes_written: int, content_length: int | None) -> str:
    downloaded_kb = _bytes_to_kb(bytes_written)
    if content_length is None:
        return f"Downloaded {downloaded_kb:,} KB"

    total_kb = _bytes_to_kb(content_length)
    percent = bytes_written / content_length * 100 if content_length else 100.0
    return f"Downloaded {downloaded_kb:,} KB / {total_kb:,} KB ({percent:.1f}%)"


def _bytes_to_kb(value: int) -> int:
    return value // 1024
