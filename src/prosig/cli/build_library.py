from pathlib import Path
from typing import Annotated

import typer

from prosig.cli.logging import get_logger
from prosig.go.build import MF_ROOT, build_go_pkl, format_log_number


def build_library(
    go_obo: Annotated[
        Path,
        typer.Option("--go-obo", help="Path to the GO OBO file."),
    ] = Path("go-basic.obo"),
    swissprot: Annotated[
        Path,
        typer.Option(
            "--swissprot",
            help="Path to the Swiss-Prot flat-file annotation source.",
        ),
    ] = Path("uniprot_sprot.dat.gz"),
    go_out: Annotated[
        Path,
        typer.Option(
            "--go-out",
            help="Path to write the compact GO graph and IC artifact.",
        ),
    ] = Path("go_graph.pkl"),
    write_report: Annotated[
        Path | None,
        typer.Option(
            "--write-report",
            help="Optional path to write a build validation report.",
        ),
    ] = None,
) -> None:
    """Build the compact GO graph and IC artifact."""
    artifact = build_go_pkl(
        go_obo=go_obo,
        swissprot=swissprot,
        go_out=go_out,
        report_out=write_report,
    )
    meta = artifact["meta"]
    logger = get_logger()
    logger.info(
        "Wrote GO graph and IC artifact to %s with %s MF terms",
        go_out,
        format_log_number(meta["n_terms"]),
    )
    logger.info(
        "Used %s of %s provided accessions for IC",
        format_log_number(artifact["terms"][MF_ROOT]["count"]),
        format_log_number(meta["n_accessions_provided"]),
    )
    if write_report is not None:
        logger.info("Wrote build report to %s", write_report)
