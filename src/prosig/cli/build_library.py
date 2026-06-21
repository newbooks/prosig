from pathlib import Path
from typing import Annotated

import typer

from prosig.cli.logging import get_logger
from prosig.go.build import MF_ROOT, build_go_pkl, format_log_number
from prosig.motifs.prosite import write_prosig_motif_library


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
    prosite_dat: Annotated[
        Path,
        typer.Option(
            "--prosite-dat",
            help="Path to the PROSITE dat source file for motif translation.",
        ),
    ] = Path("prosite.dat"),
    motif_out: Annotated[
        Path,
        typer.Option(
            "--motif-out",
            help="Path to write the translated ProSig motif library TSV.",
        ),
    ] = Path("prosig_motifs.tsv"),
    write_report: Annotated[
        Path | None,
        typer.Option(
            "--write-report",
            help="Optional path to write a build validation report.",
        ),
    ] = None,
) -> None:
    """Build the compact GO graph, IC artifact, and ProSig motif library."""
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

    motif_result = write_prosig_motif_library(
        prosite_file=prosite_dat,
        output_file=motif_out,
        logger=logger,
    )
    motif_stats = motif_result.stats
    logger.info("Wrote ProSig motif library to %s", motif_result.output_file)
    logger.info(
        "Translated %s of %s PROSITE PATTERN entries into ProSig motifs",
        format_log_number(motif_stats.translated_entries),
        format_log_number(motif_stats.pattern_entries),
    )
    logger.info(
        "Skipped %s PATTERN entries without PA lines; "
        "omitted %s unsupported translations",
        format_log_number(motif_stats.skipped_pattern_entries_without_pa),
        format_log_number(motif_stats.failed_entries),
    )
    logger.info(
        "Converted %s entries to ProSig macros; translated %s ambiguous residue codes",
        format_log_number(motif_stats.macro_converted_entries),
        format_log_number(motif_stats.ambiguous_codes_translated),
    )
