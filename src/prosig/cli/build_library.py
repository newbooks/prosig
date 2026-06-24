import pickle
from importlib.resources import files
from pathlib import Path
from typing import Annotated

import typer

from prosig.cli.logging import get_logger
from prosig.go.build import (
    MF_ROOT,
    build_go_pkl,
    ensure_role_map_from_template,
    format_go_report,
    format_log_number,
    write_accession_mf_go_tsv,
    write_excluded_mf_annotation_diagnostics,
    write_go_graph_json,
)
from prosig.go.clustering import (
    cluster_accessions_by_go,
    parse_cluster_config,
)
from prosig.io.freshness import artifact_is_stale
from prosig.motifs.prosite import write_prosig_motif_library

ROLE_MAP_TEMPLATE = files("prosig.data").joinpath("role_map.yaml.template")
CLUSTER_CONFIG_TEMPLATE = files("prosig.data").joinpath(
    "cluster_config.yaml.template"
)


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
    role_map: Annotated[
        Path,
        typer.Option(
            "--role-map",
            help=(
                "Path to GO semantic role map. Created from the starter template "
                "when missing."
            ),
        ),
    ] = Path("role_map.yaml"),
    cluster_out: Annotated[
        Path,
        typer.Option(
            "--cluster-out",
            help="Path to write GO accession clusters.",
        ),
    ] = Path("go_clusters.tsv"),
    cluster_config: Annotated[
        Path,
        typer.Option(
            "--cluster-config",
            help=(
                "Path to GO clustering config. Created from the starter template "
                "when missing."
            ),
        ),
    ] = Path("cluster_config.yaml"),
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            "-f",
            help="Rebuild derived artifacts even when outputs are newer than inputs.",
        ),
    ] = False,
) -> None:
    """Build the compact GO graph, IC artifact, and ProSig motif library."""
    logger = get_logger()
    if ensure_role_map_from_template(role_map, ROLE_MAP_TEMPLATE):
        logger.info("Created starter GO semantic role map: %s", role_map)
    else:
        logger.info("Using GO semantic role map: %s", role_map)
    if ensure_role_map_from_template(cluster_config, CLUSTER_CONFIG_TEMPLATE):
        logger.info("Created starter GO clustering config: %s", cluster_config)
    else:
        logger.info("Using GO clustering config: %s", cluster_config)
    parsed_cluster_config = parse_cluster_config(cluster_config)

    go_dependencies = [go_obo, swissprot, role_map]
    if artifact_is_stale(go_out, go_dependencies, force=force):
        logger.info(
            "Building GO graph because %s is missing, stale, or rebuild was forced",
            go_out,
        )
        artifact = build_go_pkl(
            go_obo=go_obo,
            swissprot=swissprot,
            go_out=go_out,
            report_out=write_report,
            role_map=role_map,
        )
    else:
        logger.info(
            "Skipping GO graph build: %s is current with %s",
            go_out,
            _format_dependencies(go_dependencies),
        )
        artifact = _load_go_artifact(go_out)
    meta = artifact["meta"]
    logger.info(
        "GO graph and IC artifact available at %s with %s MF terms",
        go_out,
        format_log_number(meta["n_terms"]),
    )
    logger.info(
        "Used %s of %s provided accessions for IC",
        format_log_number(artifact["terms"][MF_ROOT]["count"]),
        format_log_number(meta["n_accessions_provided"]),
    )
    _refresh_go_side_artifacts(
        artifact=artifact,
        go_out=go_out,
        swissprot=swissprot,
        write_report=write_report,
        force=force,
        logger=logger,
    )

    if artifact_is_stale(motif_out, [prosite_dat], force=force):
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
            "Converted %s entries to ProSig macros; translated %s ambiguous "
            "residue codes",
            format_log_number(motif_stats.macro_converted_entries),
            format_log_number(motif_stats.ambiguous_codes_translated),
        )
    else:
        logger.info(
            "Skipping ProSig motif library build: %s is current with %s",
            motif_out,
            prosite_dat,
        )

    accession_mf_go = go_out.parent / "accession_mf_go.tsv"
    cluster_dependencies = [go_out, accession_mf_go, cluster_config]
    cluster_stats_out = Path(parsed_cluster_config.stats_file)
    cluster_meta_out = Path(parsed_cluster_config.meta_file)
    cluster_outputs = [cluster_out, cluster_stats_out, cluster_meta_out]
    if any(
        artifact_is_stale(output, cluster_dependencies, force=force)
        for output in cluster_outputs
    ):
        cluster_result = cluster_accessions_by_go(
            accession_mf_go,
            go_artifact=artifact,
            go_graph_file=go_out,
            output_file=cluster_out,
            stats_file=cluster_stats_out,
            meta_file=cluster_meta_out,
            resolution=parsed_cluster_config.resolution,
            neighbors=parsed_cluster_config.neighbors,
            term_cache_size_mb=parsed_cluster_config.term_cache_size_mb,
            profile_cache_size_mb=parsed_cluster_config.profile_cache_size_mb,
            min_informative_ic=parsed_cluster_config.min_informative_ic,
            max_posting_fraction=parsed_cluster_config.max_posting_fraction,
            max_posting_size=parsed_cluster_config.max_posting_size,
            progress_interval_seconds=parsed_cluster_config.progress_interval_seconds,
        )
        logger.info(
            "Wrote GO clusters to %s with %s clustered accessions in %s clusters",
            cluster_result.output_file,
            format_log_number(cluster_result.clustered_accessions),
            format_log_number(cluster_result.clusters),
        )
    else:
        logger.info(
            "Skipping GO clustering: %s is current with %s",
            _format_dependencies(cluster_outputs),
            _format_dependencies(cluster_dependencies),
        )


def _load_go_artifact(path: Path) -> dict:
    with path.open("rb") as handle:
        artifact = pickle.load(handle)
    if not isinstance(artifact, dict):
        raise ValueError(f"GO graph pickle must contain a dictionary artifact: {path}")
    return artifact


def _refresh_go_side_artifacts(
    *,
    artifact: dict,
    go_out: Path,
    swissprot: Path,
    write_report: Path | None,
    force: bool,
    logger,
) -> None:
    if write_report is not None:
        if artifact_is_stale(write_report, [go_out], force=force):
            logger.info("Writing GO build validation report: %s", write_report)
            write_report.parent.mkdir(parents=True, exist_ok=True)
            write_report.write_text(format_go_report(artifact), encoding="utf-8")
            logger.info("Wrote GO build validation report: %s", write_report)
        else:
            logger.info(
                "Skipping GO build validation report: %s is current with %s",
                write_report,
                go_out,
            )

    go_json_out = go_out.parent / "go_graph.json"
    if artifact_is_stale(go_json_out, [go_out], force=force):
        logger.info("Writing diagnostic GO graph JSON: %s", go_json_out)
        write_go_graph_json(go_json_out, artifact)
        logger.info("Wrote diagnostic GO graph JSON: %s", go_json_out)
    else:
        logger.info(
            "Skipping diagnostic GO graph JSON: %s is current with %s",
            go_json_out,
            go_out,
        )

    excluded_mf_annotations_out = go_out.parent / "excluded_mf_annotations.tsv"
    if artifact_is_stale(excluded_mf_annotations_out, [swissprot], force=force):
        logger.info(
            "Writing excluded MF annotation diagnostics: %s",
            excluded_mf_annotations_out,
        )
        write_excluded_mf_annotation_diagnostics(
            excluded_mf_annotations_out,
            swissprot,
        )
        logger.info(
            "Wrote excluded MF annotation diagnostics: %s",
            excluded_mf_annotations_out,
        )
    else:
        logger.info(
            "Skipping excluded MF annotation diagnostics: %s is current with %s",
            excluded_mf_annotations_out,
            swissprot,
        )

    accession_mf_go_out = go_out.parent / "accession_mf_go.tsv"
    if artifact_is_stale(accession_mf_go_out, [swissprot], force=force):
        logger.info("Writing accession MF GO terms: %s", accession_mf_go_out)
        write_accession_mf_go_tsv(accession_mf_go_out, swissprot)
        logger.info("Wrote accession MF GO terms: %s", accession_mf_go_out)
    else:
        logger.info(
            "Skipping accession MF GO terms: %s is current with %s",
            accession_mf_go_out,
            swissprot,
        )


def _format_dependencies(dependencies: list[Path]) -> str:
    return ", ".join(str(dependency) for dependency in dependencies)
