from __future__ import annotations

import csv
import json
import pickle
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Any

import typer

from prosig.go.describe import describe_go_function
from prosig.go.similarity import GoSimilarity
from prosig.library import resolve_core_library
from prosig.motifs.scanning import motif_present, read_prosig_motif_library


@dataclass(frozen=True)
class _QuerySequence:
    query: str
    sequence: str


@dataclass(frozen=True)
class _ClusterMeta:
    cluster_id: str
    composed_go: tuple[str, ...]
    composed_description: str


def scan(
    seq: Annotated[
        str | None,
        typer.Option("--seq", help="Protein sequence string to scan."),
    ] = None,
    fasta: Annotated[
        Path | None,
        typer.Option("--fasta", help="FASTA file containing query sequence(s)."),
    ] = None,
    library_dir: Annotated[
        Path | None,
        typer.Option(
            "--library-dir",
            help=(
                "Directory containing the complete ProSig runtime library. "
                "If omitted, scan uses all core files from the current "
                "directory when any are present, otherwise packaged defaults."
            ),
        ),
    ] = None,
    min_weight: Annotated[
        float,
        typer.Option(
            "--min-weight",
            help="Minimum motif-cluster weight retained as an inferred GO set.",
        ),
    ] = 2.0,
    top_n: Annotated[
        int,
        typer.Option(
            "--top-n",
            help="Maximum inferred GO sets to report; use 0 to report all.",
        ),
    ] = 5,
    json_out: Annotated[
        Path | None,
        typer.Option("--json-out", help="Optional path to write JSON output."),
    ] = None,
) -> None:
    """Scan query sequence(s), infer motif-supported GO sets, and print results."""
    if top_n < 0:
        raise typer.BadParameter("must be at least 0", param_hint="--top-n")
    queries = _resolve_queries(
        seq=seq,
        fasta=fasta,
    )
    try:
        library = resolve_core_library(library_dir)
    except FileNotFoundError as exc:
        raise typer.BadParameter(str(exc)) from exc
    motifs = read_prosig_motif_library(library.path("prosig_motifs.tsv"))
    motifs_by_name = {motif.name: motif for motif in motifs}
    scoreboard = _load_scoreboard(library.path("motif_cluster_scoreboard.pkl"))
    clusters = _load_cluster_meta(library.path("clusters_meta.tsv"))
    calibration = _load_calibration(library.path("motif_cluster_scoreboard_meta.json"))
    similarity = _load_similarity_if_available(library.path("go_graph.pkl"))

    reports = [
        _scan_one_query(
            query,
            motifs=motifs,
            motifs_by_name=motifs_by_name,
            scoreboard=scoreboard,
            clusters=clusters,
            calibration=calibration,
            similarity=similarity,
            min_weight=min_weight,
            top_n=top_n,
        )
        for query in queries
    ]
    payload = {
        "queries": reports,
        "min_weight": min_weight,
        "top_n": top_n,
        "library": {
            "source": library.source,
            "directory": str(library.directory),
        },
    }

    if json_out is not None:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return

    typer.echo(
        _format_scan_reports(
            reports,
            min_weight=min_weight,
            top_n=top_n,
            library_source=library.source,
            library_dir=library.directory,
        )
    )


def _resolve_queries(
    *,
    seq: str | None,
    fasta: Path | None,
) -> list[_QuerySequence]:
    inputs = [seq is not None, fasta is not None]
    if sum(inputs) != 1:
        raise typer.BadParameter("provide exactly one of --seq, --fasta")
    if seq is not None:
        sequence = _normalize_sequence(seq)
        if not sequence:
            raise typer.BadParameter("--seq cannot be empty")
        return [_QuerySequence("sequence", sequence)]
    if fasta is not None:
        return list(_read_fasta_queries(fasta))
    raise AssertionError("unreachable query input state")


def _read_fasta_queries(path: Path):
    try:
        handle = path.open("r", encoding="utf-8")
    except FileNotFoundError as exc:
        raise typer.BadParameter(f"FASTA file not found: {path}") from exc
    with handle:
        current_id: str | None = None
        parts: list[str] = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield _QuerySequence(
                        current_id,
                        _normalize_sequence("".join(parts)),
                    )
                current_id = line[1:].split(None, 1)[0] or "query"
                parts = []
                continue
            parts.append(line)
        if current_id is not None:
            yield _QuerySequence(current_id, _normalize_sequence("".join(parts)))
    if current_id is None:
        raise typer.BadParameter(f"FASTA file contains no records: {path}")


def _normalize_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper()


def _load_scoreboard(path: Path) -> dict[str, dict[str, dict[str, Any]]]:
    try:
        with path.open("rb") as handle:
            artifact = pickle.load(handle)
    except FileNotFoundError as exc:
        raise typer.BadParameter(f"motif scoreboard not found: {path}") from exc
    if not isinstance(artifact, dict) or not isinstance(artifact.get("weights"), dict):
        raise typer.BadParameter(f"invalid motif scoreboard artifact: {path}")
    return artifact["weights"]


def _load_cluster_meta(path: Path) -> dict[str, _ClusterMeta]:
    clusters: dict[str, _ClusterMeta] = {}
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except FileNotFoundError as exc:
        raise typer.BadParameter(f"cluster metadata file not found: {path}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        required = {"cluster_id", "composed_go"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise typer.BadParameter(f"{path} missing required column(s): {missing}")
        for row in reader:
            cluster_id = str(row.get("cluster_id", "")).strip()
            if not cluster_id:
                continue
            clusters[cluster_id] = _ClusterMeta(
                cluster_id=cluster_id,
                composed_go=_parse_go_terms(row.get("composed_go", "")),
                composed_description=str(
                    row.get("composed_description", "") or ""
                ).strip(),
            )
    return clusters


def _parse_go_terms(value: str | None) -> tuple[str, ...]:
    if value is None:
        return ()
    return tuple(
        dict.fromkeys(
            token.strip()
            for token in value.replace(",", ";").split(";")
            if token.strip()
        )
    )


def _load_calibration(path: Path) -> list[dict[str, Any]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return []
    except json.JSONDecodeError as exc:
        raise typer.BadParameter(f"invalid scoreboard metadata JSON: {path}") from exc
    calibration = payload.get("stats", {}).get("calibration", [])
    return calibration if isinstance(calibration, list) else []


def _load_similarity_if_available(path: Path) -> GoSimilarity | None:
    try:
        return GoSimilarity.from_pickle(path)
    except (FileNotFoundError, ValueError):
        return None


def _scan_one_query(
    query: _QuerySequence,
    *,
    motifs,
    motifs_by_name,
    scoreboard: dict[str, dict[str, dict[str, Any]]],
    clusters: dict[str, _ClusterMeta],
    calibration: list[dict[str, Any]],
    similarity: GoSimilarity | None,
    min_weight: float,
    top_n: int,
) -> dict[str, Any]:
    hit_motifs = [
        motif.name
        for motif in motifs
        if motif_present(query.sequence, motif)
    ]
    cluster_predictions: dict[str, dict[str, Any]] = {}
    for motif_id in hit_motifs:
        for cluster_id, record in scoreboard.get(motif_id, {}).items():
            weight = float(record.get("weight", 0.0))
            if weight < min_weight:
                continue
            current = cluster_predictions.get(cluster_id)
            if current is None or weight > current["weight"]:
                cluster_predictions[cluster_id] = {
                    "cluster_id": cluster_id,
                    "weight": weight,
                    "motif_id": motif_id,
                    "signature": motifs_by_name[motif_id].pattern,
                }

    go_set_predictions = _cluster_predictions_to_go_sets(
        query.query,
        cluster_predictions,
        clusters=clusters,
        calibration=calibration,
        similarity=similarity,
    )
    if top_n > 0:
        go_set_predictions = go_set_predictions[:top_n]
    return {
        "query": query.query,
        "sequence_length": len(query.sequence),
        "matched_motifs": hit_motifs,
        "inferred_go_sets": go_set_predictions,
    }


def _cluster_predictions_to_go_sets(
    query_name: str,
    cluster_predictions: dict[str, dict[str, Any]],
    *,
    clusters: dict[str, _ClusterMeta],
    calibration: list[dict[str, Any]],
    similarity: GoSimilarity | None,
) -> list[dict[str, Any]]:
    predictions_by_go: dict[tuple[str, ...], dict[str, Any]] = {}
    for cluster_id, prediction in cluster_predictions.items():
        cluster = clusters.get(cluster_id)
        go_terms = cluster.composed_go if cluster is not None else ()
        key = go_terms or (cluster_id,)
        description = _cluster_description(query_name, cluster, similarity)
        current = predictions_by_go.get(key)
        if current is None:
            predictions_by_go[key] = {
                "cluster_ids": [cluster_id],
                "go_terms": list(go_terms),
                "composed_description": description,
                "weight": prediction["weight"],
                "motif_id": prediction["motif_id"],
                "signature": prediction["signature"],
                "calibrated_confidence": _calibrated_confidence(
                    prediction["weight"],
                    calibration,
                ),
            }
            continue
        current["cluster_ids"].append(cluster_id)
        if prediction["weight"] > current["weight"]:
            current["weight"] = prediction["weight"]
            current["composed_description"] = description
            current["motif_id"] = prediction["motif_id"]
            current["signature"] = prediction["signature"]
            current["calibrated_confidence"] = _calibrated_confidence(
                prediction["weight"],
                calibration,
            )

    return sorted(
        predictions_by_go.values(),
        key=lambda prediction: (
            -prediction["weight"],
            ";".join(prediction["go_terms"]),
        ),
    )


def _cluster_description(
    query_name: str,
    cluster: _ClusterMeta | None,
    similarity: GoSimilarity | None,
) -> str:
    if cluster is None:
        return "NA"
    if cluster.composed_description:
        return cluster.composed_description
    if similarity is None or not cluster.composed_go:
        return "NA"
    return describe_go_function(
        query_name,
        cluster.composed_go,
        similarity.terms,
    ).summary


def _calibrated_confidence(
    weight: float,
    calibration: list[dict[str, Any]],
) -> dict[str, Any] | None:
    eligible = [
        point
        for point in calibration
        if float(point.get("weight_threshold", 0.0)) <= weight
    ]
    if not eligible:
        return None
    point = max(eligible, key=lambda item: float(item.get("weight_threshold", 0.0)))
    return {
        "weight_threshold": point.get("weight_threshold"),
        "set_accuracy": point.get("set_accuracy"),
        "top1_accuracy": point.get("top1_accuracy"),
        "top3_accuracy": point.get("top3_accuracy"),
        "coverage": point.get("coverage"),
    }


def _format_scan_reports(
    reports: list[dict[str, Any]],
    *,
    min_weight: float,
    top_n: int,
    library_source: str,
    library_dir: Path,
) -> str:
    blocks = [f"Library:       {library_source} ({library_dir})", ""]
    blocks.extend(
        _format_one_scan_report(report, min_weight=min_weight, top_n=top_n)
        for report in reports
    )
    return "\n\n".join(blocks)


def _format_one_scan_report(
    report: dict[str, Any],
    *,
    min_weight: float,
    top_n: int,
) -> str:
    limit_text = "all" if top_n == 0 else f"top {top_n}"
    lines = [
        f"Query:          {report['query']}",
        f"Sequence size:  {report['sequence_length']} aa",
        f"Matched motifs: {len(report['matched_motifs'])}",
        f"Inferred GO sets ({limit_text}, weight >= {min_weight:g}):",
    ]
    predictions = report["inferred_go_sets"]
    if not predictions:
        lines.append("None")
        return "\n".join(lines)
    for index, prediction in enumerate(predictions, start=1):
        lines.extend(_format_prediction(index, prediction))
    return "\n".join(lines)


def _format_prediction(index: int, prediction: dict[str, Any]) -> list[str]:
    confidence = prediction.get("calibrated_confidence")
    confidence_text = "NA"
    if confidence is not None and confidence.get("set_accuracy") is not None:
        confidence_text = (
            f"{float(confidence['set_accuracy']):.4f}".rstrip("0").rstrip(".")
            + f" (set_acc @ >= {float(confidence['weight_threshold']):g})"
        )
    description_lines = _wrapped_value(
        "Description:",
        prediction["composed_description"],
    )
    return [
        "",
        f"{index}. {';'.join(prediction['go_terms']) or 'NA'}",
        f"Signature:     {prediction['signature']}",
        f"Clusters:       {','.join(prediction['cluster_ids'])}",
        *description_lines,
        f"GO terms:       {';'.join(prediction['go_terms']) or 'NA'}",
        f"Weight:         {prediction['weight']:.4f}".rstrip("0").rstrip("."),
        f"Confidence:     {confidence_text}",
    ]


def _wrapped_value(label: str, value: str, *, width: int = 79) -> list[str]:
    label_width = 15
    text = value or "NA"
    wrapped = textwrap.wrap(
        text,
        width=max(20, width - label_width),
        break_long_words=False,
        break_on_hyphens=False,
    ) or ["NA"]
    lines = [f"{label:<{label_width}}{wrapped[0]}"]
    lines.extend(f"{'':<{label_width}}{line}" for line in wrapped[1:])
    return lines
