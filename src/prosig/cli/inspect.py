from __future__ import annotations

import csv
import json
import pickle
import re
import textwrap
from collections import deque
from dataclasses import asdict
from pathlib import Path
from typing import Annotated, Any

import typer

from prosig.go.build import MF_ROOT
from prosig.go.describe import describe_go_function
from prosig.go.similarity import (
    GoBestMatch,
    GoSetSimilarityResult,
    GoSimilarity,
    is_go_term_set_input,
    load_accession_mf_go_terms,
    resolve_go_set_query,
)
from prosig.library import resolve_core_library

CLUSTER_ID_PATTERN = re.compile(r"^cluster_\d+$")

inspect_app = typer.Typer(
    help="Inspect ProSig data artifacts and diagnostic calculations.",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)


def _load_go_similarity(go_graph: Path) -> GoSimilarity:
    try:
        return GoSimilarity.from_pickle(go_graph)
    except FileNotFoundError as exc:
        raise typer.BadParameter(f"GO graph file not found: {go_graph}") from exc
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


def _resolve_runtime_library(library_dir: Path | None):
    try:
        return resolve_core_library(library_dir)
    except FileNotFoundError as exc:
        raise typer.BadParameter(str(exc), param_hint="--library-dir") from exc


def _resolve_go_graph_path(
    go_graph: Path | None,
    library_dir: Path | None,
) -> Path:
    if go_graph is not None:
        return go_graph
    return _resolve_runtime_library(library_dir).path("go_graph.pkl")


@inspect_app.command(name="go-summary")
def go_summary(
    go_graph: Annotated[
        Path | None,
        typer.Option(
            "--go-graph",
            help=(
                "Optional path to a compact GO graph pickle. If omitted, "
                "the resolved runtime library is used."
            ),
        ),
    ] = None,
    library_dir: Annotated[
        Path | None,
        typer.Option(
            "--library-dir",
            help=(
                "Directory containing the complete ProSig runtime library. "
                "If omitted, inspect uses all core files from the current "
                "directory when any are present, otherwise packaged defaults."
            ),
        ),
    ] = None,
) -> None:
    """Summarize the GO graph and IC artifact."""
    resolved_go_graph = _resolve_go_graph_path(go_graph, library_dir)
    similarity = _load_go_similarity(resolved_go_graph)
    meta = similarity.meta
    typer.echo(f"path\t{resolved_go_graph}")
    typer.echo(f"namespace\t{meta.get('namespace', 'unknown')}")
    typer.echo(f"schema_version\t{meta.get('schema_version', 'unknown')}")
    typer.echo(f"terms\t{len(similarity.terms)}")
    ic_terms = sum(
        1 for term in similarity.terms.values() if term.get("ic") is not None
    )
    typer.echo(f"ic_terms\t{ic_terms}")
    typer.echo(f"created_at\t{meta.get('created_at', 'unknown')}")


@inspect_app.command(name="go-term")
def go_term(
    go_id: Annotated[str, typer.Argument(help="Molecular Function GO term ID.")],
    go_graph: Annotated[
        Path | None,
        typer.Option(
            "--go-graph",
            help=(
                "Optional path to a compact GO graph pickle. If omitted, "
                "the resolved runtime library is used."
            ),
        ),
    ] = None,
    library_dir: Annotated[
        Path | None,
        typer.Option(
            "--library-dir",
            help=(
                "Directory containing the complete ProSig runtime library. "
                "If omitted, inspect uses all core files from the current "
                "directory when any are present, otherwise packaged defaults."
            ),
        ),
    ] = None,
    show_ancestors: Annotated[
        bool,
        typer.Option("--ancestors", help="Print ancestors including the term itself."),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Write diagnostic output as JSON."),
    ] = False,
) -> None:
    """Inspect one Molecular Function GO term."""
    similarity = _load_go_similarity(_resolve_go_graph_path(go_graph, library_dir))
    term = similarity.term(go_id)
    if term is None:
        raise typer.BadParameter(f"GO term not found in MF graph: {go_id}")

    payload = asdict(term)
    if show_ancestors:
        payload["ancestors_including_self"] = sorted(
            similarity.ancestors_including_self(go_id)
        )

    if json_output:
        typer.echo(json.dumps(payload, indent=2, sort_keys=True))
        return

    for key, value in payload.items():
        if isinstance(value, tuple | list):
            value = ",".join(value)
        typer.echo(f"{key}\t{value}")


@inspect_app.command(name="go-sim")
def go_sim(
    go1: Annotated[str, typer.Argument(help="First Molecular Function GO term ID.")],
    go2: Annotated[str, typer.Argument(help="Second Molecular Function GO term ID.")],
    go_graph: Annotated[
        Path | None,
        typer.Option(
            "--go-graph",
            help=(
                "Optional path to a compact GO graph pickle. If omitted, "
                "the resolved runtime library is used."
            ),
        ),
    ] = None,
    library_dir: Annotated[
        Path | None,
        typer.Option(
            "--library-dir",
            help=(
                "Directory containing the complete ProSig runtime library. "
                "If omitted, inspect uses all core files from the current "
                "directory when any are present, otherwise packaged defaults."
            ),
        ),
    ] = None,
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose",
            "-v",
            help="Explain the terms, common ancestors, MICA, and Lin formula.",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Write diagnostic output as JSON."),
    ] = False,
    tree_style: Annotated[
        str,
        typer.Option(
            "--tree-style",
            case_sensitive=False,
            help="Tree connector style for verbose output: unicode or ascii.",
        ),
    ] = "unicode",
) -> None:
    """Inspect Lin similarity between two Molecular Function GO terms."""
    tree_style = tree_style.lower()
    if tree_style not in {"unicode", "ascii"}:
        raise typer.BadParameter("choose one of: unicode, ascii")

    similarity = _load_go_similarity(_resolve_go_graph_path(go_graph, library_dir))
    if not verbose and not json_output:
        typer.echo(_format_score(similarity.lin(go1, go2)))
        return

    result = similarity.lin_with_details(go1, go2)
    if json_output:
        payload = asdict(result)
        typer.echo(json.dumps(payload, indent=2, sort_keys=True))
        return

    if verbose:
        typer.echo(_format_go_sim_verbose(similarity, result, tree_style=tree_style))
        return


@inspect_app.command(name="go-set-sim")
def go_set_sim(
    set1: Annotated[
        str,
        typer.Argument(
            help=(
                "First GO set or accession. Quote GO set expressions in the "
                "shell, e.g. '(GO:0003677;GO:0004386)'."
            )
        ),
    ],
    set2: Annotated[
        str,
        typer.Argument(
            help=(
                "Second GO set or accession. Quote GO set expressions in the "
                "shell, e.g. 'GO:0005524;GO:0046872'."
            )
        ),
    ],
    library_dir: Annotated[
        Path | None,
        typer.Option(
            "--library-dir",
            help=(
                "Directory containing the complete ProSig runtime library. "
                "If omitted, inspect uses all core files from the current "
                "directory when any are present, otherwise packaged defaults."
            ),
        ),
    ] = None,
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose",
            "-v",
            help="Explain AMB directional best matches and final score.",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Write diagnostic output as JSON."),
    ] = False,
) -> None:
    """Inspect AMB Lin similarity between two MF GO term sets."""
    library = _resolve_runtime_library(library_dir)
    similarity = _load_go_similarity(library.path("go_graph.pkl"))
    terms1, terms2 = _resolve_go_set_queries(
        set1,
        set2,
        library.path("accession_mf_go.tsv"),
    )

    if not verbose and not json_output:
        typer.echo(_format_score(similarity.set_lin_amb(terms1, terms2)))
        return

    result = similarity.set_lin_amb_with_details(
        terms1,
        terms2,
        query1=set1,
        query2=set2,
    )
    if json_output:
        typer.echo(json.dumps(asdict(result), indent=2, sort_keys=True))
        return

    typer.echo(_format_go_set_sim_verbose(similarity, result))


@inspect_app.command(name="function")
def function(
    query: Annotated[
        str,
        typer.Argument(
            help=(
                "Accession, cluster ID, or MF GO set. Quote GO set expressions "
                "in the shell, e.g. 'GO:0004672;GO:0005524'."
            )
        ),
    ],
    library_dir: Annotated[
        Path | None,
        typer.Option(
            "--library-dir",
            help=(
                "Directory containing the complete ProSig runtime library. "
                "If omitted, inspect uses all core files from the current "
                "directory when any are present, otherwise packaged defaults."
            ),
        ),
    ] = None,
    max_modifiers: Annotated[
        int,
        typer.Option(
            "--max-modifiers",
            min=0,
            help="Maximum number of binding modifiers in the composed summary.",
        ),
    ] = 3,
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose",
            "-v",
            help="Show resolved GO terms and composition diagnostics.",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Write function description as JSON."),
    ] = False,
) -> None:
    """Describe function from an accession, cluster ID, or MF GO term set."""
    library = _resolve_runtime_library(library_dir)
    similarity = _load_go_similarity(library.path("go_graph.pkl"))
    terms = _resolve_function_query(
        query,
        accession_go=library.path("accession_mf_go.tsv"),
        cluster_meta=library.path("clusters_meta.tsv"),
    )
    result = describe_go_function(
        query,
        terms,
        similarity.terms,
        max_modifiers=max_modifiers,
    )

    if json_output:
        typer.echo(json.dumps(result.asdict(), indent=2, sort_keys=True))
        return

    if not verbose:
        typer.echo(result.summary)
        return

    typer.echo(f"Query: {query}")
    if not is_go_term_set_input(query):
        typer.echo(f"Resolved terms: {';'.join(terms)}")
    typer.echo("")
    typer.echo("GO terms:")
    for term in result.terms:
        suffix_parts = []
        if term.role:
            suffix_parts.append(f"role={term.role}")
        if term.dropped:
            suffix_parts.append("dropped=ancestor")
        if term.missing:
            suffix_parts.append("missing")
        suffix = f" ({', '.join(suffix_parts)})" if suffix_parts else ""
        typer.echo(f"- {term.go_id} {term.name}{suffix}")
    typer.echo("")
    typer.echo(f"Function: {result.summary}")


@inspect_app.command(name="cluster")
def cluster(
    cluster_id: Annotated[
        str,
        typer.Argument(help="Function cluster ID, e.g. cluster_0008."),
    ],
    library_dir: Annotated[
        Path | None,
        typer.Option(
            "--library-dir",
            help=(
                "Directory containing the complete ProSig runtime library. "
                "If omitted, inspect uses all core files from the current "
                "directory when any are present, otherwise packaged defaults."
            ),
        ),
    ] = None,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Write cluster inspection report as JSON."),
    ] = False,
) -> None:
    """Inspect one function cluster and its positive motif-cluster weights."""
    if not _is_cluster_id_input(cluster_id):
        raise typer.BadParameter("cluster ID must look like cluster_0008")
    library = _resolve_runtime_library(library_dir)
    try:
        cluster_record = _load_cluster_record(
            library.path("clusters_meta.tsv"),
            cluster_id,
        )
        motif_descriptions = _load_motif_descriptions(
            library.path("prosig_motifs.tsv"),
        )
        motif_records = _load_cluster_motif_records(
            library.path("motif_cluster_scoreboard.pkl"),
            cluster_id,
        )
    except FileNotFoundError as exc:
        raise typer.BadParameter(f"required file not found: {exc.filename}") from exc
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    composed_description = cluster_record.get("composed_description", "").strip()
    composed_go = _parse_cluster_composed_go(cluster_record.get("composed_go", ""))
    if not composed_description:
        composed_description = _describe_cluster_go_terms(
            cluster_id,
            composed_go,
            library.path("go_graph.pkl"),
        )

    motifs = [
        _format_cluster_motif_record(record, motif_descriptions)
        for record in motif_records
    ]
    payload = {
        "cluster_id": cluster_id,
        "cluster_size": _parse_cluster_size(cluster_record),
        "composed_go": list(composed_go),
        "composed_description": composed_description,
        "motifs": motifs,
    }

    if json_output:
        typer.echo(json.dumps(payload, indent=2, sort_keys=True))
        return

    typer.echo(_format_cluster_report(payload))


def _resolve_go_set_queries(
    set1: str,
    set2: str,
    accession_go: Path,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    accession_terms: dict[str, tuple[str, ...]] = {}
    if not is_go_term_set_input(set1) or not is_go_term_set_input(set2):
        try:
            accession_terms = load_accession_mf_go_terms(accession_go)
        except FileNotFoundError as exc:
            raise typer.BadParameter(
                f"accession GO file not found: {accession_go}"
            ) from exc
        except ValueError as exc:
            raise typer.BadParameter(str(exc)) from exc
    try:
        return (
            resolve_go_set_query(set1, accession_terms),
            resolve_go_set_query(set2, accession_terms),
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


def _resolve_go_set_query(query: str, accession_go: Path) -> tuple[str, ...]:
    accession_terms: dict[str, tuple[str, ...]] = {}
    if not is_go_term_set_input(query):
        try:
            accession_terms = load_accession_mf_go_terms(accession_go)
        except FileNotFoundError as exc:
            raise typer.BadParameter(
                f"accession GO file not found: {accession_go}"
            ) from exc
        except ValueError as exc:
            raise typer.BadParameter(str(exc)) from exc
    try:
        return resolve_go_set_query(query, accession_terms)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


def _load_cluster_record(cluster_meta: Path, cluster_id: str) -> dict[str, str]:
    with cluster_meta.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        required = {"cluster_id", "size", "composed_go"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise ValueError(
                f"{cluster_meta} missing required column(s): {missing}"
            )
        for row in reader:
            if str(row.get("cluster_id", "")).strip() != cluster_id:
                continue
            return {key: value or "" for key, value in row.items()}
    raise ValueError(f"cluster ID not found in {cluster_meta}: {cluster_id}")


def _parse_cluster_size(cluster_record: dict[str, str]) -> int:
    raw_size = cluster_record.get("size", "").strip()
    try:
        return int(raw_size)
    except ValueError as exc:
        raise ValueError(f"Invalid cluster size: {raw_size}") from exc


def _describe_cluster_go_terms(
    cluster_id: str,
    composed_go: tuple[str, ...],
    go_graph: Path,
) -> str:
    if not composed_go:
        return "NA"
    try:
        similarity = _load_go_similarity(go_graph)
    except typer.BadParameter:
        return "NA"
    result = describe_go_function(cluster_id, composed_go, similarity.terms)
    return result.summary


def _load_motif_descriptions(motif_library: Path) -> dict[str, str]:
    descriptions: dict[str, str] = {}
    with motif_library.open("r", encoding="utf-8", newline="") as handle:
        data_lines = [
            line for line in handle
            if not line.startswith("#")
        ]
    if not data_lines:
        raise ValueError(f"{motif_library} does not contain a motif TSV header")
    reader = csv.DictReader(data_lines, delimiter="\t")
    fieldnames = set(reader.fieldnames or ())
    required = {"name", "description"}
    if not required.issubset(fieldnames):
        missing = ", ".join(sorted(required - fieldnames))
        raise ValueError(f"{motif_library} missing required column(s): {missing}")
    for row in reader:
        motif_id = str(row.get("name", "")).strip()
        if motif_id:
            descriptions[motif_id] = str(row.get("description", "")).strip()
    return descriptions


def _load_cluster_motif_records(
    motif_scoreboard: Path,
    cluster_id: str,
) -> list[dict[str, Any]]:
    with motif_scoreboard.open("rb") as handle:
        artifact = pickle.load(handle)
    if not isinstance(artifact, dict):
        raise ValueError(
            f"Motif scoreboard must contain a dictionary: {motif_scoreboard}"
        )
    weights = artifact.get("weights")
    if not isinstance(weights, dict):
        raise ValueError(f"Motif scoreboard missing weights: {motif_scoreboard}")

    records: list[dict[str, Any]] = []
    for motif_id, cluster_weights in weights.items():
        if not isinstance(cluster_weights, dict):
            continue
        raw_record = cluster_weights.get(cluster_id)
        if not isinstance(raw_record, dict):
            continue
        record = dict(raw_record)
        record.setdefault("motif_id", motif_id)
        record.setdefault("cluster_id", cluster_id)
        records.append(record)
    records.sort(
        key=lambda record: (
            -_as_float(record.get("weight")),
            str(record.get("motif_id", "")),
        )
    )
    return records


def _format_cluster_motif_record(
    record: dict[str, Any],
    motif_descriptions: dict[str, str],
) -> dict[str, Any]:
    motif_id = str(record.get("motif_id", "")).strip()
    tp = _as_int(record.get("TP"))
    fp = _as_int(record.get("FP"))
    fn = _as_int(record.get("FN"))
    tn = _as_int(record.get("TN"))
    weight = _as_float(record.get("weight"))
    return {
        "motif_id": motif_id,
        "description": motif_descriptions.get(motif_id, ""),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "odds_ratio": _format_odds_ratio(tp, fp, fn, tn),
        "weight": _format_float(weight),
    }


def _format_cluster_report(payload: dict[str, Any]) -> str:
    cluster_id = str(payload["cluster_id"])
    cluster_size = int(payload["cluster_size"])
    lines = [
        *_format_wrapped_label("Cluster ID:", cluster_id),
        *_format_wrapped_label("Cluster Size:", str(cluster_size)),
        *_format_wrapped_label("Composed GO:", ";".join(payload["composed_go"])),
        *_format_wrapped_label("Description:", payload["composed_description"]),
        "",
        "Motif Hits:",
    ]
    motifs = payload["motifs"]
    if not motifs:
        lines.append("None")
        return "\n".join(lines)

    for index, motif in enumerate(motifs, start=1):
        lines.extend(
            _format_cluster_motif_section(
                index,
                motif,
                cluster_id=cluster_id,
                cluster_size=cluster_size,
            )
        )
    return "\n".join(lines)


def _format_wrapped_label(
    label: str,
    value: Any,
    *,
    label_width: int = 17,
    width: int = 79,
) -> list[str]:
    text = str(value) if value else "NA"
    wrapped = textwrap.wrap(
        text,
        width=max(20, width - label_width),
        break_long_words=False,
        break_on_hyphens=False,
    ) or ["NA"]
    lines = [f"{label:<{label_width}}{wrapped[0]}"]
    lines.extend(f"{'':<{label_width}}{line}" for line in wrapped[1:])
    return lines


def _format_cluster_motif_section(
    index: int,
    motif: dict[str, Any],
    *,
    cluster_id: str,
    cluster_size: int,
) -> list[str]:
    inside_label = f"In {cluster_id}"
    outside_label = f"Outside {cluster_id}"
    row_label_width = 15
    count_width = max(len(inside_label), len(outside_label), 12) + 2
    separator = "-" * (row_label_width + count_width * 2)
    lines = [
        "",
        f"{index}. {motif['motif_id']}",
    ]
    if motif["description"]:
        lines.append(str(motif["description"]))
    lines.extend(
        [
            (
                f"{'':<{row_label_width}}"
                f"{inside_label:>{count_width}}"
                f"{outside_label:>{count_width}}"
            ),
            separator,
            (
                f"{'Motif present':<{row_label_width}}"
                f"{motif['TP']:>{count_width}}"
                f"{motif['FP']:>{count_width}}"
            ),
            (
                f"{'Motif absent':<{row_label_width}}"
                f"{motif['FN']:>{count_width}}"
                f"{motif['TN']:>{count_width}}"
            ),
            separator,
            f"Odds ratio:  {motif['odds_ratio']}",
            f"Weight:      {motif['weight']}",
        ]
    )
    if motif["TP"] + motif["FN"] != cluster_size:
        lines.append(
            "Note: TP + FN does not equal the reported cluster size; "
            "the scoreboard may be stale relative to clusters_meta.tsv."
        )
    return lines


def _as_int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid integer in motif scoreboard: {value}") from exc


def _as_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid float in motif scoreboard: {value}") from exc


def _format_odds_ratio(tp: int, fp: int, fn: int, tn: int) -> str:
    numerator = tp * tn
    denominator = fp * fn
    if denominator == 0:
        return "inf" if numerator > 0 else "NA"
    return _format_float(numerator / denominator)


def _format_float(value: float) -> str:
    if value == float("inf"):
        return "inf"
    if value == float("-inf"):
        return "-inf"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _resolve_function_query(
    query: str,
    *,
    accession_go: Path,
    cluster_meta: Path,
) -> tuple[str, ...]:
    if _is_cluster_id_input(query):
        try:
            return _load_cluster_composed_go(cluster_meta, query)
        except FileNotFoundError as exc:
            raise typer.BadParameter(
                f"cluster metadata file not found: {cluster_meta}"
            ) from exc
        except ValueError as exc:
            raise typer.BadParameter(str(exc)) from exc
    return _resolve_go_set_query(query, accession_go)


def _is_cluster_id_input(query: str) -> bool:
    return CLUSTER_ID_PATTERN.fullmatch(query.strip()) is not None


def _load_cluster_composed_go(cluster_meta: Path, cluster_id: str) -> tuple[str, ...]:
    with cluster_meta.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        required = {"cluster_id", "composed_go"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise ValueError(
                f"{cluster_meta} missing required column(s): {missing}"
            )
        for row in reader:
            if str(row.get("cluster_id", "")).strip() != cluster_id:
                continue
            return _parse_cluster_composed_go(row.get("composed_go", ""))
    raise ValueError(f"cluster ID not found in {cluster_meta}: {cluster_id}")


def _parse_cluster_composed_go(value: str | None) -> tuple[str, ...]:
    if value is None:
        return ()
    return tuple(
        dict.fromkeys(
            token.strip()
            for token in value.replace(",", ";").split(";")
            if token.strip()
        )
    )


def _format_score(score: float | None) -> str:
    if score is None:
        return "NA"
    return f"{score:.4f}".rstrip("0").rstrip(".")


def _format_go_set_sim_verbose(
    similarity: GoSimilarity,
    result: GoSetSimilarityResult,
) -> str:
    lines = [
        f"Score: {_format_score(result.similarity)}",
        f"Status: {result.status}",
    ]
    if result.reason:
        lines.append(f"Reason: {result.reason}")
    lines.extend(
        [
            "",
            _format_go_set_query("A", result.query1, result.terms1),
            _format_go_set_query("B", result.query2, result.terms2),
        ]
    )
    if result.missing_terms1:
        lines.append(f"Ignored missing A terms: {';'.join(result.missing_terms1)}")
    if result.missing_terms2:
        lines.append(f"Ignored missing B terms: {';'.join(result.missing_terms2)}")

    mean_1_to_2 = _mean_best_match_score(result.best_matches_1_to_2)
    mean_2_to_1 = _mean_best_match_score(result.best_matches_2_to_1)
    lines.extend(
        [
            "",
            "GO term descriptions:",
            *_format_go_set_term_descriptions(similarity, result),
            "",
            "A -> B best matches:",
            *_format_best_matches(
                result.best_matches_1_to_2,
            ),
            f"A -> B average max: {_format_score(mean_1_to_2)}",
            "",
            "B -> A best matches:",
            *_format_best_matches(
                result.best_matches_2_to_1,
            ),
            f"B -> A average max: {_format_score(mean_2_to_1)}",
            "",
            "Formula: AMB(A, B) = (mean(A -> B) + mean(B -> A)) / 2",
            (
                "                   = "
                f"({_format_score(mean_1_to_2)} + {_format_score(mean_2_to_1)}) / 2"
            ),
            f"                   = {_format_score(result.similarity)}",
        ]
    )
    return "\n".join(lines)


def _format_go_set_query(label: str, query: str, terms: tuple[str, ...]) -> str:
    if is_go_term_set_input(query):
        return f"{label} query: {query}"
    return f"{label} query: {query} ({';'.join(terms)})"


def _format_go_set_term_descriptions(
    similarity: GoSimilarity,
    result: GoSetSimilarityResult,
) -> list[str]:
    lines: list[str] = []
    seen: set[str] = set()
    for label, terms in (("A", result.valid_terms1), ("B", result.valid_terms2)):
        for go_id in terms:
            if go_id in seen:
                continue
            seen.add(go_id)
            description = _format_go_id_and_name(similarity, go_id)
            lines.append(f"- {label}: {description}")
    if not lines:
        return ["- none"]
    return lines


def _format_best_matches(
    matches: tuple[GoBestMatch, ...],
) -> list[str]:
    if not matches:
        return ["- none"]
    return [
        (
            f"- {match.source} --{_format_fixed_score(match.score)}--> "
            f"{match.target}"
        )
        for match in matches
    ]


def _format_fixed_score(score: float | None) -> str:
    if score is None:
        return "NA"
    return f"{score:.4f}"


def _format_go_id_and_name(similarity: GoSimilarity, go_id: str) -> str:
    term = similarity.term(go_id)
    if term is None:
        return go_id
    return f"{term.go_id} {term.name}"


def _mean_best_match_score(matches: tuple[GoBestMatch, ...]) -> float | None:
    if not matches:
        return None
    return sum(match.score for match in matches) / len(matches)


def _format_go_sim_verbose(
    similarity: GoSimilarity,
    result: Any,
    *,
    tree_style: str = "unicode",
) -> str:
    term1 = similarity.term(result.go1)
    term2 = similarity.term(result.go2)
    lines = [
        f"Score: {_format_score(result.similarity)}",
        f"Status: {result.status}",
    ]
    if result.reason:
        lines.append(f"Reason: {result.reason}")
    lines.extend(
        [
            "",
            "Input terms:",
            f"- A: {_format_term_description(term1, result.go1)}",
            f"- B: {_format_term_description(term2, result.go2)}",
            "",
            "Common ancestors:",
        ]
    )

    if result.common_ancestors:
        for go_id in sorted(
            result.common_ancestors,
            key=lambda ancestor: _common_ancestor_sort_key(similarity, ancestor),
        ):
            description = _format_term_description(similarity.term(go_id), go_id)
            lines.append(f"- {description}")
    else:
        lines.append("- none")

    formula_prefix = "Formula: Lin(A, B)"
    lines.extend(
        [
            "",
            "Compact GO path:",
            *_format_go_path_graph(similarity, result, tree_style=tree_style),
            "",
            f"MICA: {_format_mica_description(similarity, result.mica)}",
            (
                f"{formula_prefix} = "
                "2 * IC(MICA) / (IC(A) + IC(B))"
            ),
            (
                f"{'':<{len(formula_prefix)}} = "
                f"2 * {_format_number(result.ic_mica)} / "
                f"({_format_number(result.ic_go1)} + {_format_number(result.ic_go2)}) "
            ),
            f"{'':<{len(formula_prefix)}} = {_format_score(result.similarity)}",
        ]
    )
    return "\n".join(lines)


def _format_go_path_graph(
    similarity: GoSimilarity,
    result: Any,
    *,
    tree_style: str,
) -> list[str]:
    path_a = _canonical_root_path(similarity, result.go1, through=result.mica)
    path_b = _canonical_root_path(similarity, result.go2, through=result.mica)
    if not path_a and not path_b:
        return ["- unavailable"]

    tree: dict[str, Any] = {}
    path_sets: dict[str, set[str]] = {}
    for label, path in (("A", path_a), ("B", path_b)):
        if not path:
            continue
        _insert_path(tree, path)
        for go_id in path:
            path_sets.setdefault(go_id, set()).add(label)

    lines: list[str] = []
    _render_go_tree(
        similarity,
        tree,
        lines,
        path_sets=path_sets,
        mica=result.mica,
        go1=result.go1,
        go2=result.go2,
        tree_style=tree_style,
    )
    return lines


def _canonical_root_path(
    similarity: GoSimilarity,
    go_id: str,
    *,
    through: str | None,
) -> tuple[str, ...]:
    if similarity.term(go_id) is None:
        return ()
    if through and similarity.term(through) is not None:
        root_to_mica = _path_from_ancestor_to_term(similarity, MF_ROOT, through)
        mica_to_term = _path_from_ancestor_to_term(similarity, through, go_id)
        if root_to_mica and mica_to_term:
            return (*root_to_mica, *mica_to_term[1:])
    return _path_from_ancestor_to_term(similarity, MF_ROOT, go_id)


def _path_from_ancestor_to_term(
    similarity: GoSimilarity,
    ancestor: str,
    go_id: str,
) -> tuple[str, ...]:
    if ancestor == go_id:
        return (go_id,)
    if ancestor not in similarity.ancestors_including_self(go_id):
        return ()

    queue: deque[tuple[str, tuple[str, ...]]] = deque([(go_id, (go_id,))])
    seen = {go_id}
    while queue:
        current, path = queue.popleft()
        term = similarity.term(current)
        if term is None:
            continue
        for parent in sorted(term.parents):
            if parent in seen:
                continue
            parent_path = (*path, parent)
            if parent == ancestor:
                return tuple(reversed(parent_path))
            seen.add(parent)
            queue.append((parent, parent_path))
    return ()


def _insert_path(tree: dict[str, Any], path: tuple[str, ...]) -> None:
    cursor = tree
    for go_id in path:
        cursor = cursor.setdefault(go_id, {})


def _render_go_tree(
    similarity: GoSimilarity,
    tree: dict[str, Any],
    lines: list[str],
    *,
    path_sets: dict[str, set[str]],
    mica: str | None,
    go1: str,
    go2: str,
    tree_style: str,
    prefix: str = "",
) -> None:
    items = sorted(
        tree.items(),
        key=lambda item: _tree_node_sort_key(similarity, item[0]),
    )
    connector_middle, connector_last, prefix_middle, prefix_last = (
        _tree_style_tokens(tree_style)
    )
    for index, (go_id, children) in enumerate(items):
        is_last = index == len(items) - 1
        connector = connector_last if is_last else connector_middle
        child_prefix = prefix_last if is_last else prefix_middle
        node = _format_graph_node(similarity, go_id, path_sets, mica, go1, go2)
        lines.append(f"{prefix}{connector}{node}")
        _render_go_tree(
            similarity,
            children,
            lines,
            path_sets=path_sets,
            mica=mica,
            go1=go1,
            go2=go2,
            tree_style=tree_style,
            prefix=f"{prefix}{child_prefix}",
        )


def _tree_style_tokens(tree_style: str) -> tuple[str, str, str, str]:
    if tree_style == "ascii":
        return "|- ", "`- ", "|  ", "   "
    return "├── ", "└── ", "│   ", "    "


def _tree_node_sort_key(similarity: GoSimilarity, go_id: str) -> tuple[int, str]:
    term = similarity.term(go_id)
    return (term.depth if term and term.depth is not None else 10**9, go_id)


def _format_graph_node(
    similarity: GoSimilarity,
    go_id: str,
    path_sets: dict[str, set[str]],
    mica: str | None,
    go1: str,
    go2: str,
) -> str:
    term = similarity.term(go_id)
    name = term.name if term is not None else "missing"
    markers: list[str] = []
    if go_id == mica:
        markers.append("MICA")
    elif path_sets.get(go_id) == {"A", "B"}:
        markers.append("common")
    if go_id == go1:
        markers.append("A")
    if go_id == go2:
        markers.append("B")
    suffix = f" [{' '.join(markers)}]" if markers else ""
    return f"{go_id} {name}{suffix}"


def _common_ancestor_sort_key(
    similarity: GoSimilarity,
    go_id: str,
) -> tuple[float, str]:
    term = similarity.term(go_id)
    return (term.ic if term and term.ic is not None else float("inf"), go_id)


def _format_mica_description(similarity: GoSimilarity, mica: str | None) -> str:
    if mica is None:
        return "NA"
    term = similarity.term(mica)
    if term is None:
        return mica
    return f"{term.go_id} {term.name}"


def _format_term_description(term: Any, go_id: str) -> str:
    if term is None:
        return f"{go_id} missing"
    return (
        f"{term.go_id} {term.name} "
        f"(IC={_format_number(term.ic)}, freq={_format_number(term.freq)}, "
        f"depth={term.depth if term.depth is not None else 'NA'})"
    )


def _format_number(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.4f}".rstrip("0").rstrip(".")
