#!/usr/bin/env python3
"""Plot motif enrichment across function clusters as an interactive treemap."""

from __future__ import annotations

import argparse
import csv
import math
import pickle
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from plot_function_label_treemap import (  # noqa: E402
    _build_cluster_groups,
    _load_cluster_function_descriptions,
    _load_cluster_rows,
    _load_go_terms,
)

LIGHT_GREY = "#d3d3d3"
DARK_BLUE = "#08306b"
COLOR_MAX = 10.0


@dataclass(frozen=True)
class EnrichmentGroup:
    cluster_id: str
    label: str
    count: int
    hover_label: str
    description: str
    weight: float


def main() -> None:
    args = _parse_args()
    go_terms = _load_go_terms(args.go_graph)
    cluster_rows = _load_cluster_rows(args.clusters)
    cluster_functions = _load_cluster_function_descriptions(
        args.clusters_meta,
        go_terms=go_terms,
        max_modifiers=args.max_modifiers,
    )
    cluster_groups = _build_cluster_groups(
        cluster_ids=(row["cluster_id"] for row in cluster_rows),
        cluster_functions=cluster_functions,
        top_groups=args.top_cluster_groups,
    )
    motif_weights = _load_motif_weights(args.scoreboard, args.motif)
    enrichment_groups = [
        EnrichmentGroup(
            cluster_id=_raw_cluster_id(group.group_id),
            label=group.label,
            count=group.count,
            hover_label=group.hover_label,
            description=group.details,
            weight=motif_weights.get(_raw_cluster_id(group.group_id), 0.0),
        )
        for group in cluster_groups
    ]

    output = args.output or Path(
        f"motif_enrichment_{_safe_filename_component(args.motif)}.html"
    )
    summary_tsv = args.summary_tsv or output.with_suffix(".tsv")
    _write_summary(summary_tsv, args.motif, enrichment_groups)
    if not args.summary_only:
        _write_treemap(
            output=output,
            motif=args.motif,
            groups=enrichment_groups,
            title=args.title or f"Motif enrichment: {args.motif}",
            width=args.width,
            height=args.height,
        )
        print(f"Wrote motif enrichment treemap to {output}")
    print(f"Wrote motif enrichment table to {summary_tsv}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot one motif's enrichment weights across function clusters. "
            "Tile area represents cluster population and color represents weight."
        )
    )
    parser.add_argument("motif", help="Motif ID in motif_cluster_scoreboard.pkl.")
    parser.add_argument("--clusters", type=Path, default=Path("clusters.tsv"))
    parser.add_argument(
        "--clusters-meta", type=Path, default=Path("clusters_meta.tsv")
    )
    parser.add_argument("--go-graph", type=Path, default=Path("go_graph.pkl"))
    parser.add_argument(
        "--scoreboard",
        type=Path,
        default=Path("motif_cluster_scoreboard.pkl"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output HTML or static image path. Default includes the motif ID.",
    )
    parser.add_argument(
        "--summary-tsv",
        type=Path,
        help="Auditable TSV path. Defaults alongside --output.",
    )
    parser.add_argument(
        "--top-cluster-groups",
        type=int,
        default=0,
        help="Number of largest clusters to plot. Default 0 plots all clusters.",
    )
    parser.add_argument("--max-modifiers", type=int, default=3)
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Write only the TSV summary without importing Plotly.",
    )
    parser.add_argument("--title")
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--height", type=int, default=900)
    args = parser.parse_args()
    if args.top_cluster_groups < 0:
        parser.error("--top-cluster-groups must be at least 0")
    if args.width < 400 or args.height < 300:
        parser.error("--width and --height are too small for a readable treemap")
    if args.max_modifiers < 0:
        parser.error("--max-modifiers must be at least 0")
    return args


def _load_motif_weights(path: Path, motif: str) -> dict[str, float]:
    if not path.is_file():
        raise SystemExit(f"Motif-cluster scoreboard not found: {path}")
    with path.open("rb") as handle:
        artifact = pickle.load(handle)
    if not isinstance(artifact, dict) or not isinstance(artifact.get("weights"), dict):
        raise ValueError(f"Expected motif-cluster scoreboard with weights: {path}")
    motif_records = artifact["weights"].get(motif)
    if not isinstance(motif_records, dict):
        raise SystemExit(f"Motif not found in {path}: {motif}")

    weights: dict[str, float] = {}
    for cluster_id, record in motif_records.items():
        if not isinstance(record, dict) or "weight" not in record:
            raise ValueError(f"Invalid scoreboard record for {motif}/{cluster_id}")
        weight = float(record["weight"])
        if not math.isfinite(weight) or weight < 0.0:
            raise ValueError(
                f"Invalid enrichment weight for {motif}/{cluster_id}: {weight}"
            )
        weights[str(cluster_id)] = weight
    return weights


def _safe_filename_component(value: str) -> str:
    component = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return component or "motif"


def _raw_cluster_id(group_id: str) -> str:
    prefix = "cluster:"
    if not group_id.startswith(prefix):
        raise ValueError(f"Unexpected cluster treemap group ID: {group_id}")
    return group_id.removeprefix(prefix)


def _write_summary(
    path: Path,
    motif: str,
    groups: list[EnrichmentGroup],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "rank",
                "motif_id",
                "cluster_id",
                "cluster_label",
                "member_count",
                "weight",
            ]
        )
        for rank, group in enumerate(groups, start=1):
            writer.writerow(
                [
                    rank,
                    motif,
                    group.cluster_id,
                    group.description,
                    group.count,
                    group.weight,
                ]
            )


def _write_treemap(
    *,
    output: Path,
    motif: str,
    groups: list[EnrichmentGroup],
    title: str,
    width: int,
    height: int,
) -> None:
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise SystemExit(
            "Plotly is required to draw the treemap. Install it with: "
            "python -m pip install plotly"
        ) from exc
    if not groups:
        raise SystemExit("No clusters available for the motif enrichment treemap.")

    ids = ["all", *(f"cluster:{group.cluster_id}" for group in groups)]
    labels = ["Clusters", *(group.label for group in groups)]
    parents = ["", *("all" for _group in groups)]
    values = [sum(group.count for group in groups), *(group.count for group in groups)]
    colors = [0.0, *(min(group.weight, COLOR_MAX) for group in groups)]
    customdata = [
        [motif, "", 0.0],
        *(
            [motif, group.hover_label, group.weight]
            for group in groups
        ),
    ]
    hovertemplates = [
        "<b>Clusters</b><br>Members shown: %{value:,}<extra></extra>",
        *(
            "<b>%{customdata[1]}</b><br>Members: %{value:,}<br>"
            "Motif: %{customdata[0]}<br>Enrichment weight: "
            "%{customdata[2]:.4f}<extra></extra>"
            for _group in groups
        ),
    ]
    fig = go.Figure(
        go.Treemap(
            ids=ids,
            labels=labels,
            parents=parents,
            values=values,
            customdata=customdata,
            hovertemplate=hovertemplates,
            branchvalues="total",
            marker={
                "colors": colors,
                "colorscale": [[0.0, LIGHT_GREY], [1.0, DARK_BLUE]],
                "cmin": 0.0,
                "cmax": COLOR_MAX,
                "colorbar": {"title": "Weight", "tickvals": [0, 2, 4, 6, 8, 10]},
                "line": {"width": 1, "color": "white"},
            },
            tiling={"packing": "squarify"},
        )
    )
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        margin={"t": 60, "l": 10, "r": 10, "b": 10},
        font={"size": 14},
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() == ".html":
        fig.write_html(output, include_plotlyjs="cdn")
        return
    try:
        fig.write_image(output)
    except ValueError as exc:
        raise SystemExit(
            "Static image export requires kaleido. Install it with: "
            "python -m pip install kaleido"
        ) from exc


if __name__ == "__main__":
    main()
