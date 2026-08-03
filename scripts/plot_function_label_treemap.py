#!/usr/bin/env python3
"""Plot accession label populations as treemaps."""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if REPO_SRC.is_dir():
    sys.path.insert(0, str(REPO_SRC))

from prosig.go.describe import describe_go_function  # noqa: E402


@dataclass(frozen=True)
class PopulationGroup:
    group_id: str
    label: str
    count: int
    hover_label: str
    details: str


def main() -> None:
    args = _parse_args()
    go_terms = _load_go_terms(args.go_graph)
    cluster_rows = _load_cluster_rows(args.clusters)
    accession_terms = _load_accession_terms(args.accession_go)
    cluster_functions = _load_cluster_function_descriptions(
        args.clusters_meta,
        go_terms=go_terms,
        max_modifiers=args.max_modifiers,
    )

    go_groups = _build_go_set_groups(
        member_ids=(row["member_id"] for row in cluster_rows),
        accession_terms=accession_terms,
        go_terms=go_terms,
        top_groups=args.top_go_groups,
        include_unannotated=args.include_unannotated,
    )
    cluster_groups = _build_cluster_groups(
        cluster_ids=(row["cluster_id"] for row in cluster_rows),
        cluster_functions=cluster_functions,
        top_groups=args.top_cluster_groups,
    )

    _write_population_tsv(args.go_summary_tsv, go_groups)
    _write_population_tsv(args.cluster_summary_tsv, cluster_groups)

    if args.summary_only:
        print(f"Wrote GO label population table to {args.go_summary_tsv}")
        print(f"Wrote cluster label population table to {args.cluster_summary_tsv}")
        return

    _write_treemap(
        output=args.go_output,
        groups=go_groups,
        title=args.go_title,
        root_label="GO label sets",
        width=args.width,
        height=args.height,
    )
    _write_treemap(
        output=args.cluster_output,
        groups=cluster_groups,
        title=args.cluster_title,
        root_label="Clusters",
        width=args.width,
        height=args.height,
    )
    print(f"Wrote GO label population treemap to {args.go_output}")
    print(f"Wrote GO label population table to {args.go_summary_tsv}")
    print(f"Wrote cluster label population treemap to {args.cluster_output}")
    print(f"Wrote cluster label population table to {args.cluster_summary_tsv}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot two treemaps from cluster membership and accession MF GO labels: "
            "one by exact GO label set and one by cluster_id."
        )
    )
    parser.add_argument(
        "--clusters",
        type=Path,
        default=Path("clusters.tsv"),
        help=(
            "Cluster membership file in the current directory. Accepts TSV or CSV "
            "with member_id and cluster_id columns."
        ),
    )
    parser.add_argument(
        "--accession-go",
        type=Path,
        default=Path("accession_mf_go.tsv"),
        help="Headerless TSV in the current directory: accession<TAB>GO:...;GO:...",
    )
    parser.add_argument(
        "--go-graph",
        type=Path,
        default=Path("go_graph.pkl"),
        help="GO graph pickle used to resolve GO IDs to names.",
    )
    parser.add_argument(
        "--clusters-meta",
        type=Path,
        default=Path("clusters_meta.tsv"),
        help=(
            "Cluster metadata TSV with cluster_id and composed_go columns, used "
            "to describe cluster functions."
        ),
    )
    parser.add_argument(
        "--go-output",
        type=Path,
        default=Path("go_label_population_treemap.html"),
        help="Output path for the GO label-set population treemap.",
    )
    parser.add_argument(
        "--cluster-output",
        type=Path,
        default=Path("cluster_label_population_treemap.html"),
        help="Output path for the cluster population treemap.",
    )
    parser.add_argument(
        "--go-summary-tsv",
        type=Path,
        default=Path("go_label_population_treemap.tsv"),
        help="Auditable table for the GO label-set plot.",
    )
    parser.add_argument(
        "--cluster-summary-tsv",
        type=Path,
        default=Path("cluster_label_population_treemap.tsv"),
        help="Auditable table for the cluster label plot.",
    )
    parser.add_argument(
        "--top-go-groups",
        type=int,
        default=80,
        help="Number of largest GO label-set groups to plot. Use 0 for all groups.",
    )
    parser.add_argument(
        "--top-cluster-groups",
        type=int,
        default=120,
        help="Number of largest cluster_id groups to plot. Use 0 for all groups.",
    )
    parser.add_argument(
        "--include-unannotated",
        action="store_true",
        help="Include accessions without MF GO annotations as their own GO group.",
    )
    parser.add_argument(
        "--max-modifiers",
        type=int,
        default=3,
        help="Maximum number of binding modifiers in cluster function summaries.",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Write only TSV summaries, without importing Plotly or drawing.",
    )
    parser.add_argument(
        "--go-title",
        default="GO label population",
    )
    parser.add_argument(
        "--cluster-title",
        default="Cluster label population",
    )
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--height", type=int, default=900)
    args = parser.parse_args()

    if args.top_go_groups < 0:
        parser.error("--top-go-groups must be at least 0")
    if args.top_cluster_groups < 0:
        parser.error("--top-cluster-groups must be at least 0")
    if args.width < 400 or args.height < 300:
        parser.error("--width and --height are too small for a readable treemap")
    if args.max_modifiers < 0:
        parser.error("--max-modifiers must be at least 0")
    return args


def _load_go_terms(path: Path) -> dict[str, dict[str, Any]]:
    if not path.is_file():
        raise SystemExit(f"GO graph file not found: {path}")

    with path.open("rb") as handle:
        artifact = pickle.load(handle)
    if not isinstance(artifact, dict) or not isinstance(artifact.get("terms"), dict):
        raise ValueError(f"Expected GO graph artifact with a terms dictionary: {path}")
    return artifact["terms"]


def _load_cluster_rows(path: Path) -> list[dict[str, str]]:
    delimiter = _delimiter_for(path)
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None or not {"member_id", "cluster_id"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                f"{path} must contain member_id and cluster_id columns; "
                f"found {reader.fieldnames}"
            )
        for row in reader:
            member_id = row["member_id"].strip()
            cluster_id = row["cluster_id"].strip()
            if member_id and cluster_id:
                rows.append({"member_id": member_id, "cluster_id": cluster_id})
    return rows


def _delimiter_for(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    return "\t"


def _load_accession_terms(path: Path) -> dict[str, tuple[str, ...]]:
    accession_terms: dict[str, tuple[str, ...]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            accession = row[0].strip()
            terms = tuple(sorted({term for term in row[1].split(";") if term}))
            if accession and terms:
                accession_terms[accession] = terms
    return accession_terms


def _load_cluster_function_descriptions(
    path: Path,
    *,
    go_terms: dict[str, dict[str, Any]],
    max_modifiers: int,
) -> dict[str, str]:
    descriptions: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        required = {"cluster_id", "composed_go"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise ValueError(f"{path} missing required column(s): {missing}")
        for row in reader:
            cluster_id = row["cluster_id"].strip()
            terms = _parse_go_set(row.get("composed_go", ""))
            if not cluster_id:
                continue
            result = describe_go_function(
                cluster_id,
                terms,
                go_terms,
                max_modifiers=max_modifiers,
            )
            descriptions[cluster_id] = result.summary
    return descriptions


def _parse_go_set(value: str | None) -> tuple[str, ...]:
    if value is None:
        return ()
    return tuple(
        dict.fromkeys(
            token.strip()
            for token in value.replace(",", ";").split(";")
            if token.strip()
        )
    )


def _build_go_set_groups(
    *,
    member_ids: Any,
    accession_terms: dict[str, tuple[str, ...]],
    go_terms: dict[str, dict[str, Any]],
    top_groups: int,
    include_unannotated: bool,
) -> list[PopulationGroup]:
    counts: Counter[tuple[str, ...]] = Counter()
    for member_id in member_ids:
        terms = accession_terms.get(member_id, ())
        if not terms and not include_unannotated:
            continue
        counts[terms] += 1

    groups: list[PopulationGroup] = []
    for terms, count in _most_common(counts, top_groups):
        if terms:
            go_ids = ";".join(terms)
            names = tuple(_go_name(go_id, go_terms) for go_id in terms)
            label = _format_go_set_label(terms, names)
            hover_label = "<br>".join(
                f"{go_id}: {name}" for go_id, name in zip(terms, names, strict=True)
            )
            details = " | ".join(
                f"{go_id}: {name}" for go_id, name in zip(terms, names, strict=True)
            )
            group_id = f"go_set:{go_ids}"
        else:
            label = "No MF GO annotation"
            hover_label = label
            details = label
            group_id = "go_set:unannotated"

        groups.append(
            PopulationGroup(
                group_id=group_id,
                label=label,
                count=count,
                hover_label=hover_label,
                details=details,
            )
        )
    return groups


def _build_cluster_groups(
    *,
    cluster_ids: Any,
    cluster_functions: dict[str, str],
    top_groups: int,
) -> list[PopulationGroup]:
    counts = Counter(cluster_ids)
    groups: list[PopulationGroup] = []
    for cluster_id, count in _most_common(counts, top_groups):
        description = cluster_functions.get(
            cluster_id,
            f"{cluster_id} has no resolved function description.",
        )
        label_phrase = _cluster_label_phrase(cluster_id, description)
        groups.append(
            PopulationGroup(
                group_id=f"cluster:{cluster_id}",
                label=f"{cluster_id}<br>{label_phrase}",
                count=count,
                hover_label=f"{cluster_id}<br>{description}",
                details=description,
            )
        )
    return groups


def _cluster_label_phrase(cluster_id: str, description: str) -> str:
    prefix = f"{cluster_id} is annotated as "
    if description.startswith(prefix):
        phrase = description.removeprefix(prefix).rstrip(".")
        if len(phrase) > 90:
            return f"{phrase[:87]}..."
        return phrase
    if len(description) > 110:
        return f"{description[:107]}..."
    return description


def _most_common(counter: Counter[Any], limit: int) -> list[tuple[Any, int]]:
    if limit == 0:
        return counter.most_common()
    return counter.most_common(limit)


def _go_name(go_id: str, go_terms: dict[str, dict[str, Any]]) -> str:
    return str(go_terms.get(go_id, {}).get("name") or go_id)


def _format_go_set_label(go_ids: tuple[str, ...], names: tuple[str, ...]) -> str:
    visible_names = "<br>".join(names[:3])
    if len(names) > 3:
        visible_names = f"{visible_names}<br>+{len(names) - 3} more"
    return f"{visible_names}<br>{';'.join(go_ids)}"


def _write_population_tsv(path: Path, groups: list[PopulationGroup]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "group_id", "label", "member_count", "details"])
        for rank, group in enumerate(groups, start=1):
            writer.writerow(
                [rank, group.group_id, group.label, group.count, group.details]
            )


def _write_treemap(
    *,
    output: Path,
    groups: list[PopulationGroup],
    title: str,
    root_label: str,
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
        raise SystemExit(f"No groups available for {title!r}.")

    ids = ["all", *(group.group_id for group in groups)]
    labels = [root_label, *(group.label for group in groups)]
    parents = ["", *("all" for _group in groups)]
    values = [sum(group.count for group in groups), *(group.count for group in groups)]
    customdata = [
        [""],
        *([group.hover_label] for group in groups),
    ]
    hovertemplates = [
        "<b>%{label}</b><br>Members shown: %{value:,}<extra></extra>",
        *(
            "<b>%{customdata[0]}</b><br>Members: %{value:,}<extra></extra>"
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
            marker={"line": {"width": 1, "color": "white"}},
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
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
