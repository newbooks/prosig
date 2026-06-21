from __future__ import annotations

import json
from collections import deque
from dataclasses import asdict
from pathlib import Path
from typing import Annotated, Any

import typer

from prosig.go.build import MF_ROOT
from prosig.go.similarity import GoSimilarity

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


@inspect_app.command(name="go-summary")
def go_summary(
    go_graph: Annotated[
        Path,
        typer.Option("--go-graph", help="Path to the compact GO graph pickle."),
    ] = Path("go_graph.pkl"),
) -> None:
    """Summarize the GO graph and IC artifact."""
    similarity = _load_go_similarity(go_graph)
    meta = similarity.meta
    typer.echo(f"path\t{go_graph}")
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
        Path,
        typer.Option("--go-graph", help="Path to the compact GO graph pickle."),
    ] = Path("go_graph.pkl"),
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
    similarity = _load_go_similarity(go_graph)
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
        Path,
        typer.Option("--go-graph", help="Path to the compact GO graph pickle."),
    ] = Path("go_graph.pkl"),
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

    similarity = _load_go_similarity(go_graph)
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


def _format_score(score: float | None) -> str:
    if score is None:
        return "NA"
    return f"{score:.4f}".rstrip("0").rstrip(".")


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
