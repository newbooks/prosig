#!/usr/bin/env python3
"""Plot a 2D metric MDS embedding from a cluster similarity matrix."""

from __future__ import annotations

import argparse
import os
import sys
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


def main() -> None:
    args = _parse_args()
    similarity = _read_similarity_matrix(args.input)
    distance = _similarity_to_distance(similarity)
    _validate_distance_matrix(distance, atol=args.atol)

    coordinates = _embed_mds(
        distance,
        random_seed=args.seed,
        max_iter=args.max_iter,
        n_init=args.n_init,
    )
    coordinates_path = args.out_prefix.with_name(
        f"{args.out_prefix.name}_coordinates.tsv"
    )
    _write_coordinates(coordinates_path, coordinates)
    _plot_mds(
        coordinates=coordinates,
        similarity=similarity,
        out_prefix=args.out_prefix,
        draw_lines=args.draw_lines,
        line_threshold=args.line_threshold,
        max_lines=args.max_lines,
        max_labels=args.max_labels,
        seed=args.seed,
        min_radius=args.min_radius,
        max_radius=args.max_radius,
        width=args.width,
        height=args.height,
        dpi=args.dpi,
    )
    print(f"Wrote {_append_extension(args.out_prefix, '.png')}")
    print(f"Wrote {_append_extension(args.out_prefix, '.pdf')}")
    print(f"Wrote {coordinates_path}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Embed a tab-separated functional similarity matrix into 2D using "
            "metric MDS and save PNG/PDF plots plus coordinates."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Tab-separated similarity matrix with cluster IDs in row/column headers.",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("cluster_mds"),
        help=(
            "Output prefix. Writes <prefix>.png, <prefix>.pdf, and "
            "<prefix>_coordinates.tsv."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Fixed random seed for reproducible MDS initialization.",
    )
    parser.add_argument("--n-init", type=int, default=8)
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument(
        "--atol",
        type=float,
        default=1e-8,
        help="Absolute tolerance used for diagonal and symmetry validation.",
    )
    parser.add_argument(
        "--draw-lines",
        action="store_true",
        help="Draw light gray pairwise lines with alpha/width scaled by similarity.",
    )
    parser.add_argument(
        "--line-threshold",
        type=float,
        default=0.0,
        help="Only draw lines for similarities at or above this value.",
    )
    parser.add_argument(
        "--max-lines",
        type=int,
        default=20000,
        help="Maximum lines to draw after thresholding. Use 0 for all lines.",
    )
    parser.add_argument(
        "--max-labels",
        type=int,
        default=200,
        help=(
            "Maximum point labels to draw. Use 0 for all labels. Labels are "
            "always written to the coordinates TSV."
        ),
    )
    parser.add_argument(
        "--min-radius",
        type=float,
        default=0.05,
        help="Minimum visual circle radius in data units.",
    )
    parser.add_argument(
        "--max-radius",
        type=float,
        default=0.25,
        help="Maximum visual circle radius in data units.",
    )
    parser.add_argument("--width", type=float, default=9.0, help="Figure width.")
    parser.add_argument("--height", type=float, default=7.0, help="Figure height.")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    if args.n_init < 1:
        parser.error("--n-init must be at least 1")
    if args.max_iter < 1:
        parser.error("--max-iter must be at least 1")
    if args.atol <= 0:
        parser.error("--atol must be positive")
    if not 0.0 <= args.line_threshold <= 1.0:
        parser.error("--line-threshold must be between 0 and 1")
    if args.max_lines < 0:
        parser.error("--max-lines must be at least 0")
    if args.max_labels < 0:
        parser.error("--max-labels must be at least 0")
    if args.min_radius <= 0:
        parser.error("--min-radius must be positive")
    if args.max_radius < args.min_radius:
        parser.error("--max-radius must be greater than or equal to --min-radius")
    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be positive")
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")
    return args


def _read_similarity_matrix(path: Path):
    try:
        import pandas as pd
    except ImportError as exc:
        raise SystemExit(
            "pandas is required. Install it with: python -m pip install pandas"
        ) from exc

    try:
        matrix = pd.read_csv(path, sep="\t", index_col=0, converters={0: str})
    except FileNotFoundError as exc:
        raise SystemExit(f"Input matrix not found: {path}") from exc

    if matrix.empty:
        raise SystemExit(f"Input matrix is empty: {path}")
    try:
        matrix = matrix.apply(pd.to_numeric, errors="raise")
    except Exception as exc:
        raise SystemExit(
            f"Input matrix contains non-numeric similarity values: {path}"
        ) from exc

    # Column headers are always read as strings, while pandas may infer
    # numeric-looking row IDs as integers. Keep both axes in the same string
    # domain so downstream label-based lookups work consistently.
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix


def _similarity_to_distance(similarity):
    return 1.0 - similarity


def _validate_distance_matrix(distance, *, atol: float) -> None:
    import numpy as np

    row_ids = list(distance.index.astype(str))
    col_ids = list(distance.columns.astype(str))

    if distance.shape[0] != distance.shape[1]:
        raise SystemExit(
            "Similarity matrix must be square; got "
            f"{distance.shape[0]} rows and {distance.shape[1]} columns."
        )
    if row_ids != col_ids:
        raise SystemExit(
            "Row and column cluster IDs must match in the same order."
        )

    values = distance.to_numpy(dtype=float)
    similarity_values = 1.0 - values
    if np.any(similarity_values < -atol) or np.any(similarity_values > 1.0 + atol):
        raise SystemExit("Similarity values must be between 0 and 1.")

    diagonal = np.diag(values)
    if not np.allclose(diagonal, 0.0, atol=atol):
        raise SystemExit("Diagonal distances must be zero.")
    if not np.allclose(values, values.T, atol=atol):
        max_delta = float(np.max(np.abs(values - values.T)))
        raise SystemExit(
            "Distance matrix must be symmetric; maximum absolute difference is "
            f"{max_delta:.6g}."
        )


def _embed_mds(distance, *, random_seed: int, max_iter: int, n_init: int):
    try:
        from sklearn.manifold import MDS
    except ImportError as exc:
        raise SystemExit(
            "scikit-learn is required. Install it with: "
            "python -m pip install scikit-learn"
        ) from exc

    kwargs = {
        "n_components": 2,
        "metric": True,
        "dissimilarity": "precomputed",
        "random_state": random_seed,
        "n_init": n_init,
        "max_iter": max_iter,
        "eps": 1e-9,
    }
    try:
        model = MDS(normalized_stress="auto", **kwargs)
    except TypeError:
        model = MDS(**kwargs)

    embedding = model.fit_transform(distance.to_numpy(dtype=float))
    try:
        import pandas as pd
    except ImportError as exc:
        raise SystemExit(
            "pandas is required. Install it with: python -m pip install pandas"
        ) from exc

    return pd.DataFrame(
        {
            "cluster_id": distance.index.astype(str),
            "mds1": embedding[:, 0],
            "mds2": embedding[:, 1],
        }
    )


def _write_coordinates(path: Path, coordinates) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    coordinates.to_csv(path, sep="\t", index=False, float_format="%.8f")


def _append_extension(prefix: Path, extension: str) -> Path:
    return prefix.with_name(f"{prefix.name}{extension}")


def _plot_mds(
    *,
    coordinates,
    similarity,
    out_prefix: Path,
    draw_lines: bool,
    line_threshold: float,
    max_lines: int,
    max_labels: int,
    seed: int,
    min_radius: float,
    max_radius: float,
    width: float,
    height: float,
    dpi: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required. Install it with: python -m pip install matplotlib"
        ) from exc

    _apply_plot_style()
    plot_data = coordinates.set_index("cluster_id")
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)

    if draw_lines:
        _draw_similarity_lines(
            ax,
            plot_data=plot_data,
            similarity=similarity,
            threshold=line_threshold,
            max_lines=max_lines,
        )

    _draw_cluster_circles(
        ax,
        plot_data=plot_data,
        seed=seed,
        min_radius=min_radius,
        max_radius=max_radius,
    )
    _draw_point_labels(ax, plot_data=plot_data, max_labels=max_labels)
    _center_plot_on_origin(ax, plot_data=plot_data, padding=max_radius)

    ax.set_title("2D projection of challenging set", fontsize=16, pad=18)
    ax.text(
        0.5,
        1.01,
        "Distance = 1 − functional similarity",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10,
        color="#555555",
    )
    ax.set_xlabel("Distance")
    ax.set_ylabel("Distance")
    ax.grid(True, color="#e8e8e8", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_aspect("equal", adjustable="box")

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_append_extension(out_prefix, ".png"), dpi=dpi)
    fig.savefig(_append_extension(out_prefix, ".pdf"))
    plt.close(fig)


def _apply_plot_style() -> None:
    try:
        import seaborn as sns
    except ImportError:
        return

    sns.set_theme(
        context="paper",
        style="whitegrid",
        palette="deep",
        font_scale=1.05,
        rc={
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": "#e8e8e8",
            "grid.linewidth": 0.8,
        },
    )


def _draw_cluster_circles(
    ax,
    *,
    plot_data,
    seed: int,
    min_radius: float,
    max_radius: float,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    rng = np.random.default_rng(seed)
    cmap = plt.get_cmap("turbo")
    for _cluster_id, row in plot_data.iterrows():
        radius = float(rng.uniform(min_radius, max_radius))
        base_color = cmap(float(rng.uniform(0.05, 0.95)))
        image = _radial_circle_image(base_color, pixels=96)
        ax.imshow(
            image,
            extent=(
                row["mds1"] - radius,
                row["mds1"] + radius,
                row["mds2"] - radius,
                row["mds2"] + radius,
            ),
            origin="lower",
            interpolation="bilinear",
            zorder=3,
        )


def _radial_circle_image(base_color, *, pixels: int):
    import numpy as np

    yy, xx = np.ogrid[-1.0:1.0:complex(pixels), -1.0:1.0:complex(pixels)]
    radius = np.sqrt(xx * xx + yy * yy)
    clipped = np.clip(radius, 0.0, 1.0)
    center = np.array(base_color[:3]) * 0.42
    edge = 1.0 - (1.0 - np.array(base_color[:3])) * 0.28
    rgb = center * (1.0 - clipped[..., None]) + edge * clipped[..., None]
    alpha = np.where(radius <= 1.0, 0.56 * (1.0 - 0.25 * clipped), 0.0)
    return np.dstack([rgb, alpha])


def _center_plot_on_origin(ax, *, plot_data, padding: float) -> None:
    max_abs = max(
        float(abs(plot_data["mds1"]).max()),
        float(abs(plot_data["mds2"]).max()),
        padding,
    )
    limit = max_abs + padding * 1.25
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.axhline(0.0, color="#9a9a9a", linewidth=0.9, alpha=0.55, zorder=0)
    ax.axvline(0.0, color="#9a9a9a", linewidth=0.9, alpha=0.55, zorder=0)


def _draw_similarity_lines(
    ax,
    *,
    plot_data,
    similarity,
    threshold: float,
    max_lines: int,
) -> None:
    pairs = []
    cluster_ids = list(plot_data.index)
    for left, right in combinations(cluster_ids, 2):
        score = float(similarity.loc[left, right])
        if score >= threshold:
            pairs.append((score, left, right))

    pairs.sort(reverse=True)
    if max_lines and len(pairs) > max_lines:
        print(
            f"Drawing strongest {max_lines} of {len(pairs)} eligible pairwise lines. "
            "Use --max-lines 0 to draw all.",
            file=sys.stderr,
        )
        pairs = pairs[:max_lines]

    for score, left, right in pairs:
        x_values = [plot_data.loc[left, "mds1"], plot_data.loc[right, "mds1"]]
        y_values = [plot_data.loc[left, "mds2"], plot_data.loc[right, "mds2"]]
        ax.plot(
            x_values,
            y_values,
            color="#bdbdbd",
            alpha=0.08 + 0.42 * score,
            linewidth=0.2 + 1.8 * score,
            zorder=1,
        )


def _draw_point_labels(ax, *, plot_data, max_labels: int) -> None:
    rows = list(plot_data.iterrows())
    if max_labels and len(rows) > max_labels:
        print(
            f"Skipping point labels because there are {len(rows)} clusters. "
            f"Use --max-labels 0 to label all points.",
            file=sys.stderr,
        )
        return

    for cluster_id, row in rows:
        ax.annotate(
            cluster_id,
            (row["mds1"], row["mds2"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
            color="#222222",
            zorder=4,
        )


if __name__ == "__main__":
    main()
