"""Score card rendering for feature quality metrics."""

from __future__ import annotations

import csv
import math
import os
import re
import tempfile
from pathlib import Path

METRIC_COLUMNS = ["compactness", "separation", "gradient", "specificity"]
SUPPORTED_FORMATS = {"png", "svg"}
_BAR_COLOR = "#48c879"
_BAR_COLORS = [
    "#48c879",
    "#4f83cc",
    "#f59f3a",
    "#a66dd4",
    "#e05f72",
    "#36a3a1",
]
_TEXT_COLOR = "#5f6368"
_GRID_COLOR = "#d9dee5"
_REFERENCE_TICKS = [0.0, 0.25, 0.5, 0.75, 1.0]


def write_feature_quality_score_cards(
    *,
    metrics_file: Path,
    output_dir: Path,
    image_format: str = "png",
    features: list[str] | None = None,
    output_file: Path | None = None,
) -> list[Path]:
    """Write score cards from feature metric rows."""
    image_format = image_format.lower()
    if image_format not in SUPPORTED_FORMATS:
        raise ValueError("image_format must be one of: png, svg")
    if output_file is not None and not features:
        raise ValueError("output_file can only be used with selected features")
    if output_file is not None and output_file.suffix.lower() not in {
        f".{image_format}",
        "",
    }:
        raise ValueError(f"output_file suffix must match .{image_format}")

    metrics = _load_feature_quality_metrics(metrics_file)
    selected_metrics = _select_feature_metrics(metrics, features=features)
    if not selected_metrics:
        raise ValueError("No feature metrics available to plot")

    output_dir.mkdir(parents=True, exist_ok=True)
    if features:
        card_file = (
            output_file
            or output_dir / f"selected_features_score_card.{image_format}"
        )
        if card_file.suffix == "":
            card_file = card_file.with_suffix(f".{image_format}")
        card_file.parent.mkdir(parents=True, exist_ok=True)
        _plot_score_card_rows(selected_metrics, card_file)
        return [card_file]

    card_files = []
    for row in selected_metrics:
        card_file = output_dir / (
            f"{_slugify_filename(str(row['feature']))}_score_card.{image_format}"
        )
        _plot_score_card_rows([row], card_file)
        card_files.append(card_file)
    return card_files


def _load_feature_quality_metrics(metrics_file: Path) -> list[dict[str, str]]:
    with metrics_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        required = {"feature", *METRIC_COLUMNS}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise ValueError(
                f"{metrics_file} must contain metric columns: {missing}"
            )
        return [
            {key: value or "" for key, value in row.items()}
            for row in reader
        ]


def _select_feature_metrics(
    metrics: list[dict[str, str]],
    *,
    features: list[str] | None,
) -> list[dict[str, str]]:
    if not features:
        return metrics

    metrics_by_feature = {
        str(row["feature"]): row
        for row in metrics
    }
    missing_features = [
        feature
        for feature in features
        if feature not in metrics_by_feature
    ]
    if missing_features:
        raise ValueError(f"Feature(s) not found: {', '.join(missing_features)}")
    return [metrics_by_feature[feature] for feature in features]


def _plot_score_card_rows(metrics: list[dict[str, str]], output_file: Path) -> None:
    _configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Patch

    panel_count = len(metrics)
    figure_height = 4.4 if panel_count > 1 else 3.8
    fig, ax = plt.subplots(figsize=(9.5, figure_height), dpi=160)
    fig.patch.set_facecolor("white")

    if panel_count > 1:
        _draw_multi_feature_score_card(ax, metrics, FancyBboxPatch, Patch)
        fig.subplots_adjust(top=0.86, bottom=0.16, left=0.15, right=0.96)
    else:
        _draw_score_card_panel(ax, metrics[0], FancyBboxPatch)
        fig.subplots_adjust(top=0.96, bottom=0.18, left=0.15, right=0.96)

    fig.savefig(output_file, bbox_inches="tight", dpi=300)
    plt.close(fig)


def _configure_matplotlib() -> None:
    config_dir = Path(tempfile.gettempdir()) / "prosig-matplotlib"
    config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(config_dir))
    import matplotlib

    matplotlib.use("Agg")


def _draw_score_card_panel(ax, row: dict[str, str], rounded_bar_cls) -> None:
    values = _normalized_metric_values(row)
    metric_labels = [metric.title() for metric in METRIC_COLUMNS]
    y_positions = list(reversed(range(len(METRIC_COLUMNS))))

    ax.set_xlim(-0.01, 1.12)
    ax.set_ylim(-0.6, len(METRIC_COLUMNS) - 0.4)
    ax.set_facecolor("white")
    ax.set_title(
        str(row["feature"]),
        loc="left",
        fontsize=12,
        fontweight="semibold",
        color=_TEXT_COLOR,
        pad=12,
    )
    _style_score_axis(ax)

    for y_position, value, metric in zip(
        y_positions,
        values,
        METRIC_COLUMNS,
        strict=True,
    ):
        _draw_rounded_bar(
            ax,
            rounded_bar_cls,
            y_position=y_position,
            width=value,
            height=0.74,
            color=_BAR_COLOR,
        )
        ax.text(
            value + 0.018,
            y_position,
            _format_metric_label(row[metric]),
            ha="left",
            va="center",
            fontsize=10,
            color=_TEXT_COLOR,
        )

    ax.set_xticks(_REFERENCE_TICKS)
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=10)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(metric_labels, fontsize=10.5, color=_TEXT_COLOR)
    ax.tick_params(axis="x", colors=_TEXT_COLOR, length=0, pad=8)
    ax.tick_params(axis="y", length=0, pad=10)
    _hide_spines(ax)


def _draw_multi_feature_score_card(
    ax,
    metrics: list[dict[str, str]],
    rounded_bar_cls,
    patch_cls,
) -> None:
    metric_labels = [metric.title() for metric in METRIC_COLUMNS]
    metric_y_positions = list(reversed(range(len(METRIC_COLUMNS))))
    feature_count = len(metrics)
    group_height = 0.78
    slot_height = group_height / feature_count
    bar_height = min(0.24, slot_height * 0.72)
    colors = [_BAR_COLORS[index % len(_BAR_COLORS)] for index in range(feature_count)]
    legend_handles = [
        patch_cls(facecolor=color, edgecolor="none", label=str(row["feature"]))
        for color, row in zip(colors, metrics, strict=True)
    ]

    ax.set_xlim(-0.01, 1.16)
    ax.set_ylim(-0.6, len(METRIC_COLUMNS) - 0.4)
    ax.set_facecolor("white")
    ax.set_title(
        "Feature quality score card",
        loc="left",
        fontsize=12,
        fontweight="semibold",
        color=_TEXT_COLOR,
        pad=12,
    )
    _style_score_axis(ax)

    for feature_index, row in enumerate(metrics):
        values = _normalized_metric_values(row)
        color = colors[feature_index]
        offset = ((feature_count - 1) / 2 - feature_index) * slot_height
        for metric_y_position, value, metric in zip(
            metric_y_positions,
            values,
            METRIC_COLUMNS,
            strict=True,
        ):
            y_position = metric_y_position + offset
            _draw_rounded_bar(
                ax,
                rounded_bar_cls,
                y_position=y_position,
                width=value,
                height=bar_height,
                color=color,
            )
            ax.text(
                value + 0.014,
                y_position,
                _format_metric_label(row[metric]),
                ha="left",
                va="center",
                fontsize=8.8,
                color=_TEXT_COLOR,
            )

    ax.set_xticks(_REFERENCE_TICKS)
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=10)
    ax.set_yticks(metric_y_positions)
    ax.set_yticklabels(metric_labels, fontsize=10.5, color=_TEXT_COLOR)
    ax.tick_params(axis="x", colors=_TEXT_COLOR, length=0, pad=8)
    ax.tick_params(axis="y", length=0, pad=10)
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.13),
        frameon=False,
        fontsize=9,
        ncols=1,
    )
    _hide_spines(ax)


def _style_score_axis(ax) -> None:
    ax.set_axisbelow(True)
    ax.xaxis.grid(
        True,
        color=_GRID_COLOR,
        linestyle=(0, (4, 4)),
        linewidth=0.9,
        alpha=0.7,
    )
    ax.yaxis.grid(False)


def _hide_spines(ax) -> None:
    for spine in ax.spines.values():
        spine.set_visible(False)


def _draw_rounded_bar(
    ax,
    rounded_bar_cls,
    *,
    y_position: float,
    width: float,
    height: float,
    color: str,
) -> None:
    if width <= 0:
        return
    bar = rounded_bar_cls(
        (0.0, y_position - height / 2),
        width,
        height,
        boxstyle="round,pad=0,rounding_size=0.025",
        linewidth=0,
        facecolor=color,
        antialiased=True,
        zorder=3,
    )
    ax.add_patch(bar)


def _normalized_metric_values(row: dict[str, str]) -> list[float]:
    values = []
    for metric in METRIC_COLUMNS:
        try:
            value = float(row[metric])
        except ValueError:
            value = 0.0
        if math.isnan(value):
            value = 0.0
        values.append(max(0.0, min(1.0, value)))
    return values


def _format_metric_label(value: str) -> str:
    try:
        numeric_value = float(value)
    except ValueError:
        return "nan"
    if math.isnan(numeric_value):
        return "nan"
    return f"{numeric_value:.4f}"


def _slugify_filename(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return slug.strip("._") or "feature"
