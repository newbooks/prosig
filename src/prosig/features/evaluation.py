"""Evaluate scalar feature tables against GO-derived clusters."""

from __future__ import annotations

import csv
import logging
import math
from collections.abc import Callable, Hashable
from dataclasses import dataclass
from pathlib import Path

from prosig.features.quality import (
    compute_compactness_score,
    compute_gradient_score,
    compute_separation_score,
    compute_specificity_score,
)
from prosig.go.similarity import GoSimilarity, load_accession_mf_go_terms

LOGGER = logging.getLogger(__name__)
FEATURE_SCORE_COLUMNS = [
    "feature",
    "compactness",
    "separation",
    "gradient",
    "specificity",
]
DEFAULT_MIN_CLUSTER_SIZE = 10


@dataclass(frozen=True)
class FeatureEvaluationResult:
    """Feature evaluation output and filtering summary."""

    output_file: Path
    feature_count: int
    input_rows: int
    input_clusters: int
    retained_rows: int
    retained_clusters: int
    skipped_clusters: int


@dataclass(frozen=True)
class _FeatureTable:
    rows: list[dict[str, str]]
    feature_columns: list[str]


def evaluate_feature_file(
    feature_file: str | Path,
    *,
    output_file: str | Path = "feature_scores.tsv",
    go_graph_file: str | Path = "go_graph.pkl",
    accession_go_file: str | Path = "accession_mf_go.tsv",
    min_cluster_size: int = DEFAULT_MIN_CLUSTER_SIZE,
) -> FeatureEvaluationResult:
    """Evaluate numeric feature columns and write a feature score TSV."""
    if min_cluster_size < 1:
        raise ValueError("minimum cluster size must be at least 1")

    feature_path = Path(feature_file)
    output_path = Path(output_file)
    table = _load_feature_table(feature_path)
    input_clusters = len({row["cluster_id"] for row in table.rows})
    LOGGER.info(
        "Loaded %s feature rows across %s clusters with %s feature columns from %s",
        f"{len(table.rows):,}",
        f"{input_clusters:,}",
        f"{len(table.feature_columns):,}",
        feature_path,
    )

    retained_rows, skipped_sizes = _filter_rows_by_unique_cluster_size(
        table.rows,
        min_cluster_size=min_cluster_size,
    )
    _log_cluster_filtering(skipped_sizes, min_cluster_size=min_cluster_size)
    retained_clusters = len({row["cluster_id"] for row in retained_rows})
    if retained_clusters < 2:
        raise ValueError(
            "At least two clusters must remain after filtering by minimum "
            f"cluster size ({min_cluster_size})"
        )
    LOGGER.info(
        "Using %s clusters and %s feature rows for feature evaluation after filtering",
        f"{retained_clusters:,}",
        f"{len(retained_rows):,}",
    )

    similarity = GoSimilarity.from_pickle(go_graph_file)
    accession_terms = load_accession_mf_go_terms(accession_go_file)
    centroid_by_cluster = _cluster_centroid_accessions(
        retained_rows,
        accession_terms=accession_terms,
        similarity=similarity,
    )
    functional_similarity_fn = _go_centroid_similarity_fn(
        centroid_by_cluster=centroid_by_cluster,
        accession_terms=accession_terms,
        similarity=similarity,
    )

    grouped_values = _feature_values_by_cluster(retained_rows, table.feature_columns)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FEATURE_SCORE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for feature_column in table.feature_columns:
            LOGGER.info("Evaluating feature: %s", feature_column)
            clusters = grouped_values[feature_column]
            row = {
                "feature": feature_column,
                **compute_compactness_score(clusters),
                **compute_separation_score(clusters),
                **compute_gradient_score(
                    clusters,
                    functional_similarity_fn=functional_similarity_fn,
                ),
                **compute_specificity_score(clusters),
            }
            writer.writerow(_format_score_row(row))

    LOGGER.info(
        "Wrote feature evaluation metrics for %s features to %s",
        f"{len(table.feature_columns):,}",
        output_path,
    )
    return FeatureEvaluationResult(
        output_file=output_path,
        feature_count=len(table.feature_columns),
        input_rows=len(table.rows),
        input_clusters=input_clusters,
        retained_rows=len(retained_rows),
        retained_clusters=retained_clusters,
        skipped_clusters=len(skipped_sizes),
    )


def _load_feature_table(path: Path) -> _FeatureTable:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required_columns = {"member_id", "cluster_id"}
        missing_columns = required_columns - set(fieldnames)
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"{path} must contain columns: {missing}")
        feature_columns = [
            column for column in fieldnames if column not in required_columns
        ]
        if not feature_columns:
            raise ValueError(f"{path} must contain at least one feature column")

        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            member_id = str(row.get("member_id", "")).strip()
            cluster_id = str(row.get("cluster_id", "")).strip()
            if not member_id:
                raise ValueError(f"{path} line {line_number} has empty member_id")
            if not cluster_id:
                raise ValueError(f"{path} line {line_number} has empty cluster_id")
            parsed_row = {"member_id": member_id, "cluster_id": cluster_id}
            for feature_column in feature_columns:
                raw_value = str(row.get(feature_column, "")).strip()
                try:
                    numeric_value = float(raw_value)
                except ValueError as exc:
                    raise ValueError(
                        f"Feature column is not numeric: {feature_column}"
                    ) from exc
                parsed_row[feature_column] = str(numeric_value)
            rows.append(parsed_row)

    return _FeatureTable(rows=rows, feature_columns=feature_columns)


def _filter_rows_by_unique_cluster_size(
    rows: list[dict[str, str]],
    *,
    min_cluster_size: int,
) -> tuple[list[dict[str, str]], dict[str, int]]:
    members_by_cluster: dict[str, set[str]] = {}
    for row in rows:
        members_by_cluster.setdefault(row["cluster_id"], set()).add(row["member_id"])

    retained_clusters = {
        cluster_id
        for cluster_id, members in members_by_cluster.items()
        if len(members) >= min_cluster_size
    }
    skipped_sizes = {
        cluster_id: len(members)
        for cluster_id, members in members_by_cluster.items()
        if cluster_id not in retained_clusters
    }
    retained_rows = [
        row for row in rows
        if row["cluster_id"] in retained_clusters
    ]
    return retained_rows, skipped_sizes


def _log_cluster_filtering(
    skipped_sizes: dict[str, int],
    *,
    min_cluster_size: int,
) -> None:
    if skipped_sizes:
        LOGGER.info(
            "Skipped %s clusters with fewer than %s unique members",
            f"{len(skipped_sizes):,}",
            f"{min_cluster_size:,}",
        )
        for cluster_id, size in sorted(skipped_sizes.items()):
            LOGGER.debug(
                "Skipped cluster %s with %s unique members",
                cluster_id,
                f"{size:,}",
            )
        return
    LOGGER.info(
        "Skipped 0 clusters with fewer than %s unique members",
        f"{min_cluster_size:,}",
    )


def _cluster_centroid_accessions(
    rows: list[dict[str, str]],
    *,
    accession_terms: dict[str, tuple[str, ...]],
    similarity: GoSimilarity,
) -> dict[str, str]:
    missing_accessions = sorted(
        {
            row["member_id"]
            for row in rows
            if not accession_terms.get(row["member_id"])
        }
    )
    if missing_accessions:
        raise ValueError(
            "No MF GO terms found for member(s): "
            + ", ".join(missing_accessions[:20])
        )

    members_by_cluster: dict[str, set[str]] = {}
    for row in rows:
        members_by_cluster.setdefault(row["cluster_id"], set()).add(row["member_id"])

    centroid_by_cluster: dict[str, str] = {}
    for cluster_id, members in sorted(members_by_cluster.items()):
        accessions = sorted(members)
        if len(accessions) == 1:
            centroid_by_cluster[cluster_id] = accessions[0]
            continue
        mean_similarities = []
        for accession in accessions:
            scores = [
                _go_set_similarity_value(
                    similarity,
                    accession_terms[accession],
                    accession_terms[other_accession],
                )
                for other_accession in accessions
                if other_accession != accession
            ]
            mean_score = sum(scores) / len(scores) if scores else 0.0
            mean_similarities.append((mean_score, accession))
        best_score = max(score for score, _ in mean_similarities)
        centroid_by_cluster[cluster_id] = min(
            accession
            for score, accession in mean_similarities
            if score == best_score
        )

    LOGGER.info(
        "Selected GO centroid accessions for %s clusters",
        f"{len(centroid_by_cluster):,}",
    )
    return centroid_by_cluster


def _go_centroid_similarity_fn(
    *,
    centroid_by_cluster: dict[str, str],
    accession_terms: dict[str, tuple[str, ...]],
    similarity: GoSimilarity,
) -> Callable[[Hashable, Hashable], float]:
    cache: dict[tuple[str, str], float] = {}

    def cluster_similarity(cluster_a: Hashable, cluster_b: Hashable) -> float:
        cluster_id_a = str(cluster_a)
        cluster_id_b = str(cluster_b)
        cache_key = tuple(sorted((cluster_id_a, cluster_id_b)))
        if cache_key not in cache:
            accession_a = centroid_by_cluster[cluster_id_a]
            accession_b = centroid_by_cluster[cluster_id_b]
            cache[cache_key] = _go_set_similarity_value(
                similarity,
                accession_terms[accession_a],
                accession_terms[accession_b],
            )
        return cache[cache_key]

    return cluster_similarity


def _go_set_similarity_value(
    similarity: GoSimilarity,
    terms_a: tuple[str, ...],
    terms_b: tuple[str, ...],
) -> float:
    score = similarity.set_lin_amb(terms_a, terms_b)
    if score is None:
        return 0.0
    score = float(score)
    return max(0.0, min(1.0, score))


def _feature_values_by_cluster(
    rows: list[dict[str, str]],
    feature_columns: list[str],
) -> dict[str, dict[str, list[float]]]:
    grouped_values = {
        feature_column: {}
        for feature_column in feature_columns
    }
    for row in rows:
        cluster_id = row["cluster_id"]
        for feature_column in feature_columns:
            value = float(row[feature_column])
            if math.isfinite(value):
                grouped_values[feature_column].setdefault(cluster_id, []).append(value)
    return grouped_values


def _format_score_row(row: dict[str, float | str]) -> dict[str, str]:
    return {
        "feature": str(row["feature"]),
        "compactness": _format_metric(row["compactness"]),
        "separation": _format_metric(row["separation"]),
        "gradient": _format_metric(row["gradient"]),
        "specificity": _format_metric(row["specificity"]),
    }


def _format_metric(value: float | str) -> str:
    numeric_value = float(value)
    if math.isnan(numeric_value):
        return "nan"
    return f"{numeric_value:.4f}"
