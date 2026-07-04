"""Scalar feature quality metrics for GO-clustered protein data."""

from __future__ import annotations

import math
import statistics
from collections.abc import Callable, Hashable, Iterable, Mapping

from scipy.stats import spearmanr

ClusterFeatureValues = Mapping[Hashable, Iterable[float]]
FunctionalSimilarityFn = Callable[[Hashable, Hashable], float]


def _valid_feature_clusters(
    clusters: ClusterFeatureValues,
    *,
    min_cluster_size: int,
) -> dict[Hashable, list[float]]:
    valid_clusters: dict[Hashable, list[float]] = {}
    for cluster_id, values in clusters.items():
        if cluster_id is None:
            raise ValueError("cluster IDs must not be missing")
        numeric_values = [
            float(value)
            for value in values
            if value is not None and math.isfinite(float(value))
        ]
        if len(numeric_values) >= min_cluster_size:
            valid_clusters[cluster_id] = numeric_values
    return valid_clusters


def _sample_variance(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return float(statistics.variance(values))


def _sample_std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return float(statistics.stdev(values))


def compute_compactness_score(
    clusters: ClusterFeatureValues,
    *,
    min_cluster_size: int = 10,
    epsilon: float = 1e-12,
) -> dict[str, float]:
    """Return the median relative within-cluster tightness for one feature."""
    valid_clusters = _valid_feature_clusters(
        clusters,
        min_cluster_size=min_cluster_size,
    )
    if not valid_clusters:
        return {"compactness": float("nan")}

    values = [
        value
        for cluster_values in valid_clusters.values()
        for value in cluster_values
    ]
    total_variance = _sample_variance(values)
    if total_variance <= epsilon:
        return {"compactness": 0.0}

    cluster_scores = []
    for cluster_values in valid_clusters.values():
        cluster_variance = _sample_variance(cluster_values)
        relative_dispersion = cluster_variance / (total_variance + epsilon)
        cluster_scores.append(max(0.0, min(1.0, 1.0 - relative_dispersion)))

    return {"compactness": float(statistics.median(cluster_scores))}


def compute_separation_score(
    clusters: ClusterFeatureValues,
    *,
    k: float = 1.0,
    min_cluster_size: int = 2,
) -> dict[str, float]:
    """Return the fraction of cluster pairs with separated feature means."""
    valid_clusters = _valid_feature_clusters(
        clusters,
        min_cluster_size=min_cluster_size,
    )
    if len(valid_clusters) < 2:
        return {"separation": float("nan")}

    cluster_stats = [
        (
            cluster_id,
            float(statistics.fmean(values)),
            _sample_std(values),
        )
        for cluster_id, values in valid_clusters.items()
    ]
    separable_pairs = 0
    total_pairs = 0
    for index, (_, mean_i, std_i) in enumerate(cluster_stats):
        for _, mean_j, std_j in cluster_stats[index + 1:]:
            total_pairs += 1
            if abs(mean_i - mean_j) > k * (std_i + std_j):
                separable_pairs += 1

    if total_pairs == 0:
        return {"separation": float("nan")}
    return {"separation": float(separable_pairs / total_pairs)}


def compute_gradient_score(
    clusters: ClusterFeatureValues,
    *,
    functional_similarity_fn: FunctionalSimilarityFn,
    min_cluster_size: int = 2,
) -> dict[str, float]:
    """Return Spearman association between functional and feature distance."""
    valid_clusters = _valid_feature_clusters(
        clusters,
        min_cluster_size=min_cluster_size,
    )
    if len(valid_clusters) < 2:
        return {"gradient": float("nan")}

    cluster_means = [
        (cluster_id, float(statistics.fmean(values)))
        for cluster_id, values in valid_clusters.items()
    ]
    functional_distances: list[float] = []
    feature_distances: list[float] = []
    for index, (cluster_i, mean_i) in enumerate(cluster_means):
        for cluster_j, mean_j in cluster_means[index + 1:]:
            similarity = float(functional_similarity_fn(cluster_i, cluster_j))
            if not 0.0 <= similarity <= 1.0:
                raise ValueError(
                    "functional_similarity_fn must return values in [0, 1]"
                )
            functional_distances.append(1.0 - similarity)
            feature_distances.append(abs(mean_i - mean_j))

    if not functional_distances:
        return {"gradient": float("nan")}
    if len(set(functional_distances)) == 1:
        return {"gradient": float("nan")}
    if len(set(feature_distances)) == 1:
        return {"gradient": 0.0}

    rho = float(spearmanr(functional_distances, feature_distances).statistic)
    if math.isnan(rho):
        return {"gradient": float("nan")}
    return {"gradient": float(max(0.0, rho))}


def compute_specificity_score(
    clusters: ClusterFeatureValues,
    *,
    min_cluster_size: int = 2,
    neighbor_fraction: float = 0.05,
    neighbor_k: int | None = None,
    epsilon: float = 1e-12,
) -> dict[str, float]:
    """Return a normalized best local inter-cluster feature gap."""
    valid_clusters = _valid_feature_clusters(
        clusters,
        min_cluster_size=min_cluster_size,
    )
    if len(valid_clusters) < 2:
        return {"specificity": 0.0}

    values = [
        value
        for cluster_values in valid_clusters.values()
        for value in cluster_values
    ]
    global_std = _sample_std(values)
    if global_std <= epsilon:
        return {"specificity": 0.0}

    cluster_stats = [
        (
            cluster_id,
            float(statistics.fmean(cluster_values)),
            _sample_std(cluster_values),
        )
        for cluster_id, cluster_values in valid_clusters.items()
    ]
    cluster_count = len(cluster_stats)
    local_neighbor_count = (
        round(neighbor_fraction * cluster_count)
        if neighbor_k is None
        else neighbor_k
    )
    local_neighbor_count = max(1, local_neighbor_count)
    local_neighbor_count = min(local_neighbor_count, cluster_count - 1)

    local_gaps: list[float] = []
    for cluster_i, mean_i, std_i in cluster_stats:
        gaps = []
        for cluster_j, mean_j, std_j in cluster_stats:
            if cluster_i == cluster_j:
                continue
            raw_gap = abs(mean_i - mean_j) - (std_i + std_j)
            gaps.append(max(0.0, raw_gap))
        local_gaps.append(float(max(sorted(gaps)[:local_neighbor_count])))

    specificity_raw = max(local_gaps)
    specificity = specificity_raw / (specificity_raw + global_std + epsilon)
    return {"specificity": float(specificity)}
