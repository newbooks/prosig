"""Motif-cluster score board construction for function prediction."""

from __future__ import annotations

import csv
import json
import math
import pickle
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

JEFFREYS_PRIOR = 0.5
CALIBRATION_WEIGHT_THRESHOLDS = (2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)


@dataclass(frozen=True)
class MotifClusterWeight:
    """Stored statistics for one positive motif-cluster association."""

    motif_id: str
    cluster_id: str
    TP: int
    FP: int
    FN: int
    TN: int
    support: int
    cluster_frequency: float
    background_frequency: float
    weight: float


@dataclass(frozen=True)
class MotifScoreboardStats:
    """Summary counts for motif-cluster scoreboard construction."""

    accessions_in_clusters: int
    clusters: int
    clusters_below_min_size: int
    eligible_clusters: int
    motifs: int
    motif_hit_rows: int
    duplicate_motif_hits: int
    motif_hits_outside_clusters: int
    potential_motif_cluster_combos: int
    ignored_cluster_size: int
    ignored_low_support: int
    ignored_background_unavailable: int
    ignored_non_positive_weight: int
    total_nonzero_weights: int
    stored_weights: int
    min_cluster_size: int
    min_support: int
    calibration: tuple[MotifScoreboardCalibration, ...]


@dataclass(frozen=True)
class MotifScoreboardCalibration:
    """Internal motif-inference calibration at one weight threshold."""

    weight_threshold: float
    eligible_accessions: int
    covered_accessions: int
    top1_correct_accessions: int
    top3_correct_accessions: int
    set_correct_accessions: int
    coverage: float
    top1_accuracy: float | None
    top3_accuracy: float | None
    set_accuracy: float | None
    avg_predictions: float


def build_motif_cluster_scoreboard(
    *,
    cluster_file: str | Path = "clusters.tsv",
    motif_hits_file: str | Path = "motif_features.tsv",
    output_file: str | Path = "motif_cluster_scoreboard.pkl",
    meta_file: str | Path | None = "motif_cluster_scoreboard_meta.json",
    min_cluster_size: int = 10,
    min_support: int = 5,
) -> MotifScoreboardStats:
    """Build and write positive motif-cluster log-enrichment weights.

    The motif hit input is the sparse ProSig motif feature TSV. Each
    accession-motif row is treated as binary motif presence, regardless of
    match count. Only positive weights are stored in the pickle artifact.
    """
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be at least 1")
    if min_support < 1:
        raise ValueError("min_support must be at least 1")

    cluster_path = Path(cluster_file)
    motif_hits_path = Path(motif_hits_file)
    output_path = Path(output_file)
    meta_path = Path(meta_file) if meta_file is not None else None

    cluster_by_accession = _load_cluster_membership(cluster_path)
    members_by_cluster = _members_by_cluster(cluster_by_accession)
    cluster_sizes = {
        cluster_id: len(members)
        for cluster_id, members in members_by_cluster.items()
    }
    motif_accessions, motif_hit_rows, duplicate_hits, outside_hits = _load_motif_hits(
        motif_hits_path,
        clustered_accessions=set(cluster_by_accession),
    )

    weights: dict[str, dict[str, dict[str, Any]]] = {}
    total_accessions = len(cluster_by_accession)
    small_clusters = {
        cluster_id
        for cluster_id, size in cluster_sizes.items()
        if size < min_cluster_size
    }
    eligible_clusters = [
        cluster_id
        for cluster_id in sorted(members_by_cluster)
        if cluster_id not in small_clusters
    ]

    ignored_low_support = 0
    ignored_background_unavailable = 0
    ignored_non_positive_weight = 0
    stored_weights = 0

    for motif_id in sorted(motif_accessions):
        accessions_with_motif = motif_accessions[motif_id]
        motif_total = len(accessions_with_motif)
        for cluster_id in eligible_clusters:
            cluster_members = members_by_cluster[cluster_id]
            cluster_size = cluster_sizes[cluster_id]
            tp = len(accessions_with_motif & cluster_members)
            if tp < min_support:
                ignored_low_support += 1
                continue
            fn = cluster_size - tp
            fp = motif_total - tp
            outside_count = total_accessions - cluster_size
            if outside_count == 0:
                ignored_background_unavailable += 1
                continue
            tn = outside_count - fp
            cluster_frequency = _smoothed_frequency(
                tp,
                cluster_size,
                pseudocount=JEFFREYS_PRIOR,
            )
            background_frequency = _smoothed_frequency(
                fp,
                outside_count,
                pseudocount=JEFFREYS_PRIOR,
            )
            weight = _log2_enrichment(cluster_frequency, background_frequency)
            if weight <= 0.0:
                ignored_non_positive_weight += 1
                continue

            record = MotifClusterWeight(
                motif_id=motif_id,
                cluster_id=cluster_id,
                TP=tp,
                FP=fp,
                FN=fn,
                TN=tn,
                support=tp,
                cluster_frequency=cluster_frequency,
                background_frequency=background_frequency,
                weight=weight,
            )
            weights.setdefault(motif_id, {})[cluster_id] = asdict(record)
            stored_weights += 1

    calibration = _calibrate_motif_predictions(
        cluster_by_accession=cluster_by_accession,
        cluster_sizes=cluster_sizes,
        motif_accessions=motif_accessions,
        weights=weights,
        min_cluster_size=min_cluster_size,
    )
    stats = MotifScoreboardStats(
        accessions_in_clusters=total_accessions,
        clusters=len(members_by_cluster),
        clusters_below_min_size=len(small_clusters),
        eligible_clusters=len(eligible_clusters),
        motifs=len(motif_accessions),
        motif_hit_rows=motif_hit_rows,
        duplicate_motif_hits=duplicate_hits,
        motif_hits_outside_clusters=outside_hits,
        potential_motif_cluster_combos=len(members_by_cluster) * len(motif_accessions),
        ignored_cluster_size=len(small_clusters) * len(motif_accessions),
        ignored_low_support=ignored_low_support,
        ignored_background_unavailable=ignored_background_unavailable,
        ignored_non_positive_weight=ignored_non_positive_weight,
        total_nonzero_weights=stored_weights,
        stored_weights=stored_weights,
        min_cluster_size=min_cluster_size,
        min_support=min_support,
        calibration=calibration,
    )

    artifact = {
        "schema_version": "1.0",
        "kind": "motif_cluster_scoreboard",
        "parameters": {
            "min_cluster_size": min_cluster_size,
            "min_support": min_support,
            "pseudocount": JEFFREYS_PRIOR,
            "smoothing": "Jeffreys prior",
            "weight": (
                "log2(((TP + 0.5) / (TP + FN + 1)) / "
                "((FP + 0.5) / (FP + TN + 1)))"
            ),
            "motif_presence": "binary",
        },
        "weights": weights,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as handle:
        pickle.dump(artifact, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if meta_path is not None:
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        meta = {
            "schema_version": "1.0",
            "kind": "motif_cluster_scoreboard_meta",
            "inputs": {
                "clusters": str(cluster_path),
                "motif_hits": str(motif_hits_path),
            },
            "output": str(output_path),
            "stats": asdict(stats),
        }
        meta_path.write_text(
            json.dumps(meta, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    return stats


def _load_cluster_membership(path: Path) -> dict[str, str]:
    cluster_by_accession: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        required = {"member_id", "cluster_id"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise ValueError(f"{path} missing required column(s): {missing}")
        for row in reader:
            accession = str(row.get("member_id", "")).strip()
            cluster_id = str(row.get("cluster_id", "")).strip()
            if not accession or not cluster_id:
                continue
            if accession in cluster_by_accession:
                raise ValueError(f"Cluster member appears more than once: {accession}")
            cluster_by_accession[accession] = cluster_id
    return cluster_by_accession


def _members_by_cluster(cluster_by_accession: dict[str, str]) -> dict[str, set[str]]:
    members: dict[str, set[str]] = {}
    for accession, cluster_id in cluster_by_accession.items():
        members.setdefault(cluster_id, set()).add(accession)
    return members


def _load_motif_hits(
    path: Path,
    *,
    clustered_accessions: set[str],
) -> tuple[dict[str, set[str]], int, int, int]:
    motif_accessions: dict[str, set[str]] = {}
    motif_hit_rows = 0
    duplicate_hits = 0
    outside_hits = 0
    seen_hits: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        required = {"accession", "motif_id"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise ValueError(f"{path} missing required column(s): {missing}")
        for row in reader:
            accession = str(row.get("accession", "")).strip()
            if accession.startswith("#"):
                continue
            motif_id = str(row.get("motif_id", "")).strip()
            if not accession or not motif_id:
                continue
            if _hit_count(row) <= 0:
                continue
            motif_hit_rows += 1
            if accession not in clustered_accessions:
                outside_hits += 1
                continue
            hit_key = (accession, motif_id)
            if hit_key in seen_hits:
                duplicate_hits += 1
                continue
            seen_hits.add(hit_key)
            motif_accessions.setdefault(motif_id, set()).add(accession)
    return motif_accessions, motif_hit_rows, duplicate_hits, outside_hits


def _calibrate_motif_predictions(
    *,
    cluster_by_accession: dict[str, str],
    cluster_sizes: dict[str, int],
    motif_accessions: dict[str, set[str]],
    weights: dict[str, dict[str, dict[str, Any]]],
    min_cluster_size: int,
) -> tuple[MotifScoreboardCalibration, ...]:
    eligible_accession_clusters = {
        accession: cluster_id
        for accession, cluster_id in cluster_by_accession.items()
        if cluster_sizes[cluster_id] >= min_cluster_size
    }
    motif_hits_by_accession = _motif_hits_by_accession(motif_accessions)
    calibration: list[MotifScoreboardCalibration] = []
    for threshold in CALIBRATION_WEIGHT_THRESHOLDS:
        motif_predictions = _motif_predictions_at_threshold(weights, threshold)
        covered_accessions = 0
        top1_correct_accessions = 0
        top3_correct_accessions = 0
        set_correct_accessions = 0
        total_predictions = 0
        for accession, true_cluster_id in eligible_accession_clusters.items():
            prediction_scores: dict[str, float] = {}
            for motif_id in motif_hits_by_accession.get(accession, ()):
                for cluster_id, weight in motif_predictions.get(motif_id, {}).items():
                    prediction_scores[cluster_id] = max(
                        prediction_scores.get(cluster_id, -math.inf),
                        weight,
                    )
            if not prediction_scores:
                continue
            covered_accessions += 1
            total_predictions += len(prediction_scores)
            ranked_predictions = sorted(
                prediction_scores,
                key=lambda cluster_id: (-prediction_scores[cluster_id], cluster_id),
            )
            if ranked_predictions[0] == true_cluster_id:
                top1_correct_accessions += 1
            if true_cluster_id in ranked_predictions[:3]:
                top3_correct_accessions += 1
            if true_cluster_id in prediction_scores:
                set_correct_accessions += 1

        eligible_accessions = len(eligible_accession_clusters)
        coverage = (
            covered_accessions / eligible_accessions
            if eligible_accessions
            else 0.0
        )
        top1_accuracy = (
            top1_correct_accessions / covered_accessions
            if covered_accessions
            else None
        )
        top3_accuracy = (
            top3_correct_accessions / covered_accessions
            if covered_accessions
            else None
        )
        set_accuracy = (
            set_correct_accessions / covered_accessions
            if covered_accessions
            else None
        )
        avg_predictions = (
            total_predictions / covered_accessions
            if covered_accessions
            else 0.0
        )
        calibration.append(
            MotifScoreboardCalibration(
                weight_threshold=threshold,
                eligible_accessions=eligible_accessions,
                covered_accessions=covered_accessions,
                top1_correct_accessions=top1_correct_accessions,
                top3_correct_accessions=top3_correct_accessions,
                set_correct_accessions=set_correct_accessions,
                coverage=coverage,
                top1_accuracy=top1_accuracy,
                top3_accuracy=top3_accuracy,
                set_accuracy=set_accuracy,
                avg_predictions=avg_predictions,
            )
        )
    return tuple(calibration)


def _motif_hits_by_accession(
    motif_accessions: dict[str, set[str]],
) -> dict[str, set[str]]:
    hits_by_accession: dict[str, set[str]] = {}
    for motif_id, accessions in motif_accessions.items():
        for accession in accessions:
            hits_by_accession.setdefault(accession, set()).add(motif_id)
    return hits_by_accession


def _motif_predictions_at_threshold(
    weights: dict[str, dict[str, dict[str, Any]]],
    threshold: float,
) -> dict[str, dict[str, float]]:
    predictions: dict[str, dict[str, float]] = {}
    for motif_id, cluster_weights in weights.items():
        for cluster_id, record in cluster_weights.items():
            weight = float(record.get("weight", 0.0))
            if weight >= threshold:
                predictions.setdefault(motif_id, {})[cluster_id] = weight
    return predictions


def _hit_count(row: dict[str, Any]) -> int:
    raw_present = row.get("motif_present")
    if raw_present is not None and str(raw_present).strip() != "":
        return 1 if _parse_bool(str(raw_present).strip()) else 0
    raw_count = row.get("count")
    if raw_count is None or str(raw_count).strip() == "":
        return 1
    try:
        return int(str(raw_count).strip())
    except ValueError as exc:
        raise ValueError(f"Invalid motif hit count: {raw_count}") from exc


def _parse_bool(value: str) -> bool:
    normalized = value.lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n"}:
        return False
    raise ValueError(f"Invalid motif_present value: {value}")


def _smoothed_frequency(
    successes: int,
    total: int,
    *,
    pseudocount: float,
) -> float:
    return (successes + pseudocount) / (total + 2 * pseudocount)


def _log2_enrichment(
    cluster_frequency: float,
    background_frequency: float,
) -> float:
    if cluster_frequency <= 0.0:
        return -math.inf
    if background_frequency <= 0.0:
        return math.inf
    return math.log2(cluster_frequency / background_frequency)
