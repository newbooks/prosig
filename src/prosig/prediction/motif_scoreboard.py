"""Motif-cluster score board construction for function prediction."""

from __future__ import annotations

import csv
import json
import math
import pickle
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


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
            cluster_frequency = tp / cluster_size
            background_frequency = fp / outside_count
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
    )

    artifact = {
        "schema_version": "1.0",
        "kind": "motif_cluster_scoreboard",
        "parameters": {
            "min_cluster_size": min_cluster_size,
            "min_support": min_support,
            "weight": "log2((TP / (TP + FN)) / (FP / (FP + TN)))",
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


def _log2_enrichment(
    cluster_frequency: float,
    background_frequency: float,
) -> float:
    if cluster_frequency <= 0.0:
        return -math.inf
    if background_frequency <= 0.0:
        return math.inf
    return math.log2(cluster_frequency / background_frequency)
