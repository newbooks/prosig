"""GO similarity based accession clustering utilities."""

from __future__ import annotations

import json
import logging
import math
import statistics
import time
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from prosig.go.build import LOGGER_NAME
from prosig.go.similarity import (
    BoundedProfilePairCache,
    CacheStats,
    FastGoSimilarityIndex,
    ProfilePairCache,
    build_fast_go_similarity_index,
    build_lin_similarity_matrix,
    load_accession_mf_go_terms,
    set_lin_amb_fast_for_valid_profiles,
    valid_go_profile,
)

PROGRESS_LOG_INTERVAL_SECONDS = 60.0
DEFAULT_TERM_CACHE_SIZE_MB = 256
DEFAULT_PROFILE_CACHE_SIZE_MB = 128
DEFAULT_MIN_INFORMATIVE_IC = 0.5
DEFAULT_MIN_SIMILARITY = 0.5
DEFAULT_MIN_CLUSTER_SIMILARITY = 0.25
DEFAULT_MAX_POSTING_FRACTION = 0.05
DEFAULT_NEIGHBORS = 10
DEFAULT_RESOLUTION = 2.0
BYTES_PER_MEGABYTE = 1_048_576
LEIDEN_SEED = 0


@dataclass(frozen=True)
class CandidateIndex:
    """Filtered candidate index and construction diagnostics."""

    postings_by_term: dict[str, list[int]]
    terms_by_accession_index: list[tuple[str, ...]]
    informative_terms_before_filtering: int
    informative_terms_after_filtering: int
    posting_cap: int
    fallback_accessions_after_filtering: int


@dataclass(frozen=True)
class GoClusteringResult:
    """Written cluster outputs and summary counts."""

    output_file: Path
    stats_file: Path | None
    meta_file: Path | None
    input_accessions: int
    clustered_accessions: int
    excluded_accessions: int
    edges: int
    clusters: int


@dataclass(frozen=True)
class CompleteLinkageRefinementResult:
    """Written complete-linkage outputs and refinement summary counts."""

    output_file: Path
    stats_file: Path | None
    meta_file: Path
    clustered_accessions: int
    leiden_clusters: int
    refined_clusters: int
    leiden_singletons: int
    refined_singletons: int
    leiden_clusters_split: int
    refinement_pairs_scored: int


@dataclass(frozen=True)
class _RefinedCluster:
    members: tuple[str, ...]
    sim_ave: float | None
    sim_min: float | None
    sim_max: float | None


@dataclass(frozen=True)
class _GoCandidateScore:
    go_id: str
    support: float
    score: float


@dataclass(frozen=True)
class GoClusteringConfig:
    """User-editable GO clustering tuning values."""

    stats_file: str = "leiden_clusters_stats.json"
    meta_file: str = "leiden_clusters_meta.tsv"
    neighbors: int = DEFAULT_NEIGHBORS
    resolution: float = DEFAULT_RESOLUTION
    progress_interval_seconds: float = PROGRESS_LOG_INTERVAL_SECONDS
    term_cache_size_mb: int = DEFAULT_TERM_CACHE_SIZE_MB
    profile_cache_size_mb: int = DEFAULT_PROFILE_CACHE_SIZE_MB
    min_informative_ic: float = DEFAULT_MIN_INFORMATIVE_IC
    min_similarity: float = DEFAULT_MIN_SIMILARITY
    max_posting_fraction: float = DEFAULT_MAX_POSTING_FRACTION
    max_posting_size: int = 0


def cluster_accessions_by_go(
    accession_go_file: str | Path,
    *,
    go_artifact: dict[str, Any] | None = None,
    go_graph_file: str | Path = "go_graph.pkl",
    output_file: str | Path = "leiden_clusters.tsv",
    stats_file: str | Path | None = "leiden_clusters_stats.json",
    meta_file: str | Path | None = "leiden_clusters_meta.tsv",
    resolution: float = DEFAULT_RESOLUTION,
    neighbors: int = DEFAULT_NEIGHBORS,
    term_cache_size_mb: int = DEFAULT_TERM_CACHE_SIZE_MB,
    profile_cache_size_mb: int = DEFAULT_PROFILE_CACHE_SIZE_MB,
    min_informative_ic: float = DEFAULT_MIN_INFORMATIVE_IC,
    min_similarity: float = DEFAULT_MIN_SIMILARITY,
    max_posting_fraction: float = DEFAULT_MAX_POSTING_FRACTION,
    max_posting_size: int = 0,
    progress_interval_seconds: float = PROGRESS_LOG_INTERVAL_SECONDS,
) -> GoClusteringResult:
    """Cluster accessions with Leiden over a sparse GO-similarity kNN graph."""
    _validate_parameters(
        resolution=resolution,
        neighbors=neighbors,
        term_cache_size_mb=term_cache_size_mb,
        profile_cache_size_mb=profile_cache_size_mb,
        min_informative_ic=min_informative_ic,
        min_similarity=min_similarity,
        max_posting_fraction=max_posting_fraction,
        max_posting_size=max_posting_size,
        progress_interval_seconds=progress_interval_seconds,
    )
    logger = logging.getLogger(LOGGER_NAME)
    accession_go_path = Path(accession_go_file)
    go_graph_path = Path(go_graph_file)
    output_path = Path(output_file)
    stats_path = Path(stats_file) if stats_file is not None else None
    meta_path = Path(meta_file) if meta_file is not None else None
    artifact = (
        go_artifact
        if go_artifact is not None
        else _load_go_artifact(go_graph_file)
    )
    go_index = build_fast_go_similarity_index(artifact)
    accession_terms_raw = load_accession_mf_go_terms(accession_go_path)
    input_accessions = len(accession_terms_raw)
    accession_terms = _valid_accession_terms(go_index, accession_terms_raw)
    accessions = sorted(accession_terms)
    logger.info(
        "Loaded %s accessions for GO clustering; %s retained after GO/IC filtering",
        f"{input_accessions:,}",
        f"{len(accessions):,}",
    )

    if term_cache_size_mb:
        logger.info(
            "Ignoring GO term-pair cache size because clustering uses the "
            "precomputed Lin matrix"
        )
    lin_similarity_matrix = build_lin_similarity_matrix(
        go_index,
        logger=logger,
        progress_interval_seconds=progress_interval_seconds,
    )
    profile_pair_cache = _make_profile_cache(profile_cache_size_mb, logger)
    candidate_index = build_candidate_index(
        go_index=go_index,
        accessions=accessions,
        accession_terms=accession_terms,
        min_informative_ic=min_informative_ic,
        max_posting_fraction=max_posting_fraction,
        max_posting_size=max_posting_size,
        progress_interval_seconds=progress_interval_seconds,
    )
    logger.info(
        "Filtered GO candidate index: informative terms %s -> %s; "
        "posting cap=%s; fallback accessions=%s",
        f"{candidate_index.informative_terms_before_filtering:,}",
        f"{candidate_index.informative_terms_after_filtering:,}",
        f"{candidate_index.posting_cap:,}",
        f"{candidate_index.fallback_accessions_after_filtering:,}",
    )
    edges = knn_edges_from_go_similarity(
        go_index=go_index,
        accessions=accessions,
        accession_terms=accession_terms,
        candidate_index=candidate_index,
        neighbors=neighbors,
        min_similarity=min_similarity,
        profile_pair_cache=profile_pair_cache,
        lin_similarity_matrix=lin_similarity_matrix,
        progress_interval_seconds=progress_interval_seconds,
    )
    active_indices = _active_accession_indices(edges)
    active_accessions = [accessions[index] for index in active_indices]
    active_edges = _remap_edges_to_active_indices(edges, active_indices)
    cluster_by_accession = _run_leiden(
        active_accessions,
        active_edges,
        resolution=resolution,
    )
    _write_cluster_tsv(output_path, active_accessions, cluster_by_accession)
    if meta_path is not None:
        _write_cluster_meta_tsv(
            meta_path,
            go_index=go_index,
            accession_terms=accession_terms,
            active_accessions=active_accessions,
            cluster_by_accession=cluster_by_accession,
            profile_pair_cache=profile_pair_cache,
            lin_similarity_matrix=lin_similarity_matrix,
        )
        logger.info("Saved GO cluster metadata to %s", meta_path)
    excluded_accession_count = input_accessions - len(active_accessions)
    cluster_count = len(set(cluster_by_accession.values()))
    if excluded_accession_count:
        logger.info(
            "Excluded %s accessions during GO/IC filtering or because they had "
            "no GO-similarity edges meeting the minimum similarity",
            f"{excluded_accession_count:,}",
        )
    logger.info(
        "Clustered %s accessions into %s GO clusters with k=%s, "
        "minimum similarity=%.3f, and resolution=%.3f",
        f"{len(active_accessions):,}",
        f"{cluster_count:,}",
        f"{neighbors:,}",
        min_similarity,
        resolution,
    )
    logger.info(
        "Leiden singleton clusters: %s",
        f"{_cluster_singleton_count(cluster_by_accession):,}",
    )
    _log_cluster_size_summary(logger, cluster_by_accession)
    if stats_path is not None:
        _write_stats_json(
            stats_path,
            algorithm="go_set_similarity_knn_leiden",
            similarity="lin_amb",
            partition="RBConfigurationVertexPartition",
            resolution=resolution,
            neighbors=neighbors,
            min_informative_ic=min_informative_ic,
            min_similarity=min_similarity,
            max_posting_fraction=max_posting_fraction,
            max_posting_size=max_posting_size,
            input_accessions=input_accessions,
            cleaned_accessions=len(accessions),
            clustered_accessions=len(active_accessions),
            excluded_accessions=excluded_accession_count,
            edges=len(edges),
            clusters=cluster_count,
            cluster_by_accession=cluster_by_accession,
            candidate_index=candidate_index,
            profile_cache_stats=_cache_stats(
                profile_pair_cache,
                budget_mb=profile_cache_size_mb,
            ),
            accession_go_file=accession_go_path,
            go_graph_file=go_graph_path,
            output_file=output_path,
            meta_file=meta_path,
        )
        logger.info("Saved GO cluster statistics to %s", stats_path)
    logger.info("Saved GO clusters to %s", output_path)
    return GoClusteringResult(
        output_file=output_path,
        stats_file=stats_path,
        meta_file=meta_path,
        input_accessions=input_accessions,
        clustered_accessions=len(active_accessions),
        excluded_accessions=excluded_accession_count,
        edges=len(edges),
        clusters=cluster_count,
    )


def refine_go_clusters_complete_linkage(
    accession_go_file: str | Path,
    leiden_cluster_file: str | Path,
    *,
    go_artifact: dict[str, Any] | None = None,
    go_graph_file: str | Path = "go_graph.pkl",
    output_file: str | Path = "clusters.tsv",
    meta_file: str | Path = "clusters_meta.tsv",
    stats_file: str | Path | None = "clusters_stats.json",
    min_cluster_similarity: float = DEFAULT_MIN_CLUSTER_SIMILARITY,
    profile_cache_size_mb: int = DEFAULT_PROFILE_CACHE_SIZE_MB,
    progress_interval_seconds: float = PROGRESS_LOG_INTERVAL_SECONDS,
) -> CompleteLinkageRefinementResult:
    """Refine Leiden communities with an all-pairs complete-linkage threshold."""
    _validate_refinement_parameters(
        min_cluster_similarity=min_cluster_similarity,
        profile_cache_size_mb=profile_cache_size_mb,
        progress_interval_seconds=progress_interval_seconds,
    )
    logger = logging.getLogger(LOGGER_NAME)
    accession_go_path = Path(accession_go_file)
    leiden_cluster_path = Path(leiden_cluster_file)
    go_graph_path = Path(go_graph_file)
    output_path = Path(output_file)
    meta_path = Path(meta_file)
    stats_path = Path(stats_file) if stats_file is not None else None
    artifact = (
        go_artifact
        if go_artifact is not None
        else _load_go_artifact(go_graph_path)
    )
    go_index = build_fast_go_similarity_index(artifact)
    accession_terms = _valid_accession_terms(
        go_index,
        load_accession_mf_go_terms(accession_go_path),
    )
    leiden_cluster_by_accession = _load_cluster_tsv(leiden_cluster_path)
    missing_profiles = sorted(set(leiden_cluster_by_accession) - set(accession_terms))
    if missing_profiles:
        preview = ", ".join(missing_profiles[:5])
        raise ValueError(
            "Leiden cluster members are missing valid GO profiles: "
            f"{preview}"
        )

    lin_similarity_matrix = build_lin_similarity_matrix(
        go_index,
        logger=logger,
        progress_interval_seconds=progress_interval_seconds,
    )
    profile_pair_cache = _make_profile_cache(profile_cache_size_mb, logger)
    members_by_leiden_cluster = _cluster_members(
        sorted(leiden_cluster_by_accession),
        leiden_cluster_by_accession,
    )
    leiden_singletons = sum(
        len(members) == 1 for members in members_by_leiden_cluster.values()
    )
    logger.info(
        "Starting complete-linkage refinement of %s Leiden clusters with "
        "minimum cluster similarity=%.3f; Leiden singletons=%s",
        f"{len(members_by_leiden_cluster):,}",
        min_cluster_similarity,
        f"{leiden_singletons:,}",
    )

    refined_clusters: list[_RefinedCluster] = []
    split_clusters = 0
    refinement_pairs_scored = 0
    last_log_time = time.monotonic()
    for cluster_index, cluster_id in enumerate(
        sorted(members_by_leiden_cluster),
        start=1,
    ):
        members = sorted(members_by_leiden_cluster[cluster_id])
        community_refined, pair_count = _refine_one_leiden_cluster(
            go_index=go_index,
            members=members,
            accession_terms=accession_terms,
            min_cluster_similarity=min_cluster_similarity,
            profile_pair_cache=profile_pair_cache,
            lin_similarity_matrix=lin_similarity_matrix,
        )
        refined_clusters.extend(community_refined)
        refinement_pairs_scored += pair_count
        if len(community_refined) > 1:
            split_clusters += 1
        if _should_log_progress(
            last_log_time,
            interval_seconds=progress_interval_seconds,
        ):
            last_log_time = time.monotonic()
            logger.info(
                "Refined %s/%s Leiden clusters; final clusters so far=%s",
                f"{cluster_index:,}",
                f"{len(members_by_leiden_cluster):,}",
                f"{len(refined_clusters):,}",
            )

    refined_clusters.sort(
        key=lambda cluster: (
            cluster.members[0],
            len(cluster.members),
            cluster.members,
        )
    )
    final_cluster_by_accession: dict[str, str] = {}
    final_summary_by_cluster: dict[
        str, tuple[float | None, float | None, float | None]
    ] = {}
    for index, cluster in enumerate(refined_clusters, start=1):
        cluster_id = f"cluster_{index:04d}"
        for accession in cluster.members:
            final_cluster_by_accession[accession] = cluster_id
        final_summary_by_cluster[cluster_id] = (
            cluster.sim_ave,
            cluster.sim_min,
            cluster.sim_max,
        )

    active_accessions = sorted(final_cluster_by_accession)
    _write_cluster_tsv(output_path, active_accessions, final_cluster_by_accession)
    _write_cluster_meta_from_summaries(
        meta_path,
        active_accessions=active_accessions,
        cluster_by_accession=final_cluster_by_accession,
        summary_by_cluster=final_summary_by_cluster,
        go_index=go_index,
        accession_terms=accession_terms,
    )
    refined_singletons = sum(
        len(cluster.members) == 1 for cluster in refined_clusters
    )
    logger.info(
        "Complete-linkage refinement produced %s clusters from %s Leiden "
        "clusters; split=%s; final singletons=%s",
        f"{len(refined_clusters):,}",
        f"{len(members_by_leiden_cluster):,}",
        f"{split_clusters:,}",
        f"{refined_singletons:,}",
    )
    _log_cluster_size_summary(logger, final_cluster_by_accession)
    if stats_path is not None:
        _write_refinement_stats_json(
            stats_path,
            min_cluster_similarity=min_cluster_similarity,
            clustered_accessions=len(active_accessions),
            leiden_clusters=len(members_by_leiden_cluster),
            refined_clusters=len(refined_clusters),
            leiden_singletons=leiden_singletons,
            refined_singletons=refined_singletons,
            leiden_clusters_split=split_clusters,
            refinement_pairs_scored=refinement_pairs_scored,
            leiden_cluster_by_accession=leiden_cluster_by_accession,
            final_cluster_by_accession=final_cluster_by_accession,
            accession_go_file=accession_go_path,
            go_graph_file=go_graph_path,
            leiden_cluster_file=leiden_cluster_path,
            output_file=output_path,
            meta_file=meta_path,
        )
    logger.info(
        "Saved complete-linkage GO clusters to %s and metadata to %s",
        output_path,
        meta_path,
    )
    return CompleteLinkageRefinementResult(
        output_file=output_path,
        stats_file=stats_path,
        meta_file=meta_path,
        clustered_accessions=len(active_accessions),
        leiden_clusters=len(members_by_leiden_cluster),
        refined_clusters=len(refined_clusters),
        leiden_singletons=leiden_singletons,
        refined_singletons=refined_singletons,
        leiden_clusters_split=split_clusters,
        refinement_pairs_scored=refinement_pairs_scored,
    )


def parse_cluster_config(path: str | Path) -> GoClusteringConfig:
    """Parse the restricted flat cluster_config.yaml structure."""
    values: dict[str, str] = {}
    allowed_keys = {
        "stats_file",
        "meta_file",
        "neighbors",
        "resolution",
        "progress_interval_seconds",
        "term_cache_size_mb",
        "profile_cache_size_mb",
        "min_informative_ic",
        "min_similarity",
        "max_posting_fraction",
        "max_posting_size",
    }
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        line = _strip_yaml_inline_comment(raw_line).strip()
        if not line:
            continue
        if ":" not in line:
            raise ValueError(f"Invalid cluster config entry: {line}")
        key, value = [part.strip() for part in line.split(":", 1)]
        if key not in allowed_keys:
            raise ValueError(f"Unsupported cluster config key: {key}")
        if not value:
            raise ValueError(f"Cluster config key {key} must have a value")
        values[key] = value

    config = GoClusteringConfig(
        stats_file=values.get("stats_file", "leiden_clusters_stats.json"),
        meta_file=values.get("meta_file", "leiden_clusters_meta.tsv"),
        neighbors=_parse_int_config(values, "neighbors", DEFAULT_NEIGHBORS),
        resolution=_parse_float_config(values, "resolution", DEFAULT_RESOLUTION),
        progress_interval_seconds=_parse_float_config(
            values,
            "progress_interval_seconds",
            PROGRESS_LOG_INTERVAL_SECONDS,
        ),
        term_cache_size_mb=_parse_int_config(
            values,
            "term_cache_size_mb",
            DEFAULT_TERM_CACHE_SIZE_MB,
        ),
        profile_cache_size_mb=_parse_int_config(
            values,
            "profile_cache_size_mb",
            DEFAULT_PROFILE_CACHE_SIZE_MB,
        ),
        min_informative_ic=_parse_float_config(
            values,
            "min_informative_ic",
            DEFAULT_MIN_INFORMATIVE_IC,
        ),
        min_similarity=_parse_float_config(
            values,
            "min_similarity",
            DEFAULT_MIN_SIMILARITY,
        ),
        max_posting_fraction=_parse_float_config(
            values,
            "max_posting_fraction",
            DEFAULT_MAX_POSTING_FRACTION,
        ),
        max_posting_size=_parse_int_config(values, "max_posting_size", 0),
    )
    if not config.stats_file.strip():
        raise ValueError("Cluster config key stats_file must have a non-empty value")
    if not config.meta_file.strip():
        raise ValueError("Cluster config key meta_file must have a non-empty value")
    _validate_parameters(
        resolution=config.resolution,
        neighbors=config.neighbors,
        term_cache_size_mb=config.term_cache_size_mb,
        profile_cache_size_mb=config.profile_cache_size_mb,
        min_informative_ic=config.min_informative_ic,
        min_similarity=config.min_similarity,
        max_posting_fraction=config.max_posting_fraction,
        max_posting_size=config.max_posting_size,
        progress_interval_seconds=config.progress_interval_seconds,
    )
    return config


def build_candidate_index(
    *,
    go_index: FastGoSimilarityIndex,
    accessions: list[str],
    accession_terms: dict[str, tuple[str, ...]],
    min_informative_ic: float = DEFAULT_MIN_INFORMATIVE_IC,
    max_posting_fraction: float = DEFAULT_MAX_POSTING_FRACTION,
    max_posting_size: int = 0,
    progress_interval_seconds: float = PROGRESS_LOG_INTERVAL_SECONDS,
) -> CandidateIndex:
    """Return a broad-ancestor-filtered inverted candidate index."""
    logger = logging.getLogger(LOGGER_NAME)
    raw_terms_by_accession_index: list[tuple[str, ...]] = []
    raw_postings: dict[str, list[int]] = {}
    last_log_time = time.monotonic()
    for accession_index, accession in enumerate(accessions):
        informative_terms = tuple(
            sorted(
                _informative_terms_for_accession(
                    go_index,
                    accession_terms[accession],
                    min_informative_ic=min_informative_ic,
                )
            )
        )
        raw_terms_by_accession_index.append(informative_terms)
        for informative_term in informative_terms:
            raw_postings.setdefault(informative_term, []).append(accession_index)
        if _should_log_progress(
            last_log_time,
            interval_seconds=progress_interval_seconds,
        ):
            last_log_time = time.monotonic()
            logger.info(
                "Indexed GO candidate terms for %s/%s accessions; "
                "informative GO terms=%s",
                f"{accession_index + 1:,}",
                f"{len(accessions):,}",
                f"{len(raw_postings):,}",
            )

    posting_cap = _posting_cap(
        len(accessions),
        max_posting_fraction=max_posting_fraction,
        max_posting_size=max_posting_size,
    )
    kept_terms = {
        term for term, postings in raw_postings.items() if len(postings) <= posting_cap
    }
    fallback_accessions = 0
    filtered_terms_by_accession_index: list[tuple[str, ...]] = []
    for informative_terms in raw_terms_by_accession_index:
        filtered_terms = tuple(term for term in informative_terms if term in kept_terms)
        if not filtered_terms and informative_terms:
            filtered_terms = (
                min(
                    informative_terms,
                    key=lambda term: (
                        len(raw_postings[term]),
                        -go_index.ic_by_term.get(term, 0.0),
                        term,
                    ),
                ),
            )
            fallback_accessions += 1
        filtered_terms_by_accession_index.append(filtered_terms)

    filtered_postings: dict[str, list[int]] = {}
    for accession_index, informative_terms in enumerate(
        filtered_terms_by_accession_index
    ):
        for informative_term in informative_terms:
            filtered_postings.setdefault(informative_term, []).append(accession_index)

    return CandidateIndex(
        postings_by_term=filtered_postings,
        terms_by_accession_index=filtered_terms_by_accession_index,
        informative_terms_before_filtering=len(raw_postings),
        informative_terms_after_filtering=len(filtered_postings),
        posting_cap=posting_cap,
        fallback_accessions_after_filtering=fallback_accessions,
    )


def knn_edges_from_go_similarity(
    *,
    go_index: FastGoSimilarityIndex,
    accessions: list[str],
    accession_terms: dict[str, tuple[str, ...]],
    candidate_index: CandidateIndex,
    neighbors: int,
    min_similarity: float = DEFAULT_MIN_SIMILARITY,
    profile_pair_cache: ProfilePairCache | None = None,
    lin_similarity_matrix: np.ndarray | None = None,
    progress_interval_seconds: float = PROGRESS_LOG_INTERVAL_SECONDS,
) -> dict[tuple[int, int], float]:
    """Return sparse undirected kNN edges without storing all pairwise scores."""
    if min_similarity <= 0.0 or min_similarity > 1.0:
        raise ValueError("cluster minimum similarity must be in (0, 1]")
    if len(accessions) < 2:
        return {}

    logger = logging.getLogger(LOGGER_NAME)
    indices_by_profile = _accession_indices_by_go_profile(accessions, accession_terms)
    comparable_profile_cache = {
        profile: _go_profile_has_comparable_terms(go_index, profile)
        for profile in indices_by_profile
    }

    edges: dict[tuple[int, int], float] = {}
    seen_candidate = bytearray(len(accessions))
    last_log_time = time.monotonic()
    for accession_index, accession in enumerate(accessions):
        profile = accession_terms[accession]
        profile_indices = indices_by_profile[profile]
        if comparable_profile_cache[profile] and len(profile_indices) - 1 >= neighbors:
            for other_index in _first_profile_indices_by_accession(
                profile_indices,
                accession_index=accession_index,
                neighbors=neighbors,
            ):
                edge = tuple(sorted((accession_index, other_index)))
                edges[edge] = 1.0
            if _should_log_progress(
                last_log_time,
                interval_seconds=progress_interval_seconds,
            ):
                last_log_time = time.monotonic()
                logger.info(
                    "Built GO kNN candidates for %s/%s accessions; "
                    "retained %s edges",
                    f"{accession_index + 1:,}",
                    f"{len(accessions):,}",
                    f"{len(edges):,}",
                )
            continue

        touched_indices = [accession_index]
        seen_candidate[accession_index] = 1
        top_candidates: list[tuple[float, int]] = []
        for same_profile_index in profile_indices:
            if same_profile_index == accession_index:
                continue
            seen_candidate[same_profile_index] = 1
            touched_indices.append(same_profile_index)
            if comparable_profile_cache[profile]:
                _add_top_candidate(
                    top_candidates,
                    (1.0, same_profile_index),
                    neighbors=neighbors,
                    accessions=accessions,
                )

        informative_terms = _candidate_index_terms_for_accession(
            candidate_index,
            accession_index,
        )
        candidate_indices = _candidate_indices_for_informative_terms(
            informative_terms=informative_terms,
            postings_by_term=candidate_index.postings_by_term,
            seen_candidate=seen_candidate,
            touched_indices=touched_indices,
        )
        for other_index in candidate_indices:
            similarity = set_lin_amb_fast_for_valid_profiles(
                go_index,
                accession_terms[accession],
                accession_terms[accessions[other_index]],
                profile_pair_cache=profile_pair_cache,
                lin_similarity_matrix=lin_similarity_matrix,
            )
            if similarity is not None and similarity >= min_similarity:
                _add_top_candidate(
                    top_candidates,
                    (similarity, other_index),
                    neighbors=neighbors,
                    accessions=accessions,
                )

        for touched_index in touched_indices:
            seen_candidate[touched_index] = 0

        top_candidates.sort(key=lambda item: (-item[0], accessions[item[1]]))
        for similarity, top_candidate_index in top_candidates:
            edge = tuple(sorted((accession_index, top_candidate_index)))
            edges[edge] = max(edges.get(edge, 0.0), similarity)

        if _should_log_progress(
            last_log_time,
            interval_seconds=progress_interval_seconds,
        ):
            last_log_time = time.monotonic()
            logger.info(
                "Built GO kNN candidates for %s/%s accessions; retained %s edges",
                f"{accession_index + 1:,}",
                f"{len(accessions):,}",
                f"{len(edges):,}",
            )

    logger.info(
        "Built GO kNN candidates for %s accessions; retained %s edges",
        f"{len(accessions):,}",
        f"{len(edges):,}",
    )
    return edges


def _candidate_index_terms_for_accession(
    candidate_index: CandidateIndex,
    accession_index: int,
) -> tuple[str, ...]:
    return candidate_index.terms_by_accession_index[accession_index]


def _validate_parameters(
    *,
    resolution: float,
    neighbors: int,
    term_cache_size_mb: int,
    profile_cache_size_mb: int,
    min_informative_ic: float,
    min_similarity: float,
    max_posting_fraction: float,
    max_posting_size: int,
    progress_interval_seconds: float,
) -> None:
    if resolution <= 0.0:
        raise ValueError("cluster resolution must be greater than 0")
    if neighbors < 1:
        raise ValueError("cluster neighbors must be at least 1")
    if term_cache_size_mb < 0:
        raise ValueError("cluster term cache size must be non-negative")
    if profile_cache_size_mb < 0:
        raise ValueError("cluster profile cache size must be non-negative")
    if min_informative_ic < 0.0:
        raise ValueError("cluster minimum informative IC must be non-negative")
    if min_similarity <= 0.0 or min_similarity > 1.0:
        raise ValueError("cluster minimum similarity must be in (0, 1]")
    if max_posting_fraction <= 0.0 or max_posting_fraction > 1.0:
        raise ValueError("cluster max posting fraction must be in (0, 1]")
    if max_posting_size < 0:
        raise ValueError("cluster max posting size must be non-negative")
    if progress_interval_seconds <= 0.0:
        raise ValueError("cluster progress interval must be greater than 0")


def _validate_refinement_parameters(
    *,
    min_cluster_similarity: float,
    profile_cache_size_mb: int,
    progress_interval_seconds: float,
) -> None:
    if min_cluster_similarity <= 0.0 or min_cluster_similarity > 1.0:
        raise ValueError("minimum cluster similarity must be in (0, 1]")
    if profile_cache_size_mb < 0:
        raise ValueError("cluster profile cache size must be non-negative")
    if progress_interval_seconds <= 0.0:
        raise ValueError("cluster progress interval must be greater than 0")


def _parse_int_config(
    values: dict[str, str],
    key: str,
    default: int,
) -> int:
    try:
        return int(values.get(key, str(default)))
    except ValueError as exc:
        raise ValueError(f"Cluster config key {key} must be an integer") from exc


def _parse_float_config(
    values: dict[str, str],
    key: str,
    default: float,
) -> float:
    try:
        return float(values.get(key, str(default)))
    except ValueError as exc:
        raise ValueError(f"Cluster config key {key} must be a number") from exc


def _strip_yaml_inline_comment(value: str) -> str:
    in_single = False
    in_double = False
    result: list[str] = []
    for char in value:
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        elif char == "#" and not in_single and not in_double:
            break
        result.append(char)
    return "".join(result)


def _load_go_artifact(path: str | Path) -> dict[str, Any]:
    import pickle

    with Path(path).open("rb") as handle:
        artifact = pickle.load(handle)
    if not isinstance(artifact, dict):
        raise ValueError(f"GO graph pickle must contain a dictionary artifact: {path}")
    return artifact


def _load_cluster_tsv(path: str | Path) -> dict[str, str]:
    cluster_by_accession: dict[str, str] = {}
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "member_id\tcluster_id":
        raise ValueError(
            "Cluster TSV must start with header: member_id\\tcluster_id"
        )
    for line_number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 2 or not all(fields):
            raise ValueError(f"Invalid cluster TSV row {line_number}: {line}")
        accession, cluster_id = fields
        if accession in cluster_by_accession:
            raise ValueError(f"Duplicate cluster member: {accession}")
        cluster_by_accession[accession] = cluster_id
    return cluster_by_accession


def _valid_accession_terms(
    go_index: FastGoSimilarityIndex,
    accession_terms: dict[str, tuple[str, ...]],
) -> dict[str, tuple[str, ...]]:
    valid_terms_by_accession = {
        accession: valid_go_profile(go_index, terms)
        for accession, terms in accession_terms.items()
    }
    return {
        accession: terms
        for accession, terms in valid_terms_by_accession.items()
        if terms
    }


def _make_profile_cache(
    profile_cache_size_mb: int,
    logger: logging.Logger,
) -> ProfilePairCache | None:
    if profile_cache_size_mb == 0:
        return None
    cache = BoundedProfilePairCache(profile_cache_size_mb * BYTES_PER_MEGABYTE)
    logger.info(
        "Using GO profile-pair similarity cache: budget=%s MB, entries~%s",
        f"{profile_cache_size_mb:,}",
        f"{cache.max_entries:,}",
    )
    return cache


def _cache_stats(
    cache: ProfilePairCache | None,
    *,
    budget_mb: int,
) -> CacheStats:
    if isinstance(cache, BoundedProfilePairCache):
        return cache.stats(budget_mb=budget_mb)
    return CacheStats(
        budget_mb=budget_mb,
        max_entries=0,
        entries=len(cache) if isinstance(cache, dict) else 0,
        hits=0,
        misses=0,
        evictions=0,
    )


def _informative_terms_for_accession(
    go_index: FastGoSimilarityIndex,
    terms: tuple[str, ...],
    *,
    min_informative_ic: float,
) -> set[str]:
    informative_terms: set[str] = set()
    for term in terms:
        for candidate_term in go_index.ancestors_by_term.get(term, frozenset()):
            if go_index.ic_by_term.get(candidate_term, 0.0) >= min_informative_ic:
                informative_terms.add(candidate_term)
    return informative_terms


def _posting_cap(
    accession_count: int,
    *,
    max_posting_fraction: float,
    max_posting_size: int,
) -> int:
    if max_posting_size > 0:
        return max_posting_size
    return max(1, math.ceil(accession_count * max_posting_fraction))


def _candidate_indices_for_informative_terms(
    *,
    informative_terms: tuple[str, ...],
    postings_by_term: dict[str, list[int]],
    seen_candidate: bytearray,
    touched_indices: list[int],
) -> Iterable[int]:
    for informative_term in informative_terms:
        for candidate_index in postings_by_term.get(informative_term, ()):
            if seen_candidate[candidate_index]:
                continue
            seen_candidate[candidate_index] = 1
            touched_indices.append(candidate_index)
            yield candidate_index


def _candidate_is_better(
    candidate: tuple[float, int],
    current: tuple[float, int],
    accessions: list[str],
) -> bool:
    candidate_similarity, candidate_index = candidate
    current_similarity, current_index = current
    if candidate_similarity != current_similarity:
        return candidate_similarity > current_similarity
    return accessions[candidate_index] < accessions[current_index]


def _add_top_candidate(
    top_candidates: list[tuple[float, int]],
    candidate: tuple[float, int],
    *,
    neighbors: int,
    accessions: list[str],
) -> None:
    if len(top_candidates) < neighbors:
        top_candidates.append(candidate)
        return

    worst_position = 0
    for position, current in enumerate(top_candidates[1:], start=1):
        if _candidate_is_better(top_candidates[worst_position], current, accessions):
            worst_position = position

    if _candidate_is_better(candidate, top_candidates[worst_position], accessions):
        top_candidates[worst_position] = candidate


def _go_profile_has_comparable_terms(
    go_index: FastGoSimilarityIndex,
    profile: tuple[str, ...],
) -> bool:
    return any(go_index.ic_by_term.get(term, 0.0) > 0.0 for term in profile)


def _accession_indices_by_go_profile(
    accessions: list[str],
    accession_terms: dict[str, tuple[str, ...]],
) -> dict[tuple[str, ...], list[int]]:
    indices_by_profile: dict[tuple[str, ...], list[int]] = {}
    for accession_index, accession in enumerate(accessions):
        indices_by_profile.setdefault(accession_terms[accession], []).append(
            accession_index
        )
    for profile_indices in indices_by_profile.values():
        profile_indices.sort(key=lambda accession_index: accessions[accession_index])
    return indices_by_profile


def _first_profile_indices_by_accession(
    profile_indices: list[int],
    *,
    accession_index: int,
    neighbors: int,
) -> list[int]:
    selected_indices: list[int] = []
    for candidate_index in profile_indices:
        if candidate_index == accession_index:
            continue
        selected_indices.append(candidate_index)
        if len(selected_indices) == neighbors:
            break
    return selected_indices


def _active_accession_indices(edges: dict[tuple[int, int], float]) -> list[int]:
    active_indices: set[int] = set()
    for source_index, target_index in edges:
        active_indices.add(source_index)
        active_indices.add(target_index)
    return sorted(active_indices)


def _remap_edges_to_active_indices(
    edges: dict[tuple[int, int], float],
    active_indices: list[int],
) -> dict[tuple[int, int], float]:
    active_position_by_index = {
        accession_index: active_position
        for active_position, accession_index in enumerate(active_indices)
    }
    return {
        (
            active_position_by_index[source_index],
            active_position_by_index[target_index],
        ): weight
        for (source_index, target_index), weight in edges.items()
    }


def _run_leiden(
    active_accessions: list[str],
    active_edges: dict[tuple[int, int], float],
    *,
    resolution: float,
) -> dict[str, str]:
    if not active_accessions:
        return {}
    try:
        import igraph as ig
        import leidenalg
    except ImportError as exc:
        raise ValueError(
            "GO clustering requires igraph and leidenalg. "
            "Install with: pip install -e .[cluster]"
        ) from exc

    graph = ig.Graph(n=len(active_accessions), edges=list(active_edges))
    graph.vs["name"] = active_accessions
    graph.es["weight"] = list(active_edges.values())
    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=LEIDEN_SEED,
    )
    return _cluster_ids_from_membership(active_accessions, list(partition.membership))


def _cluster_ids_from_membership(
    accessions: list[str],
    membership: list[int],
) -> dict[str, str]:
    communities: dict[int, list[str]] = {}
    for accession, community_id in zip(accessions, membership, strict=True):
        communities.setdefault(community_id, []).append(accession)

    sorted_communities = sorted(
        (sorted(members) for members in communities.values()),
        key=lambda members: (members[0], len(members)),
    )

    cluster_by_accession: dict[str, str] = {}
    for index, members in enumerate(sorted_communities, start=1):
        cluster_id = f"cluster_{index:04d}"
        for accession in members:
            cluster_by_accession[accession] = cluster_id
    return cluster_by_accession


def _write_cluster_tsv(
    output_path: Path,
    active_accessions: list[str],
    cluster_by_accession: dict[str, str],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("member_id\tcluster_id\n")
        for accession in active_accessions:
            handle.write(f"{accession}\t{cluster_by_accession[accession]}\n")


def _write_cluster_meta_tsv(
    meta_path: Path,
    *,
    go_index: FastGoSimilarityIndex,
    accession_terms: dict[str, tuple[str, ...]],
    active_accessions: list[str],
    cluster_by_accession: dict[str, str],
    profile_pair_cache: ProfilePairCache | None,
    lin_similarity_matrix: np.ndarray | None,
) -> None:
    members_by_cluster = _cluster_members(active_accessions, cluster_by_accession)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    with meta_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "cluster_id\tsim_ave\tsim_min\tsim_max\tsize\t"
            "composed_go\n"
        )
        for cluster_id in sorted(members_by_cluster):
            members = members_by_cluster[cluster_id]
            sim_ave, sim_min, sim_max = _pairwise_cluster_similarity_summary(
                go_index,
                members,
                accession_terms=accession_terms,
                profile_pair_cache=profile_pair_cache,
                lin_similarity_matrix=lin_similarity_matrix,
            )
            similarity_texts = [
                "NA" if similarity is None else f"{similarity:7.5f}"
                for similarity in (sim_ave, sim_min, sim_max)
            ]
            similarity_columns = "\t".join(similarity_texts)
            composed_go = _format_composed_go(
                _synthesize_cluster_go_terms(go_index, members, accession_terms)
            )
            handle.write(
                f"{cluster_id}\t"
                f"{similarity_columns}\t"
                f"{len(members)}\t"
                f"{composed_go}\n"
            )


def _write_cluster_meta_from_summaries(
    meta_path: Path,
    *,
    active_accessions: list[str],
    cluster_by_accession: dict[str, str],
    summary_by_cluster: dict[
        str, tuple[float | None, float | None, float | None]
    ],
    go_index: FastGoSimilarityIndex,
    accession_terms: dict[str, tuple[str, ...]],
) -> None:
    members_by_cluster = _cluster_members(active_accessions, cluster_by_accession)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    with meta_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "cluster_id\tsim_ave\tsim_min\tsim_max\tsize\t"
            "composed_go\n"
        )
        for cluster_id in sorted(members_by_cluster):
            similarities = summary_by_cluster[cluster_id]
            similarity_columns = "\t".join(
                "NA" if similarity is None else f"{similarity:7.5f}"
                for similarity in similarities
            )
            members = members_by_cluster[cluster_id]
            composed_go = _format_composed_go(
                _synthesize_cluster_go_terms(go_index, members, accession_terms)
            )
            handle.write(
                f"{cluster_id}\t{similarity_columns}\t"
                f"{len(members)}\t"
                f"{composed_go}\n"
            )


def _synthesize_cluster_go_terms(
    go_index: FastGoSimilarityIndex,
    members: list[str],
    accession_terms: dict[str, tuple[str, ...]],
    *,
    max_terms: int = 10,
    min_support: float = 0.1,
    relative_drop_cutoff: float = 0.5,
    parent_coverage_cutoff: float = 0.8,
) -> list[_GoCandidateScore]:
    """Return representative GO terms for a cluster using equal accession votes."""
    if max_terms < 1:
        raise ValueError("max_terms must be at least 1")
    if not 0.0 <= min_support <= 1.0:
        raise ValueError("minimum support must be between 0 and 1")
    if relative_drop_cutoff < 0.0:
        raise ValueError("relative drop cutoff must be non-negative")
    if not 0.0 <= parent_coverage_cutoff <= 1.0:
        raise ValueError("parent coverage cutoff must be between 0 and 1")

    annotated_members = 0
    term_votes: Counter[str] = Counter()
    for member in members:
        terms = accession_terms.get(member, ())
        if not terms:
            continue
        annotated_members += 1
        supported_terms = _propagated_go_terms_for_member(go_index, terms)
        term_votes.update(supported_terms)

    if annotated_members == 0:
        return []

    candidates: list[_GoCandidateScore] = []
    for go_id, votes in term_votes.items():
        support = votes / annotated_members
        if support < min_support:
            continue
        ic = go_index.ic_by_term.get(go_id)
        if ic is None or ic <= 0.0:
            continue
        candidates.append(_GoCandidateScore(go_id, support, support * ic))

    ranked = sorted(
        candidates,
        key=lambda item: (-item.score, -item.support, item.go_id),
    )
    selected = _suppress_parent_go_terms(
        ranked,
        go_index,
        parent_coverage_cutoff=parent_coverage_cutoff,
    )
    selected = _apply_relative_go_score_drop(
        selected,
        relative_drop_cutoff=relative_drop_cutoff,
    )
    return selected[:max_terms]


def _propagated_go_terms_for_member(
    go_index: FastGoSimilarityIndex,
    terms: tuple[str, ...],
) -> set[str]:
    propagated_terms: set[str] = set()
    for term in terms:
        for candidate_term in go_index.ancestors_by_term.get(term, ()):
            ic = go_index.ic_by_term.get(candidate_term)
            if ic is None or ic <= 0.0:
                continue
            propagated_terms.add(candidate_term)
    return propagated_terms


def _suppress_parent_go_terms(
    ranked_terms: list[_GoCandidateScore],
    go_index: FastGoSimilarityIndex,
    *,
    parent_coverage_cutoff: float,
) -> list[_GoCandidateScore]:
    selected: list[_GoCandidateScore] = []
    for candidate in ranked_terms:
        if _go_candidate_is_suppressed(
            candidate,
            selected,
            go_index,
            parent_coverage_cutoff=parent_coverage_cutoff,
        ):
            continue
        selected = [
            current
            for current in selected
            if not _go_candidate_suppresses_selected(
                candidate,
                current,
                go_index,
                parent_coverage_cutoff=parent_coverage_cutoff,
            )
        ]
        selected.append(candidate)
    return selected


def _go_candidate_is_suppressed(
    candidate: _GoCandidateScore,
    selected: list[_GoCandidateScore],
    go_index: FastGoSimilarityIndex,
    *,
    parent_coverage_cutoff: float,
) -> bool:
    candidate_ancestors = go_index.ancestors_by_term.get(candidate.go_id, ())
    for current in selected:
        current_ancestors = go_index.ancestors_by_term.get(current.go_id, ())
        if current.go_id in candidate_ancestors:
            coverage = (
                candidate.support / current.support
                if current.support > 0
                else 0.0
            )
            if coverage >= parent_coverage_cutoff:
                continue
            return True
        if candidate.go_id in current_ancestors:
            coverage = (
                current.support / candidate.support
                if candidate.support > 0
                else 0.0
            )
            if coverage >= parent_coverage_cutoff:
                return True
    return False


def _go_candidate_suppresses_selected(
    candidate: _GoCandidateScore,
    current: _GoCandidateScore,
    go_index: FastGoSimilarityIndex,
    *,
    parent_coverage_cutoff: float,
) -> bool:
    candidate_ancestors = go_index.ancestors_by_term.get(candidate.go_id, ())
    if current.go_id in candidate_ancestors:
        coverage = candidate.support / current.support if current.support > 0 else 0.0
        return coverage >= parent_coverage_cutoff

    current_ancestors = go_index.ancestors_by_term.get(current.go_id, ())
    if candidate.go_id in current_ancestors:
        coverage = current.support / candidate.support if candidate.support > 0 else 0.0
        return coverage < parent_coverage_cutoff

    return False


def _apply_relative_go_score_drop(
    ranked_terms: list[_GoCandidateScore],
    *,
    relative_drop_cutoff: float,
) -> list[_GoCandidateScore]:
    if len(ranked_terms) < 2:
        return ranked_terms
    selected = [ranked_terms[0]]
    previous_score = ranked_terms[0].score
    for candidate in ranked_terms[1:]:
        if (
            previous_score > 0.0
            and candidate.score / previous_score < relative_drop_cutoff
        ):
            break
        selected.append(candidate)
        previous_score = candidate.score
    return selected


def _format_composed_go(terms: list[_GoCandidateScore]) -> str:
    return ";".join(term.go_id for term in terms)


def _refine_one_leiden_cluster(
    *,
    go_index: FastGoSimilarityIndex,
    members: list[str],
    accession_terms: dict[str, tuple[str, ...]],
    min_cluster_similarity: float,
    profile_pair_cache: ProfilePairCache | None,
    lin_similarity_matrix: np.ndarray | None,
) -> tuple[list[_RefinedCluster], int]:
    if len(members) == 1:
        return [_RefinedCluster((members[0],), None, None, None)], 0

    try:
        from scipy.cluster.hierarchy import fcluster, linkage
    except ImportError as exc:
        raise ValueError(
            "Complete-linkage GO cluster refinement requires scipy. "
            "Install with: pip install scipy"
        ) from exc

    pair_count = len(members) * (len(members) - 1) // 2
    distances = np.empty(pair_count, dtype=np.float64)
    offset = 0
    for member_index, accession_a in enumerate(members[:-1]):
        profile_a = accession_terms[accession_a]
        for accession_b in members[member_index + 1 :]:
            score = set_lin_amb_fast_for_valid_profiles(
                go_index,
                profile_a,
                accession_terms[accession_b],
                profile_pair_cache=profile_pair_cache,
                lin_similarity_matrix=lin_similarity_matrix,
            )
            similarity = 0.0 if score is None else min(1.0, max(0.0, score))
            distances[offset] = 1.0 - similarity
            offset += 1

    hierarchy = linkage(distances, method="complete")
    labels = fcluster(
        hierarchy,
        t=1.0 - min_cluster_similarity,
        criterion="distance",
    )
    member_indices_by_label: dict[int, list[int]] = {}
    for member_index, label in enumerate(labels):
        member_indices_by_label.setdefault(int(label), []).append(member_index)

    refined_clusters: list[_RefinedCluster] = []
    for member_indices in member_indices_by_label.values():
        cluster_members = tuple(members[index] for index in member_indices)
        sim_ave, sim_min, sim_max = _condensed_similarity_summary(
            distances,
            member_count=len(members),
            member_indices=member_indices,
        )
        if sim_min is not None and sim_min + 1e-12 < min_cluster_similarity:
            raise RuntimeError(
                "Complete-linkage refinement violated minimum cluster "
                f"similarity: {sim_min} < {min_cluster_similarity}"
            )
        refined_clusters.append(
            _RefinedCluster(cluster_members, sim_ave, sim_min, sim_max)
        )
    return refined_clusters, pair_count


def _condensed_similarity_summary(
    distances: np.ndarray,
    *,
    member_count: int,
    member_indices: list[int],
) -> tuple[float | None, float | None, float | None]:
    if len(member_indices) < 2:
        return None, None, None
    total = 0.0
    pair_count = 0
    sim_min = math.inf
    sim_max = -math.inf
    for position, member_a in enumerate(member_indices[:-1]):
        for member_b in member_indices[position + 1 :]:
            distance_index = _condensed_distance_index(
                member_count,
                member_a,
                member_b,
            )
            similarity = 1.0 - float(distances[distance_index])
            total += similarity
            pair_count += 1
            sim_min = min(sim_min, similarity)
            sim_max = max(sim_max, similarity)
    return total / pair_count, sim_min, sim_max


def _condensed_distance_index(member_count: int, member_a: int, member_b: int) -> int:
    if member_a == member_b:
        raise ValueError("Condensed distance index requires distinct members")
    if member_a > member_b:
        member_a, member_b = member_b, member_a
    return (
        member_count * member_a
        - member_a * (member_a + 1) // 2
        + member_b
        - member_a
        - 1
    )


def _cluster_members(
    active_accessions: list[str],
    cluster_by_accession: dict[str, str],
) -> dict[str, list[str]]:
    members_by_cluster: dict[str, list[str]] = {}
    for accession in active_accessions:
        cluster_id = cluster_by_accession[accession]
        members_by_cluster.setdefault(cluster_id, []).append(accession)
    return members_by_cluster


def _pairwise_cluster_similarity_summary(
    go_index: FastGoSimilarityIndex,
    members: list[str],
    *,
    accession_terms: dict[str, tuple[str, ...]],
    profile_pair_cache: ProfilePairCache | None,
    lin_similarity_matrix: np.ndarray | None,
) -> tuple[float | None, float | None, float | None]:
    if len(members) < 2:
        return None, None, None

    total = 0.0
    pair_count = 0
    sim_min = math.inf
    sim_max = -math.inf
    for index, accession_a in enumerate(members[:-1]):
        profile_a = accession_terms[accession_a]
        for accession_b in members[index + 1 :]:
            score = set_lin_amb_fast_for_valid_profiles(
                go_index,
                profile_a,
                accession_terms[accession_b],
                profile_pair_cache=profile_pair_cache,
                lin_similarity_matrix=lin_similarity_matrix,
            )
            similarity = 0.0 if score is None else score
            total += similarity
            pair_count += 1
            sim_min = min(sim_min, similarity)
            sim_max = max(sim_max, similarity)
    return total / pair_count, sim_min, sim_max


def _write_stats_json(
    stats_path: Path,
    *,
    algorithm: str,
    similarity: str,
    partition: str,
    resolution: float,
    neighbors: int,
    min_informative_ic: float,
    min_similarity: float,
    max_posting_fraction: float,
    max_posting_size: int,
    input_accessions: int,
    cleaned_accessions: int,
    clustered_accessions: int,
    excluded_accessions: int,
    edges: int,
    clusters: int,
    cluster_by_accession: dict[str, str],
    candidate_index: CandidateIndex,
    profile_cache_stats: CacheStats,
    accession_go_file: Path,
    go_graph_file: Path,
    output_file: Path,
    meta_file: Path | None,
) -> None:
    min_size, mean_size, median_size, max_size = _cluster_size_summary(
        cluster_by_accession
    )
    stats = {
        "algorithm": algorithm,
        "similarity": similarity,
        "partition": partition,
        "resolution": resolution,
        "neighbors": neighbors,
        "seed": LEIDEN_SEED,
        "min_informative_ic": min_informative_ic,
        "min_similarity": min_similarity,
        "max_posting_fraction": max_posting_fraction,
        "max_posting_size": max_posting_size,
        "input_accessions": input_accessions,
        "cleaned_accessions": cleaned_accessions,
        "clustered_accessions": clustered_accessions,
        "excluded_accessions": excluded_accessions,
        "informative_terms_before_filtering": (
            candidate_index.informative_terms_before_filtering
        ),
        "informative_terms_after_filtering": (
            candidate_index.informative_terms_after_filtering
        ),
        "posting_cap": candidate_index.posting_cap,
        "fallback_accessions_after_filtering": (
            candidate_index.fallback_accessions_after_filtering
        ),
        "edges": edges,
        "clusters": clusters,
        "singletons": _cluster_singleton_count(cluster_by_accession),
        "cluster_size_min": min_size,
        "cluster_size_mean": mean_size,
        "cluster_size_median": median_size,
        "cluster_size_max": max_size,
        "lin_matrix": {"dtype": "float32", "storage": "memory"},
        "profile_cache": asdict(profile_cache_stats),
        "dependencies": {
            "go_graph": str(go_graph_file),
            "accession_go": str(accession_go_file),
        },
        "outputs": {
            "clusters": str(output_file),
            "stats": str(stats_path),
            "meta": str(meta_file) if meta_file is not None else None,
        },
    }
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")


def _write_refinement_stats_json(
    stats_path: Path,
    *,
    min_cluster_similarity: float,
    clustered_accessions: int,
    leiden_clusters: int,
    refined_clusters: int,
    leiden_singletons: int,
    refined_singletons: int,
    leiden_clusters_split: int,
    refinement_pairs_scored: int,
    leiden_cluster_by_accession: dict[str, str],
    final_cluster_by_accession: dict[str, str],
    accession_go_file: Path,
    go_graph_file: Path,
    leiden_cluster_file: Path,
    output_file: Path,
    meta_file: Path,
) -> None:
    pre_min, pre_mean, pre_median, pre_max = _cluster_size_summary(
        leiden_cluster_by_accession
    )
    post_min, post_mean, post_median, post_max = _cluster_size_summary(
        final_cluster_by_accession
    )
    stats = {
        "algorithm": "go_set_similarity_leiden_complete_linkage",
        "similarity": "lin_amb",
        "min_cluster_similarity": min_cluster_similarity,
        "clustered_accessions": clustered_accessions,
        "leiden_clusters": leiden_clusters,
        "refined_clusters": refined_clusters,
        "clusters": refined_clusters,
        "leiden_singletons": leiden_singletons,
        "refined_singletons": refined_singletons,
        "leiden_clusters_split": leiden_clusters_split,
        "refinement_pairs_scored": refinement_pairs_scored,
        "pre_refinement_cluster_size_min": pre_min,
        "pre_refinement_cluster_size_mean": pre_mean,
        "pre_refinement_cluster_size_median": pre_median,
        "pre_refinement_cluster_size_max": pre_max,
        "cluster_size_min": post_min,
        "cluster_size_mean": post_mean,
        "cluster_size_median": post_median,
        "cluster_size_max": post_max,
        "dependencies": {
            "go_graph": str(go_graph_file),
            "accession_go": str(accession_go_file),
            "leiden_clusters": str(leiden_cluster_file),
        },
        "outputs": {
            "clusters": str(output_file),
            "stats": str(stats_path),
            "meta": str(meta_file),
        },
    }
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")


def _cluster_singleton_count(cluster_by_accession: dict[str, str]) -> int:
    cluster_sizes: dict[str, int] = {}
    for cluster_id in cluster_by_accession.values():
        cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
    return sum(size == 1 for size in cluster_sizes.values())


def _cluster_size_summary(
    cluster_by_accession: dict[str, str],
) -> tuple[int, float, float, int]:
    cluster_sizes_by_id: dict[str, int] = {}
    for cluster_id in cluster_by_accession.values():
        cluster_sizes_by_id[cluster_id] = cluster_sizes_by_id.get(cluster_id, 0) + 1
    cluster_sizes = list(cluster_sizes_by_id.values())
    if not cluster_sizes:
        return 0, 0.0, 0.0, 0
    return (
        min(cluster_sizes),
        statistics.fmean(cluster_sizes),
        statistics.median(cluster_sizes),
        max(cluster_sizes),
    )


def _log_cluster_size_summary(
    logger: logging.Logger,
    cluster_by_accession: dict[str, str],
) -> None:
    if not cluster_by_accession:
        logger.info("GO cluster sizes: no clusters reported")
        return
    min_size, average_size, median_size, max_size = _cluster_size_summary(
        cluster_by_accession
    )
    logger.info(
        "GO cluster sizes: min=%s, average=%.2f, median=%.2f, max=%s",
        f"{min_size:,}",
        average_size,
        median_size,
        f"{max_size:,}",
    )


def _should_log_progress(last_log_time: float, *, interval_seconds: float) -> bool:
    return time.monotonic() - last_log_time >= interval_seconds


def _log_edge_progress(
    logger: logging.Logger,
    *,
    accession_index: int,
    accession_count: int,
    edge_count: int,
    last_log_time_ref: list[float],
    progress_interval_seconds: float,
) -> None:
    if _should_log_progress(
        last_log_time_ref[0],
        interval_seconds=progress_interval_seconds,
    ):
        last_log_time_ref[0] = time.monotonic()
        logger.info(
            "Built GO kNN candidates for %s/%s accessions; retained %s edges",
            f"{accession_index + 1:,}",
            f"{accession_count:,}",
            f"{edge_count:,}",
        )
