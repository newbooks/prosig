#!/usr/bin/env python3
"""Select a reduced, sequence-diverse cluster panel for simulations."""

from __future__ import annotations

import argparse
import csv
import json
import pickle
import random
import sys
import time
from concurrent.futures import Executor, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if REPO_SRC.is_dir():
    sys.path.insert(0, str(REPO_SRC))

from prosig.go.describe import describe_go_function  # noqa: E402
from prosig.go.similarity import (  # noqa: E402
    FastGoSimilarityIndex,
    build_fast_go_similarity_index,
    set_lin_amb_fast_for_valid_profiles,
    valid_go_profile,
)
from prosig.library import resolve_core_library  # noqa: E402
from prosig.sequence_identity import sequence_identity  # noqa: E402
from prosig.sequences import indexed_fasta_sequence  # noqa: E402


@dataclass(frozen=True)
class ClusterCandidate:
    cluster_id: str
    members: tuple[str, ...]
    selected_accessions: tuple[str, ...]
    composed_go: tuple[str, ...]
    description: str
    within_mean_identity: float
    within_max_identity: float
    within_min_identity: float
    within_identity_matrix: tuple[tuple[float, ...], ...]
    within_go_similarity_matrix: tuple[tuple[float, ...], ...]


def main() -> None:
    args = _parse_args()
    rng = random.Random(args.seed)
    started_at = time.monotonic()

    def progress(message: str) -> None:
        elapsed = time.monotonic() - started_at
        print(f"[{elapsed:8.1f}s] {message}", file=sys.stderr, flush=True)

    with resolve_core_library() as library:
        clusters_path = _resolve_path(
            None,
            default_candidates=(
                Path("work/clusters.tsv"),
                Path("clusters.tsv"),
            ),
        )
        clusters_meta_path = _resolve_path(
            args.clusters_meta,
            default_candidates=(
                clusters_path.with_name("clusters_meta.tsv"),
                library.path("clusters_meta.tsv"),
            ),
        )
        accession_go_path = _resolve_path(
            args.accession_go,
            default_candidates=(
                clusters_path.with_name("accession_mf_go.tsv"),
                library.path("accession_mf_go.tsv"),
            ),
        )
        fasta_path = _resolve_path(
            None,
            default_candidates=(
                clusters_path.with_name("accession.fasta"),
                Path("work/accession.fasta"),
                Path("accession.fasta"),
            ),
        )
        fasta_index_path = _resolve_path(
            None,
            default_candidates=(
                fasta_path.with_suffix(fasta_path.suffix + ".idx"),
                fasta_path.with_name("accession.fasta.idx"),
            ),
        )

        progress("Loading GO graph and input tables")
        go_artifact = _load_pickle(library.path("go_graph.pkl"))
        go_index = build_fast_go_similarity_index(go_artifact)
        cluster_members = _load_cluster_members(clusters_path)
        accession_terms = _load_accession_terms(accession_go_path, go_index)
        cluster_go = _load_cluster_go(clusters_meta_path, go_index)
        progress(
            f"Loaded {len(cluster_members):,} clusters and "
            f"{len(accession_terms):,} accession GO profiles"
        )
        prefiltered_cluster_ids = _prefilter_cluster_ids(
            cluster_members=cluster_members,
            cluster_go=cluster_go,
            go_index=go_index,
            cluster_count=args.cluster_count,
            accessions_per_cluster=args.accessions_per_cluster,
            bound_min=args.bound_min,
            bound_max=args.bound_max,
            max_prefiltered_clusters=args.max_prefiltered_clusters,
            rng=rng,
        )
        progress(
            f"Sequence-screening {len(prefiltered_cluster_ids):,} clusters "
            f"(up to {args.candidate_pool_size} accessions each, "
            f"{args.alignment_workers} alignment workers)"
        )

        with ThreadPoolExecutor(
            max_workers=args.alignment_workers,
            thread_name_prefix="needle",
        ) as alignment_executor:
            candidates = _build_cluster_candidates(
                cluster_ids=prefiltered_cluster_ids,
                cluster_members=cluster_members,
                cluster_go=cluster_go,
                accession_terms=accession_terms,
                go_index=go_index,
                go_terms=go_artifact["terms"],
                fasta_path=fasta_path,
                fasta_index_path=fasta_index_path,
                accessions_per_cluster=args.accessions_per_cluster,
                max_pairwise_identity=args.max_pairwise_identity,
                min_within_cluster_go_similarity=(
                    args.min_within_cluster_go_similarity
                ),
                candidate_pool_size=args.candidate_pool_size,
                rng=rng,
                progress=progress,
                progress_every=args.progress_every,
                alignment_executor=alignment_executor,
            )

    if len(candidates) < args.cluster_count:
        raise SystemExit(
            "Not enough clusters passed the size, GO, and sequence-diversity "
            f"filters: {len(candidates)} available, need {args.cluster_count}."
        )

    progress(f"Selecting panel from {len(candidates):,} candidate clusters")
    selected, matrix = _select_cluster_panel(
        candidates,
        go_index=go_index,
        cluster_count=args.cluster_count,
        bound_min=args.bound_min,
        bound_max=args.bound_max,
        attempts=args.attempts,
        rng=rng,
    )
    if not selected:
        raise SystemExit(
            "No cluster panel found. Try relaxing --bound-min/--bound-max, "
            "--max-pairwise-identity, or increasing --attempts."
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    _write_cluster_report(args.out_dir / "selected_clusters.tsv", selected)
    _write_accession_report(
        args.out_dir / "selected_accessions.tsv",
        selected,
        accession_terms=accession_terms,
    )
    _write_combined_similarity_matrices(
        args.out_dir / "sequence_identity.txt",
        selected,
        go_similarity=False,
    )
    _write_combined_similarity_matrices(
        args.out_dir / "go_similarity.txt",
        selected,
        go_similarity=True,
    )
    _write_matrix(args.out_dir / "cluster_gosim_matrix.tsv", selected, matrix)
    _write_summary(
        args.out_dir / "summary.json",
        args=args,
        selected=selected,
        matrix=matrix,
        candidate_count=len(candidates),
    )
    print(f"Wrote selected simulation cluster panel to {args.out_dir}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Find 10 GO-similar clusters and 10 sequence-diverse accessions "
            "per cluster for simulation handoff."
        )
    )
    parser.add_argument("--clusters-meta", type=Path)
    parser.add_argument("--accession-go", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("simulation_panel"))
    parser.add_argument("--bound-min", type=float, default=0.4)
    parser.add_argument("--bound-max", type=float, default=0.6)
    parser.add_argument("--cluster-count", type=int, default=10)
    parser.add_argument("--accessions-per-cluster", type=int, default=10)
    parser.add_argument(
        "--max-pairwise-identity",
        type=float,
        default=0.9,
        help=(
            "Maximum allowed pairwise sequence identity among the selected "
            "accessions in each cluster."
        ),
    )
    parser.add_argument(
        "--min-within-cluster-go-similarity",
        type=float,
        default=0.9,
        help=(
            "Minimum allowed pairwise GO similarity among the selected "
            "accessions in each cluster."
        ),
    )
    parser.add_argument("--candidate-pool-size", type=int, default=20)
    parser.add_argument(
        "--max-prefiltered-clusters",
        type=int,
        default=2000,
        help=(
            "Maximum number of GO-compatible clusters to sequence-screen. "
            "Increase this if no panel is found."
        ),
    )
    parser.add_argument("--attempts", type=int, default=50000)
    parser.add_argument(
        "--alignment-workers",
        type=int,
        default=8,
        help=(
            "Number of concurrent EMBOSS Needle alignments (default: 8). "
            "Try the number of available CPU cores."
        ),
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1,
        help="Report sequence-screening progress every N clusters (default: 1).",
    )
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    if not 0.0 <= args.bound_min <= args.bound_max <= 1.0:
        parser.error("--bound-min and --bound-max must satisfy 0 <= min <= max <= 1")
    if args.cluster_count < 1:
        parser.error("--cluster-count must be at least 1")
    if args.accessions_per_cluster < 2:
        parser.error("--accessions-per-cluster must be at least 2")
    if not 0.0 <= args.max_pairwise_identity <= 1.0:
        parser.error("--max-pairwise-identity must be between 0 and 1")
    if not 0.0 <= args.min_within_cluster_go_similarity <= 1.0:
        parser.error(
            "--min-within-cluster-go-similarity must be between 0 and 1"
        )
    if args.candidate_pool_size < args.accessions_per_cluster:
        parser.error("--candidate-pool-size must be at least --accessions-per-cluster")
    if args.max_prefiltered_clusters < args.cluster_count:
        parser.error("--max-prefiltered-clusters must be at least --cluster-count")
    if args.attempts < 1:
        parser.error("--attempts must be at least 1")
    if args.alignment_workers < 1:
        parser.error("--alignment-workers must be at least 1")
    if args.progress_every < 1:
        parser.error("--progress-every must be at least 1")
    return args


def _resolve_path(path: Path | None, *, default_candidates: tuple[Path, ...]) -> Path:
    if path is not None:
        if not path.exists():
            raise SystemExit(f"Required path does not exist: {path}")
        return path
    for candidate in default_candidates:
        if candidate.exists():
            return candidate
    expected = ", ".join(str(candidate) for candidate in default_candidates)
    raise SystemExit(f"Could not find any default path among: {expected}")


def _load_pickle(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        artifact = pickle.load(handle)
    if not isinstance(artifact, dict):
        raise ValueError(f"Expected dictionary pickle artifact: {path}")
    return artifact


def _load_cluster_members(path: Path) -> dict[str, tuple[str, ...]]:
    members: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(path, reader, {"member_id", "cluster_id"})
        for row in reader:
            accession = row["member_id"].strip()
            cluster_id = row["cluster_id"].strip()
            if accession and cluster_id:
                members.setdefault(cluster_id, []).append(accession)
    return {cluster_id: tuple(sorted(values)) for cluster_id, values in members.items()}


def _load_accession_terms(
    path: Path,
    go_index: FastGoSimilarityIndex,
) -> dict[str, tuple[str, ...]]:
    accession_terms: dict[str, tuple[str, ...]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            accession = row[0].strip()
            terms = valid_go_profile(go_index, row[1].split(";"))
            if accession and terms:
                accession_terms[accession] = terms
    return accession_terms


def _load_cluster_go(
    path: Path,
    go_index: FastGoSimilarityIndex,
) -> dict[str, tuple[str, ...]]:
    cluster_go: dict[str, tuple[str, ...]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(path, reader, {"cluster_id", "composed_go"})
        for row in reader:
            cluster_id = row["cluster_id"].strip()
            terms = valid_go_profile(go_index, row["composed_go"].split(";"))
            if cluster_id and terms:
                cluster_go[cluster_id] = terms
    return cluster_go


def _require_columns(path: Path, reader: csv.DictReader, columns: set[str]) -> None:
    fieldnames = set(reader.fieldnames or ())
    missing = columns - fieldnames
    if missing:
        raise ValueError(f"{path} missing required column(s): {', '.join(missing)}")


def _build_cluster_candidates(
    *,
    cluster_ids: list[str],
    cluster_members: dict[str, tuple[str, ...]],
    cluster_go: dict[str, tuple[str, ...]],
    accession_terms: dict[str, tuple[str, ...]],
    go_index: FastGoSimilarityIndex,
    go_terms: dict[str, dict[str, Any]],
    fasta_path: Path,
    fasta_index_path: Path,
    accessions_per_cluster: int,
    max_pairwise_identity: float,
    min_within_cluster_go_similarity: float,
    candidate_pool_size: int,
    rng: random.Random,
    progress: Any | None = None,
    progress_every: int = 1,
    alignment_executor: Executor | None = None,
) -> list[ClusterCandidate]:
    candidates: list[ClusterCandidate] = []
    alignment_count = 0
    screen_started_at = time.monotonic()
    total = len(cluster_ids)

    def report(cluster_number: int) -> None:
        if progress is None or (
            cluster_number % progress_every != 0 and cluster_number != total
        ):
            return
        elapsed = time.monotonic() - screen_started_at
        rate = cluster_number / elapsed if elapsed else 0.0
        eta = (total - cluster_number) / rate if rate else 0.0
        progress(
            f"Screened {cluster_number:,}/{total:,} clusters "
            f"({cluster_number / total:.1%}); {len(candidates):,} passed; "
            f"{alignment_count:,} Needle alignments; ETA {eta:.0f}s"
        )

    for cluster_number, cluster_id in enumerate(cluster_ids, start=1):
        composed_go = cluster_go.get(cluster_id)
        if not composed_go:
            report(cluster_number)
            continue
        available = [
            accession
            for accession in cluster_members[cluster_id]
            if accession in accession_terms
        ]
        if len(available) < accessions_per_cluster:
            report(cluster_number)
            continue

        alignment_stats: dict[str, Any] = {"count": 0}
        sampled = _select_diverse_accessions(
            available,
            fasta_path=fasta_path,
            fasta_index_path=fasta_index_path,
            accession_terms=accession_terms,
            go_index=go_index,
            accessions_per_cluster=accessions_per_cluster,
            max_pairwise_identity=max_pairwise_identity,
            min_go_similarity=min_within_cluster_go_similarity,
            candidate_pool_size=candidate_pool_size,
            rng=rng,
            alignment_stats=alignment_stats,
            alignment_executor=alignment_executor,
        )
        alignment_count += alignment_stats["count"]
        if sampled is not None:
            identities = alignment_stats["selected_identities"]
            identity_matrix = alignment_stats["selected_identity_matrix"]
            go_similarity_matrix = alignment_stats["selected_go_similarity_matrix"]
            description = describe_go_function(cluster_id, composed_go, go_terms).summary
            candidates.append(
                ClusterCandidate(
                    cluster_id=cluster_id,
                    members=cluster_members[cluster_id],
                    selected_accessions=tuple(sampled),
                    composed_go=composed_go,
                    description=description,
                    within_mean_identity=sum(identities) / len(identities),
                    within_max_identity=max(identities),
                    within_min_identity=min(identities),
                    within_identity_matrix=identity_matrix,
                    within_go_similarity_matrix=go_similarity_matrix,
                )
            )

        report(cluster_number)
    rng.shuffle(candidates)
    return candidates


def _prefilter_cluster_ids(
    *,
    cluster_members: dict[str, tuple[str, ...]],
    cluster_go: dict[str, tuple[str, ...]],
    go_index: FastGoSimilarityIndex,
    cluster_count: int,
    accessions_per_cluster: int,
    bound_min: float,
    bound_max: float,
    max_prefiltered_clusters: int,
    rng: random.Random,
) -> list[str]:
    eligible = [
        cluster_id
        for cluster_id, members in cluster_members.items()
        if len(members) >= accessions_per_cluster and cluster_id in cluster_go
    ]
    rng.shuffle(eligible)
    eligible.sort(
        key=lambda cluster_id: (
            -len(cluster_members[cluster_id]),
            cluster_id,
        )
    )
    seeded_top = eligible[: max_prefiltered_clusters // 2]
    remaining = eligible[max_prefiltered_clusters // 2 :]
    rng.shuffle(remaining)
    return (seeded_top + remaining)[:max_prefiltered_clusters]


def _select_diverse_accessions(
    accessions: list[str],
    *,
    fasta_path: Path,
    fasta_index_path: Path,
    accession_terms: dict[str, tuple[str, ...]],
    go_index: FastGoSimilarityIndex,
    accessions_per_cluster: int,
    max_pairwise_identity: float,
    min_go_similarity: float,
    candidate_pool_size: int,
    rng: random.Random,
    alignment_stats: dict[str, Any] | None = None,
    alignment_executor: Executor | None = None,
) -> dict[str, str] | None:
    shuffled = list(accessions)
    rng.shuffle(shuffled)
    pool: dict[str, str] = {}
    for accession in shuffled:
        sequence = indexed_fasta_sequence(accession, fasta_path, fasta_index_path)
        if sequence:
            pool[accession] = sequence
        if len(pool) >= candidate_pool_size:
            break
    if len(pool) < accessions_per_cluster:
        return None

    go_similarity_cache: dict[tuple[str, str], float | None] = {}
    identity_cache: dict[tuple[str, str], float] = {}
    best: dict[str, str] | None = None
    for start in list(pool)[: min(10, len(pool))]:
        working_pool = dict(pool)
        selected = {start: pool[start]}
        while len(selected) < accessions_per_cluster:
            for candidate in tuple(working_pool):
                if candidate in selected:
                    continue
                if any(
                    (similarity := _accession_go_similarity(
                        go_index,
                        candidate,
                        selected_accession,
                        accession_terms=accession_terms,
                        similarity_cache=go_similarity_cache,
                    )) is None or similarity < min_go_similarity
                    for selected_accession in selected
                ):
                    working_pool.pop(candidate)
            candidate = _most_distant_accession(
                working_pool,
                selected,
                identity_cache=identity_cache,
                alignment_executor=alignment_executor,
            )
            if candidate is None:
                break
            identity = max(
                _cached_sequence_identity(
                    candidate,
                    working_pool[candidate],
                    selected_accession,
                    sequence,
                    identity_cache,
                )
                for selected_accession, sequence in selected.items()
            )
            if identity > max_pairwise_identity:
                break
            selected[candidate] = working_pool[candidate]
        if len(selected) == accessions_per_cluster:
            identities = _pairwise_identities(selected, identity_cache=identity_cache)
            if max(identities) <= max_pairwise_identity:
                if alignment_stats is not None:
                    alignment_stats["count"] = len(identity_cache)
                    alignment_stats["selected_identities"] = identities
                    alignment_stats["selected_identity_matrix"] = (
                        _sequence_identity_matrix(selected, identity_cache)
                    )
                    alignment_stats["selected_go_similarity_matrix"] = (
                        _accession_go_similarity_matrix(
                            selected,
                            go_index=go_index,
                            accession_terms=accession_terms,
                            similarity_cache=go_similarity_cache,
                        )
                    )
                return selected
            if best is None or max(identities) < max(
                _pairwise_identities(best, identity_cache=identity_cache)
            ):
                best = selected
    result = best if best is not None and len(best) == accessions_per_cluster else None
    if alignment_stats is not None:
        alignment_stats["count"] = len(identity_cache)
        if result is not None:
            alignment_stats["selected_identities"] = _pairwise_identities(
                result, identity_cache=identity_cache
            )
            alignment_stats["selected_identity_matrix"] = (
                _sequence_identity_matrix(result, identity_cache)
            )
            alignment_stats["selected_go_similarity_matrix"] = (
                _accession_go_similarity_matrix(
                    result,
                    go_index=go_index,
                    accession_terms=accession_terms,
                    similarity_cache=go_similarity_cache,
                )
            )
    return result


def _accession_go_similarity(
    go_index: FastGoSimilarityIndex,
    accession_a: str,
    accession_b: str,
    *,
    accession_terms: dict[str, tuple[str, ...]],
    similarity_cache: dict[tuple[str, str], float | None],
) -> float | None:
    key = tuple(sorted((accession_a, accession_b)))
    if key not in similarity_cache:
        similarity_cache[key] = set_lin_amb_fast_for_valid_profiles(
            go_index,
            accession_terms[accession_a],
            accession_terms[accession_b],
        )
    return similarity_cache[key]


def _accession_go_similarity_matrix(
    sequences: dict[str, str],
    *,
    go_index: FastGoSimilarityIndex,
    accession_terms: dict[str, tuple[str, ...]],
    similarity_cache: dict[tuple[str, str], float | None],
) -> tuple[tuple[float, ...], ...]:
    accessions = list(sequences)
    return tuple(
        tuple(
            _accession_go_similarity(
                go_index,
                accession_a,
                accession_b,
                accession_terms=accession_terms,
                similarity_cache=similarity_cache,
            )
            or 0.0
            for accession_b in accessions
        )
        for accession_a in accessions
    )


def _most_distant_accession(
    pool: dict[str, str],
    selected: dict[str, str],
    *,
    identity_cache: dict[tuple[str, str], float],
    alignment_executor: Executor | None = None,
) -> str | None:
    comparisons = [
        (accession, sequence, selected_accession, selected_sequence)
        for accession, sequence in pool.items()
        if accession not in selected
        for selected_accession, selected_sequence in selected.items()
    ]
    _populate_identity_cache(
        comparisons,
        cache=identity_cache,
        executor=alignment_executor,
    )
    best_accession: str | None = None
    best_identity = 2.0
    for accession, sequence in pool.items():
        if accession in selected:
            continue
        nearest_identity = max(
            _cached_sequence_identity(
                accession,
                sequence,
                selected_accession,
                selected_sequence,
                identity_cache,
            )
            for selected_accession, selected_sequence in selected.items()
        )
        if nearest_identity < best_identity:
            best_accession = accession
            best_identity = nearest_identity
    return best_accession


def _populate_identity_cache(
    comparisons: list[tuple[str, str, str, str]],
    *,
    cache: dict[tuple[str, str], float],
    executor: Executor | None,
) -> None:
    missing: dict[tuple[str, str], tuple[tuple[str, str], str, str]] = {}
    for accession_a, sequence_a, accession_b, sequence_b in comparisons:
        key = tuple(sorted((accession_a, accession_b)))
        if key not in cache:
            missing[key] = (key, sequence_a, sequence_b)
    if not missing:
        return
    if executor is None:
        results = map(_align_sequence_pair, missing.values())
    else:
        results = executor.map(_align_sequence_pair, missing.values())
    for key, identity in results:
        cache[key] = identity


def _align_sequence_pair(
    comparison: tuple[tuple[str, str], str, str],
) -> tuple[tuple[str, str], float]:
    key, sequence_a, sequence_b = comparison
    return key, sequence_identity(sequence_a, sequence_b)


def _cached_sequence_identity(
    accession_a: str,
    sequence_a: str,
    accession_b: str,
    sequence_b: str,
    cache: dict[tuple[str, str], float],
) -> float:
    key = tuple(sorted((accession_a, accession_b)))
    if key not in cache:
        cache[key] = sequence_identity(sequence_a, sequence_b)
    return cache[key]


def _pairwise_identities(
    sequences: dict[str, str],
    *,
    identity_cache: dict[tuple[str, str], float] | None = None,
) -> list[float]:
    cache = identity_cache if identity_cache is not None else {}
    accessions = list(sequences)
    identities: list[float] = []
    for index, accession_a in enumerate(accessions[:-1]):
        for accession_b in accessions[index + 1 :]:
            identities.append(
                _cached_sequence_identity(
                    accession_a,
                    sequences[accession_a],
                    accession_b,
                    sequences[accession_b],
                    cache,
                )
            )
    return identities


def _sequence_identity_matrix(
    sequences: dict[str, str],
    identity_cache: dict[tuple[str, str], float],
) -> tuple[tuple[float, ...], ...]:
    accessions = list(sequences)
    return tuple(
        tuple(
            1.0
            if accession_a == accession_b
            else _cached_sequence_identity(
                accession_a,
                sequences[accession_a],
                accession_b,
                sequences[accession_b],
                identity_cache,
            )
            for accession_b in accessions
        )
        for accession_a in accessions
    )


def _select_cluster_panel(
    candidates: list[ClusterCandidate],
    *,
    go_index: FastGoSimilarityIndex,
    cluster_count: int,
    bound_min: float,
    bound_max: float,
    attempts: int,
    rng: random.Random,
) -> tuple[list[ClusterCandidate], list[list[float]]] | tuple[None, None]:
    similarity_cache: dict[tuple[str, str], float | None] = {}
    for _ in range(attempts):
        seed = rng.choice(candidates)
        selected = [seed]
        remaining = list(candidates)
        rng.shuffle(remaining)
        for candidate in remaining:
            if candidate.cluster_id == seed.cluster_id:
                continue
            similarities = [
                _cluster_similarity(
                    go_index,
                    candidate,
                    current,
                    similarity_cache=similarity_cache,
                )
                for current in selected
            ]
            if all(
                score is not None and bound_min <= score <= bound_max
                for score in similarities
            ):
                selected.append(candidate)
            if len(selected) == cluster_count:
                matrix = _cluster_similarity_matrix(
                    go_index,
                    selected,
                    similarity_cache=similarity_cache,
                )
                return selected, matrix
    return None, None


def _cluster_similarity(
    go_index: FastGoSimilarityIndex,
    cluster_a: ClusterCandidate,
    cluster_b: ClusterCandidate,
    *,
    similarity_cache: dict[tuple[str, str], float | None],
) -> float | None:
    key = tuple(sorted((cluster_a.cluster_id, cluster_b.cluster_id)))
    if key not in similarity_cache:
        similarity_cache[key] = set_lin_amb_fast_for_valid_profiles(
            go_index,
            cluster_a.composed_go,
            cluster_b.composed_go,
        )
    return similarity_cache[key]


def _cluster_similarity_matrix(
    go_index: FastGoSimilarityIndex,
    selected: list[ClusterCandidate],
    *,
    similarity_cache: dict[tuple[str, str], float | None],
) -> list[list[float]]:
    matrix: list[list[float]] = []
    for row_cluster in selected:
        row: list[float] = []
        for column_cluster in selected:
            if row_cluster.cluster_id == column_cluster.cluster_id:
                row.append(1.0)
            else:
                row.append(
                    _cluster_similarity(
                        go_index,
                        row_cluster,
                        column_cluster,
                        similarity_cache=similarity_cache,
                    )
                    or 0.0
                )
        matrix.append(row)
    return matrix


def _write_cluster_report(path: Path, selected: list[ClusterCandidate]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "cluster_id",
                "cluster_size",
                "selected_accession_count",
                "composed_go",
                "description",
                "selected_accessions",
                "within_mean_sequence_identity",
                "within_max_sequence_identity",
                "within_min_sequence_identity",
            ]
        )
        for cluster in selected:
            writer.writerow(
                [
                    cluster.cluster_id,
                    len(cluster.members),
                    len(cluster.selected_accessions),
                    ";".join(cluster.composed_go),
                    cluster.description,
                    ";".join(cluster.selected_accessions),
                    f"{cluster.within_mean_identity:.5f}",
                    f"{cluster.within_max_identity:.5f}",
                    f"{cluster.within_min_identity:.5f}",
                ]
            )


def _write_accession_report(
    path: Path,
    selected: list[ClusterCandidate],
    *,
    accession_terms: dict[str, tuple[str, ...]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cluster_id", "accession", "GO_terms"])
        for cluster in selected:
            for accession in cluster.selected_accessions:
                writer.writerow(
                    [cluster.cluster_id, accession, ";".join(accession_terms[accession])]
                )


def _write_combined_similarity_matrices(
    path: Path,
    selected: list[ClusterCandidate],
    *,
    go_similarity: bool,
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for cluster_number, cluster in enumerate(selected):
            if cluster_number:
                handle.write("\n")
            handle.write(f"[{cluster.cluster_id}]\n")
            accessions = cluster.selected_accessions
            matrix = (
                cluster.within_go_similarity_matrix
                if go_similarity
                else cluster.within_identity_matrix
            )
            handle.write(f"{'':>10} " + " ".join(f"{item:>10}" for item in accessions))
            handle.write("\n")
            for accession, row in zip(accessions, matrix, strict=True):
                handle.write(
                    f"{accession:>10} "
                    + " ".join(f"{score:10.5f}" for score in row)
                    + "\n"
                )


def _write_matrix(
    path: Path,
    selected: list[ClusterCandidate],
    matrix: list[list[float]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cluster_id", *(cluster.cluster_id for cluster in selected)])
        for cluster, row in zip(selected, matrix, strict=True):
            writer.writerow([cluster.cluster_id, *(f"{score:.5f}" for score in row)])


def _write_summary(
    path: Path,
    *,
    args: argparse.Namespace,
    selected: list[ClusterCandidate],
    matrix: list[list[float]],
    candidate_count: int,
) -> None:
    off_diagonal = [
        matrix[row][column]
        for row in range(len(matrix))
        for column in range(len(matrix))
        if row < column
    ]
    summary = {
        "parameters": {
            "bound_min": args.bound_min,
            "bound_max": args.bound_max,
            "cluster_count": args.cluster_count,
            "accessions_per_cluster": args.accessions_per_cluster,
            "max_pairwise_identity": args.max_pairwise_identity,
            "min_within_cluster_go_similarity": (
                args.min_within_cluster_go_similarity
            ),
            "alignment_workers": args.alignment_workers,
            "seed": args.seed,
            "attempts": args.attempts,
        },
        "candidate_clusters_after_filters": candidate_count,
        "selected_clusters": [cluster.cluster_id for cluster in selected],
        "off_diagonal_similarity": {
            "min": min(off_diagonal),
            "max": max(off_diagonal),
            "mean": sum(off_diagonal) / len(off_diagonal),
        },
    }
    path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
