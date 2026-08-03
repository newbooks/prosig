#!/usr/bin/env python3
"""Compare background, in-cluster, and motif-scan prediction similarities."""

from __future__ import annotations

import argparse
import csv
import json
import pickle
import random
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if REPO_SRC.is_dir():
    sys.path.insert(0, str(REPO_SRC))

from prosig.cli.scan import (  # noqa: E402
    _load_calibration,
    _load_cluster_meta,
    _load_scoreboard,
    _QuerySequence,
    _scan_one_query,
)
from prosig.go.similarity import (  # noqa: E402
    FastGoSimilarityIndex,
    build_fast_go_similarity_index,
    set_lin_amb_fast_for_valid_profiles,
    valid_go_profile,
)
from prosig.library import resolve_core_library  # noqa: E402
from prosig.motifs.scanning import read_prosig_motif_library  # noqa: E402
from prosig.sequences import indexed_fasta_sequence  # noqa: E402


@dataclass(frozen=True)
class SelectedAccession:
    cluster_id: str
    accession: str
    true_go: tuple[str, ...]


@dataclass(frozen=True)
class SimilarityStats:
    metric: str
    count: int
    missing_count: int
    mean: float | None
    median: float | None
    stdev: float | None
    minimum: float | None
    q25: float | None
    q75: float | None
    maximum: float | None


def main() -> None:
    args = _parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    with resolve_core_library(args.library_dir) as library:
        go_artifact = _load_pickle(library.path("go_graph.pkl"))
        go_index = build_fast_go_similarity_index(go_artifact)
        accession_terms = _load_accession_terms(
            args.accession_go or library.path("accession_mf_go.tsv"),
            go_index,
        )
        if args.selected_accessions is None:
            selected = _sample_selected_accessions(
                args.clusters,
                accession_terms=accession_terms,
                fasta_path=args.fasta,
                fasta_index_path=args.fasta_index,
                cluster_count=args.cluster_count,
                accessions_per_cluster=args.accessions_per_cluster,
                rng=rng,
            )
        else:
            selected = _load_selected_accessions(
                args.selected_accessions,
                accession_terms=accession_terms,
            )
        if len(selected) < 2:
            raise SystemExit("Need at least two selected accessions.")

        pairwise_rows = _pairwise_similarity_rows(selected, go_index)
        prediction_rows = _prediction_similarity_rows(
            selected,
            fasta_path=args.fasta,
            fasta_index_path=args.fasta_index,
            library_dir=library.directory,
            go_index=go_index,
            min_weight=args.min_weight,
        )

    summary = [
        _stats(
            "background_similarity",
            [
                row["similarity"]
                for row in pairwise_rows
                if row["pair_scope"] == "background"
            ],
        ),
        _stats(
            "in_cluster_similarity",
            [
                row["similarity"]
                for row in pairwise_rows
                if row["pair_scope"] == "in_cluster"
            ],
        ),
        _stats(
            "prediction_power",
            [row["prediction_similarity"] for row in prediction_rows],
        ),
    ]

    _write_pairwise(args.out_dir / "pairwise_similarity.tsv", pairwise_rows)
    _write_predictions(args.out_dir / "prediction_similarity.tsv", prediction_rows)
    _write_selected_accessions(args.out_dir / "selected_accessions.tsv", selected)
    _write_summary(args.out_dir / "similarity_summary.tsv", summary)
    _write_summary_json(args.out_dir / "similarity_summary.json", summary)

    print(f"Wrote scan similarity comparison to {args.out_dir}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare GO set similarity among selected accessions, within "
            "clusters, and between accessions and their top ProSig scan match."
        )
    )
    parser.add_argument(
        "--selected-accessions",
        type=Path,
        help=(
            "Optional preselected TSV with cluster_id and accession columns. "
            "If omitted, accessions are sampled from --clusters."
        ),
    )
    parser.add_argument(
        "--clusters",
        type=Path,
        default=Path("work/clusters.tsv"),
        help="TSV with member_id and cluster_id columns used for random sampling.",
    )
    parser.add_argument(
        "--library-dir",
        type=Path,
        default=Path("work"),
        help="Directory containing ProSig runtime library files.",
    )
    parser.add_argument(
        "--accession-go",
        type=Path,
        help=(
            "Optional accession_mf_go.tsv override. Defaults to the file in "
            "--library-dir."
        ),
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path("work/accession.fasta"),
        help="Indexed FASTA containing selected accession sequences.",
    )
    parser.add_argument(
        "--fasta-index",
        type=Path,
        default=Path("work/accession.fasta.idx"),
        help="Byte-offset index for --fasta.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("work/random_cluster_scan_similarity"),
    )
    parser.add_argument(
        "--min-weight",
        type=float,
        default=2.0,
        help="Minimum motif-cluster weight retained by the scan step.",
    )
    parser.add_argument("--cluster-count", type=int, default=10)
    parser.add_argument("--accessions-per-cluster", type=int, default=10)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    if args.cluster_count < 1:
        parser.error("--cluster-count must be at least 1")
    if args.accessions_per_cluster < 2:
        parser.error("--accessions-per-cluster must be at least 2")
    return args


def _load_pickle(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        artifact = pickle.load(handle)
    if not isinstance(artifact, dict):
        raise ValueError(f"Expected dictionary pickle artifact: {path}")
    return artifact


def _load_selected_accessions(
    path: Path,
    *,
    accession_terms: dict[str, tuple[str, ...]],
) -> list[SelectedAccession]:
    selected: list[SelectedAccession] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(path, reader, {"cluster_id", "accession"})
        for row in reader:
            cluster_id = row["cluster_id"].strip()
            accession = row["accession"].strip()
            if not cluster_id or not accession or accession in seen:
                continue
            true_go = accession_terms.get(accession)
            if not true_go:
                raise ValueError(
                    f"No valid Molecular Function GO profile for {accession}"
                )
            selected.append(
                SelectedAccession(
                    cluster_id=cluster_id,
                    accession=accession,
                    true_go=true_go,
                )
            )
            seen.add(accession)
    return selected


def _sample_selected_accessions(
    path: Path,
    *,
    accession_terms: dict[str, tuple[str, ...]],
    fasta_path: Path,
    fasta_index_path: Path,
    cluster_count: int,
    accessions_per_cluster: int,
    rng: random.Random,
) -> list[SelectedAccession]:
    members_by_cluster = _load_cluster_members(path)
    eligible: list[tuple[str, list[str]]] = []
    for cluster_id, members in sorted(members_by_cluster.items()):
        accessions = [
            accession
            for accession in members
            if accession in accession_terms
            and indexed_fasta_sequence(accession, fasta_path, fasta_index_path)
            is not None
        ]
        if len(accessions) >= accessions_per_cluster:
            eligible.append((cluster_id, accessions))

    if len(eligible) < cluster_count:
        raise SystemExit(
            "Not enough clusters have the requested number of accessions with "
            "valid GO profiles and sequences: "
            f"{len(eligible)} available, need {cluster_count}."
        )

    selected: list[SelectedAccession] = []
    for cluster_id, accessions in rng.sample(eligible, cluster_count):
        for accession in rng.sample(accessions, accessions_per_cluster):
            selected.append(
                SelectedAccession(
                    cluster_id=cluster_id,
                    accession=accession,
                    true_go=accession_terms[accession],
                )
            )
    selected.sort(key=lambda row: (row.cluster_id, row.accession))
    return selected


def _load_cluster_members(path: Path) -> dict[str, list[str]]:
    members: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(path, reader, {"member_id", "cluster_id"})
        for row in reader:
            accession = row["member_id"].strip()
            cluster_id = row["cluster_id"].strip()
            if accession and cluster_id:
                members.setdefault(cluster_id, []).append(accession)
    return members


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


def _require_columns(
    path: Path,
    reader: csv.DictReader,
    required: set[str],
) -> None:
    fieldnames = set(reader.fieldnames or ())
    missing = required - fieldnames
    if missing:
        raise ValueError(f"{path} missing required column(s): {', '.join(missing)}")


def _pairwise_similarity_rows(
    selected: list[SelectedAccession],
    go_index: FastGoSimilarityIndex,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    cache: dict[tuple[tuple[int, ...], tuple[int, ...]], float | None] = {}
    for index, accession_a in enumerate(selected[:-1]):
        for accession_b in selected[index + 1 :]:
            similarity = set_lin_amb_fast_for_valid_profiles(
                go_index,
                accession_a.true_go,
                accession_b.true_go,
                profile_pair_cache=cache,
            )
            scope = (
                "in_cluster"
                if accession_a.cluster_id == accession_b.cluster_id
                else "background"
            )
            rows.append(
                {
                    "pair_scope": scope,
                    "cluster_id_a": accession_a.cluster_id,
                    "accession_a": accession_a.accession,
                    "true_go_a": ";".join(accession_a.true_go),
                    "cluster_id_b": accession_b.cluster_id,
                    "accession_b": accession_b.accession,
                    "true_go_b": ";".join(accession_b.true_go),
                    "similarity": similarity,
                }
            )
            if scope == "in_cluster":
                rows.append({**rows[-1], "pair_scope": "background"})
    return rows


def _prediction_similarity_rows(
    selected: list[SelectedAccession],
    *,
    fasta_path: Path,
    fasta_index_path: Path,
    library_dir: Path,
    go_index: FastGoSimilarityIndex,
    min_weight: float,
) -> list[dict[str, Any]]:
    motifs = read_prosig_motif_library(library_dir / "prosig_motifs.tsv")
    motifs_by_name = {motif.name: motif for motif in motifs}
    scoreboard = _load_scoreboard(library_dir / "motif_cluster_scoreboard.pkl")
    clusters = _load_cluster_meta(library_dir / "clusters_meta.tsv")
    calibration = _load_calibration(library_dir / "motif_cluster_scoreboard_meta.json")
    rows: list[dict[str, Any]] = []
    cache: dict[tuple[tuple[int, ...], tuple[int, ...]], float | None] = {}

    for selected_accession in selected:
        sequence = indexed_fasta_sequence(
            selected_accession.accession,
            fasta_path,
            fasta_index_path,
        )
        if sequence is None:
            raise ValueError(f"No sequence found for {selected_accession.accession}")
        report = _scan_one_query(
            _QuerySequence(selected_accession.accession, sequence),
            motifs=motifs,
            motifs_by_name=motifs_by_name,
            scoreboard=scoreboard,
            clusters=clusters,
            calibration=calibration,
            similarity=None,
            min_weight=min_weight,
            top_n=1,
        )
        predictions = report["inferred_go_sets"]
        top_prediction = predictions[0] if predictions else None
        prediction_go = _prediction_go_profile(top_prediction, go_index)
        prediction_similarity = (
            set_lin_amb_fast_for_valid_profiles(
                go_index,
                selected_accession.true_go,
                prediction_go,
                profile_pair_cache=cache,
            )
            if prediction_go
            else None
        )
        rows.append(
            {
                "cluster_id": selected_accession.cluster_id,
                "accession": selected_accession.accession,
                "true_go": ";".join(selected_accession.true_go),
                "matched_motif_count": len(report["matched_motifs"]),
                "prediction_go": ";".join(prediction_go),
                "prediction_similarity": prediction_similarity,
                "prediction_weight": (
                    top_prediction.get("weight") if top_prediction is not None else None
                ),
                "prediction_motif_id": (
                    top_prediction.get("motif_id")
                    if top_prediction is not None
                    else ""
                ),
                "prediction_cluster_ids": (
                    ",".join(top_prediction.get("cluster_ids", ()))
                    if top_prediction is not None
                    else ""
                ),
            }
        )
    return rows


def _prediction_go_profile(
    prediction: dict[str, Any] | None,
    go_index: FastGoSimilarityIndex,
) -> tuple[str, ...]:
    if prediction is None:
        return ()
    raw_terms = prediction.get("go_terms", ())
    if not isinstance(raw_terms, list):
        return ()
    return valid_go_profile(go_index, [str(term) for term in raw_terms])


def _stats(metric: str, values: list[float | None]) -> SimilarityStats:
    observed = sorted(value for value in values if value is not None)
    if not observed:
        return SimilarityStats(
            metric=metric,
            count=0,
            missing_count=len(values),
            mean=None,
            median=None,
            stdev=None,
            minimum=None,
            q25=None,
            q75=None,
            maximum=None,
        )
    return SimilarityStats(
        metric=metric,
        count=len(observed),
        missing_count=len(values) - len(observed),
        mean=statistics.fmean(observed),
        median=statistics.median(observed),
        stdev=statistics.stdev(observed) if len(observed) > 1 else 0.0,
        minimum=observed[0],
        q25=_quantile(observed, 0.25),
        q75=_quantile(observed, 0.75),
        maximum=observed[-1],
    )


def _quantile(sorted_values: list[float], fraction: float) -> float:
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = position - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight


def _write_pairwise(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "pair_scope",
        "cluster_id_a",
        "accession_a",
        "true_go_a",
        "cluster_id_b",
        "accession_b",
        "true_go_b",
        "similarity",
    ]
    _write_dict_rows(path, fieldnames, rows)


def _write_predictions(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "cluster_id",
        "accession",
        "true_go",
        "matched_motif_count",
        "prediction_go",
        "prediction_similarity",
        "prediction_weight",
        "prediction_motif_id",
        "prediction_cluster_ids",
    ]
    _write_dict_rows(path, fieldnames, rows)


def _write_selected_accessions(
    path: Path,
    selected: list[SelectedAccession],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cluster_id", "accession", "true_go"])
        for row in selected:
            writer.writerow([row.cluster_id, row.accession, ";".join(row.true_go)])


def _write_dict_rows(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, Any]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    fieldname: _format_cell(row.get(fieldname))
                    for fieldname in fieldnames
                }
            )


def _write_summary(path: Path, rows: list[SimilarityStats]) -> None:
    fieldnames = [
        "metric",
        "count",
        "missing_count",
        "mean",
        "median",
        "stdev",
        "min",
        "q25",
        "q75",
        "max",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fieldnames)
        for row in rows:
            writer.writerow(
                [
                    row.metric,
                    row.count,
                    row.missing_count,
                    _format_cell(row.mean),
                    _format_cell(row.median),
                    _format_cell(row.stdev),
                    _format_cell(row.minimum),
                    _format_cell(row.q25),
                    _format_cell(row.q75),
                    _format_cell(row.maximum),
                ]
            )


def _write_summary_json(path: Path, rows: list[SimilarityStats]) -> None:
    payload = [
        {
            "metric": row.metric,
            "count": row.count,
            "missing_count": row.missing_count,
            "mean": row.mean,
            "median": row.median,
            "stdev": row.stdev,
            "min": row.minimum,
            "q25": row.q25,
            "q75": row.q75,
            "max": row.maximum,
        }
        for row in rows
    ]
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _format_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


if __name__ == "__main__":
    main()
