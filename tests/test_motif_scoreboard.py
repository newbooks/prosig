import json
import math
import pickle
from pathlib import Path

import pytest

from prosig.prediction.motif_scoreboard import build_motif_cluster_scoreboard


def test_build_motif_cluster_scoreboard_writes_positive_weights_and_meta(
    tmp_path: Path,
) -> None:
    clusters = tmp_path / "clusters.tsv"
    motif_hits = tmp_path / "motif_features.tsv"
    output = tmp_path / "motif_cluster_scoreboard.pkl"
    meta = tmp_path / "motif_cluster_scoreboard_meta.json"
    clusters.write_text(_cluster_rows(), encoding="utf-8")
    motif_hits.write_text(
        "accession\tmotif_id\tmotif_present\n"
        + "".join(f"A{i}\tmotif_strong\tTrue\n" for i in range(1, 7))
        + "".join(f"B{i}\tmotif_strong\tTrue\n" for i in range(1, 3))
        + "".join(f"A{i}\tmotif_low_support\tTrue\n" for i in range(1, 5))
        + "".join(f"A{i}\tmotif_non_positive\tTrue\n" for i in range(1, 6))
        + "".join(f"B{i}\tmotif_non_positive\tTrue\n" for i in range(1, 11))
        + "Z1\tmotif_strong\tTrue\n"
        + "A1\tmotif_strong\tTrue\n",
        encoding="utf-8",
    )

    stats = build_motif_cluster_scoreboard(
        cluster_file=clusters,
        motif_hits_file=motif_hits,
        output_file=output,
        meta_file=meta,
        min_cluster_size=10,
        min_support=5,
    )

    assert stats.accessions_in_clusters == 30
    assert stats.clusters == 3
    assert stats.clusters_below_min_size == 1
    assert stats.eligible_clusters == 2
    assert stats.motifs == 3
    assert stats.ignored_cluster_size == 3
    assert stats.ignored_low_support == 3
    assert stats.ignored_non_positive_weight == 1
    assert stats.stored_weights == 2
    assert stats.total_nonzero_weights == 2
    assert stats.motif_hits_outside_clusters == 1
    assert stats.duplicate_motif_hits == 1

    with output.open("rb") as handle:
        artifact = pickle.load(handle)
    assert artifact["schema_version"] == "1.0"
    assert set(artifact["weights"]) == {"motif_non_positive", "motif_strong"}
    record = artifact["weights"]["motif_strong"]["cluster_0001"]
    assert record["TP"] == 6
    assert record["FP"] == 2
    assert record["FN"] == 4
    assert record["TN"] == 18
    assert record["support"] == 6
    assert record["cluster_frequency"] == pytest.approx(6.5 / 11)
    assert record["background_frequency"] == pytest.approx(2.5 / 21)
    assert record["weight"] == pytest.approx(math.log2((6.5 / 11) / (2.5 / 21)))
    assert "cluster_0001" not in artifact["weights"]["motif_non_positive"]
    assert "cluster_0002" in artifact["weights"]["motif_non_positive"]

    meta_payload = json.loads(meta.read_text(encoding="utf-8"))
    assert meta_payload["stats"]["stored_weights"] == 2
    assert meta_payload["stats"]["ignored_cluster_size"] == 3
    assert meta_payload["stats"]["ignored_low_support"] == 3
    assert meta_payload["stats"]["ignored_non_positive_weight"] == 1


def test_build_motif_cluster_scoreboard_smooths_zero_false_positives(
    tmp_path: Path,
) -> None:
    clusters = tmp_path / "clusters.tsv"
    motif_hits = tmp_path / "motif_features.tsv"
    output = tmp_path / "motif_cluster_scoreboard.pkl"
    clusters.write_text(
        "member_id\tcluster_id\n"
        + "".join(f"A{i}\tcluster_0001\n" for i in range(1, 11))
        + "".join(f"B{i}\tcluster_0002\n" for i in range(1, 11)),
        encoding="utf-8",
    )
    motif_hits.write_text("accession\tmotif_id\nA1\tmotif_specific\n", encoding="utf-8")

    stats = build_motif_cluster_scoreboard(
        cluster_file=clusters,
        motif_hits_file=motif_hits,
        output_file=output,
        meta_file=None,
        min_cluster_size=10,
        min_support=1,
    )

    assert stats.stored_weights == 1
    with output.open("rb") as handle:
        artifact = pickle.load(handle)
    record = artifact["weights"]["motif_specific"]["cluster_0001"]
    assert record["TP"] == 1
    assert record["FP"] == 0
    assert record["cluster_frequency"] == pytest.approx(1.5 / 11)
    assert record["background_frequency"] == pytest.approx(0.5 / 11)
    assert record["weight"] == pytest.approx(math.log2(3))
    assert math.isfinite(record["weight"])
    assert artifact["parameters"]["smoothing"] == "Jeffreys prior"
    assert artifact["parameters"]["pseudocount"] == 0.5


def test_build_motif_cluster_scoreboard_reports_calibration(
    tmp_path: Path,
) -> None:
    clusters = tmp_path / "clusters.tsv"
    motif_hits = tmp_path / "motif_features.tsv"
    output = tmp_path / "motif_cluster_scoreboard.pkl"
    meta = tmp_path / "motif_cluster_scoreboard_meta.json"
    clusters.write_text(
        "member_id\tcluster_id\n"
        + "".join(f"A{i}\tcluster_0001\n" for i in range(1, 11))
        + "".join(f"B{i}\tcluster_0002\n" for i in range(1, 11))
        + "".join(f"S{i}\tcluster_0003\n" for i in range(1, 5)),
        encoding="utf-8",
    )
    motif_hits.write_text(
        "accession\tmotif_id\n"
        + "".join(f"A{i}\tmotif_specific\n" for i in range(1, 6)),
        encoding="utf-8",
    )

    stats = build_motif_cluster_scoreboard(
        cluster_file=clusters,
        motif_hits_file=motif_hits,
        output_file=output,
        meta_file=meta,
        min_cluster_size=10,
        min_support=5,
    )

    calibration_by_threshold = {
        point.weight_threshold: point
        for point in stats.calibration
    }
    assert set(calibration_by_threshold) == {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}
    assert calibration_by_threshold[2.0].eligible_accessions == 20
    assert calibration_by_threshold[2.0].covered_accessions == 5
    assert calibration_by_threshold[2.0].top1_correct_accessions == 5
    assert calibration_by_threshold[2.0].top3_correct_accessions == 5
    assert calibration_by_threshold[2.0].set_correct_accessions == 5
    assert calibration_by_threshold[2.0].coverage == pytest.approx(0.25)
    assert calibration_by_threshold[2.0].top1_accuracy == pytest.approx(1.0)
    assert calibration_by_threshold[2.0].top3_accuracy == pytest.approx(1.0)
    assert calibration_by_threshold[2.0].set_accuracy == pytest.approx(1.0)
    assert calibration_by_threshold[2.0].avg_predictions == pytest.approx(1.0)
    assert calibration_by_threshold[3.0].covered_accessions == 5
    assert calibration_by_threshold[4.0].covered_accessions == 0
    assert calibration_by_threshold[4.0].top1_accuracy is None
    assert calibration_by_threshold[4.0].top3_accuracy is None
    assert calibration_by_threshold[4.0].set_accuracy is None
    assert calibration_by_threshold[4.0].avg_predictions == 0.0

    meta_payload = json.loads(meta.read_text(encoding="utf-8"))
    assert meta_payload["stats"]["calibration"][0] == {
        "avg_predictions": 1.0,
        "coverage": 0.25,
        "covered_accessions": 5,
        "eligible_accessions": 20,
        "set_accuracy": 1.0,
        "set_correct_accessions": 5,
        "top1_accuracy": 1.0,
        "top1_correct_accessions": 5,
        "top3_accuracy": 1.0,
        "top3_correct_accessions": 5,
        "weight_threshold": 2.0,
    }


def test_build_motif_cluster_scoreboard_accepts_sparse_two_column_hits(
    tmp_path: Path,
) -> None:
    clusters = tmp_path / "clusters.tsv"
    motif_hits = tmp_path / "motif_features.tsv"
    output = tmp_path / "motif_cluster_scoreboard.pkl"
    clusters.write_text(_cluster_rows(), encoding="utf-8")
    motif_hits.write_text(
        "accession\tmotif_id\n"
        + "".join(f"A{i}\tmotif_strong\n" for i in range(1, 7))
        + "".join(f"B{i}\tmotif_strong\n" for i in range(1, 3)),
        encoding="utf-8",
    )

    stats = build_motif_cluster_scoreboard(
        cluster_file=clusters,
        motif_hits_file=motif_hits,
        output_file=output,
        meta_file=None,
        min_cluster_size=10,
        min_support=5,
    )

    assert stats.motif_hit_rows == 8
    with output.open("rb") as handle:
        artifact = pickle.load(handle)
    assert "cluster_0001" in artifact["weights"]["motif_strong"]


def test_build_motif_cluster_scoreboard_skips_completion_marker(
    tmp_path: Path,
) -> None:
    clusters = tmp_path / "clusters.tsv"
    motif_hits = tmp_path / "motif_features.tsv"
    output = tmp_path / "motif_cluster_scoreboard.pkl"
    clusters.write_text(_cluster_rows(), encoding="utf-8")
    motif_hits.write_text(
        "accession\tmotif_id\n"
        + "".join(f"A{i}\tmotif_strong\n" for i in range(1, 7))
        + "".join(f"B{i}\tmotif_strong\n" for i in range(1, 3))
        + "# completed\ttrue\n",
        encoding="utf-8",
    )

    stats = build_motif_cluster_scoreboard(
        cluster_file=clusters,
        motif_hits_file=motif_hits,
        output_file=output,
        meta_file=None,
        min_cluster_size=10,
        min_support=5,
    )

    assert stats.motif_hit_rows == 8
    assert stats.motif_hits_outside_clusters == 0
    assert stats.motifs == 1


def test_build_motif_cluster_scoreboard_validates_required_columns(
    tmp_path: Path,
) -> None:
    clusters = tmp_path / "clusters.tsv"
    motif_hits = tmp_path / "motif_features.tsv"
    clusters.write_text("member_id\tcluster_id\nA1\tcluster_0001\n", encoding="utf-8")
    motif_hits.write_text("accession\tcount\nA1\t1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing required column"):
        build_motif_cluster_scoreboard(
            cluster_file=clusters,
            motif_hits_file=motif_hits,
            output_file=tmp_path / "scoreboard.pkl",
            meta_file=None,
            min_cluster_size=1,
            min_support=1,
        )


def _cluster_rows() -> str:
    rows = ["member_id\tcluster_id\n"]
    rows.extend(f"A{i}\tcluster_0001\n" for i in range(1, 11))
    rows.extend(f"B{i}\tcluster_0002\n" for i in range(1, 16))
    rows.extend(f"C{i}\tcluster_0003\n" for i in range(1, 6))
    return "".join(rows)
