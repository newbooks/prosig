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
    assert record["cluster_frequency"] == pytest.approx(0.6)
    assert record["background_frequency"] == pytest.approx(2 / 20)
    assert record["weight"] == pytest.approx(math.log2(0.6 / (2 / 20)))
    assert "cluster_0001" not in artifact["weights"]["motif_non_positive"]
    assert "cluster_0002" in artifact["weights"]["motif_non_positive"]

    meta_payload = json.loads(meta.read_text(encoding="utf-8"))
    assert meta_payload["stats"]["stored_weights"] == 2
    assert meta_payload["stats"]["ignored_cluster_size"] == 3
    assert meta_payload["stats"]["ignored_low_support"] == 3
    assert meta_payload["stats"]["ignored_non_positive_weight"] == 1


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
