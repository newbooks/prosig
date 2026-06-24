import json
import math
from pathlib import Path

import pytest

from prosig.go.clustering import (
    CandidateIndex,
    build_candidate_index,
    cluster_accessions_by_go,
    knn_edges_from_go_similarity,
    parse_cluster_config,
)
from prosig.go.similarity import build_fast_go_similarity_index


def test_build_candidate_index_filters_broad_ancestors_with_fallback() -> None:
    artifact = _small_artifact()
    index = build_fast_go_similarity_index(artifact)
    accessions = ["P1", "P2", "P3"]
    accession_terms = {
        "P1": ("GO:0000002",),
        "P2": ("GO:0000003",),
        "P3": ("GO:0000001",),
    }

    candidate_index = build_candidate_index(
        go_index=index,
        accessions=accessions,
        accession_terms=accession_terms,
        min_informative_ic=0.5,
        max_posting_fraction=0.34,
        max_posting_size=0,
    )

    assert candidate_index.posting_cap == 2
    assert candidate_index.informative_terms_before_filtering == 3
    assert candidate_index.postings_by_term["GO:0000001"] == [2]
    assert candidate_index.fallback_accessions_after_filtering == 1
    assert candidate_index.terms_by_accession_index[2] == ("GO:0000001",)


def test_cluster_accessions_by_go_writes_clusters_and_stats(tmp_path: Path) -> None:
    accession_go = tmp_path / "accession_mf_go.tsv"
    cluster_out = tmp_path / "go_clusters.tsv"
    stats_out = tmp_path / "go_clusters_stats.json"
    meta_out = tmp_path / "go_clusters_meta.tsv"
    accession_go.write_text(
        "P1\tGO:0000002\n"
        "P2\tGO:0000002\n"
        "P3\tGO:0000003\n"
        "P4\tGO:9999999\n",
        encoding="utf-8",
    )

    result = cluster_accessions_by_go(
        accession_go,
        go_artifact=_small_artifact(),
        go_graph_file=tmp_path / "go_graph.pkl",
        output_file=cluster_out,
        stats_file=stats_out,
        meta_file=meta_out,
        neighbors=1,
        term_cache_size_mb=1,
        profile_cache_size_mb=1,
        max_posting_fraction=1.0,
    )

    assert result.clustered_accessions == 3
    assert result.input_accessions == 4
    assert result.excluded_accessions == 1
    assert result.input_accessions == (
        result.clustered_accessions + result.excluded_accessions
    )
    assert result.meta_file == meta_out
    assert result.clusters >= 1
    assert cluster_out.read_text(encoding="utf-8").splitlines()[0] == (
        "member_id\tcluster_id"
    )
    assert "P1\tcluster_" in cluster_out.read_text(encoding="utf-8")
    stats = json.loads(stats_out.read_text(encoding="utf-8"))
    assert stats["algorithm"] == "go_set_similarity_knn_leiden"
    assert stats["clustered_accessions"] == 3
    assert stats["input_accessions"] == 4
    assert stats["cleaned_accessions"] == 3
    assert stats["excluded_accessions"] == 1
    assert stats["input_accessions"] == (
        stats["clustered_accessions"] + stats["excluded_accessions"]
    )
    assert stats["outputs"]["meta"] == str(meta_out)
    assert stats["lin_matrix"]["dtype"] == "float32"
    assert stats["lin_matrix"]["storage"] == "memory"
    assert stats["profile_cache"]["budget_mb"] == 1
    meta_lines = meta_out.read_text(encoding="utf-8").splitlines()
    assert meta_lines[0] == "cluster_id\tincluster_sim\tcomposed_description"
    assert any(line.startswith("cluster_0001\t") for line in meta_lines[1:])


def test_parse_cluster_config_reads_flat_yaml_and_validates(tmp_path: Path) -> None:
    config_path = tmp_path / "cluster_config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "stats_file: custom_stats.json",
                "meta_file: custom_meta.tsv",
                "neighbors: 7",
                "resolution: 0.75",
                "progress_interval_seconds: 12.5",
                "term_cache_size_mb: 16",
                "profile_cache_size_mb: 8",
                "min_informative_ic: 0.75",
                "max_posting_fraction: 0.2",
                "max_posting_size: 50",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    config = parse_cluster_config(config_path)

    assert config.stats_file == "custom_stats.json"
    assert config.meta_file == "custom_meta.tsv"
    assert config.neighbors == 7
    assert config.resolution == 0.75
    assert config.progress_interval_seconds == 12.5
    assert config.term_cache_size_mb == 16
    assert config.profile_cache_size_mb == 8
    assert config.min_informative_ic == 0.75
    assert config.max_posting_fraction == 0.2
    assert config.max_posting_size == 50

    config_path.write_text("max_posting_fraction: 2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="max posting fraction"):
        parse_cluster_config(config_path)

    config_path.write_text("neighbors: 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cluster neighbors"):
        parse_cluster_config(config_path)

    config_path.write_text("resolution: 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cluster resolution"):
        parse_cluster_config(config_path)


def test_identical_profile_ties_use_accession_order() -> None:
    artifact = _small_artifact()
    index = build_fast_go_similarity_index(artifact)
    accessions = ["P4", "P5", "P3", "P2", "P1"]
    accession_terms = {
        accession: ("GO:0000002",)
        for accession in accessions
    }
    candidate_index = CandidateIndex(
        postings_by_term={},
        terms_by_accession_index=[() for _ in accessions],
        informative_terms_before_filtering=0,
        informative_terms_after_filtering=0,
        posting_cap=0,
        fallback_accessions_after_filtering=0,
    )

    edges = knn_edges_from_go_similarity(
        go_index=index,
        accessions=accessions,
        accession_terms=accession_terms,
        candidate_index=candidate_index,
        neighbors=2,
        progress_interval_seconds=60.0,
    )

    p5_index = accessions.index("P5")
    p1_index = accessions.index("P1")
    p2_index = accessions.index("P2")
    assert (min(p5_index, p1_index), max(p5_index, p1_index)) in edges
    assert (min(p5_index, p2_index), max(p5_index, p2_index)) in edges
    assert all(weight == 1.0 for weight in edges.values())


def _small_artifact() -> dict:
    return {
        "meta": {
            "schema_version": "1.0",
            "namespace": "molecular_function",
        },
        "terms": {
            "GO:0003674": {
                "name": "molecular_function",
                "parents": [],
                "children": ["GO:0000001"],
                "ancestors": set(),
                "depth": 0,
                "freq": 1.0,
                "ic": 0.0,
            },
            "GO:0000001": {
                "name": "parent activity",
                "parents": ["GO:0003674"],
                "children": ["GO:0000002", "GO:0000003"],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": math.exp(-1.0),
                "ic": 1.0,
            },
            "GO:0000002": {
                "name": "child A activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": math.exp(-2.0),
                "ic": 2.0,
            },
            "GO:0000003": {
                "name": "child B activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": math.exp(-3.0),
                "ic": 3.0,
            },
        },
    }
