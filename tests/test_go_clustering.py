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
    refine_go_clusters_complete_linkage,
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

    assert result.clustered_accessions == 2
    assert result.input_accessions == 4
    assert result.excluded_accessions == 2
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
    assert stats["clustered_accessions"] == 2
    assert stats["input_accessions"] == 4
    assert stats["cleaned_accessions"] == 3
    assert stats["excluded_accessions"] == 2
    assert stats["min_similarity"] == 0.5
    assert stats["input_accessions"] == (
        stats["clustered_accessions"] + stats["excluded_accessions"]
    )
    assert stats["outputs"]["meta"] == str(meta_out)
    assert stats["lin_matrix"]["dtype"] == "float32"
    assert stats["lin_matrix"]["storage"] == "memory"
    assert stats["profile_cache"]["budget_mb"] == 1
    meta_lines = meta_out.read_text(encoding="utf-8").splitlines()
    assert meta_lines[0] == (
        "cluster_id\tsim_ave\tsim_min\tsim_max\tsize\tcomposed_go"
    )
    assert any(line.startswith("cluster_0001\t") for line in meta_lines[1:])
    meta_rows = [line.split("\t") for line in meta_lines[1:]]
    assert sum(int(row[4]) for row in meta_rows) == result.clustered_accessions
    for row in meta_rows:
        sim_ave, sim_min, sim_max = row[1:4]
        assert row[5] == "GO:0000002"
        if int(row[4]) == 1:
            assert (sim_ave, sim_min, sim_max) == ("NA", "NA", "NA")
        else:
            assert all(len(value) == 7 and value.count(".") == 1 for value in row[1:4])
            assert float(sim_min) <= float(sim_ave) <= float(sim_max)


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
                "min_similarity: 0.65",
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
    assert config.min_similarity == 0.65
    assert config.max_posting_fraction == 0.2
    assert config.max_posting_size == 50

    default_config_path = tmp_path / "default_cluster_config.yaml"
    default_config_path.write_text("", encoding="utf-8")
    assert parse_cluster_config(default_config_path).resolution == 2.0
    assert parse_cluster_config(default_config_path).min_similarity == 0.5
    assert parse_cluster_config(default_config_path).stats_file == (
        "leiden_clusters_stats.json"
    )
    assert parse_cluster_config(default_config_path).meta_file == (
        "leiden_clusters_meta.tsv"
    )

    config_path.write_text("max_posting_fraction: 2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="max posting fraction"):
        parse_cluster_config(config_path)

    config_path.write_text("neighbors: 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cluster neighbors"):
        parse_cluster_config(config_path)

    config_path.write_text("resolution: 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cluster resolution"):
        parse_cluster_config(config_path)

    config_path.write_text("min_similarity: 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="minimum similarity"):
        parse_cluster_config(config_path)

    config_path.write_text("min_similarity: 1.1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="minimum similarity"):
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


def test_knn_edges_enforce_inclusive_minimum_similarity() -> None:
    index = build_fast_go_similarity_index(_small_artifact())
    accessions = ["P1", "P2"]
    accession_terms = {
        "P1": ("GO:0000002",),
        "P2": ("GO:0000003",),
    }
    candidate_index = build_candidate_index(
        go_index=index,
        accessions=accessions,
        accession_terms=accession_terms,
        max_posting_fraction=1.0,
    )

    included_edges = knn_edges_from_go_similarity(
        go_index=index,
        accessions=accessions,
        accession_terms=accession_terms,
        candidate_index=candidate_index,
        neighbors=1,
        min_similarity=0.4,
    )
    excluded_edges = knn_edges_from_go_similarity(
        go_index=index,
        accessions=accessions,
        accession_terms=accession_terms,
        candidate_index=candidate_index,
        neighbors=1,
        min_similarity=0.400001,
    )

    assert included_edges == {(0, 1): pytest.approx(0.4)}
    assert excluded_edges == {}


def test_complete_linkage_refines_leiden_clusters_and_keeps_singletons(
    tmp_path: Path,
) -> None:
    accession_go = tmp_path / "accession_mf_go.tsv"
    leiden_clusters = tmp_path / "leiden_clusters.tsv"
    clusters_out = tmp_path / "clusters.tsv"
    meta_out = tmp_path / "clusters_meta.tsv"
    stats_out = tmp_path / "clusters_stats.json"
    accession_go.write_text(
        "P1\tGO:0000002\n"
        "P2\tGO:0000003\n"
        "P3\tGO:0000001\n"
        "P4\tGO:0000002\n",
        encoding="utf-8",
    )
    leiden_clusters.write_text(
        "member_id\tcluster_id\n"
        "P1\tcluster_0001\n"
        "P2\tcluster_0001\n"
        "P3\tcluster_0001\n"
        "P4\tcluster_0002\n",
        encoding="utf-8",
    )

    result = refine_go_clusters_complete_linkage(
        accession_go,
        leiden_clusters,
        go_artifact=_small_artifact(),
        go_graph_file=tmp_path / "go_graph.pkl",
        output_file=clusters_out,
        meta_file=meta_out,
        stats_file=stats_out,
        min_cluster_similarity=0.5,
        profile_cache_size_mb=1,
    )

    assert result.clustered_accessions == 4
    assert result.leiden_clusters == 2
    assert result.refined_clusters == 3
    assert result.leiden_singletons == 1
    assert result.refined_singletons == 2
    assert result.leiden_clusters_split == 1
    assert result.refinement_pairs_scored == 3
    assert len(clusters_out.read_text(encoding="utf-8").splitlines()) == 5
    meta_rows = [
        line.split("\t")
        for line in meta_out.read_text(encoding="utf-8").splitlines()[1:]
    ]
    assert sum(int(row[4]) == 1 for row in meta_rows) == 2
    assert all(
        row[2] == "NA" or float(row[2]) >= 0.5
        for row in meta_rows
    )
    assert all(row[5].startswith("GO:") for row in meta_rows)
    stats = json.loads(stats_out.read_text(encoding="utf-8"))
    assert stats["algorithm"] == "go_set_similarity_leiden_complete_linkage"
    assert stats["min_cluster_similarity"] == 0.5
    assert stats["leiden_singletons"] == 1
    assert stats["refined_singletons"] == 2


def test_complete_linkage_validates_minimum_similarity(tmp_path: Path) -> None:
    accession_go = tmp_path / "accession_mf_go.tsv"
    leiden_clusters = tmp_path / "leiden_clusters.tsv"
    accession_go.write_text("P1\tGO:0000002\n", encoding="utf-8")
    leiden_clusters.write_text(
        "member_id\tcluster_id\nP1\tcluster_0001\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="minimum cluster similarity"):
        refine_go_clusters_complete_linkage(
            accession_go,
            leiden_clusters,
            go_artifact=_small_artifact(),
            min_cluster_similarity=0.0,
        )


def test_complete_linkage_metadata_limits_composed_go_to_ten_terms(
    tmp_path: Path,
) -> None:
    accession_go = tmp_path / "accession_mf_go.tsv"
    leiden_clusters = tmp_path / "leiden_clusters.tsv"
    meta_out = tmp_path / "clusters_meta.tsv"
    go_terms = tuple(f"GO:{term_index:07d}" for term_index in range(1, 12))
    accession_go.write_text(
        f"P1\t{';'.join(go_terms)}\n",
        encoding="utf-8",
    )
    leiden_clusters.write_text(
        "member_id\tcluster_id\nP1\tcluster_0001\n",
        encoding="utf-8",
    )

    refine_go_clusters_complete_linkage(
        accession_go,
        leiden_clusters,
        go_artifact=_many_term_artifact(go_terms),
        output_file=tmp_path / "clusters.tsv",
        meta_file=meta_out,
        stats_file=None,
    )

    row = meta_out.read_text(encoding="utf-8").splitlines()[1].split("\t")
    composed_go = row[5].split(";")
    assert len(composed_go) == 10
    assert composed_go == list(reversed(go_terms[1:]))


def test_composed_go_suppresses_low_coverage_child_when_parent_ranks_later(
    tmp_path: Path,
) -> None:
    accession_go = tmp_path / "accession_mf_go.tsv"
    leiden_clusters = tmp_path / "leiden_clusters.tsv"
    meta_out = tmp_path / "clusters_meta.tsv"
    accession_go.write_text(
        "P1\tGO:0000002\n"
        "P2\tGO:0000001\n"
        "P3\tGO:0000001\n"
        "P4\tGO:0000001\n",
        encoding="utf-8",
    )
    leiden_clusters.write_text(
        "member_id\tcluster_id\n"
        "P1\tcluster_0001\n"
        "P2\tcluster_0001\n"
        "P3\tcluster_0001\n"
        "P4\tcluster_0001\n",
        encoding="utf-8",
    )

    refine_go_clusters_complete_linkage(
        accession_go,
        leiden_clusters,
        go_artifact=_low_coverage_child_artifact(),
        output_file=tmp_path / "clusters.tsv",
        meta_file=meta_out,
        stats_file=None,
        min_cluster_similarity=0.1,
    )

    row = meta_out.read_text(encoding="utf-8").splitlines()[1].split("\t")
    assert row[5] == "GO:0000001"


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


def _low_coverage_child_artifact() -> dict:
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
                "children": ["GO:0000002"],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": math.exp(-1.0),
                "ic": 1.0,
            },
            "GO:0000002": {
                "name": "specific child activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": math.exp(-5.0),
                "ic": 5.0,
            },
        },
    }


def _many_term_artifact(go_terms: tuple[str, ...]) -> dict:
    terms = {
        "GO:0003674": {
            "name": "molecular_function",
            "parents": [],
            "children": list(go_terms),
            "ancestors": set(),
            "depth": 0,
            "freq": 1.0,
            "ic": 0.0,
        }
    }
    for term_index, go_id in enumerate(go_terms, start=1):
        terms[go_id] = {
            "name": f"activity {term_index}",
            "parents": ["GO:0003674"],
            "children": [],
            "ancestors": {"GO:0003674"},
            "depth": 1,
            "freq": math.exp(-float(term_index)),
            "ic": float(term_index),
        }
    return {
        "meta": {
            "schema_version": "1.0",
            "namespace": "molecular_function",
        },
        "terms": terms,
    }
