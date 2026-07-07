from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

from prosig.go.similarity import build_fast_go_similarity_index


def _load_script_module():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / (
        "compare_scan_similarity.py"
    )
    spec = importlib.util.spec_from_file_location(
        "compare_scan_similarity",
        script_path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _small_artifact() -> dict:
    return {
        "meta": {"schema_version": "1.0", "namespace": "molecular_function"},
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


def test_pairwise_rows_include_same_cluster_pairs_in_background() -> None:
    module = _load_script_module()
    go_index = build_fast_go_similarity_index(_small_artifact())
    selected = [
        module.SelectedAccession("cluster_a", "P1", ("GO:0000002",)),
        module.SelectedAccession("cluster_a", "P2", ("GO:0000002",)),
        module.SelectedAccession("cluster_b", "P3", ("GO:0000003",)),
    ]

    rows = module._pairwise_similarity_rows(selected, go_index)

    background = [row for row in rows if row["pair_scope"] == "background"]
    in_cluster = [row for row in rows if row["pair_scope"] == "in_cluster"]
    assert len(background) == 3
    assert len(in_cluster) == 1
    assert in_cluster[0]["similarity"] == 1.0


def test_stats_reports_missing_values_and_interpolated_quartiles() -> None:
    module = _load_script_module()

    stats = module._stats("example", [0.2, None, 0.4, 1.0])

    assert stats.metric == "example"
    assert stats.count == 3
    assert stats.missing_count == 1
    assert stats.mean == 0.5333333333333333
    assert stats.median == 0.4
    assert stats.q25 == 0.30000000000000004
    assert stats.q75 == 0.7


def test_sample_selected_accessions_uses_clusters_without_similarity_constraints(
    tmp_path: Path,
    monkeypatch,
) -> None:
    module = _load_script_module()
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text(
        "member_id\tcluster_id\n"
        "P1\tcluster_a\n"
        "P2\tcluster_a\n"
        "P3\tcluster_a\n"
        "P4\tcluster_b\n"
        "P5\tcluster_b\n"
        "P6\tcluster_b\n"
        "P7\tcluster_c\n"
        "P8\tcluster_c\n",
        encoding="utf-8",
    )
    accession_terms = {
        accession: ("GO:0000002",)
        for accession in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")
    }
    monkeypatch.setattr(
        module,
        "indexed_fasta_sequence",
        lambda accession, fasta_path, fasta_index_path: (
            None if accession == "P6" else "AAAA"
        ),
    )

    selected = module._sample_selected_accessions(
        clusters,
        accession_terms=accession_terms,
        fasta_path=tmp_path / "accession.fasta",
        fasta_index_path=tmp_path / "accession.fasta.idx",
        cluster_count=2,
        accessions_per_cluster=2,
        rng=module.random.Random(1),
    )

    selected_by_cluster: dict[str, list[str]] = {}
    for row in selected:
        selected_by_cluster.setdefault(row.cluster_id, []).append(row.accession)
    assert set(selected_by_cluster) == {"cluster_a", "cluster_c"}
    assert all(len(accessions) == 2 for accessions in selected_by_cluster.values())
