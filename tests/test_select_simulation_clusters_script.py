from __future__ import annotations

import importlib.util
import random
import sys
import threading
from pathlib import Path

import pytest


def _load_script_module():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / (
        "select_simulation_clusters.py"
    )
    spec = importlib.util.spec_from_file_location(
        "select_simulation_clusters",
        script_path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_sequence_identity_uses_needle_alignment(monkeypatch) -> None:
    module = _load_script_module()

    calls: list[tuple[str, str]] = []

    def fake_identity(first: str, second: str) -> float:
        calls.append((first, second))
        return 0.375

    monkeypatch.setattr(module, "sequence_identity", fake_identity)

    assert module._pairwise_identities({"A": "ABCD", "B": "ABXY"}) == [0.375]
    assert calls == [("ABCD", "ABXY")]


def test_removed_input_overrides_are_not_command_options(monkeypatch, capsys) -> None:
    module = _load_script_module()

    for option in ("--library-dir", "--clusters", "--fasta", "--fasta-index"):
        monkeypatch.setattr(sys, "argv", ["select_simulation_clusters.py", option, "x"])
        with pytest.raises(SystemExit) as exc_info:
            module._parse_args()
        assert exc_info.value.code == 2
        assert "unrecognized arguments" in capsys.readouterr().err


def test_default_output_directory_is_in_current_directory(monkeypatch) -> None:
    module = _load_script_module()
    monkeypatch.setattr(sys, "argv", ["select_simulation_clusters.py"])

    assert module._parse_args().out_dir == Path("simulation_panel")


def test_default_selection_parameters_match_standard_command(monkeypatch) -> None:
    module = _load_script_module()
    monkeypatch.setattr(sys, "argv", ["select_simulation_clusters.py"])

    args = module._parse_args()
    assert args.cluster_count == 10
    assert args.accessions_per_cluster == 10
    assert args.min_within_cluster_go_similarity == 0.9
    assert args.bound_min == 0.4
    assert args.bound_max == 0.6


def test_within_cluster_go_similarity_is_enforced(monkeypatch, tmp_path) -> None:
    module = _load_script_module()
    terms = {"A": ("A",), "B": ("B",), "C": ("C",)}

    monkeypatch.setattr(
        module,
        "indexed_fasta_sequence",
        lambda accession, _fasta, _index: accession * 4,
    )
    monkeypatch.setattr(module, "sequence_identity", lambda _a, _b: 0.5)

    def go_similarity(_index, first, second):
        return 0.89 if {first[0], second[0]} == {"B", "C"} else 0.9

    monkeypatch.setattr(module, "set_lin_amb_fast_for_valid_profiles", go_similarity)

    selected = module._select_diverse_accessions(
        ["A", "B", "C"],
        fasta_path=tmp_path / "accession.fasta",
        fasta_index_path=tmp_path / "accession.fasta.idx",
        accession_terms=terms,
        go_index=object(),
        accessions_per_cluster=3,
        max_pairwise_identity=0.9,
        min_go_similarity=0.9,
        candidate_pool_size=3,
        rng=random.Random(1),
    )

    assert selected is None

    monkeypatch.setattr(
        module,
        "set_lin_amb_fast_for_valid_profiles",
        lambda _index, _first, _second: 0.9,
    )
    selected = module._select_diverse_accessions(
        ["A", "B", "C"],
        fasta_path=tmp_path / "accession.fasta",
        fasta_index_path=tmp_path / "accession.fasta.idx",
        accession_terms=terms,
        go_index=object(),
        accessions_per_cluster=3,
        max_pairwise_identity=0.9,
        min_go_similarity=0.9,
        candidate_pool_size=3,
        rng=random.Random(1),
    )

    assert selected is not None
    assert set(selected) == {"A", "B", "C"}


def test_sequence_alignments_can_run_concurrently(monkeypatch) -> None:
    module = _load_script_module()
    barrier = threading.Barrier(3)
    worker_names: set[str] = set()

    def fake_identity(_first: str, _second: str) -> float:
        worker_names.add(threading.current_thread().name)
        barrier.wait(timeout=2)
        return 0.5

    monkeypatch.setattr(module, "sequence_identity", fake_identity)
    cache: dict[tuple[str, str], float] = {}
    with module.ThreadPoolExecutor(max_workers=3) as executor:
        selected = module._most_distant_accession(
            {"A": "AAAA", "B": "BBBB", "C": "CCCC", "D": "DDDD"},
            {"A": "AAAA"},
            identity_cache=cache,
            alignment_executor=executor,
        )

    assert selected == "B"
    assert len(cache) == 3
    assert len(worker_names) == 3


def test_accession_report_and_within_cluster_matrix(tmp_path) -> None:
    module = _load_script_module()
    cluster = module.ClusterCandidate(
        cluster_id="cluster_0001",
        members=("A", "B"),
        selected_accessions=("A", "B"),
        composed_go=("GO:0001",),
        description="example",
        within_mean_identity=0.25,
        within_max_identity=0.25,
        within_min_identity=0.25,
        within_identity_matrix=((1.0, 0.25), (0.25, 1.0)),
        within_go_similarity_matrix=((1.0, 0.95), (0.95, 1.0)),
    )

    module._write_accession_report(
        tmp_path / "selected_accessions.tsv",
        [cluster],
        accession_terms={"A": ("GO:0001", "GO:0002"), "B": ("GO:0003",)},
    )
    module._write_combined_similarity_matrices(
        tmp_path / "sequence_identity.txt", [cluster], go_similarity=False
    )
    module._write_combined_similarity_matrices(
        tmp_path / "go_similarity.txt", [cluster], go_similarity=True
    )

    assert (tmp_path / "selected_accessions.tsv").read_text() == (
        "cluster_id\taccession\tGO_terms\n"
        "cluster_0001\tA\tGO:0001;GO:0002\n"
        "cluster_0001\tB\tGO:0003\n"
    )
    assert (tmp_path / "sequence_identity.txt").read_text() == (
        "[cluster_0001]\n"
        f"{'':>10} {'A':>10} {'B':>10}\n"
        f"{'A':>10} {1.0:10.5f} {0.25:10.5f}\n"
        f"{'B':>10} {0.25:10.5f} {1.0:10.5f}\n"
    )
    assert (tmp_path / "go_similarity.txt").read_text() == (
        "[cluster_0001]\n"
        f"{'':>10} {'A':>10} {'B':>10}\n"
        f"{'A':>10} {1.0:10.5f} {0.95:10.5f}\n"
        f"{'B':>10} {0.95:10.5f} {1.0:10.5f}\n"
    )


def test_cluster_report_labels_sequence_identity_columns(tmp_path) -> None:
    module = _load_script_module()
    cluster = module.ClusterCandidate(
        cluster_id="cluster_0001",
        members=("A", "B"),
        selected_accessions=("A", "B"),
        composed_go=("GO:0001",),
        description="example",
        within_mean_identity=0.25,
        within_max_identity=0.3,
        within_min_identity=0.2,
        within_identity_matrix=((1.0, 0.25), (0.25, 1.0)),
        within_go_similarity_matrix=((1.0, 0.95), (0.95, 1.0)),
    )

    module._write_cluster_report(tmp_path / "selected_clusters.tsv", [cluster])

    header = (tmp_path / "selected_clusters.tsv").read_text().splitlines()[0]
    assert header.endswith(
        "within_mean_sequence_identity\t"
        "within_max_sequence_identity\t"
        "within_min_sequence_identity"
    )
