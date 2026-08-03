from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_script_module():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / (
        "plot_cluster_mds.py"
    )
    spec = importlib.util.spec_from_file_location("plot_cluster_mds", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_numeric_cluster_ids_remain_strings_for_line_lookup(tmp_path: Path) -> None:
    module = _load_script_module()
    matrix_path = tmp_path / "similarity.tsv"
    matrix_path.write_text(
        "cluster_id\t101\t102\n"
        "101\t1.0\t0.5\n"
        "102\t0.5\t1.0\n",
        encoding="utf-8",
    )

    similarity = module._read_similarity_matrix(matrix_path)

    assert list(similarity.index) == ["101", "102"]
    assert list(similarity.columns) == ["101", "102"]
    assert similarity.loc["101", "102"] == 0.5


def test_leading_zeros_in_cluster_ids_are_preserved(tmp_path: Path) -> None:
    module = _load_script_module()
    matrix_path = tmp_path / "similarity.tsv"
    matrix_path.write_text(
        "cluster_id\t001\t002\n"
        "001\t1.0\t0.5\n"
        "002\t0.5\t1.0\n",
        encoding="utf-8",
    )

    similarity = module._read_similarity_matrix(matrix_path)
    module._validate_distance_matrix(
        module._similarity_to_distance(similarity),
        atol=1e-8,
    )

    assert list(similarity.index) == ["001", "002"]
    assert list(similarity.columns) == ["001", "002"]
    assert similarity.loc["001", "002"] == 0.5


def test_output_extensions_preserve_dotted_prefix() -> None:
    module = _load_script_module()
    prefix = Path("results/mds.v1")

    assert module._append_extension(prefix, ".png") == Path("results/mds.v1.png")
    assert module._append_extension(prefix, ".pdf") == Path("results/mds.v1.pdf")
