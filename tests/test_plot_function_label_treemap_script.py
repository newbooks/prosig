from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


def _load_script_module():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / (
        "plot_function_label_treemap.py"
    )
    spec = importlib.util.spec_from_file_location(
        "plot_function_label_treemap",
        script_path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_missing_go_graph_is_not_replaced_with_packaged_graph(tmp_path: Path) -> None:
    module = _load_script_module()
    missing_path = tmp_path / "misspelled-go-graph.pkl"

    with pytest.raises(SystemExit, match=f"GO graph file not found: {missing_path}"):
        module._load_go_terms(missing_path)
