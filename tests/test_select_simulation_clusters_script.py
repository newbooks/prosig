from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


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


def test_sequence_identity_uses_longer_sequence_denominator() -> None:
    module = _load_script_module()

    assert module.sequence_identity("ABCD", "ABXY") == 0.5
    assert module.sequence_identity("ABCD", "AB") == 0.5
    assert module.sequence_identity("", "") == 0.0
