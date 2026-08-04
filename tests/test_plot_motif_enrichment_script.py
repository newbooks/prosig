from __future__ import annotations

import importlib.util
import pickle
import sys
from pathlib import Path

import pytest


def _load_script_module():
    script_path = (
        Path(__file__).resolve().parents[1] / "scripts" / "plot_motif_enrichment.py"
    )
    spec = importlib.util.spec_from_file_location("plot_motif_enrichment", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_motif_weights_and_missing_cluster_defaults_to_zero(tmp_path) -> None:
    module = _load_script_module()
    scoreboard = tmp_path / "scoreboard.pkl"
    scoreboard.write_bytes(
        pickle.dumps(
            {
                "weights": {
                    "MOTIF_A": {
                        "cluster_0001": {"weight": 3.25},
                        "cluster_0002": {"weight": 12.5},
                    }
                }
            }
        )
    )

    weights = module._load_motif_weights(scoreboard, "MOTIF_A")

    assert weights == {"cluster_0001": 3.25, "cluster_0002": 12.5}
    assert weights.get("cluster_0003", 0.0) == 0.0


def test_unknown_motif_is_rejected(tmp_path) -> None:
    module = _load_script_module()
    scoreboard = tmp_path / "scoreboard.pkl"
    scoreboard.write_bytes(pickle.dumps({"weights": {"MOTIF_A": {}}}))

    with pytest.raises(SystemExit, match="Motif not found"):
        module._load_motif_weights(scoreboard, "MISSING")


def test_load_motif_signature(tmp_path) -> None:
    module = _load_script_module()
    motifs = tmp_path / "prosig_motifs.tsv"
    motifs.write_text(
        "# ProSig motif library\n"
        "name\tdescription\tprosig_pattern\tstatus\n"
        "MOTIF_A\tExample\t[AG]?C\tprosig\n",
        encoding="utf-8",
    )

    assert module._load_motif_signature(motifs, "MOTIF_A") == "[AG]?C"


def test_summary_preserves_actual_weight_above_color_cap(tmp_path) -> None:
    module = _load_script_module()
    output = tmp_path / "summary.tsv"
    group = module.EnrichmentGroup(
        cluster_id="cluster_0001",
        label="cluster_0001<br>example",
        count=42,
        hover_label="cluster_0001 — example",
        description="example",
        weight=12.5,
    )

    module._write_summary(output, "MOTIF_A", [group])

    assert output.read_text() == (
        "rank\tmotif_id\tcluster_id\tcluster_label\tmember_count\tweight\n"
        "1\tMOTIF_A\tcluster_0001\texample\t42\t12.5\n"
    )


def test_safe_output_filename_component() -> None:
    module = _load_script_module()

    assert module._safe_filename_component("MOTIF/A B") == "MOTIF_A_B"
    assert module._safe_filename_component("...") == "motif"


def test_default_color_scale_is_plasma() -> None:
    module = _load_script_module()

    assert module.DEFAULT_COLOR_SCALE == "Plasma"


def test_template_group_id_is_mapped_to_scoreboard_cluster_id() -> None:
    module = _load_script_module()

    assert module._raw_cluster_id("cluster:cluster_0008") == "cluster_0008"


def test_click_script_pins_details_inside_tile() -> None:
    module = _load_script_module()

    script = module._click_to_pin_script()

    assert "plotly_click" in script
    assert "Enrichment weight" in script
    assert "Signature: %{customdata[3]}" in script
    assert "texttemplate" in script
    assert "Plotly.restyle" in script
    assert "annotations" not in script
