import json
import pickle

from typer.testing import CliRunner

from prosig.cli.app import app


def test_version_command() -> None:
    result = CliRunner().invoke(app, ["version"])

    assert result.exit_code == 0
    assert "ProSig version:" in result.stdout
    assert "Developer: Junjun Mao <junjun.mao@gmail.com>" in result.stdout


def test_short_help_option() -> None:
    result = CliRunner().invoke(app, ["-h"])

    assert result.exit_code == 0
    assert "Usage:" in result.stdout
    assert "version" in result.stdout
    assert "build-library" in result.stdout
    assert "scan" in result.stdout
    assert "inspect" in result.stdout


def test_scan_sequence_reports_inferred_go_sets(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    _write_scan_artifacts(tmp_path)

    result = CliRunner().invoke(app, ["scan", "--seq", "XXAAXX"])

    assert result.exit_code == 0
    assert "Library:       current-directory" in result.stdout
    assert "Query:          sequence" in result.stdout
    assert "Matched motifs: 1" in result.stdout
    assert "Inferred GO sets (top 5, weight >= 2):" in result.stdout
    assert "1. GO:0004672;GO:0005524" in result.stdout
    assert "Signature:     AA" in result.stdout
    assert "Clusters:       cluster_0001" in result.stdout
    assert "Description:" in result.stdout
    assert "ATP-binding protein kinase" in result.stdout
    assert "GO terms:       GO:0004672;GO:0005524" in result.stdout
    assert "Weight:         6.5" in result.stdout
    assert "Confidence:     0.91 (set_acc @ >= 5)" in result.stdout
    assert "5. GO:000005" in result.stdout
    assert "6. GO:000006" not in result.stdout


def test_scan_writes_json_output(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    _write_scan_artifacts(tmp_path)
    json_out = tmp_path / "scan.json"

    result = CliRunner().invoke(
        app,
        ["scan", "--seq", "XXAAXX", "--json-out", str(json_out)],
    )

    assert result.exit_code == 0
    assert result.stdout == ""
    payload = json.loads(json_out.read_text(encoding="utf-8"))
    prediction = payload["queries"][0]["inferred_go_sets"][0]
    assert prediction["go_terms"] == ["GO:0004672", "GO:0005524"]
    assert prediction["signature"] == "AA"
    assert prediction["motif_id"] == "MOTIF_A"
    assert prediction["weight"] == 6.5
    assert prediction["calibrated_confidence"]["set_accuracy"] == 0.91
    assert payload["top_n"] == 5
    assert payload["library"]["source"] == "current-directory"


def test_scan_top_n_zero_reports_all_inferences(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    _write_scan_artifacts(tmp_path)

    result = CliRunner().invoke(app, ["scan", "--seq", "XXAAXX", "--top-n", "0"])

    assert result.exit_code == 0
    assert "Inferred GO sets (all, weight >= 2):" in result.stdout
    assert "6. GO:000006" in result.stdout


def test_scan_reads_fasta_queries(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    _write_scan_artifacts(tmp_path)
    query_fasta = tmp_path / "queries.fasta"
    query_fasta.write_text(">query_1\nXXBBXX\n", encoding="ascii")

    result = CliRunner().invoke(app, ["scan", "--fasta", str(query_fasta)])

    assert result.exit_code == 0
    assert "Query:          query_1" in result.stdout
    assert "Matched motifs: 1" in result.stdout
    assert "Weight:         6" in result.stdout


def test_scan_accepts_complete_library_dir(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    library_dir = tmp_path / "library"
    library_dir.mkdir()
    _write_scan_artifacts(library_dir)

    result = CliRunner().invoke(
        app,
        ["scan", "--seq", "XXAAXX", "--library-dir", str(library_dir)],
    )

    assert result.exit_code == 0
    assert "Library:       user-specified" in result.stdout
    assert "1. GO:0004672;GO:0005524" in result.stdout


def test_scan_rejects_partial_current_directory_library(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "prosig_motifs.tsv").write_text(
        "name\tdescription\tprosig_pattern\tstatus\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(app, ["scan", "--seq", "XXAAXX"])

    assert result.exit_code != 0
    assert "motif_cluster_scoreboard.pkl" in result.output


def test_scan_requires_exactly_one_query_input() -> None:
    result = CliRunner().invoke(app, ["scan"])

    assert result.exit_code != 0
    assert "provide exactly one of --seq, --fasta" in result.output


def test_build_library_help_includes_options() -> None:
    result = CliRunner().invoke(app, ["build-library", "-h"])

    assert result.exit_code == 0
    assert "Build the compact GO graph" in result.stdout
    assert "--go-obo" in result.stdout
    assert "--swissprot" in result.stdout
    assert "--go-out" in result.stdout
    assert "--prosite-dat" in result.stdout
    assert "--motif-out" in result.stdout
    assert "--write-report" in result.stdout
    assert "--role-map" in result.stdout
    assert "--leiden-cluster-out" in result.stdout
    assert "--cluster-out" in result.stdout
    assert "--cluster-config" in result.stdout
    assert "--min-cluster-similarity" in result.stdout
    assert "--package" in result.stdout
    assert "--package-dir" in result.stdout
    assert "--cluster-neighbors" not in result.stdout
    assert "--cluster-resolution" not in result.stdout
    assert "--cluster-stats-out" not in result.stdout
    assert "--cluster-progress-interval" not in result.stdout
    assert "--cluster-term-cache" not in result.stdout
    assert "--cluster-profile-cache" not in result.stdout
    assert "--cluster-min-informative-ic" not in result.stdout
    assert "--cluster-max-posting-fraction" not in result.stdout
    assert "--cluster-max-posting-size" not in result.stdout
    assert "--force" in result.stdout
    assert "-f" in result.stdout
    assert "--namespace" not in result.stdout
    assert "--include-part-of" not in result.stdout
    assert "--ic-log-base" not in result.stdout
    assert "--min-count" not in result.stdout


def test_log_level_option_suppresses_info_logs() -> None:
    result = CliRunner().invoke(
        app, ["--log-level", "WARNING", "setup-data", "--dry-run"]
    )

    assert result.exit_code == 0
    assert "[INFO]:" not in result.output


def test_build_library_rejects_invalid_min_cluster_similarity() -> None:
    result = CliRunner().invoke(
        app,
        ["build-library", "--min-cluster-similarity", "0"],
    )

    assert result.exit_code != 0
    assert "--min-cluster-similarity" in result.output


def _write_scan_artifacts(tmp_path) -> None:
    (tmp_path / "prosig_motifs.tsv").write_text(
        "name\tdescription\tprosig_pattern\tstatus\n"
        "MOTIF_A\tAA motif\tAA\tprosig\n"
        "MOTIF_B\tBB motif\tBB\tprosig\n",
        encoding="utf-8",
    )
    cluster_rows = [
        "cluster_id\tsim_ave\tsim_min\tsim_max\tsize\tcomposed_go\t"
        "composed_description\n",
        "cluster_0001\tNA\tNA\tNA\t10\tGO:0004672;GO:0005524\t"
        "ATP-binding protein kinase\n",
    ]
    cluster_rows.extend(
        f"cluster_000{i}\tNA\tNA\tNA\t10\tGO:00000{i}\tFunction {i}\n"
        for i in range(2, 7)
    )
    (tmp_path / "clusters_meta.tsv").write_text(
        "".join(cluster_rows),
        encoding="utf-8",
    )
    scoreboard = {
        "schema_version": "1.0",
        "kind": "motif_cluster_scoreboard",
        "parameters": {},
        "weights": {
            "MOTIF_A": {
                f"cluster_000{i}": {
                    "motif_id": "MOTIF_A",
                    "cluster_id": f"cluster_000{i}",
                    "weight": 6.6 - i / 10,
                }
                for i in range(1, 7)
            },
            "MOTIF_B": {
                "cluster_0001": {
                    "motif_id": "MOTIF_B",
                    "cluster_id": "cluster_0001",
                    "weight": 6.0,
                }
            },
        },
    }
    with (tmp_path / "motif_cluster_scoreboard.pkl").open("wb") as handle:
        pickle.dump(scoreboard, handle)
    (tmp_path / "motif_cluster_scoreboard_meta.json").write_text(
        json.dumps(
            {
                "stats": {
                    "calibration": [
                        {
                            "weight_threshold": 5.0,
                            "set_accuracy": 0.91,
                            "top1_accuracy": 0.8,
                            "top3_accuracy": 0.9,
                            "coverage": 0.4,
                        }
                    ]
                }
            }
        ),
        encoding="utf-8",
    )
    with (tmp_path / "go_graph.pkl").open("wb") as handle:
        pickle.dump({"meta": {"namespace": "molecular_function"}, "terms": {}}, handle)
    (tmp_path / "accession_mf_go.tsv").write_text("", encoding="utf-8")
