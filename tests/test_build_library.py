import gzip
import json
import logging
import math
import os
import pickle
from importlib.resources import files
from pathlib import Path

from typer.testing import CliRunner

from prosig.cli.app import app
from prosig.go.build import (
    assign_semantic_roles,
    build_go_pkl,
    parse_role_map,
    parse_swissprot_entry,
    write_accession_mf_go_tsv,
)


def test_parse_swissprot_entry_keeps_primary_accession_and_high_quality_mf() -> None:
    accession, terms = parse_swissprot_entry(
        [
            "AC   P31946; Q53XZ2; Q96QU6;",
            "DR   GO; GO:0005634; C:nucleus; IDA:UniProtKB.",
            "DR   GO; GO:0005515; F:protein binding; IPI:UniProtKB.",
            "DR   GO; GO:0005524; F:ATP binding; IEA:InterPro.",
            "DR   GO; GO:0004672; F:protein kinase activity; IKR:UniProtKB.",
            "DR   GO; GO:0003674; F:molecular_function; ND:UniProtKB.",
            "DR   GO; GO:0005488; F:binding; NAS:UniProtKB.",
        ]
    )

    assert accession == "P31946"
    assert terms == {"GO:0005515", "GO:0005524", "GO:0004672"}


def test_build_go_pkl_keeps_mf_terms_and_computes_propagated_ic(
    tmp_path: Path, caplog
) -> None:
    go_obo = tmp_path / "go-basic.obo"
    swissprot = tmp_path / "uniprot_sprot.dat.gz"
    go_out = tmp_path / "go_graph.pkl"
    report_out = tmp_path / "go_report.txt"
    go_json_out = tmp_path / "go_graph.json"
    excluded_mf_annotations_out = tmp_path / "excluded_mf_annotations.tsv"
    accession_mf_go_out = tmp_path / "accession_mf_go.tsv"
    go_obo.write_text(_small_obo(), encoding="utf-8")
    _write_gzip(swissprot, _small_swissprot())
    caplog.set_level(logging.INFO, logger="prosig")

    artifact = build_go_pkl(
        go_obo=go_obo,
        swissprot=swissprot,
        go_out=go_out,
        report_out=report_out,
    )

    assert go_out.exists()
    with go_out.open("rb") as handle:
        saved = pickle.load(handle)
    assert saved == artifact

    terms = artifact["terms"]
    assert set(terms) == {
        "GO:0003674",
        "GO:0000001",
        "GO:0000002",
        "GO:0000003",
        "GO:0000005",
    }
    assert terms["GO:0003674"]["ic"] == 0.0
    assert terms["GO:0003674"]["freq"] == 1.0
    assert terms["GO:0003674"]["count"] == 3
    assert terms["GO:0000001"]["count"] == 3
    assert terms["GO:0000002"]["count"] == 1
    assert terms["GO:0000003"]["count"] == 2
    assert terms["GO:0000005"]["count"] == 0
    assert terms["GO:0000005"]["freq"] == 0.0
    assert terms["GO:0000005"]["ic"] is None
    assert math.isclose(terms["GO:0000002"]["freq"], 1 / 3)
    assert math.isclose(terms["GO:0000002"]["ic"], -math.log(1 / 3))
    assert terms["GO:0000002"]["ancestors"] == {"GO:0000001", "GO:0003674"}
    assert terms["GO:0000002"]["depth"] == 2
    assert terms["GO:0000001"]["children"] == [
        "GO:0000002",
        "GO:0000003",
        "GO:0000005",
    ]

    meta = artifact["meta"]
    assert meta["namespace"] == "molecular_function"
    assert meta["n_terms"] == 5
    assert meta["n_accessions_provided"] == 6
    assert meta["n_accessions_with_hq_mf_go"] == 5
    assert meta["n_accessions_with_any_mf_go"] == 5
    assert meta["n_hq_mf_go_assignments_not_in_graph"] == 2
    assert meta["n_hq_mf_go_assignments_obsolete"] == 1
    assert meta["mf_frequency_min"] == 1
    assert meta["mf_frequency_median"] == 1
    assert meta["mf_frequency_mean"] == 1.25
    assert meta["mf_frequency_max"] == 2
    assert meta["mf_frequency_status"] == "OK"
    assert "accession_to_terms" not in artifact
    removed_meta_fields = {
        "n_accessions_with_any_bp_go",
        "n_accessions_with_hq_bp_go",
        "n_accessions_with_any_cc_go",
        "n_accessions_with_hq_cc_go",
        "n_accessions_used_for_ic",
        "total_unique_accessions",
        "mf_accessions_with_hq_direct_go",
        "mf_accessions_without_hq_direct_go",
        "mf_unique_direct_go_terms",
        "mf_total_direct_go_assignments",
        "n_accessions_skipped_no_valid_mf",
        "source_checksums",
    }
    assert not removed_meta_fields & set(meta)
    report_text = report_out.read_text(encoding="utf-8")
    assert "number of accessions used for IC: 3" in report_text
    assert (
        "number of HQ MF GO assignments skipped because the GO term is obsolete: 1"
        in report_text
    )
    go_json_text = go_json_out.read_text(encoding="utf-8")
    assert go_json_text.startswith('{\n  "_comment": "Diagnostic only.')
    go_json = json.loads(go_json_text)
    assert go_json["_comment"] == (
        "Diagnostic only. Use go_graph.pkl as the runtime artifact."
    )
    assert go_json["meta"] == meta
    assert go_json["terms"]["GO:0000002"]["ancestors"] == [
        "GO:0000001",
        "GO:0003674",
    ]
    assert not (tmp_path / "go_frequency_metadata.tsv").exists()
    assert excluded_mf_annotations_out.read_text(encoding="utf-8") == (
        "accession\tgo_term\tevidence\n"
        "P00001\tGO:0000005\tNAS\n"
    )
    assert accession_mf_go_out.read_text(encoding="utf-8") == (
        "P00001\tGO:0000002\n"
        "P00002\tGO:0000003\n"
        "P00004\tGO:9999999\n"
        "P00005\tGO:0000003\n"
        "P00006\tGO:0000004\n"
    )
    expected_log_fragments = [
        "Parsing GO OBO file",
        "Parsed 5 connected Molecular Function GO terms",
        "Parsing Swiss-Prot MF GO annotations",
        "Parsed 6 primary accessions",
        "GO annotation accession summary: total=6; MF=5; MF high-quality=5; "
        "BP=2; BP high-quality=2; CC=2; CC high-quality=2",
        "Propagating MF GO annotations",
        "Calculated IC values using 3 accessions; skipped 2 HQ MF GO assignments",
        "Skipped 1 HQ MF GO assignments because the GO term is obsolete",
        "1 GO terms did not receive valid IC because no accession matched them",
        "Top 10 most frequent MF GO terms",
        "1. GO:0003674 molecular_function count=3 freq=1 ic=0",
        "3. GO:0000003 sibling activity count=2 freq=0.6667 ic=0.4055",
        "Writing GO graph and IC artifact",
        "Wrote GO graph and IC artifact",
        "Writing GO build validation report",
        "Wrote GO build validation report",
        "Writing diagnostic GO graph JSON",
        "Wrote diagnostic GO graph JSON",
        "Writing excluded MF annotation diagnostics",
        "Wrote excluded MF annotation diagnostics",
        "Writing accession MF GO terms",
        "Wrote accession MF GO terms",
    ]
    for expected in expected_log_fragments:
        assert any(expected in message for message in caplog.messages)
    assert not [
        record for record in caplog.records if record.levelno >= logging.WARNING
    ]


def test_build_library_command_writes_go_graph_pkl(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    Path("go-basic.obo").write_text(_small_obo(), encoding="utf-8")
    _write_gzip(Path("uniprot_sprot.dat.gz"), _small_swissprot())
    Path("prosite.dat").write_text(_small_prosite(), encoding="utf-8")

    result = CliRunner().invoke(app, ["build-library", "--write-report", "report.txt"])

    assert result.exit_code == 0
    assert Path("go_graph.pkl").exists()
    assert Path("report.txt").exists()
    assert Path("go_graph.json").exists()
    assert not Path("go_frequency_metadata.tsv").exists()
    assert Path("excluded_mf_annotations.tsv").exists()
    assert Path("accession_mf_go.tsv").exists()
    assert Path("go_clusters.tsv").exists()
    assert Path("go_clusters_stats.json").exists()
    assert Path("go_clusters_meta.tsv").exists()
    assert Path("cluster_config.yaml").exists()
    assert Path("role_map.yaml").exists()
    assert Path("go_terms_unknown_role.txt").exists()
    with Path("go_graph.pkl").open("rb") as handle:
        artifact = pickle.load(handle)
    assert artifact["terms"]["GO:0000002"]["semantic_role"] == {
        "role": "unknown",
        "priority": 0,
        "source": "unknown",
        "matched": None,
    }
    assert artifact["meta"]["semantic_role_assignment"] == {
        "role_map": "role_map.yaml",
        "unknown_role_report": "go_terms_unknown_role.txt",
        "n_processed": 4,
        "n_anchor": 0,
        "n_keyword": 0,
        "n_unknown": 4,
        "role_counts": {"unknown": 4},
    }
    assert Path("prosig_motifs.tsv").read_text(encoding="utf-8") == (
        "# ProSig motif library\n"
        "name\tprosite_ac\tdescription\tprosite_pattern\tprosig_pattern\tstatus\n"
        "N_GLYCOSYLATION\tPS00001\tN-glycosylation site\t"
        "N-{P}-[ST]-{P}\tN!P[ST]!P\tprosite\n"
    )


def test_build_library_skips_current_derived_artifacts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    Path("go-basic.obo").write_text(_small_obo(), encoding="utf-8")
    _write_gzip(Path("uniprot_sprot.dat.gz"), _small_swissprot())
    Path("prosite.dat").write_text(_small_prosite(), encoding="utf-8")

    first_result = CliRunner().invoke(
        app,
        ["build-library", "--write-report", "report.txt"],
    )
    assert first_result.exit_code == 0
    mtimes = {
        path: Path(path).stat().st_mtime_ns
        for path in [
            "go_graph.pkl",
            "accession_mf_go.tsv",
            "go_graph.json",
            "excluded_mf_annotations.tsv",
            "prosig_motifs.tsv",
            "report.txt",
            "go_clusters.tsv",
            "go_clusters_stats.json",
            "go_clusters_meta.tsv",
            "cluster_config.yaml",
        ]
    }

    second_result = CliRunner().invoke(
        app,
        ["build-library", "--write-report", "report.txt"],
    )

    assert second_result.exit_code == 0
    assert {
        path: Path(path).stat().st_mtime_ns for path in mtimes
    } == mtimes
    expected_log_fragments = [
        "Skipping GO graph build",
        "Skipping accession MF GO terms",
        "Skipping ProSig motif library build",
        "Skipping GO clustering",
    ]
    for expected in expected_log_fragments:
        assert expected in second_result.output


def test_build_library_force_rebuilds_current_artifacts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    Path("go-basic.obo").write_text(_small_obo(), encoding="utf-8")
    _write_gzip(Path("uniprot_sprot.dat.gz"), _small_swissprot())
    Path("prosite.dat").write_text(_small_prosite(), encoding="utf-8")

    first_result = CliRunner().invoke(app, ["build-library"])
    assert first_result.exit_code == 0
    os.utime("go-basic.obo", (100, 100))
    os.utime("uniprot_sprot.dat.gz", (100, 100))
    os.utime("prosite.dat", (100, 100))
    os.utime("role_map.yaml", (100, 100))
    os.utime("cluster_config.yaml", (100, 100))
    os.utime("go_graph.pkl", (200, 200))
    os.utime("accession_mf_go.tsv", (200, 200))
    os.utime("prosig_motifs.tsv", (200, 200))
    os.utime("go_clusters.tsv", (200, 200))
    os.utime("go_clusters_stats.json", (200, 200))
    os.utime("go_clusters_meta.tsv", (200, 200))

    forced_result = CliRunner().invoke(app, ["build-library", "-f"])

    assert forced_result.exit_code == 0
    assert Path("go_graph.pkl").stat().st_mtime > 200
    assert Path("accession_mf_go.tsv").stat().st_mtime > 200
    assert Path("prosig_motifs.tsv").stat().st_mtime > 200
    assert Path("go_clusters.tsv").stat().st_mtime > 200
    assert Path("go_clusters_stats.json").stat().st_mtime > 200
    assert Path("go_clusters_meta.tsv").stat().st_mtime > 200
    assert Path("cluster_config.yaml").stat().st_mtime == 100


def test_build_go_pkl_logs_semantic_role_assignment(
    tmp_path: Path,
    caplog,
) -> None:
    go_obo = tmp_path / "go-basic.obo"
    swissprot = tmp_path / "uniprot_sprot.dat.gz"
    go_out = tmp_path / "go_graph.pkl"
    role_map = tmp_path / "role_map.yaml"
    go_obo.write_text(_small_obo(), encoding="utf-8")
    _write_gzip(swissprot, _small_swissprot())
    role_map.write_text(
        "\n".join(
            [
                "roles:",
                "  catalytic:",
                "    priority: 100",
                "    anchors:",
                "      - GO:0000002",
                "role_rules:",
                "  binding_generic:",
                "    priority: 20",
                "    keywords:",
                '      - "binding"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    caplog.set_level(logging.INFO, logger="prosig")

    artifact = build_go_pkl(
        go_obo=go_obo,
        swissprot=swissprot,
        go_out=go_out,
        role_map=role_map,
    )

    assert artifact["meta"]["semantic_role_assignment"]["n_processed"] == 4
    assert artifact["meta"]["semantic_role_assignment"]["role_counts"] == {
        "catalytic": 1,
        "unknown": 3,
    }
    expected_log_fragments = [
        f"Loading GO semantic role map: {role_map}",
        "Assigning GO semantic roles to 4 non-root GO terms",
        "Applying Layer 1 GO anchor/ancestor role matching",
        "Applying Layer 2 keyword role matching to remaining terms",
        "Processed 4 GO terms for semantic role assignment",
        "GO semantic role layer summary:\n"
        "  total non-root terms =       4\n"
        "  anchor assigned      =       1\n"
        "  keyword assigned     =       0\n"
        "  unknown              =       3",
        "GO semantic role stats:\n"
        "  catalytic = 1\n"
        "  unknown = 3",
    ]
    for expected in expected_log_fragments:
        assert any(expected in message for message in caplog.messages)


def test_assign_semantic_roles_uses_anchors_then_keywords(tmp_path: Path) -> None:
    terms = {
        "GO:0003674": {
            "name": "molecular_function",
            "ancestors": set(),
            "ic": 0.0,
        },
        "GO:0000001": {
            "name": "catalytic activity",
            "ancestors": {"GO:0003674"},
            "ic": 1.0,
        },
        "GO:0000002": {
            "name": "child activity",
            "ancestors": {"GO:0000001", "GO:0003674"},
            "ic": 2.0,
        },
        "GO:0000003": {
            "name": "protein binding",
            "ancestors": {"GO:0003674"},
            "ic": 3.0,
        },
    }
    role_map = tmp_path / "role_map.yaml"
    role_map.write_text(
        "\n".join(
            [
                "roles:",
                "  catalytic:",
                "    priority: 100",
                "    anchors:",
                "      - GO:0000001",
                "role_rules:",
                "  binding_generic:",
                "    priority: 20",
                "    keywords:",
                '      - "binding"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    stats = assign_semantic_roles(
        terms,
        parse_role_map(role_map),
        unknown_role_out=tmp_path / "unknown.txt",
    )

    assert stats["processed"] == 3
    assert stats["anchor"] == 2
    assert stats["keyword"] == 1
    assert stats["unknown"] == 0
    assert stats["role_counts"] == {"catalytic": 2, "binding_generic": 1}
    assert terms["GO:0000002"]["semantic_role"] == {
        "role": "catalytic",
        "priority": 100,
        "source": "anchor",
        "matched": "GO:0000001",
    }
    assert terms["GO:0000003"]["semantic_role"] == {
        "role": "binding_generic",
        "priority": 20,
        "source": "keyword",
        "matched": "binding",
    }
    assert (tmp_path / "unknown.txt").read_text(encoding="utf-8") == ""


def test_assign_semantic_roles_refines_broad_binding_anchor(
    tmp_path: Path,
) -> None:
    terms = {
        "GO:0003674": {
            "name": "molecular_function",
            "ancestors": set(),
            "ic": 0.0,
        },
        "GO:0005488": {
            "name": "binding",
            "ancestors": {"GO:0003674"},
            "ic": 1.0,
        },
        "GO:0000001": {
            "name": "heme binding",
            "ancestors": {"GO:0005488", "GO:0003674"},
            "ic": 2.0,
        },
        "GO:0000002": {
            "name": "protein binding",
            "ancestors": {"GO:0005488", "GO:0003674"},
            "ic": 3.0,
        },
    }
    role_map = tmp_path / "role_map.yaml"
    role_map.write_text(
        "\n".join(
            [
                "roles:",
                "  binding:",
                "    priority: 20",
                "    anchors:",
                "      - GO:0005488",
                "  binding_cofactor:",
                "    priority: 40",
                "    keywords:",
                '      - "heme binding"',
                "  binding_generic:",
                "    priority: 20",
                "    keywords:",
                '      - "protein binding"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    stats = assign_semantic_roles(
        terms,
        parse_role_map(role_map),
        unknown_role_out=tmp_path / "unknown.txt",
    )

    assert stats["anchor"] == 1
    assert stats["keyword"] == 2
    assert terms["GO:0005488"]["semantic_role"] == {
        "role": "binding",
        "priority": 20,
        "source": "anchor",
        "matched": "GO:0005488",
    }
    assert terms["GO:0000001"]["semantic_role"] == {
        "role": "binding_cofactor",
        "priority": 40,
        "source": "keyword",
        "matched": "heme binding",
    }
    assert terms["GO:0000002"]["semantic_role"] == {
        "role": "binding_generic",
        "priority": 20,
        "source": "keyword",
        "matched": "protein binding",
    }


def test_assign_semantic_roles_accepts_grouped_role_map(tmp_path: Path) -> None:
    terms = {
        "GO:0003674": {
            "name": "molecular_function",
            "ancestors": set(),
            "ic": 0.0,
        },
        "GO:0000001": {
            "name": "receptor activity",
            "ancestors": {"GO:0003674"},
            "ic": 1.0,
        },
        "GO:0000002": {
            "name": "signaling receptor activity",
            "ancestors": {"GO:0003674"},
            "ic": 2.0,
        },
    }
    role_map = tmp_path / "role_map.yaml"
    role_map.write_text(
        "\n".join(
            [
                "roles:",
                "  receptor:",
                "    priority: 80",
                "    anchors:",
                "      - GO:0000001",
                "    keywords:",
                '      - "signaling receptor activity"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    stats = assign_semantic_roles(
        terms,
        parse_role_map(role_map),
        unknown_role_out=tmp_path / "unknown.txt",
    )

    assert stats["processed"] == 2
    assert stats["anchor"] == 1
    assert stats["keyword"] == 1
    assert stats["unknown"] == 0
    assert terms["GO:0000001"]["semantic_role"]["source"] == "anchor"
    assert terms["GO:0000002"]["semantic_role"] == {
        "role": "receptor",
        "priority": 80,
        "source": "keyword",
        "matched": "signaling receptor activity",
    }


def test_parse_role_map_strips_inline_comments_from_list_items(tmp_path: Path) -> None:
    role_map = tmp_path / "role_map.yaml"
    role_map.write_text(
        "\n".join(
            [
                "roles:",
                "  tagging:",
                "    priority: 68",
                "    anchors:",
                "      - GO:0141047   # molecular tag activity",
                "    keywords:",
                '      - "protein # tag activity"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    parsed = parse_role_map(role_map)

    assert parsed["roles"]["tagging"]["anchors"] == ["GO:0141047"]
    assert parsed["roles"]["tagging"]["keywords"] == ["protein # tag activity"]


def test_parse_role_map_strips_inline_comments_from_mapping_fields(
    tmp_path: Path,
) -> None:
    role_map = tmp_path / "role_map.yaml"
    role_map.write_text(
        "\n".join(
            [
                "roles:",
                "  catalytic:",
                "    priority: 100  # highest priority",
                "    anchors:  # anchor GO terms",
                "      - GO:0003824",
                "    keywords:  # keyword fallback",
                '      - "catalytic # activity"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    parsed = parse_role_map(role_map)

    assert parsed["roles"]["catalytic"] == {
        "priority": 100,
        "anchors": ["GO:0003824"],
        "keywords": ["catalytic # activity"],
    }


def test_packaged_role_map_template_is_available() -> None:
    template = files("prosig.data").joinpath("role_map.yaml.template")

    assert template.is_file()
    assert "roles:" in template.read_text(encoding="utf-8")


def test_write_accession_mf_go_tsv_uses_primary_accessions_and_hq_mf_terms(
    tmp_path: Path,
) -> None:
    swissprot = tmp_path / "uniprot_sprot.dat.gz"
    output = tmp_path / "accession_mf_go.tsv"
    _write_gzip(swissprot, _small_swissprot())

    write_accession_mf_go_tsv(output, swissprot)

    assert output.read_text(encoding="utf-8") == (
        "P00001\tGO:0000002\n"
        "P00002\tGO:0000003\n"
        "P00004\tGO:9999999\n"
        "P00005\tGO:0000003\n"
        "P00006\tGO:0000004\n"
    )


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _small_obo() -> str:
    return """format-version: 1.2

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0000001
name: parent activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000002
name: child activity
namespace: molecular_function
is_a: GO:0000001 ! parent activity

[Term]
id: GO:0000003
name: sibling activity
namespace: molecular_function
is_a: GO:0000001 ! parent activity

[Term]
id: GO:0000005
name: unused activity
namespace: molecular_function
is_a: GO:0000001 ! parent activity

[Term]
id: GO:0000004
name: obsolete activity
namespace: molecular_function
is_obsolete: true
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process
"""


def _small_prosite() -> str:
    return """ID   N_GLYCOSYLATION; PATTERN.
AC   PS00001;
DE   N-glycosylation site.
PA   N-{P}-[ST]-{P}.
//
"""


def _small_swissprot() -> str:
    return """ID   TEST1                  Reviewed;         10 AA.
AC   P00001; P00001-2;
DR   GO; GO:0000002; F:child activity; EXP:UniProtKB.
DR   GO; GO:0000005; F:unused activity; NAS:UniProtKB.
DR   GO; GO:0008150; P:biological_process; NAS:UniProtKB.
DR   GO; GO:0005575; C:cellular_component; NAS:UniProtKB.
DR   GO; GO:0008150; P:biological_process; EXP:UniProtKB.
DR   GO; GO:0005575; C:cellular_component; IDA:UniProtKB.
//
ID   TEST2                  Reviewed;         10 AA.
AC   P00002;
DR   GO; GO:0000003; F:sibling activity; IEA:InterPro.
DR   GO; GO:0005575; C:cellular_component; IEA:UniProtKB.
//
ID   TEST3                  Reviewed;         10 AA.
AC   P00003;
DR   GO; GO:0008150; P:biological_process; EXP:UniProtKB.
//
ID   TEST4                  Reviewed;         10 AA.
AC   P00004;
DR   GO; GO:9999999; F:missing activity; IDA:UniProtKB.
//
ID   TEST5                  Reviewed;         10 AA.
AC   P00005;
DR   GO; GO:0000003; F:sibling activity; IMP:UniProtKB.
//
ID   TEST6                  Reviewed;         10 AA.
AC   P00006;
DR   GO; GO:0000004; F:obsolete activity; IDA:UniProtKB.
//
"""
