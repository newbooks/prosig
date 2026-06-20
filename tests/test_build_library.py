import gzip
import json
import logging
import math
import pickle
from pathlib import Path

from typer.testing import CliRunner

from prosig.cli.app import app
from prosig.go.build import build_go_pkl, parse_swissprot_entry


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

    result = CliRunner().invoke(app, ["build-library", "--write-report", "report.txt"])

    assert result.exit_code == 0
    assert Path("go_graph.pkl").exists()
    assert Path("report.txt").exists()
    assert Path("go_graph.json").exists()
    assert not Path("go_frequency_metadata.tsv").exists()
    assert Path("excluded_mf_annotations.tsv").exists()


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
