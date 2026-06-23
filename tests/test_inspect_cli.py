import json
import pickle

from typer.testing import CliRunner

from prosig.cli.app import app
from prosig.go.similarity import GoSimilarity


def _write_go_graph(path) -> None:
    artifact = {
        "meta": {
            "schema_version": "1.0",
            "namespace": "molecular_function",
            "created_at": "2026-06-21",
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
                "freq": 0.5,
                "ic": 1.0,
            },
            "GO:0000002": {
                "name": "child A activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": 0.25,
                "ic": 2.0,
            },
            "GO:0000003": {
                "name": "child B activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": 0.125,
                "ic": 3.0,
            },
        },
    }
    with path.open("wb") as handle:
        pickle.dump(artifact, handle)


def _write_function_go_graph(path) -> None:
    artifact = {
        "meta": {
            "schema_version": "1.0",
            "namespace": "molecular_function",
            "created_at": "2026-06-21",
        },
        "terms": {
            "GO:0003674": {
                "name": "molecular_function",
                "parents": [],
                "children": [
                    "GO:0016740",
                    "GO:0005488",
                    "GO:0044183",
                    "GO:0003700",
                    "GO:0004497",
                    "GO:0016705",
                ],
                "ancestors": set(),
                "depth": 0,
                "freq": 1.0,
                "ic": 0.0,
            },
            "GO:0016740": {
                "name": "transferase activity",
                "parents": ["GO:0003674"],
                "children": ["GO:0004672"],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": 0.4,
                "ic": 1.0,
                "semantic_role": {
                    "role": "catalytic",
                    "priority": 100,
                    "source": "anchor",
                    "matched": "GO:0003824",
                },
            },
            "GO:0004672": {
                "name": "protein kinase activity",
                "parents": ["GO:0016740"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0016740"},
                "depth": 2,
                "freq": 0.1,
                "ic": 4.0,
                "semantic_role": {
                    "role": "catalytic",
                    "priority": 100,
                    "source": "anchor",
                    "matched": "GO:0003824",
                },
            },
            "GO:0005488": {
                "name": "binding",
                "parents": ["GO:0003674"],
                "children": ["GO:0005524", "GO:0000287"],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": 0.5,
                "ic": 0.5,
                "semantic_role": {
                    "role": "binding",
                    "priority": 10,
                    "source": "anchor",
                    "matched": "GO:0005488",
                },
            },
            "GO:0005524": {
                "name": "ATP binding",
                "parents": ["GO:0005488"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0005488"},
                "depth": 2,
                "freq": 0.2,
                "ic": 3.0,
                "semantic_role": {
                    "role": "binding_cofactor",
                    "priority": 40,
                    "source": "keyword",
                    "matched": "ATP",
                },
            },
            "GO:0000287": {
                "name": "magnesium ion binding",
                "parents": ["GO:0005488"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0005488"},
                "depth": 2,
                "freq": 0.25,
                "ic": 2.5,
                "semantic_role": {
                    "role": "binding_cofactor",
                    "priority": 40,
                    "source": "keyword",
                    "matched": "magnesium",
                },
            },
            "GO:0044183": {
                "name": "protein folding chaperone",
                "parents": ["GO:0003674"],
                "children": [],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": 0.15,
                "ic": 2.0,
                "semantic_role": {
                    "role": "chaperone",
                    "priority": 75,
                    "source": "anchor",
                    "matched": "GO:0044183",
                },
            },
            "GO:0003700": {
                "name": "DNA-binding transcription factor activity",
                "parents": ["GO:0003674"],
                "children": ["GO:0043565"],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": 0.12,
                "ic": 3.0,
                "semantic_role": {
                    "role": "transcription_factor",
                    "priority": 85,
                    "source": "anchor",
                    "matched": "GO:0003700",
                },
            },
            "GO:0043565": {
                "name": "sequence-specific DNA binding",
                "parents": ["GO:0003700"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0003700"},
                "depth": 2,
                "freq": 0.08,
                "ic": 3.5,
                "semantic_role": {
                    "role": "binding_nucleic_acid",
                    "priority": 45,
                    "source": "keyword",
                    "matched": "DNA binding",
                },
            },
            "GO:0004497": {
                "name": "monooxygenase activity",
                "parents": ["GO:0003674"],
                "children": [],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": 0.07,
                "ic": 5.0,
                "semantic_role": {
                    "role": "catalytic",
                    "priority": 100,
                    "source": "keyword",
                    "matched": "monooxygenase activity",
                },
            },
            "GO:0016705": {
                "name": "oxidoreductase activity",
                "parents": ["GO:0003674"],
                "children": [],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": 0.09,
                "ic": 4.0,
                "semantic_role": {
                    "role": "catalytic",
                    "priority": 100,
                    "source": "keyword",
                    "matched": "oxidoreductase activity",
                },
            },
        },
    }
    with path.open("wb") as handle:
        pickle.dump(artifact, handle)


def test_inspect_help_lists_diagnostic_commands() -> None:
    result = CliRunner().invoke(app, ["inspect", "-h"])

    assert result.exit_code == 0
    assert "go-summary" in result.stdout
    assert "go-term" in result.stdout
    assert "go-sim" in result.stdout
    assert "go-set-sim" in result.stdout
    assert "function" in result.stdout
    assert "go-similarity" not in result.stdout


def test_inspect_go_term_outputs_term_details(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-term",
            "GO:0000002",
            "--go-graph",
            str(go_graph),
            "--ancestors",
        ],
    )

    assert result.exit_code == 0
    assert "name\tchild A activity" in result.stdout
    assert "ic\t2.0" in result.stdout
    assert "GO:0000001,GO:0000002,GO:0003674" in result.stdout


def test_inspect_go_sim_outputs_simple_score(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-sim",
            "GO:0000002",
            "GO:0000003",
            "--go-graph",
            str(go_graph),
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == "0.4\n"


def test_inspect_go_sim_simple_output_does_not_build_details(
    tmp_path,
    monkeypatch,
) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("simple go-sim should use scalar lin()")

    monkeypatch.setattr(GoSimilarity, "lin_with_details", fail_if_called)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-sim",
            "GO:0000002",
            "GO:0000003",
            "--go-graph",
            str(go_graph),
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == "0.4\n"


def test_inspect_go_sim_verbose_explains_lin_score(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-sim",
            "GO:0000002",
            "GO:0000003",
            "--go-graph",
            str(go_graph),
            "--verbose",
        ],
    )

    assert result.exit_code == 0
    assert "Score: 0.4" in result.stdout
    assert "GO:0000002 child A activity" in result.stdout
    assert "GO:0000003 child B activity" in result.stdout
    assert "GO:0000001 parent activity" in result.stdout
    assert "Common ancestors:" in result.stdout
    assert "GO:0000001 parent activity (IC=1, freq=0.5, depth=1)" in result.stdout
    assert result.stdout.index("GO:0003674 molecular_function") < result.stdout.index(
        "GO:0000001 parent activity"
    )
    assert "Compact GO path:" in result.stdout
    assert "└── GO:0003674 molecular_function [common]" in result.stdout
    assert "    └── GO:0000001 parent activity [MICA]" in result.stdout
    assert "        ├── GO:0000002 child A activity [A]" in result.stdout
    assert "        └── GO:0000003 child B activity [B]" in result.stdout
    assert "MICA: GO:0000001 parent activity" in result.stdout
    assert "Formula: Lin(A, B) = 2 * IC(MICA) / (IC(A) + IC(B))" in result.stdout
    assert "                   = 2 * 1 / (2 + 3) " in result.stdout
    assert "                   = 0.4" in result.stdout


def test_inspect_go_sim_short_verbose_option_works(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-sim",
            "GO:0000002",
            "GO:0000003",
            "--go-graph",
            str(go_graph),
            "-v",
        ],
    )

    assert result.exit_code == 0
    assert "Score: 0.4" in result.stdout


def test_inspect_go_sim_ascii_tree_style_works(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-sim",
            "GO:0000002",
            "GO:0000003",
            "--go-graph",
            str(go_graph),
            "-v",
            "--tree-style",
            "ascii",
        ],
    )

    assert result.exit_code == 0
    assert "`- GO:0003674 molecular_function [common]" in result.stdout
    assert "   `- GO:0000001 parent activity [MICA]" in result.stdout
    assert "      |- GO:0000002 child A activity [A]" in result.stdout
    assert "      `- GO:0000003 child B activity [B]" in result.stdout


def test_inspect_go_set_sim_help_mentions_quoting() -> None:
    result = CliRunner().invoke(app, ["inspect", "go-set-sim", "-h"])

    assert result.exit_code == 0
    assert "Quote GO set expressions" in result.stdout


def test_inspect_go_set_sim_outputs_simple_score_for_direct_sets(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-set-sim",
            "(GO:0000002;GO:0000003)",
            "(GO:0000002)",
            "--go-graph",
            str(go_graph),
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == "0.85\n"


def test_inspect_go_set_sim_allows_mixed_set_and_accession(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    accession_go = tmp_path / "accession_mf_go.tsv"
    _write_go_graph(go_graph)
    accession_go.write_text("P00001\tGO:0000002\n", encoding="utf-8")

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-set-sim",
            "(GO:0000002;GO:0000003)",
            "P00001",
            "--go-graph",
            str(go_graph),
            "--accession-go",
            str(accession_go),
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == "0.85\n"


def test_inspect_go_set_sim_accepts_quoted_set_without_parentheses(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-set-sim",
            "GO:0000002;GO:0000003",
            "GO:0000002",
            "--go-graph",
            str(go_graph),
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == "0.85\n"


def test_inspect_go_set_sim_simple_output_does_not_build_details(
    tmp_path,
    monkeypatch,
) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("simple go-set-sim should use scalar set_lin_amb()")

    monkeypatch.setattr(GoSimilarity, "set_lin_amb_with_details", fail_if_called)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-set-sim",
            "(GO:0000002;GO:0000003)",
            "(GO:0000002)",
            "--go-graph",
            str(go_graph),
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == "0.85\n"


def test_inspect_go_set_sim_verbose_explains_amb_score(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-set-sim",
            "(GO:0000002;GO:0000003)",
            "(GO:0000002)",
            "--go-graph",
            str(go_graph),
            "--verbose",
        ],
    )

    assert result.exit_code == 0
    assert "Score: 0.85" in result.stdout
    assert "A query: (GO:0000002;GO:0000003)" in result.stdout
    assert "B query: (GO:0000002)" in result.stdout
    duplicated_query = "A query: (GO:0000002;GO:0000003) (GO:0000002;GO:0000003)"
    assert duplicated_query not in result.stdout
    assert "GO term descriptions:" in result.stdout
    assert "- A: GO:0000002 child A activity" in result.stdout
    assert "- A: GO:0000003 child B activity" in result.stdout
    assert "IC=2" not in result.stdout
    assert "A -> B best matches:" in result.stdout
    assert "GO:0000003 --0.4000--> GO:0000002" in result.stdout
    assert "GO:0000003 child B activity --0.4-->" not in result.stdout
    assert "A -> B average max: 0.7" in result.stdout
    assert "B -> A average max: 1" in result.stdout
    assert "Formula: AMB(A, B) = (mean(A -> B) + mean(B -> A)) / 2" in result.stdout
    assert "                   = (0.7 + 1) / 2" in result.stdout
    assert "                   = 0.85" in result.stdout


def test_inspect_go_set_sim_verbose_expands_accession_query(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    accession_go = tmp_path / "accession_mf_go.tsv"
    _write_go_graph(go_graph)
    accession_go.write_text("P00001\tGO:0000002;GO:0000003\n", encoding="utf-8")

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "go-set-sim",
            "P00001",
            "(GO:0000002)",
            "--go-graph",
            str(go_graph),
            "--accession-go",
            str(accession_go),
            "--verbose",
        ],
    )

    assert result.exit_code == 0
    assert "A query: P00001 (GO:0000002;GO:0000003)" in result.stdout
    assert "A expanded terms:" not in result.stdout
    assert "B query: (GO:0000002)" in result.stdout


def test_inspect_function_describes_direct_go_set(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_function_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "GO:0004672;GO:0005524;GO:0000287;GO:0016740",
            "--go-graph",
            str(go_graph),
        ],
    )

    assert result.exit_code == 0
    assert (
        result.stdout
        == "GO:0004672;GO:0005524;GO:0000287;GO:0016740 is annotated "
        "as an ATP- and magnesium-binding protein kinase.\n"
    )


def test_inspect_function_describes_accession(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    accession_go = tmp_path / "accession_mf_go.tsv"
    _write_function_go_graph(go_graph)
    accession_go.write_text(
        "P00001\tGO:0004672;GO:0005524;GO:0000287\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "P00001",
            "--go-graph",
            str(go_graph),
            "--accession-go",
            str(accession_go),
        ],
    )

    assert result.exit_code == 0
    assert (
        result.stdout
        == "P00001 is annotated as an ATP- and magnesium-binding protein kinase.\n"
    )


def test_inspect_function_honors_zero_max_modifiers(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_function_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "GO:0004672;GO:0005524;GO:0000287",
            "--go-graph",
            str(go_graph),
            "--max-modifiers",
            "0",
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == (
        "GO:0004672;GO:0005524;GO:0000287 is annotated as "
        "a protein kinase.\n"
    )


def test_inspect_function_excludes_dropped_ancestors_from_support(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_function_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "GO:0004672;GO:0016740",
            "--go-graph",
            str(go_graph),
            "--json",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["head"] == "GO:0004672"
    assert payload["dropped_terms"] == ["GO:0016740"]
    assert payload["supporting_terms"] == []
    assert payload["summary"] == (
        "GO:0004672;GO:0016740 is annotated as a protein kinase."
    )


def test_inspect_function_uses_compiled_role_priority_for_head(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_function_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "GO:0044183;GO:0005524",
            "--go-graph",
            str(go_graph),
            "--json",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["head"] == "GO:0044183"
    assert payload["modifiers"] == ["ATP-binding"]
    assert payload["summary"] == (
        "GO:0044183;GO:0005524 is annotated as an ATP-binding "
        "protein folding chaperone."
    )


def test_inspect_function_keeps_non_head_activity_as_support(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_function_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "GO:0004497;GO:0016705",
            "--go-graph",
            str(go_graph),
            "--json",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["head"] == "GO:0004497"
    assert payload["supporting_terms"] == ["oxidoreductase"]
    assert payload["summary"] == (
        "GO:0004497;GO:0016705 is annotated as a monooxygenase "
        "with oxidoreductase activity."
    )


def test_inspect_function_replaces_binding_prefix_in_head(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_function_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "GO:0003700;GO:0043565",
            "--go-graph",
            str(go_graph),
            "--json",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["head"] == "GO:0003700"
    assert payload["modifiers"] == ["sequence-specific DNA-binding"]
    assert payload["summary"] == (
        "GO:0003700;GO:0043565 is annotated as a sequence-specific "
        "DNA-binding transcription factor."
    )
    assert "DNA-binding DNA-binding" not in payload["summary"]


def test_inspect_function_verbose_shows_resolved_terms(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    accession_go = tmp_path / "accession_mf_go.tsv"
    _write_function_go_graph(go_graph)
    accession_go.write_text(
        "P00001\tGO:0004672;GO:0005524;GO:0000287;GO:0016740\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "P00001",
            "--go-graph",
            str(go_graph),
            "--accession-go",
            str(accession_go),
            "--verbose",
        ],
    )

    assert result.exit_code == 0
    assert "Query: P00001" in result.stdout
    assert (
        "Resolved terms: GO:0004672;GO:0005524;GO:0000287;GO:0016740"
        in result.stdout
    )
    assert "GO terms:" in result.stdout
    assert "- GO:0004672 protein kinase activity (role=catalytic)" in result.stdout
    assert "- GO:0005524 ATP binding (role=binding_cofactor)" in result.stdout
    assert "- GO:0000287 magnesium ion binding (role=binding_cofactor)" in result.stdout
    assert (
        "- GO:0016740 transferase activity (role=catalytic, dropped=ancestor)"
        in result.stdout
    )
    assert (
        "Function: P00001 is annotated as an ATP- and magnesium-binding "
        "protein kinase."
    ) in result.stdout


def test_inspect_function_json_outputs_structured_description(tmp_path) -> None:
    go_graph = tmp_path / "go_graph.pkl"
    _write_function_go_graph(go_graph)

    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "function",
            "GO:0004672;GO:0005524",
            "--go-graph",
            str(go_graph),
            "--json",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["head"] == "GO:0004672"
    assert payload["modifiers"] == ["ATP-binding"]
    assert payload["summary"] == (
        "GO:0004672;GO:0005524 is annotated as an ATP-binding protein kinase."
    )
