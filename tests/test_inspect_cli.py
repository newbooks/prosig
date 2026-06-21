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


def test_inspect_help_lists_diagnostic_commands() -> None:
    result = CliRunner().invoke(app, ["inspect", "-h"])

    assert result.exit_code == 0
    assert "go-summary" in result.stdout
    assert "go-term" in result.stdout
    assert "go-sim" in result.stdout
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
