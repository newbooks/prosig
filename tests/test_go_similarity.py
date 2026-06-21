import math

from prosig.go.similarity import (
    GoSimilarity,
    load_accession_mf_go_terms,
    parse_go_term_set,
    resolve_go_set_query,
)


def _small_artifact() -> dict:
    return {
        "meta": {
            "schema_version": "1.0",
            "namespace": "molecular_function",
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
                "children": ["GO:0000002", "GO:0000003", "GO:0000004"],
                "ancestors": {"GO:0003674"},
                "depth": 1,
                "freq": math.exp(-1.0),
                "ic": 1.0,
            },
            "GO:0000002": {
                "name": "child A activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": math.exp(-2.0),
                "ic": 2.0,
            },
            "GO:0000003": {
                "name": "child B activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": math.exp(-3.0),
                "ic": 3.0,
            },
            "GO:0000004": {
                "name": "missing IC activity",
                "parents": ["GO:0000001"],
                "children": [],
                "ancestors": {"GO:0003674", "GO:0000001"},
                "depth": 2,
                "freq": 0.0,
                "ic": None,
            },
        },
    }


def test_lin_similarity_identical_terms_return_one() -> None:
    result = GoSimilarity(_small_artifact()).lin_with_details(
        "GO:0000002",
        "GO:0000002",
    )

    assert result.status == "ok"
    assert result.mica == "GO:0000002"
    assert result.similarity == 1.0


def test_lin_similarity_siblings_use_mica() -> None:
    result = GoSimilarity(_small_artifact()).lin_with_details(
        "GO:0000002",
        "GO:0000003",
    )

    assert result.status == "ok"
    assert result.mica == "GO:0000001"
    assert result.similarity == 2 * 1.0 / (2.0 + 3.0)


def test_scalar_lin_similarity_matches_detailed_score() -> None:
    similarity = GoSimilarity(_small_artifact())

    assert similarity.lin("GO:0000002", "GO:0000003") == (
        similarity.lin_with_details("GO:0000002", "GO:0000003").similarity
    )


def test_set_lin_amb_uses_directional_best_matches() -> None:
    result = GoSimilarity(_small_artifact()).set_lin_amb_with_details(
        ("GO:0000002", "GO:0000003"),
        ("GO:0000002",),
    )

    assert result.status == "ok"
    assert result.similarity == 0.85
    best_1_to_2 = [
        (match.source, match.target, match.score)
        for match in result.best_matches_1_to_2
    ]
    best_2_to_1 = [
        (match.source, match.target, match.score)
        for match in result.best_matches_2_to_1
    ]
    assert best_1_to_2 == [
        ("GO:0000002", "GO:0000002", 1.0),
        ("GO:0000003", "GO:0000002", 0.4),
    ]
    assert best_2_to_1 == [
        ("GO:0000002", "GO:0000002", 1.0),
    ]


def test_scalar_set_lin_amb_matches_detailed_score() -> None:
    similarity = GoSimilarity(_small_artifact())

    assert similarity.set_lin_amb(("GO:0000002", "GO:0000003"), ("GO:0000002",)) == (
        similarity.set_lin_amb_with_details(
            ("GO:0000002", "GO:0000003"),
            ("GO:0000002",),
        ).similarity
    )


def test_set_lin_amb_reports_missing_terms_and_empty_cleaned_set() -> None:
    result = GoSimilarity(_small_artifact()).set_lin_amb_with_details(
        ("GO:9999999",),
        ("GO:0000002",),
    )

    assert result.status == "unavailable"
    assert result.reason == "empty_cleaned_set"
    assert result.missing_terms1 == ("GO:9999999",)
    assert result.similarity is None


def test_parse_go_term_set_accepts_semicolons_and_commas() -> None:
    assert parse_go_term_set("(GO:0000002; GO:0000003,GO:0000002)") == (
        "GO:0000002",
        "GO:0000003",
    )
    assert parse_go_term_set("GO:0000002;GO:0000003") == (
        "GO:0000002",
        "GO:0000003",
    )


def test_parse_go_term_set_rejects_bad_input() -> None:
    try:
        parse_go_term_set("(GO:00000012)")
    except ValueError as exc:
        assert "Malformed GO term" in str(exc)
    else:
        raise AssertionError("Expected malformed GO term to be rejected")


def test_accession_mf_go_loader_and_resolver(tmp_path) -> None:
    accession_go = tmp_path / "accession_mf_go.tsv"
    accession_go.write_text(
        "P00001\tGO:0000002;GO:0000003\n"
        "P00002\tGO:0000001\n",
        encoding="utf-8",
    )

    accession_terms = load_accession_mf_go_terms(accession_go)

    assert accession_terms["P00001"] == ("GO:0000002", "GO:0000003")
    assert resolve_go_set_query("P00002", accession_terms) == ("GO:0000001",)


def test_lin_similarity_missing_ic_is_unavailable() -> None:
    result = GoSimilarity(_small_artifact()).lin_with_details(
        "GO:0000004",
        "GO:0000003",
    )

    assert result.status == "unavailable"
    assert result.reason == "missing_ic_go1"
    assert result.similarity is None


def test_lin_similarity_root_zero_denominator_is_unavailable() -> None:
    result = GoSimilarity(_small_artifact()).lin_with_details(
        "GO:0003674",
        "GO:0003674",
    )

    assert result.status == "unavailable"
    assert result.reason == "zero_ic_denominator"
    assert result.mica == "GO:0003674"
    assert result.ic_mica == 0.0
