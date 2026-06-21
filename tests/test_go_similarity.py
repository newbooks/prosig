import math

from prosig.go.similarity import GoSimilarity


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
