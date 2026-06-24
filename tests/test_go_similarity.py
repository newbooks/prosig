import math

import numpy as np

from prosig.go.similarity import (
    BoundedProfilePairCache,
    BoundedTermPairCache,
    GoSimilarity,
    build_fast_go_similarity_index,
    build_lin_similarity_matrix,
    lin_fast,
    load_accession_mf_go_terms,
    parse_go_term_set,
    resolve_go_set_query,
    set_lin_amb_fast,
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


def test_fast_lin_similarity_matches_diagnostic_score() -> None:
    artifact = _small_artifact()
    diagnostic = GoSimilarity(artifact)
    fast_index = build_fast_go_similarity_index(artifact)

    assert lin_fast(fast_index, "GO:0000002", "GO:0000003") == diagnostic.lin(
        "GO:0000002",
        "GO:0000003",
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


def test_fast_set_lin_amb_matches_diagnostic_score() -> None:
    artifact = _small_artifact()
    diagnostic = GoSimilarity(artifact)
    fast_index = build_fast_go_similarity_index(artifact)

    assert set_lin_amb_fast(
        fast_index,
        ("GO:0000002", "GO:0000003"),
        ("GO:0000002",),
    ) == diagnostic.set_lin_amb(("GO:0000002", "GO:0000003"), ("GO:0000002",))


def test_lin_similarity_matrix_matches_fast_scores() -> None:
    artifact = _small_artifact()
    fast_index = build_fast_go_similarity_index(artifact)
    matrix = build_lin_similarity_matrix(fast_index)
    term_id2 = fast_index.term_id_by_term["GO:0000002"]
    term_id3 = fast_index.term_id_by_term["GO:0000003"]
    root_id = fast_index.term_id_by_term["GO:0003674"]

    assert matrix.dtype == np.float32
    assert matrix.shape == (4, 4)
    assert matrix[term_id2, term_id3] == np.float32(
        lin_fast(fast_index, "GO:0000002", "GO:0000003")
    )
    assert matrix[term_id3, term_id2] == matrix[term_id2, term_id3]
    assert np.isnan(matrix[root_id, root_id])


def test_fast_set_lin_amb_uses_lin_similarity_matrix() -> None:
    artifact = _small_artifact()
    diagnostic = GoSimilarity(artifact)
    fast_index = build_fast_go_similarity_index(artifact)
    matrix = build_lin_similarity_matrix(fast_index)

    assert math.isclose(
        set_lin_amb_fast(
            fast_index,
            ("GO:0000002", "GO:0000003"),
            ("GO:0000002",),
            lin_similarity_matrix=matrix,
        )
        or 0.0,
        diagnostic.set_lin_amb(("GO:0000002", "GO:0000003"), ("GO:0000002",)) or 0.0,
        rel_tol=1e-7,
    )


def test_bounded_similarity_caches_track_hits_and_evictions() -> None:
    term_cache = BoundedTermPairCache(TERM_ENTRY_BYTES_FOR_TEST := 256)
    term_cache.put(1, 0.1)
    assert term_cache.get(1) == (True, 0.1)
    assert term_cache.get(2) == (False, None)
    term_cache.put(2, 0.2)
    stats = term_cache.stats(budget_mb=0)

    assert TERM_ENTRY_BYTES_FOR_TEST == 256
    assert stats.max_entries == 1
    assert stats.entries == 1
    assert stats.hits == 1
    assert stats.misses == 1
    assert stats.evictions == 1

    profile_cache = BoundedProfilePairCache(512)
    key = ((1,), (2,))
    profile_cache.put(key, None)
    assert profile_cache.get(key) == (True, None)


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
