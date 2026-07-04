import math

from prosig.features.quality import compute_specificity_score


def test_specificity_selects_neighbors_by_feature_distance_not_gap() -> None:
    clusters = {
        "A": [0.0, 0.0],
        "B": [5.0, 5.0],
        "C": [0.0, 20.0],
    }

    result = compute_specificity_score(clusters, neighbor_k=1)

    global_std = math.sqrt(60.0)
    assert math.isclose(
        result["specificity"],
        5.0 / (5.0 + global_std),
    )
