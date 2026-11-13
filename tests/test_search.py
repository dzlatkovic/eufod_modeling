from pathlib import Path

import numpy as np
import pytest

from eufod.comprehensions import solve as comprehensions
from eufod.io import load_input
from eufod.reference import solve as reference
from eufod.search import predicted_shifts, solve_pearson, solve_r_factor
from eufod.search import candidate_positions, search_mapping
from eufod.io import parse_mapping_text
from eufod.model import pearson_coefficient, r_factor


def get_example():
    return load_input(
        Path(__file__).parents[1] / "data/example_input.txt"
    )


def test_reference_comprehensions_numpy_agree():
    _labels, exp, oxygen, hydrogens = get_example()
    args = (exp, oxygen, hydrogens, 300, 50)

    reference_result = reference(*args)
    comprehension_result = comprehensions(*args)
    numpy_result = solve_r_factor(*args)

    assert reference_result == comprehension_result == numpy_result


def test_pearson_returns_valid_result():
    _labels, exp, oxygen, hydrogens = get_example()

    score, position = solve_pearson(
        exp,
        oxygen,
        hydrogens,
        300,
        50,
    )

    assert 0 <= score <= 1
    assert len(position) == 3


def test_predicted_shifts_are_normalized_to_first_proton():
    _labels, exp, oxygen, hydrogens = get_example()
    _score, position = solve_r_factor(exp, oxygen, hydrogens, 300, 50)

    predicted = predicted_shifts(position, oxygen, hydrogens)

    assert predicted.shape == exp.shape
    assert predicted[0] == 1.0


@pytest.mark.parametrize("metric", ["r_factor", "pearson"])
def test_chunked_mapping_matches_independent_grid_scores(metric):
    _, exp, oxygen, hydrogens = get_example()
    oxygen = oxygen + 0.25
    hydrogens = hydrogens + 0.25
    progress = []
    result = search_mapping(exp.tolist(), oxygen.tolist(), hydrogens.tolist(),
                            step=100, metric=metric, chunk_size=1,
                            progress_callback=progress.append)
    expected_positions, expected_scores = [], []
    for position in candidate_positions(oxygen, step=100):
        if np.linalg.norm(position - oxygen) < 200:
            continue
        if np.any(np.linalg.norm(hydrogens - position, axis=1) < 200):
            continue
        shifts = predicted_shifts(position, oxygen, hydrogens)
        score = r_factor(shifts, exp) if metric == "r_factor" else pearson_coefficient(shifts, exp)
        if np.isfinite(score):
            expected_positions.append(position)
            expected_scores.append(score)
    _, positions, scores = parse_mapping_text(result.mapping_text)
    np.testing.assert_array_equal(positions, expected_positions)
    np.testing.assert_allclose(scores, expected_scores, atol=5e-9)
    best = np.argmin(expected_scores) if metric == "r_factor" else np.argmax(expected_scores)
    assert result.score == pytest.approx(expected_scores[best])
    np.testing.assert_array_equal(result.position, expected_positions[best])
    assert progress[-1] == 1
    assert all(a < b for a, b in zip(progress, progress[1:]))


def test_fractional_grid_agrees_across_implementations():
    _, exp, oxygen, hydrogens = get_example()
    args = (exp, oxygen + 0.25, hydrogens + 0.25, 300, 50)
    results = [fn(*args) for fn in (reference, comprehensions, solve_r_factor)]
    for score, position in results:
        assert score == pytest.approx(results[0][0])
        assert position == results[0][1]
        assert all(value % 1 == 0.25 for value in position)


@pytest.mark.parametrize("solver", [reference, comprehensions, solve_r_factor, solve_pearson, search_mapping])
def test_no_valid_candidate_raises(solver):
    _, exp, oxygen, hydrogens = get_example()
    with pytest.raises(ValueError, match="No valid Eu positions"):
        solver(exp, oxygen, hydrogens, offset=1, step=5)


@pytest.mark.parametrize("kwargs, message", [
    ({"step": 0}, "step"), ({"step": -1}, "step"),
    ({"offset": 0}, "offset"), ({"offset": np.nan}, "offset"),
    ({"chunk_size": 0}, "chunk_size"), ({"chunk_size": 1.5}, "chunk_size"),
    ({"metric": "unknown"}, "Unknown metric"),
])
def test_invalid_search_settings(kwargs, message):
    _, exp, oxygen, hydrogens = get_example()
    with pytest.raises(ValueError, match=message):
        search_mapping(exp, oxygen, hydrogens, **kwargs)


@pytest.mark.parametrize("experimental, metric, message", [
    ([1, 0], "r_factor", "nonzero"),
    ([1, 1], "pearson", "nonconstant"),
    ([1, np.nan], "r_factor", "finite"),
    ([1], "r_factor", "one experimental shift"),
])
def test_invalid_experimental_data(experimental, metric, message):
    with pytest.raises(ValueError, match=message):
        search_mapping(experimental, [0, 0, 0], [[100, 0, 0], [0, 100, 0]], metric=metric)


def test_undefined_calculated_pearson_is_not_exported():
    with pytest.raises(ValueError, match="No valid Eu positions"):
        search_mapping([1, 2], [0, 0, 0], [[100, 0, 0], [100, 0, 0]],
                       metric="pearson", step=100)


def test_prediction_rejects_singular_position():
    with pytest.raises(ValueError, match="cannot be normalized"):
        predicted_shifts([0, 0, 0], [0, 0, 0], [[100, 0, 0]])
