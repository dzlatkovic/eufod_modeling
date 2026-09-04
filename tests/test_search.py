from pathlib import Path

from eufod.comprehensions import solve as comprehensions
from eufod.io import load_input
from eufod.reference import solve as reference
from eufod.search import predicted_shifts, solve_pearson, solve_r_factor


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
