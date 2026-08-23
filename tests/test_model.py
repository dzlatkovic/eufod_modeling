import numpy as np

from eufod.model import pearson_coefficient, r_factor


def test_r_factor_perfect_match():
    x = np.array([1.0, 2.0, 3.0])
    assert r_factor(x, x) == 0


def test_pearson_perfect_linear_pattern():
    x = np.array([1.0, 2.0, 3.0])
    assert np.isclose(pearson_coefficient(x, x), 1.0)
