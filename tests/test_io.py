import numpy as np
import pytest

from eufod.io import parse_input_text, parse_mapping_text


def test_parse_labeled_assignments_with_comments():
    labels, shifts, oxygen, hydrogens = parse_input_text(
        """
        # Oxygen Center
        O 1.254 2.108 -0.543  # coordinating atom

        # Label, Exp_Shift, X, Y, Z
        H-1 1.000 2.110 3.450 0.120
        H-3a 0.432 3.015 0.420 1.890
        """
    )

    assert labels == ["H-1", "H-3a"]
    np.testing.assert_allclose(shifts, [1.000, 0.432])
    np.testing.assert_allclose(oxygen, [1.254, 2.108, -0.543])
    np.testing.assert_allclose(
        hydrogens,
        [[2.110, 3.450, 0.120], [3.015, 0.420, 1.890]],
    )


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("H-1 1.000 2.110 3.450 0.120", "oxygen record"),
        ("O 0 0 0\nO 1 1 1\nH-1 1 2 3 4", "exactly one oxygen"),
        ("O 0 0 0\nH-1 1 2 3", "proton assignment"),
        ("O 0 0 0\nH-1 1 2 3 4\nH-1 2 3 4 5", "labels must be unique"),
    ],
)
def test_parse_input_rejects_invalid_records(text, message):
    with pytest.raises(ValueError, match=message):
        parse_input_text(text)


def test_parse_r_factor_mapping():
    metric, coordinates, scores = parse_mapping_text(
        """\ufeff
        # Eu(fod)3 coordinate mapping
        # coordinate unit: cÅ
        # metric: r_factor
        # x y z score
        100 200 -50 5.25
        110 200 -50 5.75
        """
    )

    assert metric == "r_factor"
    np.testing.assert_allclose(coordinates, [[100, 200, -50], [110, 200, -50]])
    np.testing.assert_allclose(scores, [5.25, 5.75])


def test_parse_mapping_rejects_unknown_metric():
    with pytest.raises(ValueError, match="declare metric"):
        parse_mapping_text("# metric: unknown\n1 2 3 4")
