import numpy as np
import pytest

from eufod.io import parse_input_text


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
