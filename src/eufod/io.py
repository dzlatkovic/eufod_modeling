"""Input/output helpers."""

from pathlib import Path

import numpy as np


def parse_input_text(text: str):
    """Parse experimental shifts, oxygen coordinates, and hydrogen coordinates."""
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

    if len(lines) < 3:
        raise ValueError(
            "Input must contain shifts, oxygen coordinates, and at least one "
            "hydrogen coordinate line."
        )

    try:
        exp_delta = np.asarray([float(x) for x in lines[0].split()], dtype=float)
        oxygen = np.asarray([float(x) for x in lines[1].split()], dtype=float)
        hydrogens = np.asarray(
            [[float(x) for x in line.split()] for line in lines[2:]],
            dtype=float,
        )
    except ValueError as exc:
        raise ValueError("Input contains a non-numeric value.") from exc

    if oxygen.shape != (3,):
        raise ValueError("Oxygen coordinates must contain exactly 3 values.")

    if hydrogens.ndim != 2 or hydrogens.shape[1] != 3:
        raise ValueError("Each hydrogen coordinate must contain exactly 3 values.")

    if len(exp_delta) != len(hydrogens):
        raise ValueError(
            "Number of experimental shifts must equal the number of hydrogen coordinates."
        )

    return exp_delta, oxygen, hydrogens


def load_input(path: str | Path):
    """Load an input text file."""
    return parse_input_text(Path(path).read_text(encoding="utf-8"))


def load_input_bytes(data: bytes):
    """Load input from uploaded file bytes."""
    return parse_input_text(data.decode("utf-8"))
