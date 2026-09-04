"""Input/output helpers."""

from pathlib import Path

import numpy as np


def parse_input_text(text: str):
    """Parse labeled proton assignments and the coordinating oxygen."""
    lines = [
        line.split("#", 1)[0].strip()
        for line in text.splitlines()
        if line.split("#", 1)[0].strip()
    ]

    oxygen = None
    labels = []
    exp_delta = []
    hydrogens = []

    for line in lines:
        fields = line.split()

        if fields[0] == "O":
            if len(fields) != 4:
                raise ValueError("Oxygen record must be: O X Y Z.")
            if oxygen is not None:
                raise ValueError("Input must contain exactly one oxygen record.")
            try:
                oxygen = [float(value) for value in fields[1:]]
            except ValueError as exc:
                raise ValueError("Oxygen coordinates must be numeric.") from exc
            continue

        if len(fields) != 5:
            raise ValueError(
                "Each proton assignment must be: Label Exp_Shift X Y Z."
            )

        try:
            shift, x, y, z = (float(value) for value in fields[1:])
        except ValueError as exc:
            raise ValueError(
                "Proton assignment shifts and coordinates must be numeric."
            ) from exc

        labels.append(fields[0])
        exp_delta.append(shift)
        hydrogens.append([x, y, z])

    if oxygen is None:
        raise ValueError("Input must contain exactly one oxygen record: O X Y Z.")
    if not hydrogens:
        raise ValueError("Input must contain at least one proton assignment.")
    if len(set(labels)) != len(labels):
        raise ValueError("Proton assignment labels must be unique.")

    return (
        labels,
        np.asarray(exp_delta, dtype=float),
        np.asarray(oxygen, dtype=float),
        np.asarray(hydrogens, dtype=float),
    )


def load_input(path: str | Path):
    """Load an input text file."""
    return parse_input_text(Path(path).read_text(encoding="utf-8"))


def load_input_bytes(data: bytes):
    """Load input from uploaded file bytes."""
    return parse_input_text(data.decode("utf-8"))
