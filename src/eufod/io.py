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


def parse_mapping_text(text: str):
    """Parse a coordinate mapping exported by the EuFOD application."""
    metric = None
    rows = []

    for line in text.lstrip("\ufeff").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if stripped.lower().startswith("# metric:"):
                metric = stripped.split(":", 1)[1].strip().lower()
            continue

        fields = stripped.split()
        if len(fields) != 4:
            raise ValueError("Mapping rows must contain X Y Z and score.")
        try:
            rows.append([float(value) for value in fields])
        except ValueError as exc:
            raise ValueError("Mapping coordinates and scores must be numeric.") from exc

    if metric not in {"r_factor", "pearson_r"}:
        raise ValueError("Mapping must declare metric r_factor or pearson_r.")
    if not rows:
        raise ValueError("Mapping must contain at least one candidate position.")

    values = np.asarray(rows, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Mapping coordinates and scores must be finite.")

    return metric, values[:, :3], values[:, 3]
