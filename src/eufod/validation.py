"""Shared validation and grid definition for all search implementations."""

import numpy as np


def validate_geometry(oxygen, hydrogens):
    oxygen = np.asarray(oxygen, dtype=float)
    hydrogens = np.asarray(hydrogens, dtype=float)
    if oxygen.shape != (3,):
        raise ValueError("Oxygen must contain three coordinates.")
    if hydrogens.ndim != 2 or hydrogens.shape[1] != 3 or not len(hydrogens):
        raise ValueError("Hydrogens must be a nonempty N x 3 array.")
    if not np.isfinite(oxygen).all() or not np.isfinite(hydrogens).all():
        raise ValueError("Coordinates must be finite.")
    return oxygen, hydrogens


def validate_inputs(exp_delta, oxygen, hydrogens, metric="r_factor"):
    oxygen, hydrogens = validate_geometry(oxygen, hydrogens)
    experimental = np.asarray(exp_delta, dtype=float)
    if experimental.shape != (len(hydrogens),):
        raise ValueError("Provide one experimental shift per hydrogen.")
    if not np.isfinite(experimental).all():
        raise ValueError("Experimental shifts must be finite.")
    if metric not in {"r_factor", "pearson"}:
        raise ValueError(f"Unknown metric: {metric}")
    if metric == "r_factor" and np.any(experimental == 0):
        raise ValueError("R-factor requires nonzero experimental shifts.")
    if metric == "pearson" and (
        len(experimental) < 2 or np.all(experimental == experimental[0])
    ):
        raise ValueError("Pearson requires at least two nonconstant experimental shifts.")
    return experimental, oxygen, hydrogens


def grid_axes(oxygen, offset=300, step=5):
    """Axes include oxygen - offset and exclude oxygen + offset."""
    oxygen = np.asarray(oxygen, dtype=float)
    if oxygen.shape != (3,) or not np.isfinite(oxygen).all():
        raise ValueError("Oxygen must contain three finite coordinates.")
    for name, value in (("offset", offset), ("step", step)):
        if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be a positive finite number.")
    return [np.arange(c - offset, c + offset, step) for c in oxygen]
