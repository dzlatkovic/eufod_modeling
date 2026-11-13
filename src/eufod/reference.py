"""Straightforward Python reference implementation."""

from math import dist, isfinite

from .model import MIN_DISTANCE, delta_value, r_factor
from .validation import grid_axes, validate_inputs


def solve(exp_delta, oxygen, hydrogens, offset=300, step=5):
    """Search the grid using explicit Python loops."""
    exp_delta, oxygen, hydrogens = validate_inputs(exp_delta, oxygen, hydrogens)
    axes = grid_axes(oxygen, offset, step)
    oh_distances = [dist(oxygen, h) for h in hydrogens]
    best_r = float("inf")
    best_pos = None

    for x in axes[0]:
        for y in axes[1]:
            for z in axes[2]:
                eu = [x, y, z]

                eu_oxygen = dist(eu, oxygen)
                if eu_oxygen < MIN_DISTANCE:
                    continue

                eu_h = [dist(eu, h) for h in hydrogens]
                if any(d < MIN_DISTANCE for d in eu_h):
                    continue

                raw = [
                    delta_value(a, eu_oxygen, c)
                    for a, c in zip(eu_h, oh_distances)
                ]
                if raw[0] == 0 or not all(isfinite(v) for v in raw):
                    continue
                corrected = [v / raw[0] for v in raw]
                score = r_factor(corrected, exp_delta)

                if score < best_r:
                    best_r, best_pos = score, eu

    if best_pos is None:
        raise ValueError("No valid Eu positions were found.")
    return best_r, best_pos
