"""Straightforward Python reference implementation."""

from math import dist

from .model import MIN_DISTANCE, delta_value, r_factor


def solve(exp_delta, oxygen, hydrogens, offset=300, step=5):
    """Search the grid using explicit Python loops."""
    oh_distances = [dist(oxygen, h) for h in hydrogens]
    best_r = float("inf")
    best_pos = None

    ox, oy, oz = oxygen

    for x in range(int(ox - offset), int(ox + offset), step):
        for y in range(int(oy - offset), int(oy + offset), step):
            for z in range(int(oz - offset), int(oz + offset), step):
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
                corrected = [v / raw[0] for v in raw]
                score = r_factor(corrected, exp_delta)

                if score < best_r:
                    best_r, best_pos = score, eu

    return best_r, best_pos
