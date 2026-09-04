"""Vectorized Eu-position searches and mapping generation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .model import MIN_DISTANCE, pearson_coefficient


@dataclass
class SearchResult:
    """Result returned by the mapping-capable search."""

    score: float
    position: list[int]
    mapping_text: str | None = None
    preview: object | None = None


def candidate_positions(oxygen, offset=300, step=5):
    """Create the Cartesian search grid around oxygen in cÅ units."""
    ranges = [
        np.arange(c - offset, c + offset, step)
        for c in oxygen
    ]
    return (
        np.stack(np.meshgrid(*ranges, indexing="ij"), axis=-1)
        .reshape(-1, 3)
        .astype(float)
    )


def _predicted_shifts(positions, oxygen, hydrogens):
    """Calculate corrected shifts for all candidate positions."""
    oxygen_eu = np.linalg.norm(positions - oxygen, axis=1)
    eu_h = np.linalg.norm(
        positions[:, None, :] - hydrogens[None, :, :],
        axis=2,
    )
    oh = np.linalg.norm(hydrogens - oxygen, axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        cos_theta = (
            eu_h**2
            + oxygen_eu[:, None] ** 2
            - oh[None, :] ** 2
        ) / (2 * eu_h * oxygen_eu[:, None])

        raw = (3 * cos_theta**2 - 1) / eu_h**3

    return raw / raw[:, [0]], oxygen_eu, eu_h


def _valid_mask(oxygen_eu, eu_h):
    """Return candidates satisfying the minimum-distance constraints."""
    return (
        (oxygen_eu >= MIN_DISTANCE)
        & np.all(eu_h >= MIN_DISTANCE, axis=1)
    )


def predicted_shifts(position, oxygen, hydrogens):
    """Calculate normalized predicted shifts for one Eu position."""
    position = np.asarray(position, dtype=float).reshape(1, 3)
    shifts, _, _ = _predicted_shifts(position, oxygen, hydrogens)
    return shifts[0]


def solve_r_factor(exp_delta, oxygen, hydrogens, offset=300, step=5):
    """Vectorized search minimizing the R-factor."""
    positions = candidate_positions(oxygen, offset, step)
    shifts, oxygen_eu, eu_h = _predicted_shifts(
        positions, oxygen, hydrogens
    )
    valid = _valid_mask(oxygen_eu, eu_h)

    scores = np.full(len(positions), np.inf)
    scores[valid] = 100.0 * np.sqrt(
        np.mean(
            ((shifts[valid] - exp_delta) / exp_delta) ** 2,
            axis=1,
        )
    )

    index = np.argmin(scores)
    return float(scores[index]), positions[index].astype(int).tolist()


def solve_pearson(exp_delta, oxygen, hydrogens, offset=300, step=5):
    """Vectorized search maximizing Pearson correlation."""
    positions = candidate_positions(oxygen, offset, step)
    shifts, oxygen_eu, eu_h = _predicted_shifts(
        positions, oxygen, hydrogens
    )
    valid = _valid_mask(oxygen_eu, eu_h)

    scores = np.full(len(positions), np.nan)
    scores[valid] = pearson_coefficient(shifts[valid], exp_delta)

    index = np.nanargmax(scores)
    return float(scores[index]), positions[index].astype(int).tolist()


def _score_chunk(shifts, exp_delta, metric):
    """Score one chunk of predicted shifts."""
    if metric == "r_factor":
        return 100.0 * np.sqrt(
            np.mean(
                ((shifts - exp_delta) / exp_delta) ** 2,
                axis=1,
            )
        )

    if metric == "pearson":
        return pearson_coefficient(shifts, exp_delta)

    raise ValueError(f"Unknown metric: {metric}")


def search_mapping(
    exp_delta,
    oxygen,
    hydrogens,
    offset=300,
    step=5,
    metric="r_factor",
    include_mapping=True,
    chunk_size=25_000,
    progress_callback: Callable[[float], None] | None = None,
):
    """Search positions in chunks and optionally create a text mapping.

    Chunking keeps the full mapping search from requiring one giant NumPy
    allocation. The returned text contains only valid candidate positions.
    """
    positions = candidate_positions(oxygen, offset, step)
    total = len(positions)

    best_score = np.inf if metric == "r_factor" else -np.inf
    best_position = None
    mapping_lines = []
    preview_rows = []

    for start in range(0, total, chunk_size):
        stop = min(start + chunk_size, total)
        chunk_positions = positions[start:stop]

        shifts, oxygen_eu, eu_h = _predicted_shifts(
            chunk_positions,
            oxygen,
            hydrogens,
        )
        valid = _valid_mask(oxygen_eu, eu_h)

        scores = np.full(
            len(chunk_positions),
            np.inf if metric == "r_factor" else np.nan,
        )
        scores[valid] = _score_chunk(
            shifts[valid],
            exp_delta,
            metric,
        )

        if metric == "r_factor":
            local_index = np.argmin(scores)
            local_best = scores[local_index]
            better = local_best < best_score
        else:
            local_index = np.nanargmax(scores)
            local_best = scores[local_index]
            better = local_best > best_score

        if better:
            best_score = float(local_best)
            best_position = (
                chunk_positions[local_index]
                .astype(int)
                .tolist()
            )

        if include_mapping:
            valid_indices = np.flatnonzero(valid)
            for idx in valid_indices:
                x, y, z = chunk_positions[idx].astype(int)
                score = float(scores[idx])
                mapping_lines.append(
                    f"{x}\t{y}\t{z}\t{score:.8f}\n"
                )

                if len(preview_rows) < 100:
                    preview_rows.append(
                        {
                            "x": int(x),
                            "y": int(y),
                            "z": int(z),
                            "score": score,
                        }
                    )

        if progress_callback is not None:
            progress_callback(stop / total)

    if best_position is None:
        raise ValueError("No valid Eu positions were found.")

    mapping_text = None
    preview = None

    if include_mapping:
        metric_name = "r_factor" if metric == "r_factor" else "pearson_r"
        header = (
            "# Eu(fod)3 coordinate mapping\n"
            "# coordinate unit: cÅ\n"
            f"# metric: {metric_name}\n"
            f"# offset: {offset}\n"
            f"# step: {step}\n"
            "# valid positions only\n"
            "# x\ty\tz\tscore\n"
        )
        mapping_text = header + "".join(mapping_lines)

        # Import only when the UI actually requests a preview.
        import pandas as pd

        preview = pd.DataFrame(preview_rows)

    return SearchResult(
        score=float(best_score),
        position=best_position,
        mapping_text=mapping_text,
        preview=preview,
    )
