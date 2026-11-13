"""Vectorized Eu-position searches and mapping generation."""

from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import json
from math import prod
from numbers import Integral
from typing import Callable

import numpy as np

from .model import MIN_DISTANCE, pearson_coefficient
from .provenance import build_metadata, prepare_input
from .validation import grid_axes, validate_geometry, validate_inputs


@dataclass
class SearchResult:
    """Result returned by the mapping-capable search."""

    score: float
    position: list[float]
    mapping_text: str | None = None
    preview: object | None = None


def candidate_positions(oxygen, offset=300, step=5):
    """Create the Cartesian search grid around oxygen in cÅ units."""
    ranges = grid_axes(oxygen, offset, step)
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

        shifts = raw / raw[:, [0]]

    return shifts, oxygen_eu, eu_h


def _valid_mask(oxygen_eu, eu_h):
    """Return candidates satisfying the minimum-distance constraints."""
    return (
        (oxygen_eu >= MIN_DISTANCE)
        & np.all(eu_h >= MIN_DISTANCE, axis=1)
    )


def predicted_shifts(position, oxygen, hydrogens):
    """Calculate normalized predicted shifts for one Eu position."""
    oxygen, hydrogens = validate_geometry(oxygen, hydrogens)
    position = np.asarray(position, dtype=float).reshape(1, 3)
    if not np.isfinite(position).all():
        raise ValueError("Eu position must be finite.")
    shifts, _, _ = _predicted_shifts(position, oxygen, hydrogens)
    if not np.isfinite(shifts).all():
        raise ValueError("Predicted shifts cannot be normalized at this Eu position.")
    return shifts[0]


def solve_r_factor(exp_delta, oxygen, hydrogens, offset=300, step=5):
    """Vectorized, bounded-memory search minimizing the R-factor."""
    result = search_mapping(exp_delta, oxygen, hydrogens, offset, step,
                            metric="r_factor", include_mapping=False)
    return result.score, result.position


def solve_pearson(exp_delta, oxygen, hydrogens, offset=300, step=5):
    """Vectorized, bounded-memory search maximizing Pearson correlation."""
    result = search_mapping(exp_delta, oxygen, hydrogens, offset, step,
                            metric="pearson", include_mapping=False)
    return result.score, result.position


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
    *,
    labels=None,
    input_text=None,
):
    """Search positions in chunks and optionally create a text mapping.

    Chunking keeps the full mapping search from requiring one giant NumPy
    allocation. The returned text contains only valid candidate positions.
    """
    exp_delta, oxygen, hydrogens = validate_inputs(
        exp_delta, oxygen, hydrogens, metric
    )
    if include_mapping:
        labels, input_text = prepare_input(exp_delta, oxygen, hydrogens, labels, input_text)
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, Integral) or chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer.")
    axes = grid_axes(oxygen, offset, step)
    shape = tuple(len(axis) for axis in axes)
    total = prod(shape)

    best_score = np.inf if metric == "r_factor" else -np.inf
    best_position = None
    mapping_buffer = StringIO()
    preview_rows = []
    valid_count = 0

    for start in range(0, total, chunk_size):
        stop = min(start + chunk_size, total)
        indices = np.unravel_index(np.arange(start, stop), shape)
        chunk_positions = np.column_stack([
            axis[index] for axis, index in zip(axes, indices)
        ])

        shifts, oxygen_eu, eu_h = _predicted_shifts(
            chunk_positions,
            oxygen,
            hydrogens,
        )
        valid = _valid_mask(oxygen_eu, eu_h)

        valid &= np.isfinite(shifts).all(axis=1)
        valid_positions = chunk_positions[valid]
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            scores = _score_chunk(shifts[valid], exp_delta, metric)
        finite = np.isfinite(scores)
        scores = scores[finite]
        valid_positions = valid_positions[finite]
        valid_count += len(scores)

        if len(scores):
            local_index = np.argmin(scores) if metric == "r_factor" else np.argmax(scores)
            local_best = scores[local_index]
            better = local_best < best_score if metric == "r_factor" else local_best > best_score
            if better:
                best_score = float(local_best)
                best_position = valid_positions[local_index].tolist()

            if include_mapping:
                for (x, y, z), score in zip(valid_positions, scores):
                    mapping_buffer.write(f"{x:.17g}\t{y:.17g}\t{z:.17g}\t{score:.8f}\n")
                    if len(preview_rows) < 100:
                        preview_rows.append({"x": x, "y": y, "z": z, "score": float(score)})

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
        metadata = build_metadata(
            labels, input_text, exp_delta, oxygen, hydrogens,
            offset, step, metric_name, best_position, float(best_score),
            predicted_shifts(best_position, oxygen, hydrogens), total, valid_count,
        )
        metadata_comments = "\n".join(
            "# " + line for line in json.dumps(metadata, ensure_ascii=False, indent=2,
                                               allow_nan=False).splitlines()
        )
        mapping_text = (header + "# format_version: 2\n# METADATA_BEGIN\n"
                        + metadata_comments + "\n# METADATA_END\n"
                        + mapping_buffer.getvalue())

        # Import only when the UI actually requests a preview.
        import pandas as pd

        preview = pd.DataFrame(preview_rows)

    return SearchResult(
        score=float(best_score),
        position=best_position,
        mapping_text=mapping_text,
        preview=preview,
    )
