"""Reproducibility metadata for comment-compatible mapping TXT exports."""

from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import platform

import numpy as np

from . import __version__
from .io import parse_input_text
from .model import MIN_DISTANCE


def prepare_input(experimental, oxygen, hydrogens, labels=None, input_text=None):
    """Preserve supplied input text, or create a lossless equivalent for API calls."""
    if input_text is not None:
        parsed_labels, parsed_exp, parsed_oxygen, parsed_hydrogens = parse_input_text(input_text)
        if labels is not None and list(labels) != parsed_labels:
            raise ValueError("Input text labels do not match calculation labels.")
        for parsed, actual in ((parsed_exp, experimental), (parsed_oxygen, oxygen),
                               (parsed_hydrogens, hydrogens)):
            if not np.array_equal(parsed, actual):
                raise ValueError("Input text does not match calculation data.")
        return parsed_labels, input_text

    labels = list(labels) if labels is not None else [f"H-{i + 1}" for i in range(len(experimental))]
    if len(labels) != len(experimental):
        raise ValueError("Provide one label per experimental shift.")
    if any(not isinstance(label, str) or not label or label == "O"
           or "#" in label or any(char.isspace() for char in label) for label in labels):
        raise ValueError("Proton labels must be nonempty tokens other than O, without comments.")
    lines = ["# Coordinate unit: cÅ", "O " + " ".join(f"{v:.17g}" for v in oxygen)]
    for label, shift, xyz in zip(labels, experimental, hydrogens):
        lines.append(f"{label} {shift:.17g} " + " ".join(f"{v:.17g}" for v in xyz))
    text = "\n".join(lines) + "\n"
    parse_input_text(text)  # Also reject duplicate labels.
    return labels, text


def build_metadata(labels, input_text, experimental, oxygen, hydrogens,
                   offset, step, metric, position, score, predicted,
                   candidate_count, valid_count):
    software = {"python": platform.python_version(), "eufod-modeling": __version__}
    for package in ("numpy", "pandas", "streamlit", "py3Dmol"):
        try:
            software[package] = version(package)
        except PackageNotFoundError:
            software[package] = "not installed"
    # Identify the actual source as well as the release version (including local edits).
    source_hash = sha256()
    for path in sorted(Path(__file__).parent.glob("*.py")):
        source_hash.update(path.name.encode("utf-8") + b"\0" + path.read_bytes() + b"\0")
    software["eufod_source_sha256"] = source_hash.hexdigest()
    return {
        "format_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "metric": metric,
        "coordinate_unit": "cÅ",
        "settings": {
            "offset": float(offset), "step": float(step),
            "minimum_distance": MIN_DISTANCE,
            "grid_interval": "[oxygen - offset, oxygen + offset)",
            "normalization_reference": labels[0],
            "distance_constraints": "oxygen and supplied hydrogens only",
            "candidate_count": candidate_count, "valid_count": valid_count,
        },
        "input": {
            "text": input_text,
            "sha256": sha256(input_text.encode("utf-8")).hexdigest(),
            "labels": labels, "experimental_shifts": experimental.tolist(),
            "oxygen": oxygen.tolist(), "hydrogens": hydrogens.tolist(),
        },
        "result": {
            "position": position, "score": score,
            "predicted_shifts": predicted.tolist(),
            "residuals": (predicted - experimental).tolist(),
            "residual_definition": "predicted - experimental",
        },
        "software": software,
    }
