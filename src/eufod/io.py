"""Input/output helpers."""

from pathlib import Path
from hashlib import sha256
from datetime import datetime, timedelta
import json

import numpy as np


def parse_input_text(text: str):
    """Parse labeled proton assignments and the coordinating oxygen."""
    lines = [
        line.split("#", 1)[0].strip()
        for line in text.lstrip("\ufeff").splitlines()
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

    if not all(np.isfinite(values).all() for values in (oxygen, exp_delta, hydrogens)):
        raise ValueError("Input shifts and coordinates must be finite.")

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


def parse_mapping_metadata(text: str):
    """Read version-2 metadata, returning None for legacy mapping files.

    The coordinate-only parser deliberately continues to ignore comments.
    This reader validates metadata before the UI uses it or restores input.
    """
    version = None
    lines = []
    active = False
    seen = False
    complete = False
    for line in text.lstrip("\ufeff").splitlines():
        line = line.strip()
        if not line.startswith("#"):
            if active:
                raise ValueError("Metadata block is incomplete: expected commented lines or an end marker.")
            continue
        comment = line[1:].lstrip()
        if comment.startswith("format_version:"):
            if version is not None:
                raise ValueError("Duplicate mapping format version.")
            version = comment.split(":", 1)[1].strip()
        elif comment == "METADATA_BEGIN":
            if seen:
                raise ValueError("Duplicate metadata block.")
            seen = active = True
        elif comment == "METADATA_END":
            if not active:
                raise ValueError("Unexpected metadata end marker.")
            active = False
            complete = True
        elif active:
            lines.append(comment)
    if version is None and not seen:
        return None
    if version != "2":
        raise ValueError(f"Unsupported mapping format version: {version}.")
    if not complete or active:
        raise ValueError("Mapping metadata block is missing or incomplete.")
    try:
        def reject_constant(value):
            raise ValueError(f"Nonfinite JSON value: {value}")

        metadata = json.loads("\n".join(lines), parse_constant=reject_constant)
        if metadata["format_version"] != 2 or metadata["coordinate_unit"] != "cÅ":
            raise ValueError("Unsupported metadata version or coordinate unit.")
        if metadata["metric"] not in {"r_factor", "pearson_r"}:
            raise ValueError("Unknown metadata metric.")
        original = metadata["input"]
        labels, shifts, oxygen, hydrogens = parse_input_text(original["text"])
        if sha256(original["text"].encode("utf-8")).hexdigest() != original["sha256"]:
            raise ValueError("Input checksum does not match metadata.")
        if labels != original["labels"]:
            raise ValueError("Input labels do not match metadata.")
        for parsed, key in ((shifts, "experimental_shifts"), (oxygen, "oxygen"),
                            (hydrogens, "hydrogens")):
            if not np.array_equal(parsed, original[key]):
                raise ValueError("Input values do not match metadata.")
        result = metadata["result"]
        for key, shape in (("position", (3,)), ("predicted_shifts", shifts.shape),
                           ("residuals", shifts.shape), ("score", ())):
            array = np.asarray(result[key])
            if array.shape != shape or array.dtype.kind not in "iuf" or not np.isfinite(array).all():
                raise ValueError(f"Invalid result metadata: {key}.")
        if not np.allclose(np.asarray(result["predicted_shifts"]) - shifts,
                           result["residuals"], rtol=1e-12, atol=1e-12):
            raise ValueError("Residuals do not match metadata shifts.")
        settings = metadata["settings"]
        for key in ("offset", "step", "minimum_distance"):
            value = settings[key]
            if (type(value) not in (int, float) or not np.isfinite(value)
                    or value <= 0):
                raise ValueError(f"Invalid search setting: {key}.")
        for key in ("candidate_count", "valid_count"):
            if type(settings[key]) is not int or settings[key] <= 0:
                raise ValueError(f"Invalid candidate count: {key}.")
        if settings["valid_count"] > settings["candidate_count"]:
            raise ValueError("Valid count exceeds total candidate count.")
        if settings["normalization_reference"] != labels[0]:
            raise ValueError("Normalization reference must match the first proton.")
        if not isinstance(metadata["software"], dict) or not isinstance(metadata["created_utc"], str):
            raise ValueError("Invalid provenance metadata.")
        timestamp = datetime.fromisoformat(metadata["created_utc"].replace("Z", "+00:00"))
        if timestamp.utcoffset() != timedelta(0):
            raise ValueError("Calculation timestamp must include the UTC timezone.")
    except (KeyError, TypeError, AttributeError, ValueError, OverflowError) as exc:
        raise ValueError(f"Invalid mapping metadata: {exc}") from exc
    return metadata


def parse_mapping_document(text: str):
    """Read legacy or extended mappings and cross-check metadata against rows."""
    metric, coordinates, scores = parse_mapping_text(text)
    metadata = parse_mapping_metadata(text)
    if metadata is not None:
        if metadata["metric"] != metric:
            raise ValueError("Metadata metric does not match the mapping header.")
        if metadata["settings"]["valid_count"] != len(scores):
            raise ValueError("Mapping row count does not match metadata; the file may be incomplete.")
        best_index = np.argmin(scores) if metric == "r_factor" else np.argmax(scores)
        matches = np.all(coordinates == metadata["result"]["position"], axis=1)
        # Exported scores have eight decimal places; full-precision ties can
        # select a different row than argmin/argmax of the rounded table.
        if (not matches.any()
                or not np.isclose(scores[best_index], metadata["result"]["score"], rtol=0, atol=5.1e-9)
                or not np.any(np.isclose(scores[matches], metadata["result"]["score"], rtol=0, atol=5.1e-9))):
            raise ValueError("Metadata best result does not match the mapping.")
    return metric, coordinates, scores, metadata
