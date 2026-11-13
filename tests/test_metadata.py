import json
from pathlib import Path

import numpy as np
import pytest

from eufod.io import parse_input_text, parse_mapping_metadata, parse_mapping_text, parse_mapping_document
from eufod.search import search_mapping

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def exported():
    text = "\ufeff# Original comments: α\r\n" + (ROOT / "data/example_input.txt").read_text(encoding="utf-8")
    labels, exp, oxygen, hydrogens = parse_input_text(text)
    result = search_mapping(exp, oxygen, hydrogens, step=50, labels=labels, input_text=text)
    return text, result.mapping_text


@pytest.mark.parametrize("metric", ["r_factor", "pearson"])
def test_extended_export_reproduces_calculation(metric):
    original = (ROOT / "data/example_input.txt").read_text(encoding="utf-8")
    labels, exp, oxygen, hydrogens = parse_input_text(original)
    result = search_mapping(exp, oxygen, hydrogens, step=50, metric=metric,
                            labels=labels, input_text=original)
    metadata = parse_mapping_metadata(result.mapping_text)
    assert metadata["input"]["text"] == original
    assert metadata["input"]["labels"] == labels
    assert metadata["result"]["position"] == result.position
    assert metadata["result"]["score"] == result.score
    assert metadata["software"]["numpy"] == np.__version__
    assert len(metadata["software"]["eufod_source_sha256"]) == 64
    restored = parse_input_text(metadata["input"]["text"])
    repeat = search_mapping(*restored[1:], offset=metadata["settings"]["offset"],
                            step=metadata["settings"]["step"], metric=metric,
                            include_mapping=False)
    assert repeat.position == result.position
    assert repeat.score == result.score
    parsed_metric, positions, scores = parse_mapping_text(result.mapping_text)
    assert parsed_metric == metadata["metric"]
    assert len(positions) == metadata["settings"]["valid_count"]
    assert np.isfinite(scores).all()


def test_exact_input_text_preserved(exported):
    original, mapping = exported
    assert parse_mapping_metadata(mapping)["input"]["text"] == original


def test_api_export_synthesizes_reusable_input():
    result = search_mapping([1, 0.5], [0.25, 0, 0], [[100, 0, 0], [0, 100, 0]], step=100)
    metadata = parse_mapping_metadata(result.mapping_text)
    assert parse_input_text(metadata["input"]["text"])[0] == ["H-1", "H-2"]


def test_legacy_mapping_has_no_metadata():
    text = "# metric: r_factor\n0 0 300 2\n"
    assert parse_mapping_metadata(text) is None
    assert parse_mapping_text(text)[2][0] == 2


@pytest.mark.parametrize("change, message", [
    (lambda text: text.replace("# format_version: 2", "# format_version: 99"), "Unsupported"),
    (lambda text: text.replace("# METADATA_END", "# missing end"), "incomplete"),
    (lambda text: text.replace('"format_version": 2', '"format_version": NaN'), "Nonfinite"),
])
def test_invalid_metadata_rejected(exported, change, message):
    with pytest.raises(ValueError, match=message):
        parse_mapping_metadata(change(exported[1]))


def test_input_checksum_mismatch_rejected(exported):
    metadata = parse_mapping_metadata(exported[1])
    metadata["input"]["text"] += "# changed\n"
    text = "# format_version: 2\n# METADATA_BEGIN\n# " + json.dumps(metadata) + "\n# METADATA_END\n"
    with pytest.raises(ValueError, match="checksum"):
        parse_mapping_metadata(text)


def test_cannot_export_unrelated_input_text():
    with pytest.raises(ValueError, match="does not match"):
        search_mapping([1], [0, 0, 0], [[100, 0, 0]],
                       input_text="O 0 0 0\nH 1 200 0 0\n")


def rewrite_metadata(mapping, metadata):
    header = mapping.split("# METADATA_BEGIN\n")[0]
    rows = mapping.split("# METADATA_END\n")[1]
    return header + "# METADATA_BEGIN\n# " + json.dumps(metadata) + "\n# METADATA_END\n" + rows


@pytest.mark.parametrize("field, value", [
    ("score", "7.0"), ("score", True), ("position", ["1", "2", "3"]),
])
def test_result_metadata_requires_numbers(exported, field, value):
    metadata = parse_mapping_metadata(exported[1])
    metadata["result"][field] = value
    with pytest.raises(ValueError, match="Invalid result metadata"):
        parse_mapping_document(rewrite_metadata(exported[1], metadata))


@pytest.mark.parametrize("field, value", [
    ("step", "50"), ("offset", True), ("valid_count", -1),
    ("candidate_count", 1.5), ("normalization_reference", "other-proton"),
])
def test_invalid_settings_metadata(exported, field, value):
    metadata = parse_mapping_metadata(exported[1])
    metadata["settings"][field] = value
    with pytest.raises(ValueError, match="Invalid mapping metadata"):
        parse_mapping_document(rewrite_metadata(exported[1], metadata))


def test_document_rejects_truncation_even_with_best_row(exported):
    metadata = parse_mapping_metadata(exported[1])
    header, rows = exported[1].split("# METADATA_END\n")
    best_row = next(row for row in rows.splitlines()
                    if list(map(float, row.split()[:3])) == metadata["result"]["position"])
    with pytest.raises(ValueError, match="row count"):
        parse_mapping_document(header + "# METADATA_END\n" + best_row + "\n")


def test_document_rejects_best_score_mismatch(exported):
    metadata = parse_mapping_metadata(exported[1])
    metadata["result"]["score"] += 1
    with pytest.raises(ValueError, match="best result"):
        parse_mapping_document(rewrite_metadata(exported[1], metadata))


def test_complete_document_and_legacy_files(exported):
    metric, positions, scores, metadata = parse_mapping_document(exported[1])
    assert metric == metadata["metric"]
    assert len(positions) == len(scores) == metadata["settings"]["valid_count"]
    assert parse_mapping_document("# metric: r_factor\n0 0 300 2\n")[3] is None
