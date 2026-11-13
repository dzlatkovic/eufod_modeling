"""Exercise calculations and uploaded mappings through Streamlit's test runner."""

from pathlib import Path
from unittest.mock import patch

import pytest
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("criterion", ["R-factor", "Pearson correlation"])
def test_calculation_and_stale_result_invalidation(criterion):
    app = AppTest.from_file(str(ROOT / "app.py")).run()
    assert not app.exception
    app.main.radio[0].set_value("Example input").run()
    assert not app.exception
    app.sidebar.radio[0].set_value(criterion)
    app.sidebar.number_input[1].set_value(50)
    app.button[0].click().run(timeout=30)
    assert not app.exception
    assert app.session_state.last_result is not None
    from eufod.io import parse_mapping_metadata

    metadata = parse_mapping_metadata(app.session_state.last_result[0].mapping_text)
    assert metadata["input"]["text"].encode("utf-8") == (ROOT / "data/example_input.txt").read_bytes()
    app.sidebar.number_input[0].set_value(200).run()
    assert app.session_state.last_result is None


def test_theory_page():
    app = AppTest.from_file(str(ROOT / "pages/1_Theory.py")).run()
    assert not app.exception


@pytest.mark.parametrize("row_count", [1, 2])
def test_viewer_accepts_small_mapping(row_count):
    class Upload:
        def __init__(self, data):
            self.data = data
            self.size = len(data)

        def getvalue(self):
            return self.data

    mapping = "# metric: r_factor\n" + "0 0 300 1.0\n" * row_count
    geometry = (ROOT / "outputs/optimized_geometries/borneol_mmff94_optimized.sdf").read_bytes()
    with patch("streamlit.file_uploader", side_effect=[Upload(mapping.encode()), Upload(geometry)]):
        app = AppTest.from_file(str(ROOT / "pages/2_3D_Result_Viewer.py")).run()
    assert not app.exception
    assert app.metric[1].value == str(row_count)


@pytest.mark.parametrize("with_geometry", [False, True])
def test_metadata_viewer_optional_geometry(with_geometry):
    from eufod.io import load_input, parse_mapping_metadata
    from eufod.search import search_mapping

    labels, exp, oxygen, hydrogens = load_input(ROOT / "data/example_input.txt")
    result = search_mapping(exp, oxygen, hydrogens, step=100, labels=labels)

    class Upload:
        size = len(result.mapping_text.encode("utf-8"))

        def getvalue(self):
            return result.mapping_text.encode("utf-8")

    class Geometry:
        data = (ROOT / "outputs/optimized_geometries/borneol_mmff94_optimized.sdf").read_bytes()
        size = len(data)

        def getvalue(self):
            return self.data

    with patch("streamlit.file_uploader", side_effect=[Upload(), Geometry() if with_geometry else None]):
        app = AppTest.from_file(str(ROOT / "pages/2_3D_Result_Viewer.py")).run()
    assert not app.exception
    assert not app.error
    assert app.dataframe[0].value["Proton"].tolist() == labels
    if with_geometry:
        count = parse_mapping_metadata(result.mapping_text)["settings"]["valid_count"]
        assert app.metric[1].value == f"{count:,}"
    else:
        assert "Add a matching SDF" in app.info[1].value
