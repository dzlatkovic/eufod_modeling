"""In-browser presentation of the Eu(fod)3 theory notebook."""

from __future__ import annotations

import json
from pathlib import Path

import streamlit as st


ROOT = Path(__file__).resolve().parents[1]
THEORY_NOTEBOOK = ROOT / "notebooks" / "eufod_theory.ipynb"
THEORY_DIAGRAM = ROOT / "notebooks" / "assets" / "eufod_geometry.svg"


st.set_page_config(
    page_title="EuFOD Theory",
    page_icon="⚗️",
    layout="wide",
)

st.title("EuFOD Theory and Methodology")
st.caption("Mathematical framework for Eu(fod)₃ position modeling.")

if not THEORY_NOTEBOOK.exists():
    st.error("The theory notebook is not available.")
    st.stop()

notebook = json.loads(THEORY_NOTEBOOK.read_text(encoding="utf-8"))

for cell in notebook["cells"]:
    if cell["cell_type"] != "markdown":
        continue

    source = "".join(cell["source"])
    has_diagram = "eufod_geometry.svg" in source
    if has_diagram:
        source = source.replace(
            "![Eu–O–H Geometry](assets/eufod_geometry.svg)",
            "",
        )

    st.markdown(source)

    if has_diagram:
        st.image(str(THEORY_DIAGRAM))
