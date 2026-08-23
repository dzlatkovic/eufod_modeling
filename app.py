"""Streamlit interface for Eu(fod)3 position modeling."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from eufod.io import load_input_bytes
from eufod.search import search_mapping


ROOT = Path(__file__).resolve().parent
EXAMPLE_INPUT = ROOT / "data" / "example_input.txt"


st.set_page_config(
    page_title="Eu(fod)3 Modeling",
    page_icon="⚗️",
    layout="wide",
)

st.title("Eu(fod)₃ Position Modeling")
st.caption(
    "Search a three-dimensional Eu position grid using "
    "experimental proton chemical-shift data."
)

with st.sidebar:
    st.header("Search parameters")

    metric_label = st.radio(
        "Agreement criterion",
        ["R-factor", "Pearson correlation"],
        help=(
            "R-factor is minimized. Pearson correlation is maximized. "
            "They are different scientific objectives."
        ),
    )

    offset = st.number_input(
        "Search offset",
        min_value=1,
        max_value=1000,
        value=300,
        step=10,
    )

    step = st.number_input(
        "Grid step",
        min_value=5,
        max_value=100,
        value=5,
        step=5,
        help="Smaller values give a finer but more expensive search.",
    )

    include_mapping = st.checkbox(
        "Generate full coordinate mapping",
        value=True,
        help=(
            "Creates a tab-separated TXT file containing every "
            "valid candidate position."
        ),
    )

st.subheader("Input data")

input_mode = st.radio(
    "Input source",
    ["Example input", "Upload TXT file"],
    horizontal=True,
)

if input_mode == "Example input":
    input_bytes = EXAMPLE_INPUT.read_bytes()

    with st.expander("Show example input"):
        st.code(
            EXAMPLE_INPUT.read_text(encoding="utf-8"),
            language="text",
        )
else:
    uploaded = st.file_uploader(
        "Upload experimental shifts and coordinates",
        type=["txt"],
    )

    if uploaded is None:
        st.info("Upload a TXT file to continue.")
        st.stop()

    input_bytes = uploaded.getvalue()

try:
    experimental, oxygen, hydrogens = load_input_bytes(input_bytes)
except ValueError as exc:
    st.error(f"Invalid input: {exc}")
    st.stop()

col1, col2, col3 = st.columns(3)
col1.metric("Experimental shifts", len(experimental))
col2.metric("Hydrogen atoms", len(hydrogens))
col3.metric("Grid step", step)

st.write("Oxygen coordinates:", oxygen.tolist())

grid_points = ((2 * offset) // step) ** 3

st.write(
    f"Approximate search size: **{grid_points:,} candidate positions** "
    "(before distance filtering)."
)

if st.button("Run calculation", type="primary"):
    progress = st.progress(
        0.0,
        text="Searching candidate Eu positions...",
    )

    try:
        result = search_mapping(
            experimental,
            oxygen,
            hydrogens,
            offset=int(offset),
            step=int(step),
            metric=(
                "r_factor"
                if metric_label == "R-factor"
                else "pearson"
            ),
            include_mapping=include_mapping,
            progress_callback=lambda value: progress.progress(
                value,
                text=f"Searching... {value:.0%}",
            ),
        )
    except Exception as exc:
        progress.empty()
        st.exception(exc)
        st.stop()

    progress.progress(
        1.0,
        text="Calculation complete.",
    )

    score_name = (
        "R-factor"
        if metric_label == "R-factor"
        else "Pearson r"
    )

    st.success("Calculation complete.")

    result_col1, result_col2 = st.columns(2)

    result_col1.metric(
        "Best Eu position",
        str(result.position),
    )

    result_col2.metric(
        score_name,
        f"{result.score:.6f}",
    )

    if metric_label == "R-factor":
        st.caption(
            "Lower R-factor values indicate better agreement."
        )
    else:
        st.caption(
            "Higher Pearson correlation indicates stronger "
            "pattern agreement."
        )

    if include_mapping and result.mapping_text is not None:
        st.download_button(
            "Download full coordinate mapping",
            data=result.mapping_text.encode("utf-8"),
            file_name=(
                "eu_position_mapping_r_factor.txt"
                if metric_label == "R-factor"
                else "eu_position_mapping_pearson.txt"
            ),
            mime="text/plain",
        )