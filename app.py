"""Streamlit interface for Eu(fod)3 position modeling."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from eufod.io import load_input_bytes
from eufod.search import predicted_shifts, search_mapping


ROOT = Path(__file__).resolve().parent
EXAMPLE_INPUT = ROOT / "data" / "example_input.txt"
MAX_FILE_SIZE_MB = 10
MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_MB * 1024 * 1024


st.set_page_config(
    page_title="Eu(fod)3 Modeling",
    page_icon="⚗️",
    layout="wide",
)

st.title("EuFOD Position Modeling")
st.caption(
    "Search a three-dimensional Eu position grid using experimental proton chemical-shift data. "
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
        "Search offset (cÅ)",
        min_value=1,
        max_value=300,
        value=300,
        step=10,
        help="300 cÅ corresponds to 3.00 Å from the oxygen center.",
    )

    step = st.number_input(
        "Grid step (cÅ)",
        min_value=5,
        max_value=50,
        value=5,
        step=5,
        help="5 cÅ corresponds to 0.05 Å. Smaller values give a finer but more expensive search.",
    )

    include_mapping = st.checkbox(
        "Generate full coordinate mapping",
        value=True,
        help=(
            "Creates a tab-separated TXT file containing every "
            "valid candidate position."
        ),
    )

    st.divider()

    st.markdown(
        "### Scientific Reference\n"
        "This tool models paramagnetic shifts based on the problem described in:\n\n"
        "> Zlatković, D.; Đorđević Zlatković, M.; Radulović, N.  \n"
        "> *Problem-Solving with Python: Modeling of Lanthanide-Shift Reagent Complexes.*  \n"
        "> **J. Chem. Educ.** 2023, 100, 3620–3625.  \n"
        "> 🔗 [DOI: 10.1021/acs.jchemed.3c00613](https://doi.org/10.1021/acs.jchemed.3c00613)"
    )
    st.page_link("pages/1_Theory.py", label="Read the theory", icon="📖")
    st.link_button(
        "See how to prepare input data",
        "https://github.com/dzlatkovic/eufod_modeling/blob/master/notebooks/borneol_titration_to_input.ipynb",
        icon="🧪",
        help=(
            "Open the borneol case-study notebook: titration data, Δδ calculation, "
            "normalization, coordinate assignment, and EuFOD input generation."
        ),
    )

st.subheader("Input data")

input_mode = st.radio(
    "Input source",
    ["Upload TXT file", "Example input"],
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
    st.info(
        "**Input File Format:** Plain-text records containing:\n"
        "1. One oxygen record: `O X Y Z` (coordinates in cÅ)\n"
        "2. One proton-assignment record per proton: "
        "`Label Exp_Shift X Y Z` (coordinates in cÅ)\n\n"
        "Lines beginning with `#` are comments. "
        "💡 *Tip: Switch to the **Example input** option above to see a sample dataset derived from menthol.*  \n"
        "🧪 [Worked borneol example: titration data to EuFOD input]"
        "(https://github.com/dzlatkovic/eufod_modeling/blob/master/notebooks/borneol_titration_to_input.ipynb)"
    )
    uploaded = st.file_uploader(
    "Upload experimental shifts and coordinates",
    type=["txt"],
        )

    if uploaded is None:
        st.info("Upload a TXT file to continue.")
        st.stop()

    if uploaded.size > MAX_FILE_SIZE_BYTES:
        st.error(
            f"File size exceeds the {MAX_FILE_SIZE_MB} MB limit. "
            f"Uploaded file is {uploaded.size / (1024 * 1024):.2f} MB."
        )
        st.stop()

    input_bytes = uploaded.getvalue()



try:
    labels, experimental, oxygen, hydrogens = load_input_bytes(input_bytes)
except ValueError as exc:
    st.error(f"Invalid input: {exc}")
    st.stop()

col1, col2, col3 = st.columns(3)
col1.metric("Experimental shifts", len(experimental))
col2.metric("Hydrogen atoms", len(hydrogens))
col3.metric("Grid step (cÅ)", step)

st.caption("Oxygen center (cÅ)")
oxygen_x, oxygen_y, oxygen_z = st.columns(3)
oxygen_x.metric("O x", f"{oxygen[0]:g}")
oxygen_y.metric("O y", f"{oxygen[1]:g}")
oxygen_z.metric("O z", f"{oxygen[2]:g}")

assignment_column, _ = st.columns([1, 3])
with assignment_column:
    st.dataframe(
        {"Proton": labels, "Experimental shift": experimental},
        hide_index=True,
        use_container_width=True,
    )

# Correct inclusive grid point count: (2 * offset // step + 1)^3
steps_per_axis = (2 * offset) // step + 1
grid_points = steps_per_axis ** 3

st.write(
    f"Approximate search size: **{grid_points:,} candidate positions** "
    "(before distance filtering)."
)

# Initialize calculation result state
if "last_result" not in st.session_state:
    st.session_state.last_result = None

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
        st.session_state.last_result = (result, metric_label)
    except Exception as exc:
        progress.empty()
        st.exception(exc)
        st.stop()

    progress.progress(
        1.0,
        text="Calculation complete.",
    )

# Render results if available in session state
if st.session_state.last_result is not None:
    result, active_metric = st.session_state.last_result
    score_name = "R-factor" if active_metric == "R-factor" else "Pearson r"

    st.success("Calculation complete.")

    result_col1, result_col2 = st.columns(2)
    result_col1.metric("Best Eu position (cÅ)", str(result.position))
    if result.mapping_text is not None:
        result_col1.download_button(
            "Download full coordinate mapping",
            data=result.mapping_text.encode("utf-8"),
            file_name=(
                "eu_position_mapping_r_factor.txt"
                if active_metric == "R-factor"
                else "eu_position_mapping_pearson.txt"
            ),
            mime="text/plain",
        )
    result_col2.metric(score_name, f"{result.score:.6f}")

    if active_metric == "R-factor":
        st.caption("Lower R-factor values indicate better agreement.")
    else:
        st.caption("Higher Pearson correlation indicates stronger pattern agreement.")

    predicted = predicted_shifts(result.position, oxygen, hydrogens)
    chart_rows = [
        {
            "label": label,
            "experimental": float(experimental_value),
            "predicted": float(predicted_value),
        }
        for label, experimental_value, predicted_value in zip(
            labels,
            experimental,
            predicted,
        )
    ]
    shift_min = min(
        row["experimental"]
        for row in chart_rows
    )
    shift_min = min(
        shift_min,
        min(row["predicted"] for row in chart_rows),
    )
    shift_max = max(
        row["experimental"]
        for row in chart_rows
    )
    shift_max = max(
        shift_max,
        max(row["predicted"] for row in chart_rows),
    )
    padding = max((shift_max - shift_min) * 0.05, 0.05)
    axis_min = shift_min - padding
    axis_max = shift_max + padding

    chart_column, shift_table_column = st.columns([1, 1], gap="small")
    chart_column.subheader("Experimental vs. predicted shifts")
    chart_column.vega_lite_chart(
        {
            "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
            "width": 520,
            "height": 520,
            "data": {"values": chart_rows},
            "layer": [
                {
                    "data": {
                        "values": [
                            {
                                "experimental": axis_min,
                                "predicted": axis_min,
                            },
                            {
                                "experimental": axis_max,
                                "predicted": axis_max,
                            },
                        ]
                    },
                    "mark": {
                        "type": "line",
                        "color": "#888888",
                        "strokeDash": [6, 4],
                    },
                    "encoding": {
                        "x": {
                            "field": "experimental",
                            "type": "quantitative",
                            "scale": {"domain": [axis_min, axis_max]},
                        },
                        "y": {
                            "field": "predicted",
                            "type": "quantitative",
                            "scale": {"domain": [axis_min, axis_max]},
                        },
                    },
                },
                {
                    "mark": {"type": "point", "filled": True, "size": 90},
                    "encoding": {
                        "x": {
                            "field": "experimental",
                            "type": "quantitative",
                            "title": "Experimental relative shift",
                            "scale": {"domain": [axis_min, axis_max]},
                        },
                        "y": {
                            "field": "predicted",
                            "type": "quantitative",
                            "title": "Predicted relative shift",
                            "scale": {"domain": [axis_min, axis_max]},
                        },
                        "tooltip": [
                            {"field": "label", "title": "Proton"},
                            {
                                "field": "experimental",
                                "title": "Experimental",
                                "format": ".4f",
                            },
                            {
                                "field": "predicted",
                                "title": "Predicted",
                                "format": ".4f",
                            },
                        ],
                    },
                },
            ],
        },
        use_container_width=False,
    )
    shift_table_column.subheader("Shift comparison")
    shift_table_column.dataframe(
        {
            "Proton": labels,
            "Experimental": experimental,
            "Calculated": predicted,
        },
        hide_index=True,
        use_container_width=True,
    )
