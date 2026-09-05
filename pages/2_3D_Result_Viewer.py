"""Interactive visualization of an exported EuFOD coordinate mapping."""

from __future__ import annotations

import numpy as np
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from eufod.io import parse_mapping_text


def colour_for_difference(difference: float, upper_limit: float) -> str:
    """Return a continuous green-yellow-red colour for a score difference."""
    fraction = 0.0 if upper_limit == 0 else min(difference / upper_limit, 1.0)
    if fraction < 0.5:
        fraction *= 2
        start, end = (44, 160, 44), (255, 192, 0)
    else:
        fraction = (fraction - 0.5) * 2
        start, end = (255, 192, 0), (214, 39, 40)
    return "#" + "".join(
        f"{round(start[index] + fraction * (end[index] - start[index])):02x}"
        for index in range(3)
    )


st.set_page_config(page_title="EuFOD 3D Result Viewer", page_icon="🧭", layout="wide")
st.title("EuFOD 3D Result Viewer")
st.caption("Inspect an exported Eu coordinate mapping around the molecular geometry used for the calculation.")

st.info(
    "Upload the mapping TXT exported by EuFOD and the matching optimized SDF. "
    "The SDF must use the same coordinate frame as the EuFOD input. "
    "Eu positions are converted automatically from cÅ to Å for display."
)

mapping_upload = st.file_uploader("EuFOD coordinate mapping (TXT)", type="txt")
sdf_upload = st.file_uploader("Optimized molecular geometry (SDF)", type=["sdf", "mol"])

if mapping_upload is None or sdf_upload is None:
    st.stop()

try:
    metric, coordinates_cangstrom, scores = parse_mapping_text(
        mapping_upload.getvalue().decode("utf-8")
    )
    sdf_text = sdf_upload.getvalue().decode("utf-8")
except (UnicodeDecodeError, ValueError) as exc:
    st.error(f"Could not read the uploaded result files: {exc}")
    st.stop()

ascending = metric == "r_factor"
order = np.argsort(scores) if ascending else np.argsort(-scores)
coordinates_cangstrom = coordinates_cangstrom[order]
scores = scores[order]
best_score = scores[0]

if metric == "r_factor":
    differences = 100 * (scores - best_score) / max(abs(best_score), 1e-12)
    score_label = "R-factor"
else:
    differences = 100 * (best_score - scores) / max(abs(best_score), 1e-12)
    score_label = "Pearson r"

top_limit = min(1_000, len(scores))
top_count = st.slider(
    "Number of best positions to display",
    min_value=1,
    max_value=top_limit,
    value=min(300, top_limit),
)

display_coordinates = coordinates_cangstrom[:top_count] / 100
display_differences = differences[:top_count]
colour_limit = float(display_differences.max())

controls, metrics = st.columns([3, 2])
with controls:
    sphere_radius = st.slider("Eu sphere radius (Å)", 0.02, 0.20, 0.075, 0.005)
    sphere_opacity = st.slider("Eu sphere opacity", 0.10, 1.00, 0.60, 0.05)
with metrics:
    st.metric(f"Best {score_label}", f"{best_score:.6f}")
    st.metric("Valid positions in mapping", f"{len(scores):,}")
    st.caption(
        "Green is the best score; yellow and red represent progressively "
        "worse scores among the displayed positions. Cyan marks the best position."
    )

viewer = py3Dmol.view(width=900, height=650)
viewer.setBackgroundColor("white")
viewer.addModel(sdf_text, "sdf")
viewer.setStyle({"stick": {"radius": 0.09, "color": "#666666"}})
viewer.addStyle({"elem": "O"}, {"sphere": {"radius": 0.24, "color": "#d62728"}})

for position, difference in zip(display_coordinates, display_differences):
    viewer.addSphere(
        {
            "center": {"x": position[0], "y": position[1], "z": position[2]},
            "radius": sphere_radius,
            "color": colour_for_difference(float(difference), colour_limit),
            "opacity": sphere_opacity,
        }
    )

best_position = display_coordinates[0]
viewer.addSphere(
    {
        "center": {"x": best_position[0], "y": best_position[1], "z": best_position[2]},
        "radius": sphere_radius * 2,
        "color": "#00bcd4",
        "opacity": 1.0,
    }
)
viewer.zoomTo()
components.html(viewer._make_html(), height=660, scrolling=False)
