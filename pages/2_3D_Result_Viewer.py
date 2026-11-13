"""Interactive visualization of an exported EuFOD coordinate mapping."""

from __future__ import annotations

import numpy as np
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from eufod.io import parse_mapping_document


MAX_UPLOAD_MB = 10
MAX_UPLOAD_BYTES = MAX_UPLOAD_MB * 1024 * 1024


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
    "Upload a mapping TXT to inspect its calculation details. Add the matching optimized SDF for 3D visualization. "
    "The SDF must use the same coordinate frame as the EuFOD input. "
    "Eu positions are converted automatically from cÅ to Å for display."
)

mapping_upload = st.file_uploader("EuFOD coordinate mapping (TXT)", type="txt")
sdf_upload = st.file_uploader("Optimized molecular geometry (SDF, optional)", type=["sdf", "mol"])

if mapping_upload is None:
    st.stop()

for uploaded_file, file_description in (
    (mapping_upload, "Coordinate mapping"),
    (sdf_upload, "Molecular geometry"),
):
    if uploaded_file is not None and uploaded_file.size > MAX_UPLOAD_BYTES:
        st.error(
            f"{file_description} exceeds the {MAX_UPLOAD_MB} MB upload limit. "
            "Generate the mapping with a 10 cÅ or coarser grid for 3D visualization."
        )
        st.stop()

try:
    mapping_text = mapping_upload.getvalue().decode("utf-8")
    metric, coordinates_cangstrom, scores, metadata = parse_mapping_document(mapping_text)
except (UnicodeDecodeError, ValueError) as exc:
    st.error(f"Could not read the uploaded result files: {exc}")
    st.stop()

if metadata is not None:
    with st.expander("Calculation details", expanded=True):
        st.caption(f"Calculated {metadata['created_utc']} • {metric}")
        st.write({
            "Grid offset (cÅ)": metadata["settings"]["offset"],
            "Grid step (cÅ)": metadata["settings"]["step"],
            "Minimum distance (cÅ)": metadata["settings"]["minimum_distance"],
            "Best Eu position (cÅ)": metadata["result"]["position"],
            "Best score": metadata["result"]["score"],
        })
        st.dataframe({
            "Proton": metadata["input"]["labels"],
            "Experimental": metadata["input"]["experimental_shifts"],
            "Predicted": metadata["result"]["predicted_shifts"],
            "Residual": metadata["result"]["residuals"],
        }, hide_index=True)
        st.download_button("Download original input", metadata["input"]["text"].encode("utf-8"),
                           file_name="eufod_original_input.txt", mime="text/plain")
    with st.expander("Full calculation metadata"):
        st.json(metadata)
else:
    st.caption("Legacy mapping: calculation metadata is not included in this file.")

if sdf_upload is None:
    st.info("Add a matching SDF/MOL file to display the 3D view. It is not needed to inspect calculation details.")
    st.stop()
try:
    sdf_text = sdf_upload.getvalue().decode("utf-8")
except UnicodeDecodeError as exc:
    st.error(f"Could not read molecular geometry: {exc}")
    st.stop()

ascending = metric == "r_factor"
order = np.argsort(scores) if ascending else np.argsort(-scores)
if metadata is not None:
    # Preserve the full-precision winner when exported scores round to a tie.
    winner = np.flatnonzero(np.all(coordinates_cangstrom == metadata["result"]["position"], axis=1))[0]
    order = np.concatenate(([winner], order[order != winner]))
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
top_count = 1 if top_limit == 1 else st.slider(
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
viewer_html = viewer._make_html()
# Keep the button in the viewer's iframe so it captures the current camera
# without a Streamlit rerun resetting the rotation and zoom.
viewer_html += """
<style>
  #download-png {
    margin-top: 10px; padding: 8px 14px; border: 1px solid #bbb;
    border-radius: 6px; background: white; color: #222;
    font: 14px sans-serif; cursor: pointer;
  }
  #download-png:disabled { opacity: 0.5; cursor: wait; }
  #export-status { margin-left: 10px; font: 14px sans-serif; color: #444; }
</style>
<button id="download-png" type="button" disabled>Download view as PNG</button>
<span id="export-status" role="status" aria-live="polite"></span>
<script>
  const downloadButton = document.getElementById('download-png');
  const exportStatus = document.getElementById('export-status');
  $3Dmolpromise.then(function () {
    if (viewer_VIEWER_ID) downloadButton.disabled = false;
  });
  downloadButton.addEventListener('click', function () {
    try {
      const currentViewer = viewer_VIEWER_ID;
      currentViewer.render();
      const image = currentViewer.pngURI();
      if (!image.startsWith('data:image/png')) {
        throw new Error('No PNG image was produced.');
      }
      const link = document.createElement('a');
      link.href = image;
      link.download = 'eufod_3d_view.png';
      document.body.appendChild(link);
      link.click();
      link.remove();
      exportStatus.textContent = '';
    } catch (error) {
      exportStatus.textContent = 'Could not export the image. Please try again.';
      console.error('EuFOD PNG export failed:', error);
    }
  });
</script>
""".replace("VIEWER_ID", viewer.uniqueid)
st.caption("Rotate and zoom to frame your image, then select Download view as PNG below the viewer.")
components.html(viewer_html, height=720, scrolling=False)
