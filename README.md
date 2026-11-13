# EuFOD Position Modeling

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-App-FF4B4B.svg)](https://streamlit.io/)

EuFOD is a Python and Streamlit application for modelling paramagnetic $^1$H NMR shifts induced by Eu(fod)₃ and locating the best-fitting Eu position around a coordinating oxygen atom.

It expands on the educational problem described by Zlatković, Đorđević Zlatković, and Radulović in [*Journal of Chemical Education* 2023, 100, 3620–3625](https://doi.org/10.1021/acs.jchemed.3c00613).

## Quick start

```bash
git clone https://github.com/dzlatkovic/eufod_modeling.git
cd eufod_modeling
pip install -e ".[dev]"
streamlit run app.py
```

## What it does

- Searches a three-dimensional Eu grid using either minimum R-factor or maximum Pearson correlation.
- Compares experimental and predicted relative shifts with proton-aware tooltips.
- Exports the full valid-coordinate mapping.
- Visualizes an exported mapping as a colour-coded Eu-position cloud around an optimized SDF structure.
- Includes menthol input data and reproducible borneol/isoborneol titration case studies.

## Input format

The first proton is the normalization reference and normally has an experimental relative shift of `1.000`. Coordinates are in centiångströms (cÅ). Decimal coordinates are preserved.

```text
# Coordinate unit: cÅ
# Oxygen Center
O   -797  -153  -182

# Proton Assignments (Label, Exp_Shift, X, Y, Z)
H-1   1.000  -701    26  -183
H-2   0.623  -744  -133    78
H-3a  0.200  -589    35   183
```

All coordinates and shifts must be finite. R-factor requires nonzero experimental shifts; Pearson requires at least two nonconstant shifts. Candidates with undefined normalized shifts or scores are excluded. A search with no usable candidates reports an error.

The search interval on each axis is `[oxygen - offset, oxygen + offset)`: the upper boundary is excluded. All three implementations use the same grid.

Comments begin with `#`. See [`data/example_input.txt`](data/example_input.txt) for the complete menthol example.

## 3D result viewer

The **3D Result Viewer** page accepts a EuFOD mapping TXT. New files show calculation details and allow downloading the original input without a geometry file. Add the matching optimized SDF/MOL geometry for the 3D scene. Eu coordinates are converted from cÅ to Å automatically. The molecular geometry must use the same coordinate frame as the input used for the calculation; the export does not verify or embed an SDF.

The viewer displays the selected best positions as a green-to-red sphere cloud; cyan marks the optimal position.

Use **Download view as PNG** below the 3D viewer to save the current view,
including its rotation, zoom, colours, and visible Eu positions. The image
contains the molecular scene, without the surrounding app controls.

## Extended mapping TXT format

New mapping exports use format version 2. The four numeric columns (`x y z score`)
are unchanged. All metadata is stored in `#` comment lines, so older readers and
`pandas.read_csv(..., comment="#")` can still read the mapping.

After the original seven header lines, `# format_version: 2` and a commented
JSON block between `# METADATA_BEGIN` and `# METADATA_END` record:

- UTC calculation timestamp, metric, and coordinate units.
- Offset, step, grid interval, minimum distance, normalization reference, and
  candidate counts. Distance filtering applies to oxygen and supplied hydrogens.
- Original input text, its SHA-256 checksum, proton labels, experimental shifts,
  and oxygen/hydrogen coordinates.
- Best position and full-precision score, predicted shifts, and residuals
  (`predicted - experimental`). Mapping table scores retain eight decimal places.
- Python and package versions, plus a SHA-256 fingerprint of the installed
  EuFOD Python source files, to distinguish local edits from a release version.

`parse_mapping_text()` retains its existing three-value return. Use
`parse_mapping_metadata()` to read the additional dictionary; it returns `None`
for legacy files and rejects malformed or unsupported metadata. Checksums detect
accidental input changes, not authenticity or chemical correctness.
`parse_mapping_document()` returns `(metric, coordinates, scores, metadata)` and
also checks the row count and best result against the metadata; the viewer uses
this reader to reject incomplete or inconsistent exports.

Python callers can pass `labels=` and `input_text=` to `search_mapping()`;
supplied text must match the calculation arrays. Without original text, a
lossless input is reconstructed, using generated labels if none were supplied.
Metadata is included when full mapping export is enabled. No SDF is required.

## Case studies

- [`notebooks/borneol_titration_to_input.ipynb`](notebooks/borneol_titration_to_input.ipynb)
- [`notebooks/isoborneol_titration_to_input.ipynb`](notebooks/isoborneol_titration_to_input.ipynb)

Install the optional notebook dependencies before executing them:

```bash
python -m pip install -e ".[dev,notebooks]"
```

Alternatively, install the full app, test, and notebook environment from the
repository root with `python -m pip install -r requirements.txt`. This also
installs the local `eufod` package. In PyCharm, select the project's `venv`
interpreter for both Python scripts and notebook kernels.

Run notebooks from the repository root or the `notebooks` directory. The notebooks reuse the saved SDF in `outputs/optimized_geometries` to preserve the coordinate frame of existing mappings. To regenerate a geometry, move both its old SDF and mapping aside and regenerate the mapping from the new input.

Both notebooks document the path from EuFOD titration data to normalized shifts, optimized geometry and EuFOD input.

## Development

```bash
python -m pytest
python benchmarks/benchmark.py
```

The tests cover input validation, independent grid-score comparisons, fractional coordinates, empty search regions, mapping round trips, and Streamlit calculation/viewer flows. The NumPy search processes coordinates in chunks; exporting a full mapping still requires memory proportional to the exported text.

The repository retains reference loop, list-comprehension, and NumPy-vectorized implementations for comparison.

## Citation

If you use the underlying chemistry problem, cite:

> Zlatković, D.; Đorđević Zlatković, M.; Radulović, N. *Problem-Solving with Python: Modeling of Lanthanide-Shift Reagent Complexes.* **J. Chem. Educ.** 2023, 100, 3620–3625. [10.1021/acs.jchemed.3c00613](https://doi.org/10.1021/acs.jchemed.3c00613)
