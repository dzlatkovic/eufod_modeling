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

The first proton is the normalization reference and normally has an experimental relative shift of `1.000`. The coordinates in centiångström (cÅ).

```text
# Coordinate unit: cÅ
# Oxygen Center
O   -797  -153  -182

# Proton Assignments (Label, Exp_Shift, X, Y, Z)
H-1   1.000  -701    26  -183
H-2   0.623  -744  -133    78
H-3a  0.200  -589    35   183
```

Comments begin with `#`. See [`data/example_input.txt`](data/example_input.txt) for the complete menthol example.

## 3D result viewer

The **3D Result Viewer** page accepts a EuFOD mapping TXT and the matching optimized SDF/MOL geometry. Eu coordinates are converted from cÅ to Å automatically. The molecular geometry must use the same coordinate frame as the input used for the calculation.

The viewer displays the selected best positions as a green-to-red sphere cloud; cyan marks the optimal position.

## Case studies

- [`notebooks/borneol_titration_to_input.ipynb`](notebooks/borneol_titration_to_input.ipynb)
- [`notebooks/isoborneol_titration_to_input.ipynb`](notebooks/isoborneol_titration_to_input.ipynb)

Both notebooks document the path from EuFOD titration data to normalized shifts, optimized geometry and EuFOD input.

## Development

```bash
pytest
python benchmarks/benchmark.py
```

The repository retains reference loop, list-comprehension, and NumPy-vectorized implementations for comparison.

## Citation

If you use the underlying chemistry problem, cite:

> Zlatković, D.; Đorđević Zlatković, M.; Radulović, N. *Problem-Solving with Python: Modeling of Lanthanide-Shift Reagent Complexes.* **J. Chem. Educ.** 2023, 100, 3620–3625. [10.1021/acs.jchemed.3c00613](https://doi.org/10.1021/acs.jchemed.3c00613)
