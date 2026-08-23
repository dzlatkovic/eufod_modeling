# Eu(fod)₃ Position Modeling

[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-App-FF4B4B.svg)](https://streamlit.io/)

A modular Python framework and interactive web application for modeling paramagnetic shifts induced by lanthanide shift reagents (LSR) and determining the 3D spatial position of Eu in EuFOD–substrate complexes.

This repository expands upon the educational problem presented in:
> Zlatković, D.; Đorđević Zlatković, M.; Radulović, N. *Problem-Solving with Python: Modeling of Lanthanide-Shift Reagent Complexes.* **J. Chem. Educ.** 2023, 100 (9), 3620–3625. DOI: [10.1021/acs.jchemed.3c00613](https://doi.org/10.1021/acs.jchemed.3c00613)

---

## Key Features

* **Implementation Comparison:** Benchmarks identical grid searches across native Python `for` loops, list comprehensions, and vectorized NumPy array operations.
* **Dual Optimization Objectives:** Supports both absolute deviation (**R-factor**) and pattern similarity (**Pearson correlation**).
* **Interactive Web UI:** Streamlit application for uploading atomic coordinate datasets, configuring search parameters, and downloading full mapping grids.
* **Modular Engine:** Clean separation of chemical models, computational search algorithms, CLI scripts, and web UI.

---

## Quickstart

### Installation

```bash
git clone https://github.com/dzlatkovic/eufod_modeling.git
cd eufod_modeling
pip install -e ".[dev]"
```

### Interactive Streamlit App

Launch the interactive dashboard to run searches and visualize result maps:

```bash
streamlit run app.py
```

### CLI Implementation Scripts

Compare the execution across different computational paradigms:

```bash
python scripts/run_reference.py     # Native Python loops
python scripts/run_comprehensions.py # List comprehensions
python scripts/run_vectorized.py    # NumPy vectorization
```

---

## Computational Implementations & Benchmarking

The primary search evaluates candidate Eu³⁺ positions across a 3D Cartesian grid around the coordinating oxygen atom. Unphysical candidate positions violating predefined Eu–O or Eu–H distance constraints are filtered out prior to evaluation.

| Implementation | Engine | Key Characteristic |
| :--- | :--- | :--- |
| **Reference** | Pure Python (`for` loops) | Iterates candidate positions sequentially. Clear baseline for model validation. |
| **Comprehension** | Pure Python (`[x for x in ...]`) | Intermediate benchmark evaluating Python-level syntax overhead. |
| **Vectorized** | NumPy (`ndarray`) | Evaluates the entire 3D grid simultaneously using broad-array vectorization. |

Run the full comparative benchmark and test suite:

```bash
# Run execution time comparison
python benchmarks/benchmark.py

# Run test suite verifying mathematical consistency across implementations
pytest
```

---

## Agreement Criteria

1. **R-factor ($R$)** — Primary objective measuring absolute magnitude agreement:
   $$
R = 100 \sqrt{
\operatorname{mean}
\left[
\left(
\frac{\delta_{\mathrm{calc}} - \delta_{\mathrm{exp}}}
{\delta_{\mathrm{exp}}}
\right)^2
\right]
}
$$
2. **Pearson Correlation ($r$)** — Secondary objective evaluating pattern similarity across relative shift magnitudes.

> **Note:** Because $R$ minimizes absolute fractional deviation while $r$ maximizes trend correlation, the two objectives may yield distinct optimal coordinates.

---

## Input File Format

Input files are plain-text containing experimental relative shifts followed by 3D Cartesian coordinates for Oxygen and Hydrogens:

```text
1.00 0.623 0.200 0.158 0.241 0.663 0.628 0.643 0.288 0.239
-797 -153 -182  # Oxygen coordinates (X, Y, Z)
-701   26 -183  # H1 coordinates
-744 -133   78  # H2 coordinates
...
```

An example dataset derived from menthol titration is included in `data/example_input.txt`.

---

## Repository Structure

```text
eufod_modeling/
├── app.py                  # Streamlit web application
├── benchmarks/
│   └── benchmark.py        # Performance timing suite
├── data/
│   └── example_input.txt   # Menthol titration example data
├── scripts/                # CLI runners for each paradigm
├── src/
│   └── eufod/              # Core computational engine
│       ├── model.py        # Mathematical definitions & distance constraints
│       ├── io.py           # Parser & input validation
│       ├── reference.py    # Loop implementation
│       ├── search.py       # NumPy vectorized implementation
│       └── pearson.py      # Pearson-based objective search
└── tests/                  # Pytest verification suites
```

---

## Citation

If you use this codebase or the underlying chemistry problem in academic work, please cite the original J. Chem. Educ. publication:

```bibtex
@article{zlatković2023problem,
  title={Problem-Solving with Python: Modeling of Lanthanide-Shift Reagent Complexes},
  author={Zlatković, Dragan and {\DJ}or{\dj}ević Zlatković, Miljana and Radulovic, Niko},
  journal={Journal of Chemical Education},
  volume={100},
  number={9},
  pages={3620--3625},
  year={2023},
  publisher={ACS Publications}
}
```
