# Eu(fod)₃ Position Modeling

A Python implementation for modeling lanthanide-induced NMR shifts and determining the position of Eu³⁺ in an Eu(fod)₃–substrate complex.

This project is based on the computational problem described in:

> Zlatković, D.; Đorđević Zlatković, M.; Radulović, N.  
> *Problem-Solving with Python: Modeling of Lanthanide-Shift Reagent Complexes.*  
> **Journal of Chemical Education** 2023, *100*, 3620–3625.  
> DOI: [10.1021/acs.jchemed.3c00613](https://doi.org/10.1021/acs.jchemed.3c00613)

The original work was developed as a programming exercise for chemistry students, using NMR titration data obtained from the interaction of menthol with Eu(fod)₃. The exercise introduces fundamental programming concepts through a chemically meaningful complete-search problem.

This repository develops that original concept further by separating the chemical model from the search algorithms, comparing different implementations, investigating alternative agreement criteria, and providing an interactive Streamlit interface.

---

## Scientific background

Lanthanide shift reagents can produce characteristic paramagnetic shifts in the NMR signals of coordinated substrates. The magnitude of the induced shift depends on the spatial relationship between the lanthanide ion and the observed nucleus.

The original J. Chem. Educ. exercise uses experimental NMR shift data to determine the most probable position of the Eu³⁺ ion relative to the substrate.

The computational problem can be expressed as:

1. Define a three-dimensional search space around the coordinating oxygen.
2. Generate candidate Eu³⁺ positions.
3. Calculate Eu–H and Eu–O distances.
4. Calculate the expected relative shifts for each candidate position.
5. Compare the calculated shifts with the experimental values.
6. Identify the Eu³⁺ position giving the best agreement.

The original implementation uses a complete-search approach. This repository retains that approach as a reference implementation while also exploring more efficient numerical implementations.

---

## From the educational exercise to an optimized implementation

The original problem is particularly useful computationally because it provides a simple example of how the same scientific calculation can be implemented in progressively more efficient ways.

This repository therefore contains three implementations of the R-factor search:

```text
Python for-loops
       ↓
List comprehensions
       ↓
NumPy vectorization
```

The first two implementations operate largely at the Python level, whereas the NumPy implementation evaluates many candidate Eu positions simultaneously using array operations.

This allows the project to demonstrate an important computational chemistry concept:

> The mathematical model can remain unchanged while the computational implementation is substantially optimized.

The benchmark included in the repository compares these implementations using the same input data, search space, and R-factor calculation.

---

## Agreement criteria

### R-factor

The primary agreement criterion is the R-factor:

```text
R = 100 × sqrt(
    mean(
        ((δcalc - δexp) / δexp)²
    )
)
```

Lower values indicate better agreement between calculated and experimental relative shifts.

The R-factor is the primary criterion inherited from the original modeling problem.

### Pearson correlation

The repository also provides Pearson correlation as an alternative agreement criterion.

Pearson correlation evaluates the similarity of the pattern of calculated and experimental shifts rather than their absolute agreement.

Higher values indicate stronger correlation.

These two criteria should not be interpreted as equivalent:

```text
R-factor
    ↓
agreement of calculated shift magnitudes

Pearson correlation
    ↓
similarity of the shift pattern
```

Consequently, they can produce different optimal Eu positions.

This difference is intentionally preserved in the application so that the effect of the choice of objective function can be investigated.

---

# Interactive application

The project includes a Streamlit interface for running the calculation without directly interacting with the Python code.

Run the application with:

```bash
streamlit run app.py
```

The interface allows the user to:

- use the included example dataset;
- upload a custom TXT input file;
- select R-factor or Pearson correlation;
- define the search range;
- define the grid spacing;
- run the Eu-position search;
- view the optimal Eu³⁺ coordinates;
- view the corresponding agreement score;
- download the complete coordinate/agreement mapping.

The default and minimum grid spacing is **5 coordinate units**.


---

## Input format

The program accepts a plain-text file containing:

1. Experimental relative shifts.
2. Oxygen coordinates.
3. Hydrogen coordinates.

For example:

```text
1.00 0.623 0.200 0.158 0.241 0.663 0.628 0.643 0.288 0.239
-797 -153 -182
-701 26 -183
-744 -133 78
...
```

The first line contains the experimental relative shifts.

The second line contains the coordinates of the coordinating oxygen atom.

Each subsequent line contains the coordinates of one hydrogen atom.

The number of experimental shifts must correspond to the number of hydrogen coordinates.

### Input requirements

- All values must be numeric.
- Oxygen coordinates must contain exactly three values.
- Each hydrogen coordinate must contain exactly three values.
- The number of experimental shifts must equal the number of hydrogen coordinates.
- Blank lines and lines beginning with `#` are ignored.

---

## Project structure

```text
eufod_modeling/
│
├── app.py
├── README.md
├── LICENSE
├── requirements.txt
├── pyproject.toml
├── .gitignore
│
├── data/
│   └── example_input.txt
│
├── src/
│   └── eufod/
│       ├── __init__.py
│       ├── io.py
│       ├── model.py
│       ├── reference.py
│       ├── comprehensions.py
│       ├── search.py
│       └── pearson.py
│
├── scripts/
│   ├── run_reference.py
│   ├── run_comprehensions.py
│   └── run_vectorized.py
│
├── benchmarks/
│   └── benchmark.py
│
└── tests/
    ├── test_model.py
    └── test_search.py
```

### `model.py`

Contains the mathematical model and agreement calculations.

### `reference.py`

Contains the straightforward Python implementation using explicit `for` loops.

This implementation is retained as a reference against which optimized implementations can be compared.

### `comprehensions.py`

Contains an implementation using Python list comprehensions.

It demonstrates an intermediate step between explicit loops and NumPy vectorization.

### `search.py`

Contains the NumPy-vectorized search and the mapping-generation functionality.

### `pearson.py`

Provides the Pearson-based search as an alternative objective.

### `io.py`

Handles reading and validation of input data.

### `app.py`

Provides the Streamlit web interface.

### `benchmarks/`

Contains performance comparisons between the different implementations and agreement criteria.

### `tests/`

Contains automated tests verifying the mathematical functions and agreement between the reference and optimized implementations.

---

# Computational approaches

## 1. Python reference implementation

The reference implementation uses explicit Python `for` loops to evaluate candidate Eu positions one at a time.

This implementation prioritizes clarity and provides a transparent representation of the original computational problem.

It is intentionally retained even though it is slower than the vectorized implementation.

Run it with:

```bash
python scripts/run_reference.py
```

---

## 2. List-comprehension implementation

The comprehension implementation expresses several operations using Python list comprehensions.

It represents an intermediate implementation between explicit loops and NumPy vectorization.

Run it with:

```bash
python scripts/run_comprehensions.py
```

List comprehensions do not necessarily provide a significant performance improvement for numerical calculations when the expensive operations are still executed individually at the Python level.

---

## 3. NumPy-vectorized implementation

The NumPy implementation evaluates many candidate Eu positions simultaneously.

Instead of calculating distances and agreement values for each position independently in Python, the calculation is expressed using NumPy array operations.

This substantially reduces Python-level iteration and can provide a large performance improvement for the complete-search problem.

Run it with:

```bash
python scripts/run_vectorized.py
```

---

# Benchmarking

The benchmark is designed to answer two separate questions.

## Implementation performance

The first comparison is:

```text
Python for-loop
       ↓
List comprehensions
       ↓
NumPy vectorization
```

The `for` loop, comprehension, and NumPy R-factor implementations solve the same optimization problem.

The benchmark compares their execution times and verifies that the optimized R-factor implementation produces the same optimum as the reference implementations.

Run:

```bash
python benchmarks/benchmark.py
```

The benchmark reports execution times over multiple repetitions.

### Why compare comprehensions?

List comprehensions are often suggested as a faster alternative to explicit Python loops. In numerical problems, however, the actual performance improvement depends on where the computational work is performed.

Including the comprehension implementation therefore provides a useful intermediate benchmark between straightforward Python and vectorized numerical computation.

---

## Agreement-criterion comparison

The benchmark also compares:

```text
NumPy + R-factor
       vs
NumPy + Pearson correlation
```

This comparison has a different purpose.

R-factor and Pearson correlation optimize different mathematical objectives and therefore may identify different Eu³⁺ positions.

### R-factor

Minimize:

```text
R → 0
```

### Pearson

Maximize:

```text
r → 1
```

The two results should therefore not be interpreted as competing implementations of the same objective.

---

# Search grid

The Eu³⁺ position is determined by searching a three-dimensional Cartesian grid around the coordinating oxygen atom.

The search is controlled by two parameters:

### Offset

Defines the extent of the search region around the oxygen coordinates.

For example:

```text
offset = 300
```

searches a region extending approximately 300 coordinate units in each direction.

### Step

Defines the distance between neighboring candidate Eu positions.

The default and minimum step is:

```text
step = 5
```

A smaller step produces a finer spatial search but increases computational cost substantially.

Because the search is three-dimensional, reducing the step can increase the number of candidate positions very rapidly.

Therefore, the grid resolution should be selected according to the desired balance between spatial resolution and computational cost.

---

# Distance constraints

Candidate Eu positions are rejected when the Eu–O distance or any Eu–H distance falls below the minimum allowed distance defined by the model.

This prevents physically unrealistic positions from being considered during the search.

The precise distance constraint is defined in the implementation and should be kept consistent when comparing the different computational approaches.

---

# Mapping output

In addition to returning the optimal Eu position, the application can generate a complete mapping of the search space.

The mapping contains:

```text
x    y    z    score
```

For R-factor searches:

```text
x    y    z    r_factor
```

For Pearson searches:

```text
x    y    z    pearson_r
```

Only valid candidate positions are included.

The mapping is exported as a tab-separated text file, making it suitable for further analysis or visualization.

The complete mapping can be used to investigate the shape of the agreement landscape rather than only the global optimum.

---

# Installation

Clone the repository:

```bash
git clone <repository-url>
cd eufod_modeling
```

Create and activate a Python environment if desired.

Then install the project:

```bash
pip install -e .
```

For development and testing:

```bash
pip install -e ".[dev]"
```

Alternatively, install the dependencies directly:

```bash
pip install -r requirements.txt
```

---

# Running the application

Start the Streamlit application from the repository root:

```bash
streamlit run app.py
```

The application will open in your default web browser.

To stop the application, return to the terminal and press:

```text
Ctrl + C
```

---

# Running the command-line implementations

Reference implementation:

```bash
python scripts/run_reference.py
```

Comprehension implementation:

```bash
python scripts/run_comprehensions.py
```

NumPy implementation:

```bash
python scripts/run_vectorized.py
```

---

# Testing

The project includes automated tests using `pytest`.

Run:

```bash
pytest
```

The tests verify the mathematical functions and ensure that the reference and vectorized R-factor searches return consistent results.

---

# Reproducibility

The repository includes an example input dataset so that the calculation can be reproduced without additional data.

The same input and search parameters can be used with the different implementations to compare:

- calculated optimum;
- agreement score;
- execution time.

When comparing implementations, the search offset, grid step, input data, and distance constraints should be kept identical.

---

# Relationship to the J. Chem. Educ. paper

This repository is an extension of the computational problem presented in the original educational work rather than a reproduction of the paper's complete teaching environment.

The original paper presents a Python program based on a simple complete-search algorithm for determining the geometry of a lanthanide–substrate complex. The problem was designed as a practical programming project for chemistry students and used NMR titration data from the interaction of menthol with Eu(fod)₃.

The present repository places greater emphasis on software engineering and computational performance, including:

- separation of the chemical model from the search algorithm;
- reference and optimized implementations;
- list-comprehension comparison;
- NumPy vectorization;
- alternative agreement criteria;
- automated testing;
- performance benchmarking;
- reproducible input/output;
- downloadable search mappings;
- and an interactive Streamlit application.

The underlying chemical modeling problem remains based on the original work.

---

# Citation

If you use the underlying Eu(fod)₃ modeling problem or this project in academic work, please cite:

**Zlatković, D.; Đorđević Zlatković, M.; Radulović, N.**  
*Problem-Solving with Python: Modeling of Lanthanide-Shift Reagent Complexes.*  
**Journal of Chemical Education** 2023, *100*, 3620–3625.  
DOI: [10.1021/acs.jchemed.3c00613](https://doi.org/10.1021/acs.jchemed.3c00613)

---

# Disclaimer

This software is intended for computational exploration and educational use.

The calculated Eu³⁺ position depends on the underlying chemical model, experimental data, search-space definition, grid resolution, distance constraints, and agreement criterion.

The numerical optimum should therefore be interpreted in the context of the chemical model and experimental system rather than as an independent determination of molecular structure.
