"""Benchmark implementation and objective comparisons."""

from pathlib import Path
import statistics
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eufod.io import load_input
from eufod.reference import solve as solve_for
from eufod.comprehensions import solve as solve_comprehensions
from eufod.search import solve_r_factor, solve_pearson


REPEATS = 5
OFFSET = 300
STEP = 10


def timed(fn, *args):
    """Return median runtime and final result."""
    fn(*args)
    times = []
    result = None

    for _ in range(REPEATS):
        start = time.perf_counter()
        result = fn(*args)
        times.append(time.perf_counter() - start)

    return statistics.median(times), result


def main():
    """Run the benchmark."""
    _labels, exp, oxygen, hydrogens = load_input(
        ROOT / "data/example_input.txt"
    )

    rows = []

    implementations = [
        ("Python for-loop", solve_for, "R-factor"),
        ("Comprehensions", solve_comprehensions, "R-factor"),
        ("NumPy vectorized", solve_r_factor, "R-factor"),
        ("NumPy vectorized", solve_pearson, "Pearson"),
    ]

    print("Eu(fod)3 benchmark")
    print("=" * 75)
    print(
        f"Grid: offset={OFFSET}, step={STEP}; "
        f"repetitions={REPEATS}"
    )
    print()

    for name, fn, objective in implementations:
        runtime, result = timed(
            fn,
            exp,
            oxygen,
            hydrogens,
            OFFSET,
            STEP,
        )
        rows.append((name, objective, runtime, result))

    print(
        f"{'Implementation':<22}"
        f"{'Objective':<12}"
        f"{'Median time':>14}  Result"
    )
    print("-" * 90)

    for name, objective, runtime, result in rows:
        print(
            f"{name:<22}"
            f"{objective:<12}"
            f"{runtime:>11.4f} s  "
            f"{result}"
        )

    for_loop_time = rows[0][2]
    comprehension_time = rows[1][2]
    numpy_r_time = rows[2][2]
    numpy_pearson_time = rows[3][2]

    print()
    print("Performance comparison")
    print("-" * 75)
    print(
        "NumPy R-factor vs for-loop:       "
        f"{for_loop_time / numpy_r_time:.2f}x"
    )
    print(
        "NumPy R-factor vs comprehensions: "
        f"{comprehension_time / numpy_r_time:.2f}x"
    )
    print(
        "NumPy Pearson vs NumPy R-factor:   "
        f"{numpy_r_time / numpy_pearson_time:.2f}x"
    )

    print()
    print("Scientific comparison")
    print("-" * 75)
    print("R-factor: lower is better.")
    print("Pearson: higher is better.")
    print("Pearson is an alternative objective, not a replacement for R-factor.")

    if rows[0][3] == rows[2][3]:
        print("For-loop and NumPy R-factor: identical result.")
    else:
        print("WARNING: for-loop and NumPy R-factor differ.")

    if rows[1][3] == rows[2][3]:
        print("Comprehensions and NumPy R-factor: identical result.")
    else:
        print("WARNING: comprehensions and NumPy R-factor differ.")


if __name__ == "__main__":
    main()
