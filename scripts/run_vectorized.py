from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eufod.io import load_input
from eufod.search import solve_r_factor, solve_pearson

_labels, exp, oxygen, hydrogens = load_input(ROOT / "data/example_input.txt")

print("R-factor:", solve_r_factor(exp, oxygen, hydrogens))
print("Pearson:", solve_pearson(exp, oxygen, hydrogens))
