from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eufod.io import load_input
from eufod.reference import solve

_labels, exp, oxygen, hydrogens = load_input(ROOT / "data/example_input.txt")
print("R-factor:", solve(exp, oxygen, hydrogens))
