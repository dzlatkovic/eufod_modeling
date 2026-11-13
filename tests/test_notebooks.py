"""Optional case-study checks (install the notebooks extra to enable)."""

import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("rdkit")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import py3Dmol

from eufod.io import parse_mapping_text
from eufod.search import _predicted_shifts, _score_chunk

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("name", ["borneol", "isoborneol"])
@pytest.mark.parametrize("generate", [False, True])
def test_case_study(name, generate, tmp_path, monkeypatch):
    notebook = json.loads((ROOT / f"notebooks/{name}_titration_to_input.ipynb").read_text(encoding="utf-8"))
    monkeypatch.setattr(py3Dmol.view, "show", lambda self: None)
    monkeypatch.chdir(tmp_path if generate else ROOT)
    if generate:
        (tmp_path / "outputs").mkdir()
    namespace = {"display": lambda *args: None}
    try:
        for index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] == "code" and (not generate or index < 17):
                exec(compile("".join(cell["source"]), f"{name}:cell{index}", "exec"), namespace)
        molecule = namespace[name]
        assert molecule.GetAtomWithIdx(namespace["oxygen_atom_index"]).GetSymbol() == "O"
        for atom_index in namespace["proton_atom_indices"].values():
            assert molecule.GetAtomWithIdx(atom_index).GetSymbol() == "H"
        if generate:
            assert (tmp_path / f"outputs/optimized_geometries/{name}_mmff94_optimized.sdf").exists()
            return
        _, coordinates, scores = parse_mapping_text(
            (ROOT / f"outputs/{name}_r_factor_mapping.txt").read_text(encoding="utf-8")
        )
        by_atom = namespace["coordinates_by_atom"]
        oxygen = np.asarray(by_atom[namespace["oxygen_atom_index"]])
        shifts = namespace["model_input_shifts"]
        hydrogens = np.asarray([
            by_atom[namespace["proton_atom_indices"][label]] for label in shifts["Label"]
        ])
        experimental = shifts["Relative exp. shift"].round(6).to_numpy()
        for start in range(0, len(scores), 25_000):
            predicted, _, _ = _predicted_shifts(coordinates[start:start + 25_000], oxygen, hydrogens)
            calculated = _score_chunk(predicted, experimental, "r_factor")
            np.testing.assert_allclose(calculated, scores[start:start + 25_000], rtol=1e-12, atol=5.1e-9)
    finally:
        plt.close("all")
