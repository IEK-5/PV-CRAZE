from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).resolve().parent.parent / "data"


def load(name):
    return np.loadtxt(DATA / name, comments="#")


@pytest.fixture(scope="session")
def measured_curves():
    """The 1000 measured IV curves from data/IVs.csv as (n, 2) arrays."""
    curves = []
    for line in (DATA / "IVs.csv").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        values = np.array(line.rstrip(",").split(","), dtype=float)
        n = int(values[0])
        curves.append(np.column_stack([values[1 : 1 + n], values[1 + n : 1 + 2 * n]]))
    return curves


@pytest.fixture(scope="session")
def gt_data():
    """(T [degC], G [W/m2], ESP (N, 10)) of a diode model on a G-T grid."""
    d = load("refdiodeGT.dat")
    return d[:, 1], d[:, 0], d[:, 2:]


@pytest.fixture
def diode_curve():
    """A clean synthetic single-diode curve (V, I) inside the first quadrant."""
    V = np.linspace(0.0, 0.6, 200)
    I = 0.035 - 1e-10 * (np.exp(V / 0.03) - 1) - V / 700
    return V, I


def rel_dev(a, b):
    """Largest deviation per row, relative to the column means of b."""
    return np.max(np.abs((np.asarray(a) - np.asarray(b)) / np.mean(np.asarray(b), axis=0)), axis=1)
