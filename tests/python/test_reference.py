"""Regression tests against the reference data in tests/data.

These port the former tests/test_0*.py scripts to the new API. The Octave
tests (tests/test_0*.m) use the same reference files.
"""

import numpy as np
import pytest

import pvcraze as pc
from conftest import load, rel_dev

TH = 1e-2
DIODE_SPEC = pc.FitSpec("Diode", refine_tol=0.6e-3, population=200, max_iter=4000)
MODELS = {
    "PV-Bezier": pc.FitSpec("PV-Bezier"),
    "PV-spline": pc.FitSpec("PV-spline"),
    "KH": pc.FitSpec("KH"),
    "Diode": DIODE_SPEC,
}


def test_analyze_reference_esp_and_flags(measured_curves):
    res = pc.analyze(measured_curves)
    assert res.is_batch and len(res) == len(measured_curves)
    ref = load("refesp1.dat")
    assert np.sum(rel_dev(np.asarray(res.esp), ref[:, :10]) > TH) == 0
    ref_flags = load("refflags1.dat")[:, 0]
    np.testing.assert_array_equal(res.flags, ref_flags.astype(np.int64))


def test_analyze_flags_bad_curves():
    t = np.linspace(1, 100, num=100)
    curves = [
        np.ones((50, 2)),
        np.column_stack([t, np.cos(t * np.pi / 50)]),
        np.column_stack([t, np.sin(t * np.pi / 50)]),
        np.column_stack([t, 1 - np.exp(t * np.pi / 50)]),
        np.column_stack([t, np.exp(t * np.pi / 50)]),
    ]
    res = pc.analyze(curves)
    np.testing.assert_array_equal(res.flags, [256, 0, 2097264, 136356976, 136356976])
    assert list(res.ok) == [False, True, False, False, False]
    assert res.flag(0) == pc.Flag.INSUFF_I_RANGE  # a constant current has no range
    assert res.flag(0).messages == ["Relative Current range too small"]


def test_diode_gt_translation_reference(gt_data):
    T, G, _ = gt_data
    ref = load("refdiodeGT.dat")[:, 2:]
    I0, nkT, Iph, Rs, Rsh = 1e-8, 1.5 * 8.6173899875662609e-05 * 298.15, 0.035, 0.1, 700.0
    gt = pc.diode_gt.DiodeGT(
        pc.ModelParams("Diode", i0=I0, nkt=nkT, iph=Iph, rs=Rs, rsh=Rsh), alpha=1e-3, ea=0.8
    )
    esp = gt.esp(T=T, G=G)
    assert np.sum(rel_dev(np.asarray(esp), ref) > TH) == 0


@pytest.mark.parametrize("name", list(MODELS))
def test_fit_model_reference(gt_data, name):
    _, _, esp = gt_data
    fit = pc.fit_model(esp, MODELS[name])
    model_esp = fit.params.esp()
    ref = load(f"ref{name}.dat")
    # relative to the mean of the model ESP, as in the original test
    dev = np.max(np.abs((ref - np.asarray(model_esp)) / np.mean(np.asarray(model_esp), axis=0)), axis=1)
    assert np.sum(dev > TH) == 0


@pytest.mark.parametrize("name", list(MODELS))
def test_model_rmse_against_bezier_curves(gt_data, name):
    _, _, esp = gt_data
    curves = pc.fit_model(esp, "PV-Bezier").params.curve(100)
    fit = pc.fit_model(esp, MODELS[name])
    assert np.sum(fit.params.rmse(curves) > TH) == 0


def _dist(esp1, esp2):
    """Mean normalised distance of the five key points of two ESP batches."""
    a, b = np.asarray(esp1), np.asarray(esp2)
    voc0 = (a[:, 0] + b[:, 0]) / 2
    isc0 = (a[:, 1] + b[:, 1]) / 2
    e = np.zeros(len(a))
    e += ((a[:, 0] - b[:, 0]) / voc0) ** 2 + ((a[:, 1] - b[:, 1]) / isc0) ** 2
    for v, i in ((2, 3), (4, 5), (6, 7)):
        e += ((a[:, v] - b[:, v]) / voc0) ** 2 + ((a[:, i] - b[:, i]) / isc0) ** 2
    return np.mean(np.sqrt(e))


def test_diode_gt_fit_roundtrip(gt_data):
    T, G, esp = gt_data
    gt = pc.diode_gt.fit(esp, T, G)
    d = _dist(esp, gt.esp(T, G))
    assert np.isfinite(d) and d < 2e-2


@pytest.mark.parametrize("version", [2, 3])
def test_iec60891_fit_roundtrip(gt_data, version):
    T, G, esp = gt_data
    coef = pc.iec60891.fit(esp, T, G, version=version)
    assert coef.esp_ref is not None and not coef.esp_ref.is_batch
    d = _dist(esp, coef.translate(T, G))
    assert np.isfinite(d) and d < 15e-2


def test_model_table_is_stable():
    names = ["APWL", "B2", "Diode", "KH", "PC", "PV-Bezier", "PV-spline", "PWL", "Spline"]
    assert [m.value for m in sorted(pc.Model, key=lambda m: m.id)] == names
    assert [m.id for m in sorted(pc.Model, key=lambda m: m.id)] == list(range(len(names)))


@pytest.mark.parametrize("name", list(MODELS))
def test_iv_vi_consistent_with_curves(gt_data, name):
    _, _, esp = gt_data
    params = pc.fit_model(esp, MODELS[name]).params
    ei = ev = 0.0
    for p, c in zip(params, params.curve(100)):
        V, I = c[:, 0], c[:, 1]
        ei += np.sqrt(np.mean((I - p.iv(V)) ** 2))
        ev += np.sqrt(np.mean((V - p.vi(I)) ** 2))
    assert ei / len(params) < 1e-6
    assert ev / len(params) < 1e-6
