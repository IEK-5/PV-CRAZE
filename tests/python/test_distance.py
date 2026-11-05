"""The area distance between IV curves (pvcraze.curve_distance, ModelParams.error)."""

import numpy as np
import pytest

import pvcraze as pc

# reference: the line I = 1 - V, Voc = Isc = 1, so the distance is the plain area
LINE = np.column_stack([np.linspace(0, 1, 11), 1 - np.linspace(0, 1, 11)])


def test_identical_curves_have_no_distance():
    assert pc.curve_distance(LINE, LINE) == 0.0


def test_analytic_areas():
    # a rectangle of current 1 up to V = 1: the area between is the upper triangle
    rect = np.array([[0.0, 1.0], [1.0, 1.0], [1.0 + 1e-12, 0.0]])
    assert pc.curve_distance(LINE, rect) == pytest.approx(0.5, abs=1e-9)
    # a constant 0.5 crosses the line at V = 0.5: two triangles of 1/8
    half = np.array([[0.0, 0.5], [1.0, 0.5], [1.0 + 1e-12, 0.0]])
    assert pc.curve_distance(LINE, half) == pytest.approx(0.25, abs=1e-9)
    # a curve without active range: the full area under the reference
    assert pc.curve_distance(LINE, np.array([[0.0, 0.0], [1.0, 0.0]])) == pytest.approx(0.5)


def test_normalized_by_the_reference():
    scaled = LINE * [2.0, 3.0]  # Voc = 2, Isc = 3
    other = scaled * [1.0, 0.9]
    # area 0.1*3*2/2 = 0.3, divided by Voc*Isc = 6
    assert pc.curve_distance(scaled, other) == pytest.approx(0.05)
    # with the roles swapped the normalization is 0.9*6
    assert pc.curve_distance(other, scaled) == pytest.approx(0.05 / 0.9)


def test_only_the_active_range_counts():
    outside = np.vstack([[[-0.5, 1.5]], LINE, [[1.5, -0.5]]])
    assert pc.curve_distance(LINE, outside) == pytest.approx(0.0, abs=1e-12)
    # a curve that stops short of short and open circuit is extended linearly
    inner = LINE[2:-2]
    assert pc.curve_distance(LINE, inner) == pytest.approx(0.0, abs=1e-12)
    # point order does not matter
    rng = np.random.default_rng(1)
    assert pc.curve_distance(LINE, rng.permutation(LINE)) == pytest.approx(0.0, abs=1e-12)


def test_invalid_curves_give_nan():
    assert np.isnan(pc.curve_distance(LINE, LINE[:1]))
    assert np.isnan(pc.curve_distance(-LINE, LINE))  # the reference has no active range


@pytest.fixture
def diode():
    return pc.ModelParams("Diode", i0=1e-10, nkt=0.03, iph=0.035, rs=0.1, rsh=700)


def test_insensitive_to_the_sampling(diode):
    """Extra points in one region dilute a point-wise RMS distance, not the area."""
    # a measurement that deviates from the model beyond V = 0.4 only
    a = diode.curve(200)
    a[a[:, 0] > 0.4, 1] -= 3e-4
    # the same measurement with 2000 additional points near short circuit
    V = np.linspace(0.0, 0.1, 2000)
    b = np.vstack([a, np.column_stack([V, diode.iv(V)])])
    assert diode.error(a) == pytest.approx(diode.error(b), rel=1e-2)
    rms = [np.sqrt(np.mean(diode.distance(c)[:, 1] ** 2)) for c in (a, b)]
    assert rms[1] < 0.5 * rms[0]


def test_models_and_curves_agree(diode):
    other = pc.ModelParams("Diode", i0=1e-10, nkt=0.03, iph=0.0355, rs=0.1, rsh=700)
    d_mm = pc.curve_distance(diode, other)
    assert d_mm == pytest.approx(0.0005 / 0.035, rel=0.05)  # about dI / Isc
    assert pc.curve_distance(diode, other.curve(2000)) == pytest.approx(d_mm, rel=1e-3)
    assert pc.curve_distance(diode.curve(2000), other) == pytest.approx(d_mm, rel=1e-3)
    assert pc.curve_distance(diode.curve(2000), other.curve(2000)) == pytest.approx(d_mm, rel=1e-3)
    assert diode.error(other.curve(2000)) == pytest.approx(d_mm, rel=1e-3)


def test_batches_and_broadcasting(diode):
    batch = pc.ModelParams("Diode", i0=1e-10, nkt=0.03, iph=[0.035, 0.0355, 0.036], rs=0.1, rsh=700)
    d = pc.curve_distance(diode, batch)
    assert d.shape == (3,) and d[0] == pytest.approx(0.0, abs=1e-12) and d[1] < d[2]
    curves = batch.curve(500)
    np.testing.assert_allclose(pc.curve_distance(diode, curves), d, rtol=2e-2, atol=1e-5)
    np.testing.assert_allclose(pc.curve_distance(batch, curves), batch.error(curves))
    assert pc.curve_distance(curves, curves).shape == (3,)
    with pytest.raises(ValueError, match="3 reference curves but 2"):
        pc.curve_distance(curves, curves[:2])
