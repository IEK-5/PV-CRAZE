"""Behaviour of the Python API: input forms, errors, single vs batch."""

import inspect

import numpy as np
import pytest

import pvcraze as pc


# ---------------------------------------------------------------- analyze
def test_analyze_input_forms_agree(diode_curve):
    V, I = diode_curve
    iv = np.column_stack([V, I])
    a = pc.analyze(V, I)
    b = pc.analyze(iv)
    c = pc.analyze(iv.tolist())
    d = pc.analyze([iv, iv[::2]])
    e = pc.analyze([V, V[::2]], [I, I[::2]])
    f = pc.analyze((iv.astype(np.float32),))
    assert not a.is_batch and not b.is_batch and not c.is_batch
    assert a.esp == b.esp == c.esp
    assert d.is_batch and len(d) == 2
    assert d.esp[0] == a.esp
    np.testing.assert_array_equal(np.asarray(d.esp), np.asarray(e.esp))
    assert f.is_batch and abs(f.esp[0].voc - a.esp.voc) < 1e-4


def test_analyze_single_result_is_scalar(diode_curve):
    res = pc.analyze(*diode_curve, models="PV-Bezier")
    assert isinstance(res.esp.voc, float)
    assert isinstance(res.flags, int)
    assert res.ok is True
    assert res.flag() == pc.Flag(0) and res.flag().messages == []
    fit = res.fits["PV-Bezier"]
    assert isinstance(fit.error, float) and fit.error < 1e-2
    assert not fit.params.is_batch
    assert 0.55 < res.esp.voc < 0.62 and 0.7 < res.esp.ff < 0.85


def test_analyze_int_curves_are_converted():
    V = np.arange(0, 60)
    I = (100 - np.exp(V / 12.0)).astype(np.int64)
    res = pc.analyze(V, I)
    assert res.esp.voc > 50


def test_analyze_empty_batch():
    res = pc.analyze([], models="PV-Bezier")
    assert res.is_batch and len(res) == 0
    assert res.fits["PV-Bezier"].params.to_array().shape == (0, 10)


def test_fits_lookup_is_forgiving(diode_curve):
    res = pc.analyze(*diode_curve, models=[pc.Model.PV_BEZIER, "diode"])
    assert list(res.fits) == ["PV-Bezier", "Diode"]
    assert res.fits["pv_bezier"] is res.fits[pc.Model.PV_BEZIER]
    assert "DIODE" in res.fits and "KH" not in res.fits


@pytest.mark.parametrize(
    "kwargs, exc",
    [
        ({"models": "Bezier"}, ValueError),
        ({"models": [3]}, TypeError),
        ({"models": ["PV-Bezier", "pv-bezier"]}, ValueError),
        ({"mpp_method": "magic"}, ValueError),
        ({"mpp_dv": 1.5}, ValueError),
        ({"resample": 1}, ValueError),
        ({"identical_eps": -1.0}, ValueError),
    ],
)
def test_analyze_rejects_bad_options(diode_curve, kwargs, exc):
    with pytest.raises(exc):
        pc.analyze(*diode_curve, **kwargs)


@pytest.mark.parametrize(
    "curves",
    [np.ones((5, 3)), np.ones(5), "iv.csv", [np.ones((5, 2)), np.ones((5, 3))], np.ones((4, 2), dtype=complex)],
)
def test_analyze_rejects_bad_curves(curves):
    with pytest.raises((TypeError, ValueError)):
        pc.analyze(curves)


def test_analyze_hints_at_v_i_form(diode_curve):
    with pytest.raises(ValueError, match=r"analyze\(V, I\)"):
        pc.analyze(list(diode_curve))


def test_unknown_model_error_lists_models():
    with pytest.raises(ValueError, match="PV-Bezier"):
        pc.Model.lookup("bezier")


def test_mpp_methods(diode_curve):
    vocs = [pc.analyze(*diode_curve, mpp_method=m).esp.vmpp for m in ("savgol", "robust", "parametric", "fast")]
    assert np.ptp(vocs) < 5e-3


def test_signatures_are_introspectable():
    sig = inspect.signature(pc.analyze)
    assert "models" in sig.parameters and "mpp_method" in sig.parameters
    assert pc.analyze.__doc__ and "Parameters" in pc.analyze.__doc__


# ---------------------------------------------------------------- ESP
def test_esp_construction_and_fields():
    e = pc.ESP(voc=0.6, isc=0.035, vmpp=0.5, impp=0.032, vqmp_minus=0.42, iqmp_minus=0.034,
               vqmp_plus=0.56, iqmp_plus=0.022, gsc=1.7e-3, roc=0.96)
    assert not e.is_batch and e.voc == 0.6
    assert e.pmpp == pytest.approx(0.016) and e.ff == pytest.approx(0.016 / 0.021)
    assert pc.ESP.from_array(np.asarray(e)) == e
    assert pc.ESP.from_dict(e.to_dict()) == e
    with pytest.raises(AttributeError):
        e.voc = 1.0
    with pytest.raises(TypeError):
        len(e)
    b = pc.ESP.concat([e, e, e])
    assert b.is_batch and len(b) == 3 and b[1] == e and len(b[b.voc > 0]) == 3
    assert [x == e for x in b] == [True] * 3
    arr = np.asarray(b)
    assert arr.shape == (3, 10) and not arr.flags.writeable
    assert b.to_array().flags.writeable


def test_esp_broadcasts_fields():
    e = pc.ESP(voc=[0.6, 0.7], isc=0.035, vmpp=0.5, impp=0.032, vqmp_minus=0.42,
               iqmp_minus=0.034, vqmp_plus=0.56, iqmp_plus=0.022, gsc=1.7e-3, roc=0.96)
    assert e.is_batch and list(e.voc) == [0.6, 0.7] and list(e.isc) == [0.035, 0.035]


def test_esp_from_array_validates():
    with pytest.raises(ValueError):
        pc.ESP.from_array(np.ones((3, 9)))


def test_dataframes(diode_curve):
    pd = pytest.importorskip("pandas")
    res = pc.analyze([np.column_stack(diode_curve)] * 3, models="PV-Bezier")
    df = res.to_dataframe()
    assert list(df.columns[:10]) == list(pc.ESP.FIELDS)
    assert {"flags", "ok", "error_PV-Bezier"} <= set(df.columns) and len(df) == 3
    assert pc.ESP.from_dataframe(df) == res.esp
    assert list(res.fits["PV-Bezier"].params.to_dataframe().columns) == list(pc.ESP.FIELDS)
    assert isinstance(df, pd.DataFrame)


# ---------------------------------------------------------------- models
def test_model_lookup_and_params():
    assert pc.Model("pv bezier") is pc.Model.PV_BEZIER
    assert pc.Model.DIODE == "Diode" and str(pc.Model.DIODE) == "Diode"
    assert pc.Model.DIODE.param_names == ("i0", "nkt", "iph", "rs", "rsh")
    assert pc.Model.PV_BEZIER.param_names == pc.ESP.FIELDS
    assert pc.Model.PV_BEZIER.param_labels[4:8] == ("Vqmp-", "Iqmp-", "Vqmp+", "Iqmp+")


def test_model_params_named_and_evaluation():
    p = pc.ModelParams("Diode", i0=1e-10, nkt=0.03, iph=0.035, rs=0.1, rsh=700)
    assert p.rs == 0.1 and "rsh" in dir(p)
    assert isinstance(p.iv(0.3), float)
    assert p.iv(np.zeros((2, 3))).shape == (2, 3)
    assert p.iv(0.0) == pytest.approx(0.035, rel=1e-2)
    assert p.vi(p.iv(0.45)) == pytest.approx(0.45, abs=1e-6)
    c = p.curve(50)
    assert c.ndim == 2 and c.shape[1] == 2
    assert p.esp().voc == pytest.approx(p.vi(0.0), rel=1e-6)
    assert p.rmse(c) < 1e-6
    measured = c + [0.0, 1e-4]
    d = p.distance(measured)
    assert d.shape == c.shape and np.all(d[:, 1] >= 0)
    assert np.sqrt(np.mean(d[:, 1] ** 2)) == pytest.approx(p.rmse(measured), rel=1e-6)
    with pytest.raises(AttributeError, match="available"):
        p.voc
    with pytest.raises(TypeError, match="missing"):
        pc.ModelParams("Diode", i0=1e-10)
    with pytest.raises(ValueError, match="5 parameters"):
        pc.ModelParams("Diode", np.ones(10))


def test_model_params_batch_broadcasting():
    p = pc.ModelParams("Diode", i0=[1e-10, 2e-10], nkt=0.03, iph=0.035, rs=0.1, rsh=700)
    assert p.is_batch and len(p) == 2
    assert p.iv(0.3).shape == (2,)
    assert p.iv(np.linspace(0, 0.5, 7)).shape == (2, 7)
    assert p.iv(np.array([[0.1], [0.2]])).shape == (2, 1)
    with pytest.raises(ValueError):
        p.iv(np.zeros((3, 4)))
    curves = p.curve(30)
    assert isinstance(curves, list) and len(curves) == 2
    assert p.rmse(curves).shape == (2,)
    with pytest.raises(ValueError):
        p.rmse(curves[0])


def test_fit_model_single_and_batch(gt_data):
    _, _, esp = gt_data
    batch = pc.fit_model(esp, "KH")
    assert batch.params.is_batch and batch.error.shape == (len(esp),)
    single = pc.ESP.from_array(esp[3]).fit("KH")
    assert not single.params.is_batch and isinstance(single.error, float)
    assert single.params.to_array() == pytest.approx(batch.params[3].to_array())


@pytest.mark.parametrize(
    "kwargs",
    [{"population": 5}, {"max_iter": 1.5}, {"tol": 1.0}, {"refine": True}],
)
def test_fitspec_validation(kwargs):
    with pytest.raises(ValueError):
        pc.FitSpec("Diode", **kwargs)


def test_fitspec_spec_string():
    s = pc.FitSpec("diode", population=200, max_iter=4000, refine_tol=6e-4)._spec_string()
    assert s == "model(Diode) fitit(4000) npop(200) oeps(0.0006)"


# ---------------------------------------------------------------- flags
def test_flag_helpers():
    f = pc.Flag.MULTIPLEMPP | pc.Flag.NOMPP
    assert f.has_error and f.has_warning
    assert f.errors == pc.Flag.NOMPP and f.warnings == pc.Flag.MULTIPLEMPP
    assert f.messages == ["No power maximum", "more than one power maximum"]
    assert pc.describe_flags(int(f)) == f.messages
    assert pc.ERROR_MASK & pc.WARNING_MASK == 0


# ---------------------------------------------------------------- translation
def test_iec60891_v3_datasheet_coefficients(gt_data):
    T, G, esp = gt_data
    fitted = pc.iec60891.IEC60891v3.fit(esp, T, G)
    manual = pc.iec60891.IEC60891v3(
        alpha=fitted.alpha, beta=fitted.beta, b1=fitted.b1, b2=fitted.b2,
        rs=fitted.rs, kappa=fitted.kappa, voc_ref=fitted.esp_ref.voc,
    )
    a = fitted.translate(T=60.0, G=800.0)
    b = manual.translate(T=60.0, G=800.0, esp=fitted.esp_ref)
    assert not a.is_batch and a == b
    grid = fitted.translate(T=[25.0, 50.0], G=1000.0)
    assert grid.is_batch and len(grid) == 2
    assert grid[0].voc == pytest.approx(fitted.esp_ref.voc)
    with pytest.raises(ValueError):
        manual.translate(T=25.0, G=1000.0)
    with pytest.raises(TypeError):
        pc.iec60891.IEC60891v3(alpha=0, beta=0, b1=0, b2=0, rs=0, kappa=0)


def test_iec60891_translate_between_conditions(gt_data):
    T, G, esp = gt_data
    coef = pc.iec60891.fit(esp, T, G, version=2)
    stc = coef.translate(T=25.0, G=1000.0, esp=esp, T_from=T, G_from=G)
    assert stc.is_batch and len(stc) == len(esp)
    assert np.std(stc.voc) / np.mean(stc.voc) < 0.05
    with pytest.raises(ValueError):
        pc.iec60891.fit(esp, T, G, version=4)
    with pytest.raises(ValueError):
        coef.translate(T=[1.0, 2.0, 3.0], G=1000.0, esp=esp[:2], T_from=25.0, G_from=1000.0)


def test_diode_gt_at(gt_data):
    T, G, esp = gt_data
    gt = pc.diode_gt.fit(esp, T, G, max_points=10)
    p = gt.at(25.0, 1000.0)
    assert p.model is pc.Model.DIODE and not p.is_batch
    assert p.to_array() == pytest.approx(gt.params.to_array(), rel=1e-9)
    assert gt.at(T, G).is_batch
    with pytest.raises(ValueError):
        gt.at(-300.0, 1000.0)
    with pytest.raises(TypeError):
        pc.diode_gt.DiodeGT(pc.ModelParams("KH", np.ones(4)), 0.0, 0.0)


# ---------------------------------------------------------------- verbosity
def test_quiet_by_default_and_verbose_context(capfd, diode_curve):
    assert not any(pc.get_verbosity().values())
    pc.analyze(*diode_curve)
    assert capfd.readouterr().out == ""
    with pc.verbose(progress=False):
        assert pc.get_verbosity()["report"]
        pc.analyze(*diode_curve)
    assert "PV-CRAZE" in capfd.readouterr().out
    assert not any(pc.get_verbosity().values())
    prev = pc.set_verbosity(timing=True)
    assert prev == {"timing": False, "report": False, "progress": False, "header": False}
    pc.set_verbosity()
