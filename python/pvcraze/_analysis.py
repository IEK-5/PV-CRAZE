"""IV curve analysis and model fitting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np

from . import _core
from ._arrays import normalize_curves
from ._esp import ESP, _import_pandas, as_esp
from ._flags import Flag, error_mask, warning_mask
from ._models import FitSpec, Model, ModelLike, ModelParams, _params_from_core, coerce_specs

_MPP_METHODS = {"fast": "f", "robust": "r", "parametric": "p", "savgol": "s"}
_DEFAULTS = _core.default_analysis_settings()


@dataclass(frozen=True)
class FitResult:
    """A fitted curve model.

    Attributes
    ----------
    params:
        The fitted :class:`ModelParams` (single or batch, like the input).
    error:
        Fit error per curve (float or array), dimensionless: the RMS distance
        in the (V/Voc, I/Isc) plane. From :func:`analyze` it is measured
        against the analysed part of the IV curve; from :func:`fit_model`,
        where no curve is available, against the five ESP key points.
    spec:
        The :class:`FitSpec` used.
    """

    params: ModelParams
    error: Union[float, np.ndarray]
    spec: FitSpec

    @property
    def model(self) -> Model:
        """The fitted model."""
        return self.params.model

    def __repr__(self) -> str:
        if isinstance(self.error, float):
            return f"FitResult({self.params!r}, error={self.error:.3g})"
        return f"<FitResult {self.model.value} batch n={len(self.error)}>"


@dataclass(frozen=True)
class AnalysisResult:
    """Result of :func:`analyze`.

    Single-curve input gives single results (floats, a single :class:`ESP`),
    batch input gives batch results (arrays, a batch :class:`ESP`).

    Attributes
    ----------
    esp:
        The extracted :class:`ESP`.
    flags:
        Raw flag codes, an ``int`` (single) or an int64 array (batch). Use
        :meth:`flag` for a :class:`Flag` with messages, or :attr:`ok`.
    fits:
        ``{model name: FitResult}`` for the requested models, in request
        order. Keys are the library names, e.g. ``"PV-Bezier"``; a
        :class:`Model` member works as a key too.
    """

    esp: ESP
    flags: Union[int, np.ndarray]
    fits: Dict[str, FitResult]

    @property
    def is_batch(self) -> bool:
        return self.esp.is_batch

    @property
    def ok(self) -> Union[bool, np.ndarray]:
        """True where no error flag is set (the ESP are meaningful)."""
        m = ~error_mask(self.flags)
        return bool(m) if not self.is_batch else m

    @property
    def has_warning(self) -> Union[bool, np.ndarray]:
        """True where a warning flag is set."""
        m = warning_mask(self.flags)
        return bool(m) if not self.is_batch else m

    def flag(self, index: Optional[int] = None) -> Flag:
        """The :class:`Flag` of a curve (``index`` is required for a batch)."""
        if not self.is_batch:
            if index not in (None, 0):
                raise IndexError("single result: index must be None or 0")
            return Flag(int(self.flags))
        if index is None:
            raise TypeError("batch result: give the curve index, e.g. result.flag(3)")
        return Flag(int(self.flags[index]))

    def __len__(self) -> int:
        return len(self.esp)

    def to_dataframe(self) -> Any:
        """One row per curve: ESP columns, ``flags``, ``ok`` and ``error_<model>``. Requires pandas."""
        _import_pandas()  # fail early with an install hint
        df = self.esp.to_dataframe()
        df["flags"] = np.atleast_1d(self.flags)
        df["ok"] = np.atleast_1d(self.ok)
        for name, fit in self.fits.items():
            df[f"error_{name}"] = np.atleast_1d(fit.error)
        return df

    def __repr__(self) -> str:
        models = ", ".join(self.fits) or "none"
        if not self.is_batch:
            return f"AnalysisResult({self.esp!r}, flags={self.flag()!r}, fits: {models})"
        n_ok = int(np.sum(self.ok))
        return f"<AnalysisResult n={len(self)} ok={n_ok} fits: {models}>"


class _FitDict(dict):
    """dict keyed by model name that also accepts Model members and loose names."""

    def __getitem__(self, key: Any) -> FitResult:
        if not dict.__contains__(self, key):
            try:
                key = Model.lookup(key).value
            except (TypeError, ValueError):
                pass
        return dict.__getitem__(self, key)

    def __contains__(self, key: object) -> bool:
        try:
            self[key]
        except KeyError:
            return False
        return True


def _fits_from_core(specs: List[FitSpec], core_fits: Any, single: bool) -> Dict[str, FitResult]:
    types, pars, err = core_fits
    out = _FitDict()
    for k, spec in enumerate(specs):
        params = _params_from_core(types[k], pars[k], single)
        e = err[:, k]
        out[spec.model.value] = FitResult(params, float(e[0]) if single else e.copy(), spec)
    return out


def _fraction(name: str, value: Optional[float]) -> float:
    if value is None:
        return float(_DEFAULTS[name])
    if not 0 < value < 1:
        raise ValueError(f"{name} must be between 0 and 1, got {value!r}")
    return float(value)


def analyze(
    curves: Any,
    I: Any = None,  # noqa: E741 - physics notation
    /,
    *,
    models: Union[ModelLike, Sequence[ModelLike], None] = None,
    mpp_method: str = "savgol",
    spike_filter: bool = True,
    identical_eps: Optional[float] = None,
    resample: Optional[int] = None,
    mpp_dv: Optional[float] = None,
    mpp_di: Optional[float] = None,
    sc_dv: Optional[float] = None,
    oc_di: Optional[float] = None,
    min_rel_range: Optional[float] = None,
) -> AnalysisResult:
    """Extract the extended solar cell parameters (ESP) from measured IV curves.

    Call as ``analyze(curves)`` or ``analyze(V, I)``:

    * ``analyze(iv)`` with an ``(n, 2)`` array with columns (V, I): one curve.
    * ``analyze([iv0, iv1, ...])``: a batch; curves may differ in length.
    * ``analyze(V, I)`` with 1-D arrays: one curve.
    * ``analyze([V0, V1, ...], [I0, I1, ...])``: a batch.

    Curves use the generator convention: the current is positive between
    short circuit (V=0) and open circuit (I=0). Point order does not matter.

    Parameters
    ----------
    curves, I:
        The IV data, see above. Voltages in V, currents in A.
    models:
        Curve models to fit to each curve's ESP: a name, a :class:`Model`, a
        :class:`FitSpec`, or a sequence of those. See :class:`Model`.
    mpp_method:
        How the maximum power point is determined from 3rd order polynomial
        fits: ``"savgol"`` (Savitzky-Golay filter; robust and fast, default),
        ``"robust"`` (robust least squares), ``"parametric"`` (V, I and P as
        polynomials of the point index), or ``"fast"`` (assumes no noise).
    spike_filter:
        Replace the current of every interior point by the median of itself
        and its two neighbours (voltage order) before the analysis. This
        removes isolated spikes and is on by default. Independent of this
        setting, large isolated outliers are reported as ``Flag.OUTLIERS``,
        and three or more of them mark the curve as ``Flag.NOISY_IV``.
    identical_eps:
        If given, drop consecutive points whose V and I differ by less than
        this value.
    resample:
        If given, resample the curve to this many points, distributed evenly
        along the curve (piecewise linear interpolation). This can rescue
        curves measured with too coarse a voltage step near open circuit.
    mpp_dv, mpp_di:
        Fraction of Vmpp / Impp that sets the data range for the MPP fit: from
        ``(1 - mpp_dv) * Vmpp`` to ``(1 - mpp_di) * Impp`` (default 0.12 each).
    sc_dv:
        Voltage range for the short-circuit fit, as a fraction of Vqmp-
        (default 0.1).
    oc_di:
        Current range for the open-circuit fit, as a fraction of Iqmp+
        (default 0.1).
    min_rel_range:
        Minimum relative range ``1 - min/max`` the voltage and current data
        must span for any analysis to be attempted (default 0.5).

    Returns
    -------
    AnalysisResult
        With ``.esp``, ``.flags`` / ``.ok`` and ``.fits[model]``.

    Examples
    --------
    >>> res = pvcraze.analyze(V, I, models="PV-Bezier")
    >>> res.ok, res.esp.voc, res.esp.ff
    >>> res.fits["PV-Bezier"].params.iv(0.4)
    """
    curve_list, single = normalize_curves(curves, I)
    specs = coerce_specs(models)
    method = _MPP_METHODS.get(mpp_method, mpp_method if mpp_method in _MPP_METHODS.values() else None)
    if method is None:
        raise ValueError(f"mpp_method must be one of {', '.join(_MPP_METHODS)}, got {mpp_method!r}")
    if identical_eps is not None and not identical_eps > 0:
        raise ValueError(f"identical_eps must be positive, got {identical_eps!r}")
    if resample is not None and (int(resample) != resample or resample < 2):
        raise ValueError(f"resample must be an integer >= 2, got {resample!r}")

    esp, flags, core_fits = _core.analyze(
        curve_list,
        [s._spec_string() for s in specs],
        spike_filter=bool(spike_filter),
        identical_eps=float(identical_eps or 0.0),
        resample=int(resample or 0),
        mpp_method=method,
        mpp_dv=_fraction("mpp_dv", mpp_dv),
        mpp_di=_fraction("mpp_di", mpp_di),
        sc_dv=_fraction("sc_dv", sc_dv),
        oc_di=_fraction("oc_di", oc_di),
        min_rel_range=_fraction("min_rel_range", min_rel_range),
    )
    fits = _fits_from_core(specs, core_fits, single)
    if single:
        return AnalysisResult(ESP._wrap(esp[0]), int(flags[0]), fits)
    return AnalysisResult(ESP._wrap(esp), flags, fits)


def fit_model(esp: Any, model: ModelLike = "PV-Bezier", **options: Any) -> FitResult:
    """Fit a curve model to extended solar cell parameters.

    Parameters
    ----------
    esp:
        An :class:`ESP` (single or batch), or an array accepted by
        :meth:`ESP.from_array`.
    model:
        Model name, :class:`Model` or :class:`FitSpec` (default ``"PV-Bezier"``).
    **options:
        Fit options, see :class:`FitSpec` (e.g. ``population=200``). Not
        allowed together with a ``FitSpec``.

    Returns
    -------
    FitResult
        ``.params`` has the same single/batch form as ``esp``. ``.error`` is
        measured against the five ESP key points.

    Examples
    --------
    >>> fit = pvcraze.fit_model(esp, "Diode", population=200, max_iter=4000)
    >>> fit.params.rs, fit.params.rsh
    """
    e = as_esp(esp)
    if isinstance(model, FitSpec):
        if options:
            raise TypeError("pass fit options either in the FitSpec or as keywords, not both")
        spec = model
    else:
        spec = FitSpec(model, **options)
    single = not e.is_batch
    core_fits = _core.fit(np.atleast_2d(np.asarray(e)), [spec._spec_string()])
    return _fits_from_core([spec], core_fits, single)[spec.model.value]
