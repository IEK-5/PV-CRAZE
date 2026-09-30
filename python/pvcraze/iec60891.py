"""Irradiance and temperature translation of ESP following IEC 60891, procedure 2.

Two editions of the norm are implemented:

* :class:`IEC60891v2` -- IEC 60891 edition 2 (2009)
* :class:`IEC60891v3` -- IEC 60891 edition 3 (2021)

Each class holds the translation coefficients. Get them by fitting measured
ESP at several irradiances and temperatures, or construct the object directly
from known (e.g. datasheet) coefficients::

    from pvcraze import iec60891
    coef = iec60891.IEC60891v3.fit(esp, T=T, G=G)       # or iec60891.fit(..., version=3)
    esp_stc = coef.esp_ref                              # ESP at 25 degC, 1000 W/m2
    esp_hot = coef.translate(T=60, G=800)               # from reference conditions
    esp_2 = coef.translate(T=25, G=1000, esp=esp, T_from=T, G_from=G)   # any -> any

Temperatures are in degrees Celsius, irradiances in W/m2. All condition
arguments broadcast against each other and against the rows of ``esp``.

The ESP are translated through an interpolated IV curve, selected with the
``model`` argument: :attr:`Model.PV_SPLINE` (default) or
:attr:`Model.PV_BEZIER`. The IEC 60891 map of the (V, I) points of a curve is
affine, and a Bezier curve is invariant under affine maps, so the PV-Bezier
curve is translated exactly. With the PV-spline only the knots and slopes are
translated::

    coef = iec60891.IEC60891v3.fit(esp, T=T, G=G, model="pv-bezier")

.. note:: Fitting the translation model to measured data is experimental.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Tuple, Union

import numpy as np

from . import _core
from ._arrays import as_float
from ._esp import ESP, as_esp
from ._models import Model

__all__ = ["IEC60891v2", "IEC60891v3", "fit"]

_SCALAR = (int, float, np.integer, np.floating)
_MODELS = (Model.PV_SPLINE, Model.PV_BEZIER)


def _iec_model(model: Union[Model, str]) -> Model:
    m = Model.lookup(model)
    if m not in _MODELS:
        names = ", ".join(x.value for x in _MODELS)
        raise ValueError(f"IEC 60891 translation supports the models {names}, got {m.value!r}")
    return m


def _broadcast(esp: ESP, **conds: Any) -> Tuple[np.ndarray, dict, bool]:
    """Broadcast esp rows and condition vectors to a common length."""
    arrays = {k: as_float(v, k) for k, v in conds.items()}
    for k, a in arrays.items():
        if a.ndim > 1:
            raise ValueError(f"{k} must be a scalar or 1-D array, got shape {a.shape}")
    rows = np.asarray(esp)
    single = rows.ndim == 1 and all(a.ndim == 0 for a in arrays.values())
    rows2 = np.atleast_2d(rows)
    try:
        (n,) = np.broadcast_shapes((rows2.shape[0],), *[a.shape for a in arrays.values()])
    except ValueError:
        shapes = ", ".join(f"{k}: {a.shape}" for k, a in arrays.items())
        raise ValueError(f"esp ({rows2.shape[0]} rows) and conditions ({shapes}) do not broadcast") from None
    rows2 = np.ascontiguousarray(np.broadcast_to(rows2, (n, 10)))
    out = {k: np.ascontiguousarray(np.broadcast_to(a, (n,))) for k, a in arrays.items()}
    return rows2, out, single


def _fit_args(esp: Any, T: Any, G: Any, max_points: Optional[int], seed_index: int):
    e = as_esp(esp)
    rows = np.atleast_2d(np.asarray(e))
    n = rows.shape[0]
    t = np.ascontiguousarray(np.broadcast_to(as_float(T, "T"), (n,)))
    g = np.ascontiguousarray(np.broadcast_to(as_float(G, "G"), (n,)))
    if n < 2:
        raise ValueError("fitting needs ESP at two or more conditions")
    if max_points is not None and (int(max_points) != max_points or max_points < 1):
        raise ValueError(f"max_points must be a positive integer, got {max_points!r}")
    if not 0 <= seed_index < n:
        raise ValueError(f"seed_index must be between 0 and {n - 1}, got {seed_index!r}")
    return rows, t, g, int(max_points or 0), int(seed_index)


def _check_fit(esp_ref: np.ndarray) -> None:
    if not np.any(esp_ref) or not np.all(np.isfinite(esp_ref)):
        raise RuntimeError("IEC 60891 fit failed; check that T and G span a sufficient range")


class _Translation:
    t_ref: float
    g_ref: float
    esp_ref: Optional[ESP]
    model: Model

    _core_fn: Any = None  # the C translation routine (set by the subclasses)
    _core_fn1: Any = None  # its single-ESP, scalar-condition variant

    def _coef_args(self) -> Tuple[float, ...]:  # pragma: no cover - abstract
        raise NotImplementedError

    def _cargs(self) -> Tuple[Any, ...]:
        """The coefficient and model arguments of the C routine, computed once per object."""
        args = self.__dict__.get("_cargs_")
        if args is None:
            args = tuple(float(v) for v in self._coef_args()) + (self.model.id,)
            object.__setattr__(self, "_cargs_", args)
        return args

    def _check_model(self) -> None:
        object.__setattr__(self, "model", _iec_model(self.model))

    def _translate_rows(self, rows, c) -> np.ndarray:
        return type(self)._core_fn(rows, c["T_from"], c["G_from"], c["T"], c["G"], *self._cargs())

    def translate(
        self,
        T: Any,
        G: Any,
        *,
        esp: Any = None,
        T_from: Any = None,
        G_from: Any = None,
    ) -> ESP:
        """Translate ESP to temperature(s) ``T`` [degC] and irradiance(s) ``G`` [W/m2].

        Without ``esp`` the reference ESP (:attr:`esp_ref`, at :attr:`t_ref` and
        :attr:`g_ref`) are translated. With ``esp``, give the conditions they
        were measured at as ``T_from`` and ``G_from`` (default: the reference
        conditions). The result is a single :class:`ESP` if all inputs are
        single/scalar and a batch otherwise.
        """
        if esp is None:
            if self.esp_ref is None:
                raise ValueError("no reference ESP stored; pass esp=, T_from= and G_from=")
            if T_from is not None or G_from is not None:
                raise ValueError("T_from/G_from only apply together with esp=")
            esp = self.esp_ref
        T_from = self.t_ref if T_from is None else T_from
        G_from = self.g_ref if G_from is None else G_from
        # fast path for the most common call (one ESP, scalar conditions), e.g.
        # inside a Kalman filter; skips the general broadcasting
        if (isinstance(T, _SCALAR) and isinstance(G, _SCALAR)
                and isinstance(T_from, _SCALAR) and isinstance(G_from, _SCALAR)):
            rows = esp._data if isinstance(esp, ESP) else np.asarray(esp, dtype=np.float64)
            if rows.shape == (10,):
                args = self.__dict__.get("_cargs_") or self._cargs()
                return ESP._wrap_owned(type(self)._core_fn1(rows, T_from, G_from, T, G, *args))
        e = as_esp(esp)
        rows, c, single = _broadcast(e, T_from=T_from, G_from=G_from, T=T, G=G)
        if rows.shape[0] == 0:
            return ESP._wrap(np.empty((0, 10)))
        out = self._translate_rows(rows, c)
        return ESP._wrap(out[0] if single else out)


@dataclass(frozen=True)
class IEC60891v2(_Translation):
    """IEC 60891 edition 2 (2009), procedure 2.

    Attributes
    ----------
    alpha:
        Relative temperature coefficient of the short-circuit current [1/degC].
    beta:
        Relative temperature coefficient of the open-circuit voltage [1/degC].
    a:
        Irradiance correction factor of the open-circuit voltage [-].
    rs:
        Internal series resistance [Ohm].
    kappa:
        Temperature coefficient of the series resistance [Ohm/degC].
    t_ref, g_ref:
        Reference conditions [degC, W/m2] (default 25 degC, 1000 W/m2).
    esp_ref:
        ESP at the reference conditions, if known (set by :meth:`fit`).
    model:
        Interpolation model of the translation, :attr:`Model.PV_SPLINE`
        (default) or :attr:`Model.PV_BEZIER`. A name string is accepted.
    """

    alpha: float
    beta: float
    a: float
    rs: float
    kappa: float
    t_ref: float = 25.0
    g_ref: float = 1000.0
    esp_ref: Optional[ESP] = field(default=None, compare=False)
    model: Model = Model.PV_SPLINE

    def __post_init__(self) -> None:
        self._check_model()

    _core_fn = staticmethod(_core.iec60891_v2_translate)
    _core_fn1 = staticmethod(_core.iec60891_v2_translate1)

    def _coef_args(self) -> Tuple[float, ...]:
        return (self.alpha, self.beta, self.a, self.rs, self.kappa)

    @classmethod
    def fit(
        cls,
        esp: Any,
        T: Any,
        G: Any,
        *,
        t_ref: float = 25.0,
        g_ref: float = 1000.0,
        max_points: Optional[int] = None,
        seed_index: int = 0,
        model: Union[Model, str] = Model.PV_SPLINE,
    ) -> "IEC60891v2":
        """Fit the coefficients and the reference ESP to ESP measured at ``T``, ``G``.

        The fit looks for the coefficients that best predict the measured ESP
        from one reference ESP: for given coefficients the reference ESP is the
        mean of the measured ESP translated to ``t_ref``, ``g_ref``, and the fit
        minimizes the distance of the reference, translated back to each
        measured condition, to the measured curve. The optimization uses at
        most 128 evenly spaced conditions; the reference ESP uses all.

        Parameters
        ----------
        esp:
            Batch :class:`ESP` (or ``(N, 10)`` array), one row per condition.
        T, G:
            Temperatures [degC] and irradiances [W/m2], one per row.
        t_ref, g_ref:
            Reference conditions to fit the ESP for.
        max_points:
            If given, use furthest point sampling in the G-T plane to select at
            most this many conditions for the fit. This keeps the full range
            of conditions, but it prefers the extreme conditions.
        seed_index:
            Row the furthest point sampling starts from.
        model:
            Interpolation model used in the fit, stored in the result and used
            by its :meth:`translate`.
        """
        m = _iec_model(model)
        rows, t, g, fn, off = _fit_args(esp, T, G, max_points, seed_index)
        esp_ref, alpha, beta, a, rs, kappa = _core.iec60891_v2_fit(rows, t, g, float(t_ref), float(g_ref), fn, off, m.id)
        _check_fit(esp_ref)
        return cls(alpha, beta, a, rs, kappa, float(t_ref), float(g_ref), ESP._wrap(esp_ref), m)


@dataclass(frozen=True)
class IEC60891v3(_Translation):
    """IEC 60891 edition 3 (2021), procedure 2.

    Attributes
    ----------
    alpha:
        Relative temperature coefficient of the short-circuit current [1/degC].
    beta:
        Relative temperature coefficient of the open-circuit voltage [1/degC].
    b1, b2:
        Irradiance correction coefficients of the open-circuit voltage [-].
    rs:
        Internal series resistance [Ohm].
    kappa:
        Temperature coefficient of the series resistance [Ohm/degC].
    t_ref, g_ref:
        Reference conditions [degC, W/m2] (default 25 degC, 1000 W/m2).
    voc_ref:
        Open-circuit voltage at the reference conditions [V]. Taken from
        ``esp_ref`` when that is given.
    esp_ref:
        ESP at the reference conditions, if known (set by :meth:`fit`).
    model:
        Interpolation model of the translation, as :attr:`IEC60891v2.model`.
    """

    alpha: float
    beta: float
    b1: float
    b2: float
    rs: float
    kappa: float
    t_ref: float = 25.0
    g_ref: float = 1000.0
    voc_ref: Optional[float] = None
    esp_ref: Optional[ESP] = field(default=None, compare=False)
    model: Model = Model.PV_SPLINE

    def __post_init__(self) -> None:
        self._check_model()
        if self.esp_ref is not None:
            if self.esp_ref.is_batch:
                raise ValueError("esp_ref must be a single ESP")
            if self.voc_ref is None:
                object.__setattr__(self, "voc_ref", self.esp_ref.voc)
        if self.voc_ref is None:
            raise TypeError("IEC60891v3 needs voc_ref (or esp_ref)")

    _core_fn = staticmethod(_core.iec60891_v3_translate)
    _core_fn1 = staticmethod(_core.iec60891_v3_translate1)

    def _coef_args(self) -> Tuple[float, ...]:
        return (self.voc_ref, self.t_ref, self.g_ref, self.alpha, self.beta,
                self.b1, self.b2, self.rs, self.kappa)

    @classmethod
    def fit(
        cls,
        esp: Any,
        T: Any,
        G: Any,
        *,
        t_ref: float = 25.0,
        g_ref: float = 1000.0,
        max_points: Optional[int] = None,
        seed_index: int = 0,
        model: Union[Model, str] = Model.PV_SPLINE,
    ) -> "IEC60891v3":
        """Fit the coefficients and the reference ESP. Arguments as :meth:`IEC60891v2.fit`."""
        m = _iec_model(model)
        rows, t, g, fn, off = _fit_args(esp, T, G, max_points, seed_index)
        esp_ref, alpha, beta, b1, b2, rs, kappa = _core.iec60891_v3_fit(rows, t, g, float(t_ref), float(g_ref), fn, off, m.id)
        _check_fit(esp_ref)
        return cls(alpha, beta, b1, b2, rs, kappa, float(t_ref), float(g_ref), esp_ref=ESP._wrap(esp_ref), model=m)


def fit(
    esp: Any,
    T: Any,
    G: Any,
    *,
    version: int = 3,
    t_ref: float = 25.0,
    g_ref: float = 1000.0,
    max_points: Optional[int] = None,
    seed_index: int = 0,
    model: Union[Model, str] = Model.PV_SPLINE,
) -> Union[IEC60891v2, IEC60891v3]:
    """Fit IEC 60891 ``version`` 2 or 3 (default). See :meth:`IEC60891v2.fit`."""
    cls = {2: IEC60891v2, 3: IEC60891v3}.get(version)
    if cls is None:
        raise ValueError(f"version must be 2 or 3, got {version!r}")
    return cls.fit(esp, T, G, t_ref=t_ref, g_ref=g_ref, max_points=max_points, seed_index=seed_index, model=model)
