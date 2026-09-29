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

.. note:: Fitting the translation model to measured data is experimental.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Tuple, Union

import numpy as np

from . import _core
from ._arrays import as_float
from ._esp import ESP, as_esp

__all__ = ["IEC60891v2", "IEC60891v3", "fit"]


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

    def _translate_rows(self, rows, conds) -> np.ndarray:  # pragma: no cover - abstract
        raise NotImplementedError

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
        e = as_esp(esp)
        T_from = self.t_ref if T_from is None else T_from
        G_from = self.g_ref if G_from is None else G_from
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
    """

    alpha: float
    beta: float
    a: float
    rs: float
    kappa: float
    t_ref: float = 25.0
    g_ref: float = 1000.0
    esp_ref: Optional[ESP] = field(default=None, compare=False)

    def _translate_rows(self, rows, c) -> np.ndarray:
        return _core.iec60891_v2_translate(
            rows, c["T_from"], c["G_from"], c["T"], c["G"],
            float(self.alpha), float(self.beta), float(self.a), float(self.rs), float(self.kappa),
        )

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
    ) -> "IEC60891v2":
        """Fit the coefficients and the reference ESP to ESP measured at ``T``, ``G``.

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
            of conditions but limits the cost for large datasets.
        seed_index:
            Row the furthest point sampling starts from.
        """
        rows, t, g, fn, off = _fit_args(esp, T, G, max_points, seed_index)
        esp_ref, alpha, beta, a, rs, kappa = _core.iec60891_v2_fit(rows, t, g, float(t_ref), float(g_ref), fn, off)
        _check_fit(esp_ref)
        return cls(alpha, beta, a, rs, kappa, float(t_ref), float(g_ref), ESP._wrap(esp_ref))


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

    def __post_init__(self) -> None:
        if self.esp_ref is not None:
            if self.esp_ref.is_batch:
                raise ValueError("esp_ref must be a single ESP")
            if self.voc_ref is None:
                object.__setattr__(self, "voc_ref", self.esp_ref.voc)
        if self.voc_ref is None:
            raise TypeError("IEC60891v3 needs voc_ref (or esp_ref)")

    def _translate_rows(self, rows, c) -> np.ndarray:
        return _core.iec60891_v3_translate(
            rows, c["T_from"], c["G_from"], c["T"], c["G"],
            float(self.voc_ref), float(self.t_ref), float(self.g_ref),
            float(self.alpha), float(self.beta), float(self.b1), float(self.b2),
            float(self.rs), float(self.kappa),
        )

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
    ) -> "IEC60891v3":
        """Fit the coefficients and the reference ESP. Arguments as :meth:`IEC60891v2.fit`."""
        rows, t, g, fn, off = _fit_args(esp, T, G, max_points, seed_index)
        esp_ref, alpha, beta, b1, b2, rs, kappa = _core.iec60891_v3_fit(rows, t, g, float(t_ref), float(g_ref), fn, off)
        _check_fit(esp_ref)
        return cls(alpha, beta, b1, b2, rs, kappa, float(t_ref), float(g_ref), esp_ref=ESP._wrap(esp_ref))


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
) -> Union[IEC60891v2, IEC60891v3]:
    """Fit IEC 60891 ``version`` 2 or 3 (default). See :meth:`IEC60891v2.fit`."""
    cls = {2: IEC60891v2, 3: IEC60891v3}.get(version)
    if cls is None:
        raise ValueError(f"version must be 2 or 3, got {version!r}")
    return cls.fit(esp, T, G, t_ref=t_ref, g_ref=g_ref, max_points=max_points, seed_index=seed_index)
