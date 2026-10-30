"""Irradiance and temperature dependent one-diode model.

The diode parameters at a reference condition are translated to other
conditions with

* ``Iph = Iph_ref * (1 + alpha * (T - T_ref)) * G / G_ref``
* ``I0 = I0_ref * (T/T_ref)**3 * exp(Ea * (T - T_ref) / (k * T * T_ref))``
* ``nkT = nkT_ref * T / T_ref``

(temperatures in kelvin inside the formulas) and ``Rs``, ``Rsh`` constant::

    from pvcraze import diode_gt
    gt = diode_gt.DiodeGT.fit(esp, T=T, G=G)      # or diode_gt.fit(...)
    gt.params                                    # diode ModelParams at 25 degC, 1000 W/m2
    gt.at(T=60, G=800).iv(0.5)                   # model at other conditions
    gt.esp(T=T, G=G)                             # modelled ESP at the measured conditions

Temperatures are in degrees Celsius, irradiances in W/m2.

.. note:: Fitting the model to measured data is experimental.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from . import _core
from ._arrays import as_float
from ._esp import ESP
from ._models import Model, ModelParams, _params_from_core
from .iec60891 import _fit_args

__all__ = ["DiodeGT", "fit"]


@dataclass(frozen=True)
class DiodeGT:
    """One-diode model with irradiance and temperature dependence.

    Attributes
    ----------
    params:
        Single diode :class:`ModelParams` at the reference conditions.
    alpha:
        Relative temperature coefficient of the photocurrent [1/degC].
    ea:
        Activation energy of the saturation current [eV].
    t_ref, g_ref:
        Reference conditions [degC, W/m2] (default 25 degC, 1000 W/m2).
    """

    params: ModelParams
    alpha: float
    ea: float
    t_ref: float = 25.0
    g_ref: float = 1000.0

    def __post_init__(self) -> None:
        if not isinstance(self.params, ModelParams) or self.params.model is not Model.DIODE:
            raise TypeError("params must be ModelParams of the Diode model")
        if self.params.is_batch:
            raise ValueError("params must be a single parameter set")
        if self.t_ref <= -273.15:
            raise ValueError(f"t_ref must be above absolute zero, got {self.t_ref} degC")
        if self.g_ref <= 0:
            raise ValueError(f"g_ref must be positive, got {self.g_ref}")

    def at(self, T: Any, G: Any) -> ModelParams:
        """Diode parameters at temperature(s) ``T`` [degC] and irradiance(s) ``G`` [W/m2].

        Scalars give single :class:`ModelParams`, arrays (broadcast together)
        give a batch.
        """
        t, g = np.broadcast_arrays(as_float(T, "T"), as_float(G, "G"))
        if t.ndim > 1:
            raise ValueError("T and G must be scalars or 1-D arrays")
        if np.any(t <= -273.15):
            raise ValueError("T must be above absolute zero (degrees Celsius expected)")
        if np.any(g <= 0):
            raise ValueError("G must be positive")
        single = t.ndim == 0
        t1, g1 = np.ascontiguousarray(np.atleast_1d(t)), np.ascontiguousarray(np.atleast_1d(g))
        if t1.size == 0:
            return ModelParams._wrap(Model.DIODE, np.empty((0, len(Model.DIODE.param_names))))
        _, padded = self.params._core_args()
        types, pars = _core.diode_gt_translate(
            padded[0], float(self.t_ref), float(self.g_ref), t1, g1, float(self.alpha), float(self.ea)
        )
        return _params_from_core(types, pars, single)

    def esp(self, T: Any, G: Any) -> ESP:
        """Modelled ESP at the given conditions; shorthand for ``self.at(T, G).esp()``."""
        return self.at(T, G).esp()

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
    ) -> "DiodeGT":
        """Fit the model to ESP measured at temperatures ``T`` [degC] and irradiances ``G`` [W/m2].

        ``max_points`` and ``seed_index`` control furthest point sampling as in
        :meth:`pvcraze.iec60891.IEC60891v2.fit`.
        """
        rows, t, g, fn, off = _fit_args(esp, T, G, max_points, seed_index)
        pars, alpha, ea = _core.diode_gt_fit(rows, t, g, float(t_ref), float(g_ref), fn, off)
        params = ModelParams._from_core(Model.DIODE.id, pars)
        if not np.all(np.isfinite(params.to_array())) or not np.any(params.to_array()):
            raise RuntimeError("diode G-T fit failed; check that T and G span a sufficient range")
        return cls(params, float(alpha), float(ea), float(t_ref), float(g_ref))


def fit(esp: Any, T: Any, G: Any, **kwargs: Any) -> DiodeGT:
    """Fit a :class:`DiodeGT` model; see :meth:`DiodeGT.fit`."""
    return DiodeGT.fit(esp, T, G, **kwargs)
