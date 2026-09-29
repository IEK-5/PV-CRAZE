"""PV-CRAZE: Photo-Voltaic CuRve AnalyZEr.

Extract the *extended solar cell parameters* (ESP) from measured IV curves,
fit curve models to them, and translate them between irradiances and
temperatures.

Quick start::

    import pvcraze as pc

    res = pc.analyze(V, I, models="PV-Bezier")   # or pc.analyze(list_of_curves)
    res.ok                          # False if the analysis flagged an error
    res.flag().messages             # what went wrong / what to watch out for
    res.esp.voc, res.esp.ff         # the ESP (see pc.ESP for all fields)
    bez = res.fits["PV-Bezier"].params
    bez.iv([0.0, 0.3])              # evaluate the model
    bez.curve(200)                  # (n, 2) array with the modelled curve

    fit = pc.fit_model(res.esp, "Diode", population=200)
    fit.params.rs                   # named model parameters

Main objects:

* :func:`analyze`, :class:`AnalysisResult` -- IV curve analysis
* :class:`ESP` -- extended solar cell parameters (single or batch)
* :class:`Model`, :class:`FitSpec`, :func:`fit_model`, :class:`FitResult`
* :class:`ModelParams` -- model parameters; evaluate with ``iv``, ``vi``,
  ``curve``, ``esp``, ``rmse``
* :class:`Flag` -- analysis warnings and errors
* :mod:`pvcraze.iec60891`, :mod:`pvcraze.diode_gt` -- irradiance/temperature
  translation (degrees Celsius, W/m2)
* :func:`verbose`, :func:`set_verbosity` -- the library is silent by default

Units are SI (V, A, S, Ohm) unless noted; temperatures in degrees Celsius.
"""

from importlib.metadata import PackageNotFoundError, version as _version

from . import _core
from ._analysis import AnalysisResult, FitResult, analyze, fit_model
from ._esp import ESP
from ._flags import ERROR_MASK, WARNING_MASK, Flag, describe as describe_flags
from ._models import FitSpec, Model, ModelParams
from ._verbosity import get_verbosity, set_verbosity, verbose
from . import diode_gt, iec60891

try:
    __version__ = _version("pvcraze")
except PackageNotFoundError:  # pragma: no cover - running from a source tree
    __version__ = "unknown"

__all__ = [
    "AnalysisResult",
    "ERROR_MASK",
    "ESP",
    "FitResult",
    "FitSpec",
    "Flag",
    "Model",
    "ModelParams",
    "WARNING_MASK",
    "analyze",
    "describe_flags",
    "diode_gt",
    "fit_model",
    "get_verbosity",
    "iec60891",
    "set_verbosity",
    "verbose",
]

_core.set_verbosity(0)  # a library should not print unless asked to
