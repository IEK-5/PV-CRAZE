"""Curve models: the :class:`Model` enum, fit options and model parameters."""

from __future__ import annotations

import enum
import re
from dataclasses import dataclass
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np

from . import _core
from ._arrays import as_float, normalize_curves
from ._esp import ESP, _import_pandas


def _identifier(label: str) -> str:
    """C parameter label -> Python attribute name ("Vqmp-" -> "vqmp_minus")."""
    return label.replace("-", "_minus").replace("+", "_plus").lower()


def _normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


class Model(str, enum.Enum):
    """The curve models implemented by PV-CRAZE.

    Members compare equal to their library name, and lookup is forgiving about
    case and punctuation, so ``Model("pv-bezier")``, ``Model("PV_BEZIER")`` and
    ``Model.PV_BEZIER`` are the same. Anywhere a model is expected you can pass
    either a member or a name string.

    ``PV_BEZIER`` is the recommended interpolation model. ``SPLINE`` is the only
    model that can additionally be refined against the measured IV data
    (``FitSpec(..., refine=True)``); all other models are parameterised from the
    ESP alone.
    """

    #: Derivative respecting piecewise linear interpolation of the ESP.
    APWL = "APWL"
    #: Quadratic Bezier interpolation of the ESP.
    B2 = "B2"
    #: One-diode model with series and shunt resistance.
    DIODE = "Diode"
    #: Karmalkar-Haneefa model.
    KH = "KH"
    #: Pindado-Cubas model.
    PC = "PC"
    #: Cubic Bezier interpolation of the ESP (recommended).
    PV_BEZIER = "PV-Bezier"
    #: Monotone cubic spline interpolation of the ESP.
    PV_SPLINE = "PV-spline"
    #: Plain piecewise linear interpolation of the ESP key points.
    PWL = "PWL"
    #: Primitive cubic spline; may be non-monotone. Supports ``refine``.
    SPLINE = "Spline"

    @classmethod
    def _missing_(cls, value: object) -> Optional["Model"]:
        if isinstance(value, str):
            key = _normalize_name(value)
            for m in cls:
                if _normalize_name(m.value) == key or _normalize_name(m.name) == key:
                    return m
        return None

    @classmethod
    def lookup(cls, model: Union["Model", str]) -> "Model":
        """Return the member for ``model``, raising a helpful ``ValueError``."""
        if isinstance(model, Model):
            return model
        if not isinstance(model, str):
            raise TypeError(f"a model must be a Model member or a name string, got {type(model).__name__}")
        try:
            return cls(model)
        except ValueError:
            names = ", ".join(m.value for m in cls)
            raise ValueError(f"unknown model {model!r}; available models: {names}") from None

    @property
    def id(self) -> int:
        """Numeric model id used by the C library."""
        return _MODEL_INFO[self][0]

    @property
    def param_names(self) -> Tuple[str, ...]:
        """Names of this model's parameters, as attributes of :class:`ModelParams`."""
        return tuple(_identifier(label) for label in _MODEL_INFO[self][1])

    @property
    def param_labels(self) -> Tuple[str, ...]:
        """Parameter labels as used in the C library and the literature (e.g. ``"Vqmp-"``)."""
        return _MODEL_INFO[self][1]

    def __str__(self) -> str:
        return self.value


_MODEL_INFO: Dict[Model, Tuple[int, Tuple[str, ...]]] = {}
_MODEL_BY_ID: Dict[int, Model] = {}
for _id, _name, _labels in _core.model_table():
    _m = Model(_name)
    _MODEL_INFO[_m] = (_id, tuple(_labels))
    _MODEL_BY_ID[_id] = _m
_missing = [m for m in Model if m not in _MODEL_INFO]
if _missing:  # pragma: no cover - guards against Python/C getting out of sync
    raise ImportError(f"pvcraze._core does not provide the models {_missing}; rebuild pvcraze")

# fit option limits enforced by the C parser (fitmodels.c)
_MAX_ITER = 10000
_EPS_RANGE = (1e-12, 1e-1)
_POP_RANGE = (10, 1000)


@dataclass(frozen=True)
class FitSpec:
    """A model plus optional settings for fitting it.

    Use a plain model name where the defaults are fine and a ``FitSpec`` where
    they are not::

        pvcraze.analyze(curves, models=["PV-Bezier",
                                        FitSpec("Diode", population=200, max_iter=4000)])

    Options left at ``None`` use the library defaults. Not every option
    affects every model.

    Attributes
    ----------
    model:
        The model to fit (member or name).
    max_iter:
        Maximum iterations of the parameterisation routine (default 50, max 10000).
    tol:
        Stop criterion of the parameterisation routine (default 1e-5).
    population:
        Population size for population based optimisers (default 100, 10-1000).
    eval_max_iter:
        Maximum iterations when evaluating the model, i.e. solving I(V) or V(I)
        (default 50).
    eval_tol:
        Stop criterion when evaluating the model (default 1e-5).
    refine_tol:
        Error threshold above which expensive extra optimisation is attempted
        (default 1e-2).
    refine:
        ``Spline`` model only: also optimise the parameters against the
        measured IV data (only possible in :func:`pvcraze.analyze`).
    """

    model: Union[Model, str]
    max_iter: Optional[int] = None
    tol: Optional[float] = None
    population: Optional[int] = None
    eval_max_iter: Optional[int] = None
    eval_tol: Optional[float] = None
    refine_tol: Optional[float] = None
    refine: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "model", Model.lookup(self.model))
        for name in ("max_iter", "eval_max_iter"):
            v = getattr(self, name)
            if v is not None and (int(v) != v or not 0 <= v <= _MAX_ITER):
                raise ValueError(f"{name} must be an integer between 0 and {_MAX_ITER}, got {v!r}")
        for name in ("tol", "eval_tol", "refine_tol"):
            v = getattr(self, name)
            if v is not None and not _EPS_RANGE[0] <= v <= _EPS_RANGE[1]:
                raise ValueError(f"{name} must be between {_EPS_RANGE[0]:g} and {_EPS_RANGE[1]:g}, got {v!r}")
        p = self.population
        if p is not None and (int(p) != p or not _POP_RANGE[0] <= p <= _POP_RANGE[1]):
            raise ValueError(f"population must be an integer between {_POP_RANGE[0]} and {_POP_RANGE[1]}, got {p!r}")
        if self.refine and self.model is not Model.SPLINE:
            raise ValueError("refine is only supported by the Spline model")

    @classmethod
    def coerce(cls, spec: Union["FitSpec", Model, str]) -> "FitSpec":
        return spec if isinstance(spec, FitSpec) else cls(spec)

    def _spec_string(self) -> str:
        """The model specification string understood by the C library."""
        parts = [f"model({self.model.value})"]
        for attr, key, fmt in (
            ("max_iter", "fitit", "d"),
            ("tol", "fiteps", "r"),
            ("population", "npop", "d"),
            ("eval_max_iter", "mit", "d"),
            ("eval_tol", "meps", "r"),
            ("refine_tol", "oeps", "r"),
        ):
            v = getattr(self, attr)
            if v is not None:
                parts.append(f"{key}({int(v)})" if fmt == "d" else f"{key}({float(v)!r})")
        if self.refine:
            parts.append("refine")
        return " ".join(parts)


ModelLike = Union[Model, str, FitSpec]


class ModelParams:
    """Parameters of one curve model for a single curve or a batch of curves.

    The parameters are available as attributes named after
    :attr:`Model.param_names`, e.g. for the diode model::

        p.i0, p.nkt, p.iph, p.rs, p.rsh

    and for the ESP-based interpolation models (PV-Bezier, PV-spline, B2,
    Spline, PWL, APWL) the same names as :class:`ESP` (``voc``, ``isc``, ...).

    Like :class:`ESP`, a ``ModelParams`` is either single (attributes are
    floats) or a batch of ``N`` parameter sets (attributes are arrays) that
    supports ``len``, iteration and indexing.

    Parameters
    ----------
    model:
        The model (member or name).
    values:
        Array of shape ``(k,)`` or ``(N, k)``, with ``k = len(model.param_names)``.
        Alternatively give the parameters by name as keyword arguments.

    Examples
    --------
    >>> p = ModelParams("Diode", i0=1e-8, nkt=0.0385, iph=0.035, rs=0.1, rsh=700)
    >>> p.iv([0.0, 0.3, 0.5])        # current at the given voltages
    >>> p.esp()                      # the ESP of the modelled curve
    """

    __slots__ = ("_model", "_data")

    def __init__(self, model: Union[Model, str], values: Any = None, **named: Any) -> None:
        m = Model.lookup(model)
        names = m.param_names
        if values is not None and named:
            raise TypeError("give either values or named parameters, not both")
        if values is None:
            unknown = set(named) - set(names)
            missing = [n for n in names if n not in named]
            if unknown:
                raise TypeError(f"unknown {m.value} parameters: {', '.join(sorted(unknown))}; expected {names}")
            if missing:
                raise TypeError(f"missing {m.value} parameters: {', '.join(missing)}")
            cols = np.broadcast_arrays(*[as_float(named[n], n) for n in names])
            if cols[0].ndim > 1:
                raise ValueError("parameters must be scalars or 1-D arrays")
            data = np.stack(cols, axis=-1)
        else:
            data = as_float(values, "values")
            if data.ndim not in (1, 2) or data.shape[-1] != len(names):
                raise ValueError(
                    f"{m.value} has {len(names)} parameters {names}; values must have shape "
                    f"({len(names)},) or (N, {len(names)}), got {data.shape}"
                )
        self._init(m, data)

    def _init(self, model: Model, data: np.ndarray) -> None:
        data = np.array(data, dtype=np.float64, order="C")
        data.flags.writeable = False
        object.__setattr__(self, "_model", model)
        object.__setattr__(self, "_data", data)

    @classmethod
    def _wrap(cls, model: Model, data: np.ndarray) -> "ModelParams":
        obj = cls.__new__(cls)
        obj._init(model, data)
        return obj

    @classmethod
    def _from_core(cls, type_id: int, padded: np.ndarray) -> "ModelParams":
        model = _MODEL_BY_ID[int(type_id)]
        return cls._wrap(model, padded[..., : len(model.param_names)])

    # --------------------------------------------------------- properties
    @property
    def model(self) -> Model:
        """The curve model these parameters belong to."""
        return self._model

    @property
    def names(self) -> Tuple[str, ...]:
        """Parameter names, in column order."""
        return self._model.param_names

    @property
    def is_batch(self) -> bool:
        """True for a batch of parameter sets, False for a single one."""
        return self._data.ndim == 2

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            k = self._model.param_names.index(name)
        except ValueError:
            raise AttributeError(
                f"{self._model.value} parameters have no attribute {name!r}; available: {self._model.param_names}"
            ) from None
        return float(self._data[k]) if self._data.ndim == 1 else self._data[:, k]

    def __dir__(self) -> List[str]:
        return sorted(set(super().__dir__()) | set(self._model.param_names))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("ModelParams objects are immutable")

    def to_array(self) -> np.ndarray:
        """Return a writable copy of the ``(k,)`` or ``(N, k)`` parameter array."""
        return self._data.copy()

    def to_dict(self) -> Dict[str, Any]:
        """Return ``{name: value}``; values are floats (single) or arrays (batch)."""
        if self.is_batch:
            return {n: self._data[:, k].copy() for k, n in enumerate(self.names)}
        return {n: float(v) for n, v in zip(self.names, self._data)}

    def to_dataframe(self) -> Any:
        """Return a pandas DataFrame with one row per parameter set. Requires pandas."""
        pd = _import_pandas()
        return pd.DataFrame(np.atleast_2d(self._data), columns=list(self.names))

    def __array__(self, dtype: Any = None, copy: Any = None) -> np.ndarray:
        a = self._data if dtype is None else self._data.astype(dtype)
        return a.copy() if copy else a

    # ---------------------------------------------------------- container
    def __len__(self) -> int:
        if not self.is_batch:
            raise TypeError("single ModelParams have no len(); check .is_batch first")
        return self._data.shape[0]

    def __iter__(self) -> Iterator["ModelParams"]:
        if not self.is_batch:
            raise TypeError("single ModelParams are not iterable")
        return (ModelParams._wrap(self._model, row) for row in self._data)

    def __getitem__(self, index: Any) -> "ModelParams":
        if not self.is_batch:
            raise TypeError("single ModelParams cannot be indexed; use the parameter attributes")
        return ModelParams._wrap(self._model, self._data[index])

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ModelParams):
            return NotImplemented
        return (
            self._model is other._model
            and self._data.shape == other._data.shape
            and bool(np.array_equal(self._data, other._data))
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        if not self.is_batch:
            body = ", ".join(f"{n}={v:.6g}" for n, v in zip(self.names, self._data))
            return f"ModelParams({self._model.value!r}, {body})"
        with np.printoptions(precision=4, threshold=6, edgeitems=2):
            body = ", ".join(f"{n}={self._data[:, k]}" for k, n in enumerate(self.names[:3]))
        return f"<ModelParams {self._model.value} batch n={self._data.shape[0]}: {body}, ...>"

    # ------------------------------------------------------ core helpers
    def _core_args(self) -> Tuple[np.ndarray, np.ndarray]:
        rows = np.atleast_2d(self._data)
        padded = np.zeros((rows.shape[0], 10))
        padded[:, : rows.shape[1]] = rows
        types = np.full(rows.shape[0], self._model.id, dtype=np.int64)
        return types, padded

    def _evaluate(self, x: Any, name: str, direction: int) -> Any:
        x = as_float(x, name)
        types, pars = self._core_args()
        func = _core.vi if direction else _core.iv
        if not self.is_batch:
            out = func(types, pars, x.reshape(1, -1))
            return float(out[0, 0]) if x.ndim == 0 else out.reshape(x.shape)
        n = len(self)
        if x.ndim == 0:
            return func(types, pars, np.full((n, 1), float(x)))[:, 0]
        if x.ndim == 1:
            return func(types, pars, np.ascontiguousarray(np.broadcast_to(x, (n, x.shape[0]))))
        if x.ndim == 2 and x.shape[0] == n:
            return func(types, pars, x)
        raise ValueError(
            f"for a batch of {n} parameter sets {name} must be a scalar, a 1-D grid "
            f"shared by all sets, or an array of shape ({n}, m); got shape {x.shape}"
        )

    # -------------------------------------------------------- evaluation
    def iv(self, V: Any) -> Any:
        """Current [A] at voltage(s) ``V`` [V].

        For single parameters the result has the shape of ``V``. For a batch
        of ``N`` sets: a scalar ``V`` gives shape ``(N,)``, a 1-D grid of
        length ``m`` (shared by all sets) gives ``(N, m)``, and an ``(N, m)``
        array evaluates row ``i`` with parameter set ``i``. To evaluate one
        voltage per set use ``V[:, None]``.
        """
        return self._evaluate(V, "V", 0)

    def vi(self, I: Any) -> Any:  # noqa: E741 - physics notation
        """Voltage [V] at current(s) ``I`` [A]. Broadcasting as in :meth:`iv`."""
        return self._evaluate(I, "I", 1)

    def curve(self, npoints: int = 100) -> Any:
        """Sample the modelled IV curve between short circuit and open circuit.

        Returns an ``(n, 2)`` array with columns (V, I) for single parameters,
        and a list of such arrays for a batch. ``n`` is close to, but not
        necessarily exactly, ``npoints``.
        """
        if int(npoints) != npoints or npoints < 2:
            raise ValueError(f"npoints must be an integer >= 2, got {npoints!r}")
        out = _core.curves(*self._core_args(), int(npoints))
        return out if self.is_batch else out[0]

    def esp(self) -> ESP:
        """The extended solar cell parameters of the modelled curve(s)."""
        data = _core.model_esp(*self._core_args())
        return ESP._wrap(data if self.is_batch else data[0])

    def _curves_for(self, curves: Any, I: Any) -> List[np.ndarray]:
        lst, single = normalize_curves(curves, I)
        if single != (not self.is_batch):
            raise ValueError(
                "give one curve for single parameters and a sequence of curves for a batch"
            )
        if self.is_batch and len(lst) != len(self):
            raise ValueError(f"got {len(lst)} curves for {len(self)} parameter sets")
        return lst

    def rmse(self, curves: Any, I: Any = None) -> Any:
        """RMS distance between measured curve(s) and the model.

        Distances are Euclidean in the (V/Voc, I/Isc) plane, so the result is
        dimensionless. Every point of the given curve counts, including points
        outside the first quadrant. Takes a curve in the same forms as
        :func:`pvcraze.analyze` (one per parameter set for a batch).
        """
        lst = self._curves_for(curves, I)
        e = _core.rmse(*self._core_args(), lst)
        return e if self.is_batch else float(e[0])

    def distance(self, curves: Any, I: Any = None) -> Any:
        """Point-wise distance between measured curve(s) and the model.

        For every measured point this is the Euclidean distance to the closest
        point of the modelled curve (between short and open circuit) in the
        (V/Voc, I/Isc) plane. Returns an ``(n, 2)`` array with columns
        (V, distance) per curve, a list of those for a batch.
        """
        lst = self._curves_for(curves, I)
        out = _core.distance_curves(*self._core_args(), lst)
        for d in out:
            np.sqrt(d[:, 1], out=d[:, 1])  # the library returns squared distances
        return out if self.is_batch else out[0]


def _params_from_core(types: np.ndarray, pars: np.ndarray, single: bool) -> ModelParams:
    """Build ModelParams from one model's (N, 10) block of a core result."""
    model_id = int(types) if np.ndim(types) == 0 else int(types[0])
    p = ModelParams._from_core(model_id, pars)
    return p[0] if single else p


def coerce_specs(models: Union[ModelLike, Sequence[ModelLike], None]) -> List[FitSpec]:
    if models is None:
        return []
    if isinstance(models, (str, Model, FitSpec)):
        models = [models]
    specs = [FitSpec.coerce(m) for m in models]
    seen = set()
    for s in specs:
        if s.model in seen:
            raise ValueError(
                f"model {s.model.value} requested more than once; results are keyed by model, "
                "so call again with the other settings"
            )
        seen.add(s.model)
    return specs
