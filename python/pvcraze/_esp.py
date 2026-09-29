"""The extended solar cell parameters (ESP) container."""

from __future__ import annotations

from typing import Any, Dict, Iterable, Iterator, Mapping, Union

import numpy as np

from ._arrays import as_float

FIELDS = (
    "voc",
    "isc",
    "vmpp",
    "impp",
    "vqmp_minus",
    "iqmp_minus",
    "vqmp_plus",
    "iqmp_plus",
    "gsc",
    "roc",
)

Scalar = Union[float, np.floating]
FloatOrArray = Union[float, np.ndarray]


class ESP:
    """Extended solar cell parameters for one IV curve or a batch of curves.

    The ESP describe an illuminated IV curve (generator convention: positive
    current between short circuit and open circuit) by five key points and two
    slopes:

    ==============  ====  ======================================================
    field           unit  meaning
    ==============  ====  ======================================================
    ``voc``         V     open-circuit voltage
    ``isc``         A     short-circuit current
    ``vmpp``        V     voltage at the maximum power point (MPP)
    ``impp``        A     current at the MPP
    ``vqmp_minus``  V     voltage of the qmp- point, the maximum power point of
                          the curve segment between short circuit and the MPP,
                          measured relative to the point (0, Impp)
    ``iqmp_minus``  A     current of the qmp- point
    ``vqmp_plus``   V     voltage of the qmp+ point, the same construction for
                          the segment between the MPP and open circuit
    ``iqmp_plus``   A     current of the qmp+ point
    ``gsc``         S     -dI/dV at short circuit
    ``roc``         Ohm   -dV/dI at open circuit
    ==============  ====  ======================================================

    An ``ESP`` is either *single* (fields are floats) or a *batch* of ``N``
    curves (fields are 1-D arrays of length ``N``). A batch supports ``len``,
    iteration, integer indexing (giving a single ``ESP``) and slicing, boolean
    masks or index arrays (giving a batch). ``numpy.asarray(esp)`` returns the
    raw ``(10,)`` or ``(N, 10)`` array with columns in the order of
    :data:`ESP.FIELDS`. ESP objects are immutable.

    Construct from named values (scalars or arrays, broadcast together)::

        esp = ESP(voc=0.6, isc=0.035, vmpp=0.5, impp=0.032,
                  vqmp_minus=0.42, iqmp_minus=0.034, vqmp_plus=0.56,
                  iqmp_plus=0.022, gsc=1.7e-3, roc=0.96)

    or from an array with :meth:`from_array`, or a mapping / DataFrame with
    :meth:`from_dict` / :meth:`from_dataframe`.
    """

    FIELDS = FIELDS
    __slots__ = ("_data",)

    def __init__(
        self,
        voc: Any,
        isc: Any,
        vmpp: Any,
        impp: Any,
        vqmp_minus: Any,
        iqmp_minus: Any,
        vqmp_plus: Any,
        iqmp_plus: Any,
        gsc: Any,
        roc: Any,
    ) -> None:
        values = [voc, isc, vmpp, impp, vqmp_minus, iqmp_minus, vqmp_plus, iqmp_plus, gsc, roc]
        cols = np.broadcast_arrays(*[as_float(v, f) for v, f in zip(values, FIELDS)])
        if cols[0].ndim > 1:
            raise ValueError("ESP fields must be scalars or 1-D arrays")
        self._set(np.stack(cols, axis=-1))

    def _set(self, data: np.ndarray) -> None:
        data = np.array(data, dtype=np.float64, order="C")  # private copy
        data.flags.writeable = False
        object.__setattr__(self, "_data", data)

    @classmethod
    def _wrap(cls, data: np.ndarray) -> "ESP":
        obj = cls.__new__(cls)
        obj._set(data)
        return obj

    # ------------------------------------------------------------------ io
    @classmethod
    def from_array(cls, array: Any) -> "ESP":
        """Create from a ``(10,)`` (single) or ``(N, 10)`` (batch) array.

        Columns must be in the order of :data:`ESP.FIELDS`
        (Voc, Isc, Vmpp, Impp, Vqmp-, Iqmp-, Vqmp+, Iqmp+, Gsc, Roc).
        """
        a = as_float(array, "array")
        if a.ndim not in (1, 2) or a.shape[-1] != 10:
            raise ValueError(f"expected an array of shape (10,) or (N, 10), got {a.shape}")
        return cls._wrap(a)

    @classmethod
    def from_dict(cls, mapping: Mapping[str, Any]) -> "ESP":
        """Create from a mapping with the keys in :data:`ESP.FIELDS`."""
        missing = [f for f in FIELDS if f not in mapping]
        if missing:
            raise KeyError(f"missing ESP fields: {', '.join(missing)}")
        return cls(**{f: mapping[f] for f in FIELDS})

    @classmethod
    def from_dataframe(cls, df: Any) -> "ESP":
        """Create a batch from a pandas DataFrame with columns named as :data:`ESP.FIELDS`."""
        return cls(**{f: np.asarray(df[f].to_numpy(), dtype=np.float64) for f in FIELDS})

    @classmethod
    def concat(cls, items: Iterable["ESP"]) -> "ESP":
        """Concatenate single and batch ESP objects into one batch."""
        rows = [np.atleast_2d(e._data) for e in items]
        return cls._wrap(np.concatenate(rows, axis=0) if rows else np.empty((0, 10)))

    def to_array(self) -> np.ndarray:
        """Return a writable copy of the ``(10,)`` or ``(N, 10)`` array."""
        return self._data.copy()

    def to_dict(self) -> Dict[str, FloatOrArray]:
        """Return ``{field: value}``, values are floats (single) or array copies (batch)."""
        if self.is_batch:
            return {f: self._data[:, k].copy() for k, f in enumerate(FIELDS)}
        return {f: float(v) for f, v in zip(FIELDS, self._data)}

    def to_dataframe(self) -> Any:
        """Return a pandas DataFrame with one row per curve. Requires pandas."""
        pd = _import_pandas()
        return pd.DataFrame(np.atleast_2d(self._data), columns=list(FIELDS))

    def __array__(self, dtype: Any = None, copy: Any = None) -> np.ndarray:
        a = self._data if dtype is None else self._data.astype(dtype)
        return a.copy() if copy else a

    # ------------------------------------------------------------ fields
    @property
    def is_batch(self) -> bool:
        """True for a batch of curves, False for a single curve."""
        return self._data.ndim == 2

    def _col(self, k: int) -> FloatOrArray:
        return float(self._data[k]) if self._data.ndim == 1 else self._data[:, k]

    @property
    def voc(self) -> FloatOrArray:
        """Open-circuit voltage [V]."""
        return self._col(0)

    @property
    def isc(self) -> FloatOrArray:
        """Short-circuit current [A]."""
        return self._col(1)

    @property
    def vmpp(self) -> FloatOrArray:
        """Maximum power point voltage [V]."""
        return self._col(2)

    @property
    def impp(self) -> FloatOrArray:
        """Maximum power point current [A]."""
        return self._col(3)

    @property
    def vqmp_minus(self) -> FloatOrArray:
        """Voltage of the qmp- point (between short circuit and MPP) [V]."""
        return self._col(4)

    @property
    def iqmp_minus(self) -> FloatOrArray:
        """Current of the qmp- point (between short circuit and MPP) [A]."""
        return self._col(5)

    @property
    def vqmp_plus(self) -> FloatOrArray:
        """Voltage of the qmp+ point (between MPP and open circuit) [V]."""
        return self._col(6)

    @property
    def iqmp_plus(self) -> FloatOrArray:
        """Current of the qmp+ point (between MPP and open circuit) [A]."""
        return self._col(7)

    @property
    def gsc(self) -> FloatOrArray:
        """Short-circuit conductance -dI/dV at V=0 [S]."""
        return self._col(8)

    @property
    def roc(self) -> FloatOrArray:
        """Open-circuit resistance -dV/dI at I=0 [Ohm]."""
        return self._col(9)

    @property
    def pmpp(self) -> FloatOrArray:
        """Maximum power Vmpp * Impp [W]."""
        return self.vmpp * self.impp

    @property
    def ff(self) -> FloatOrArray:
        """Fill factor Pmpp / (Voc * Isc)."""
        return self.pmpp / (self.voc * self.isc)

    # ------------------------------------------------------ container api
    def __len__(self) -> int:
        if not self.is_batch:
            raise TypeError("a single ESP has no len(); check .is_batch first")
        return self._data.shape[0]

    def __iter__(self) -> Iterator["ESP"]:
        if not self.is_batch:
            raise TypeError("a single ESP is not iterable")
        return (ESP._wrap(row) for row in self._data)

    def __getitem__(self, index: Any) -> "ESP":
        if not self.is_batch:
            raise TypeError("a single ESP cannot be indexed; use the field attributes")
        return ESP._wrap(self._data[index])

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("ESP objects are immutable")

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ESP):
            return NotImplemented
        return self._data.shape == other._data.shape and bool(np.array_equal(self._data, other._data))

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        if not self.is_batch:
            body = ", ".join(f"{f}={v:.6g}" for f, v in zip(FIELDS, self._data))
            return f"ESP({body})"
        n = self._data.shape[0]
        with np.printoptions(precision=4, threshold=6, edgeitems=2):
            body = ", ".join(f"{f}={self._data[:, k]}" for k, f in enumerate(FIELDS[:4]))
        return f"<ESP batch n={n}: {body}, ...>"

    # ------------------------------------------------------ conveniences
    def fit(self, model: Any = "PV-Bezier", **options: Any) -> Any:
        """Fit a curve model to these ESP. Shorthand for :func:`pvcraze.fit_model`."""
        from ._analysis import fit_model

        return fit_model(self, model, **options)


def as_esp(x: Any, name: str = "esp") -> ESP:
    """Accept an :class:`ESP` or anything :meth:`ESP.from_array` accepts."""
    if isinstance(x, ESP):
        return x
    try:
        return ESP.from_array(x)
    except (TypeError, ValueError) as e:
        raise type(e)(f"{name}: {e}") from None


def _import_pandas() -> Any:
    try:
        import pandas as pd
    except ImportError:  # pragma: no cover - depends on environment
        raise ImportError(
            "this feature needs pandas; install it with `pip install pvcraze[pandas]`"
        ) from None
    return pd
