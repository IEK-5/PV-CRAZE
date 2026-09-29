"""Input normalisation helpers shared by the public modules."""

from __future__ import annotations

from typing import Any, List, Tuple

import numpy as np


def as_float(x: Any, name: str) -> np.ndarray:
    """Return ``x`` as a float64 array, rejecting non-numeric input.

    Integer input is converted; complex, boolean, string and object input
    raises ``TypeError`` instead of being silently reinterpreted.
    """
    a = np.asarray(x)
    if a.dtype.kind not in "iuf":
        raise TypeError(f"{name} must be real numbers, got dtype {a.dtype}")
    return np.asarray(a, dtype=np.float64, order="C")  # keeps 0-d input 0-d


def as_curve(x: Any, name: str = "curve") -> np.ndarray:
    """Return a single IV curve as a C-contiguous ``(n, 2)`` float64 array."""
    a = as_float(x, name)
    if a.ndim != 2 or a.shape[1] != 2:
        raise ValueError(
            f"{name} must have shape (n, 2) with columns (V, I), got shape {a.shape}"
        )
    return a


def _try_float(x: Any) -> np.ndarray | None:
    """``as_float`` for rectangular input, ``None`` for ragged nested sequences."""
    if isinstance(x, np.ndarray):
        return as_float(x, "curves")
    try:
        a = np.asarray(x)
    except ValueError:  # ragged nested sequence
        return None
    if a.dtype == object:
        return None
    return as_float(a, "curves")


def normalize_curves(curves: Any, current: Any = None) -> Tuple[List[np.ndarray], bool]:
    """Turn the accepted curve input forms into a list of ``(n, 2)`` arrays.

    Accepted forms:

    * ``curves`` is one ``(n, 2)`` array                      -> single
    * ``curves`` is a sequence of ``(n, 2)`` arrays           -> batch
    * ``curves`` is V and ``current`` is I, both 1-D          -> single
    * ``curves`` and ``current`` are equal-length sequences
      of 1-D V and I arrays                                     -> batch

    Returns the list and whether the input was a single curve.
    """
    if isinstance(curves, (str, bytes)) or np.isscalar(curves):
        raise TypeError(
            "curves must be an (n, 2) array, a sequence of (n, 2) arrays, or V and I arrays"
        )
    if current is not None:
        v, i = _try_float(curves), _try_float(current)
        if v is not None and v.ndim == 1:
            if i is None or i.shape != v.shape:
                raise ValueError("V and I must be 1-D arrays of the same length")
            return [np.column_stack([v, i])], True
        vs, is_ = list(curves), list(current)
        if len(vs) != len(is_):
            raise ValueError(f"got {len(vs)} voltage arrays but {len(is_)} current arrays")
        out = []
        for k, (vk, ik) in enumerate(zip(vs, is_)):
            vv, ii = as_float(vk, f"V[{k}]"), as_float(ik, f"I[{k}]")
            if vv.ndim != 1 or vv.shape != ii.shape:
                raise ValueError(f"V[{k}] and I[{k}] must be 1-D arrays of equal length")
            out.append(np.column_stack([vv, ii]))
        return out, False

    a = _try_float(curves)
    if a is not None:
        if a.ndim == 1 and a.size == 0:
            return [], False
        if a.ndim == 2:
            if a.shape[0] == 2 and a.shape[1] != 2:
                raise ValueError(
                    f"curves has shape {a.shape}, i.e. rows (V, I); pass analyze(V, I) "
                    "or an (n, 2) array with columns (V, I)"
                )
            return [as_curve(a)], True
        if a.ndim == 3:
            return [as_curve(c, f"curves[{k}]") for k, c in enumerate(a)], False
        raise ValueError(
            f"curves must be an (n, 2) array or a sequence of (n, 2) arrays, got shape {a.shape}"
        )
    return [as_curve(c, f"curves[{k}]") for k, c in enumerate(curves)], False


def broadcast_rows(n: int, name: str, value: Any) -> np.ndarray:
    """Broadcast a scalar or length-``n`` vector to a float64 vector of length ``n``."""
    a = as_float(value, name)
    if a.ndim > 1:
        raise ValueError(f"{name} must be a scalar or 1-D array, got shape {a.shape}")
    try:
        return np.ascontiguousarray(np.broadcast_to(a, (n,)))
    except ValueError:
        raise ValueError(f"{name} has length {a.shape[0]} but {n} values are needed") from None
