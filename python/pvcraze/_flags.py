"""Curve analysis flags."""

from __future__ import annotations

import enum
from typing import Any, Dict, List, Union

import numpy as np

from . import _core


class _FlagBase(enum.IntFlag):
    @property
    def messages(self) -> List[str]:
        """Human readable descriptions of the set flags, errors first."""
        codes = sorted((c for c in _MESSAGES if self & c), key=lambda c: (not c & ERROR_MASK, c))
        return [_MESSAGES[c] for c in codes]

    @property
    def has_error(self) -> bool:
        """True if any error flag is set, i.e. the ESP are not meaningful."""
        return bool(self & ERROR_MASK)

    @property
    def has_warning(self) -> bool:
        """True if any warning flag is set (results may be atypical or inaccurate)."""
        return bool(self & WARNING_MASK)

    @property
    def errors(self) -> "Flag":
        """The subset of set flags that are errors."""
        return Flag(self & ERROR_MASK)

    @property
    def warnings(self) -> "Flag":
        """The subset of set flags that are warnings."""
        return Flag(self & WARNING_MASK)


_TABLE = _core.flag_table()  # (code, name, message, is_error)
_MESSAGES: Dict[int, str] = {code: msg for code, _, msg, _ in _TABLE}

#: Bit mask of all error flags.
ERROR_MASK: int = sum(code for code, _, _, err in _TABLE if err)
#: Bit mask of all warning flags.
WARNING_MASK: int = sum(code for code, _, _, err in _TABLE if not err)

Flag = _FlagBase("Flag", [(name, code) for code, name, _, _ in _TABLE], module=__name__)
Flag.__doc__ = """Flags raised by the curve analysis, as an :class:`enum.IntFlag`.

A value of ``0`` (``Flag(0)``) means the curve was analysed without problems.
Flags are either *warnings* (results are atypical or inaccurate) or *errors*
(no meaningful parameters could be determined)::

    f = result.flag(3)        # Flag for curve 3
    if f.has_error:
        print(f.messages)
    Flag.NOMPP in f           # test a specific flag

Members, their messages and whether they are errors:

""" + "\n".join(
    f"* ``{name}`` ({'error' if err else 'warning'}): {msg}" for _, name, msg, err in _TABLE
)


def describe(flags: Union[int, "Flag"]) -> List[str]:
    """Return the messages for a raw flag code (see :attr:`Flag.messages`)."""
    return Flag(int(flags)).messages


def error_mask(flags: Any) -> np.ndarray:
    """Vectorised ``has_error`` for an array of raw flag codes."""
    return (np.asarray(flags, dtype=np.int64) & ERROR_MASK) != 0


def warning_mask(flags: Any) -> np.ndarray:
    """Vectorised ``has_warning`` for an array of raw flag codes."""
    return (np.asarray(flags, dtype=np.int64) & WARNING_MASK) != 0
