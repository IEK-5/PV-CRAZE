"""Control of the C library's console output."""

from __future__ import annotations

import contextlib
from typing import Dict, Iterator

from . import _core

_BITS = {"timing": 1, "report": 2, "progress": 4, "header": 8}


def get_verbosity() -> Dict[str, bool]:
    """Return which kinds of console output are enabled."""
    mask = _core.get_verbosity()
    return {name: bool(mask & bit) for name, bit in _BITS.items()}


def set_verbosity(
    *, report: bool = False, progress: bool = False, timing: bool = False, header: bool = False
) -> Dict[str, bool]:
    """Enable or disable console output of the C library, returning the previous state.

    The library is silent by default. Output is written by C code directly to
    the process's standard output, so it is not captured by
    ``contextlib.redirect_stdout`` and in Jupyter it appears in the terminal
    running the kernel.

    Parameters
    ----------
    report:
        Print a summary table of the analysis flags after :func:`pvcraze.analyze`.
    progress:
        Print a progress bar.
    timing:
        Print timing information.
    header:
        Print a banner naming the routine and the number of curves.
    """
    previous = get_verbosity()
    state = {"report": report, "progress": progress, "timing": timing, "header": header}
    _core.set_verbosity(sum(bit for name, bit in _BITS.items() if state[name]))
    return previous


@contextlib.contextmanager
def verbose(
    *, report: bool = True, progress: bool = True, timing: bool = True, header: bool = True
) -> Iterator[None]:
    """Context manager that enables console output (all kinds by default)::

        with pvcraze.verbose(progress=False):
            result = pvcraze.analyze(curves)   # prints header, timing and report
    """
    previous = set_verbosity(report=report, progress=progress, timing=timing, header=header)
    try:
        yield
    finally:
        set_verbosity(**previous)
