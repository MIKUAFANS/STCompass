"""Lazy imports for optional dependencies.

The pipelines in this package sit on top of a heavy scientific stack, and most
users only need part of it: a lab doing imaging-based assays never installs
Tangram, and a CPU-only machine cannot install a CUDA build of torch.  Importing
those packages at module import time would make ``import stcompass`` fail for
everyone, so every optional dependency is resolved at call time through
:func:`require`, which reports the exact extra to install.
"""

from __future__ import annotations

import importlib
from types import ModuleType

__all__ = ["MissingDependencyError", "require"]


class MissingDependencyError(ImportError):
    """An optional dependency is needed for the requested feature."""


# Import name -> command that installs it.
_INSTALL_HINTS: dict[str, str] = {
    "scanpy": "pip install 'stcompass'",
    "tangram": "pip install 'stcompass[tangram]'",
    "torch": "pip install 'stcompass[tangram]'",
    "singler": "pip install 'stcompass[singler]'",
    "watchdog": "pip install 'stcompass[watch]'",
    "openpyxl": "pip install 'stcompass[excel]'",
}


def require(module: str, *, feature: str) -> ModuleType:
    """Import ``module`` or raise a :class:`MissingDependencyError` explaining why.

    Parameters
    ----------
    module
        Import name, e.g. ``"tangram"`` (not the distribution name ``tangram-sc``).
    feature
        Human-readable description of what needed the module, used in the error
        message so the traceback points at the pipeline rather than at an import.
    """
    try:
        return importlib.import_module(module)
    except ImportError as exc:  # pragma: no cover - depends on the environment
        root = module.split(".")[0]
        hint = _INSTALL_HINTS.get(root, f"pip install {root}")
        raise MissingDependencyError(
            f"{feature} requires the '{root}' package, which is not installed.\n"
            f"Install it with:  {hint}"
        ) from exc
