"""Logging setup shared by the CLI and the pipelines.

The original scripts communicated through ``print`` plus hand-rolled
``open(log, "a").write(...)`` calls, which meant progress and errors were
interleaved on stdout with no timestamps and no way to turn detail down.  Here a
single :func:`configure` call installs a console handler and, optionally, a file
handler; pipeline code just calls ``logging.getLogger(__name__)``.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

__all__ = ["configure", "get_logger"]

_FORMAT = "%(asctime)s %(levelname)-8s %(name)s: %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

_configured = False


def configure(
    level: int | str = logging.INFO,
    log_file: str | Path | None = None,
    *,
    force: bool = False,
) -> None:
    """Install handlers on the ``stcompass`` logger.

    Args:
        level: Threshold for the console handler.
        log_file: Optional path that receives every record at ``DEBUG`` and
            above, so a failed batch run can be diagnosed after the fact.
        force: Re-configure even if :func:`configure` already ran.  Without this
            repeated calls are ignored, which keeps duplicate handlers (and
            duplicated lines) out of long-running processes.
    """
    global _configured
    logger = logging.getLogger("stcompass")
    if _configured and not force:
        return
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    logger.setLevel(logging.DEBUG)
    # Records are filtered per handler, not on the logger, so that a DEBUG log
    # file can coexist with an INFO console.
    logger.propagate = False

    formatter = logging.Formatter(_FORMAT, datefmt=_DATE_FORMAT)

    console = logging.StreamHandler(stream=sys.stderr)
    console.setLevel(level)
    console.setFormatter(formatter)
    logger.addHandler(console)

    if log_file is not None:
        path = Path(log_file).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(path, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    _configured = True


def get_logger(name: str) -> logging.Logger:
    """Return a logger under the ``stcompass`` namespace.

    Accepts a dunder-name (``__name__``) and keeps it as-is when it already sits
    in the package namespace, so log lines carry the originating module.
    """
    if name == "stcompass" or name.startswith("stcompass."):
        return logging.getLogger(name)
    return logging.getLogger(f"stcompass.{name}")
