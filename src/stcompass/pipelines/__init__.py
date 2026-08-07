"""The five pipeline stages.

Each module exposes a ``run_*(config)`` entry point that processes a whole tree of
samples and returns a :class:`~stcompass.pipelines._batch.BatchResult`, plus
smaller functions that operate on a single in-memory ``AnnData`` so the numerical
behaviour can be tested without touching the filesystem.

Stages are chained by path: ``qc`` reads ``paths.raw`` and writes ``paths.qc``,
which ``cluster`` and ``programs`` then read, and so on.  Nothing is implicit --
each stage is a separate command and can be run, re-run or skipped on its own.
"""

from __future__ import annotations

from ._batch import BatchResult, SampleOutcome, run_batch
from .annotation import run_annotation
from .clustering import run_clustering
from .programs import run_programs
from .qc import run_qc
from .visualize import run_plot

__all__ = [
    "BatchResult",
    "SampleOutcome",
    "run_annotation",
    "run_batch",
    "run_clustering",
    "run_plot",
    "run_programs",
    "run_qc",
]
