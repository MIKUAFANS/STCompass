"""Shared driver for the per-sample batch loop.

Every stage does the same bookkeeping: walk the input tree, skip samples that
already have an output, process each one, and keep going when a sample fails.
The original scripts each re-implemented that loop with a bare ``except`` around
the body and a hand-written append to a log file, so a failure was recorded but
never counted, and there was no way to tell a clean run from one where half the
atlas errored.

:func:`run_batch` centralises it and returns a :class:`BatchResult`, which the CLI
turns into an exit code -- a batch where every sample failed should not report
success.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..io import SamplePair
from ..logging_utils import get_logger

__all__ = ["BatchResult", "SampleOutcome", "run_batch"]

logger = get_logger(__name__)


@dataclass(frozen=True, slots=True)
class SampleOutcome:
    """What happened to one sample."""

    sample: Path
    status: str
    """``"processed"``, ``"skipped"`` or ``"failed"``."""

    detail: str = ""
    """Reason for a skip, or the exception text for a failure."""


@dataclass(slots=True)
class BatchResult:
    """Aggregate outcome of a batch run."""

    processed: list[Path] = field(default_factory=list)
    skipped: list[SampleOutcome] = field(default_factory=list)
    failed: list[SampleOutcome] = field(default_factory=list)

    @property
    def total(self) -> int:
        return len(self.processed) + len(self.skipped) + len(self.failed)

    @property
    def ok(self) -> bool:
        """True when nothing failed.

        A run with zero samples counts as OK: an empty input tree is a
        configuration question, reported separately, not a processing error.
        """
        return not self.failed

    def summary(self) -> str:
        """One-line human-readable summary for the end of a run."""
        return (
            f"{len(self.processed)} processed, "
            f"{len(self.skipped)} skipped, "
            f"{len(self.failed)} failed "
            f"(of {self.total})"
        )

    def merge(self, other: BatchResult) -> None:
        """Fold ``other`` into this result, for stages that run in several groups."""
        self.processed.extend(other.processed)
        self.skipped.extend(other.skipped)
        self.failed.extend(other.failed)


# A worker returns None when it processed the sample, or a reason string when it
# decided to skip (too few cells, no spatial coordinates, ...).  Raising means
# failure.
Worker = Callable[[SamplePair], str | None]


def run_batch(
    pairs: Iterable[SamplePair],
    worker: Worker,
    *,
    n_jobs: int = 1,
    description: str = "samples",
    progress: bool = True,
) -> BatchResult:
    """Apply ``worker`` to every pair, collecting outcomes instead of aborting.

    Args:
        pairs: Input/output pairs, typically from :func:`stcompass.io.iter_samples`.
        worker: Called once per pair. Returns ``None`` on success or a short
            reason string to record a skip; exceptions are caught and recorded.
        n_jobs: Worker processes. ``1`` runs in-process, which keeps tracebacks
            readable and is the right choice for GPU stages that manage their own
            device placement.
        description: Noun used in log messages.
        progress: Show a ``tqdm`` progress bar when available.

    Returns:
        A :class:`BatchResult`. Never raises for a per-sample failure; genuine
        programming errors in the driver itself still propagate.
    """
    items: Sequence[SamplePair] = list(pairs)
    result = BatchResult()
    if not items:
        logger.warning("No %s to process -- check the input path and filters", description)
        return result

    logger.info("Processing %d %s with n_jobs=%d", len(items), description, n_jobs)

    if n_jobs == 1:
        iterator: Iterable[SamplePair] = items
        if progress:
            iterator = _maybe_progress(items, description)
        for pair in iterator:
            _record(result, pair, _safe_call(worker, pair))
        logger.info("Finished: %s", result.summary())
        return result

    from joblib import Parallel, delayed

    # ``loky`` (joblib's default) re-imports the module in a fresh interpreter,
    # which is what isolates a segfaulting native library to one sample.
    outcomes = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(_safe_call)(worker, pair) for pair in items
    )
    for pair, outcome in zip(items, outcomes, strict=True):
        _record(result, pair, outcome)
    logger.info("Finished: %s", result.summary())
    return result


def _maybe_progress(items: Sequence[SamplePair], description: str) -> Iterable[SamplePair]:
    """Wrap ``items`` in a progress bar, degrading gracefully without tqdm."""
    try:
        from tqdm import tqdm
    except ImportError:  # pragma: no cover - tqdm is a hard dependency
        return items
    return tqdm(items, desc=description, unit="sample")


def _safe_call(worker: Worker, pair: SamplePair) -> tuple[str, str]:
    """Run ``worker``, converting any exception into a ``("failed", message)`` pair.

    Returning a tuple rather than raising keeps the parallel branch simple: joblib
    would otherwise cancel the remaining tasks on the first exception.
    """
    try:
        reason = worker(pair)
    except KeyboardInterrupt:  # pragma: no cover - interactive
        raise
    except Exception as exc:
        logger.error("Failed: %s -- %s: %s", pair.source.name, type(exc).__name__, exc)
        logger.debug("Traceback for %s", pair.source.name, exc_info=True)
        return ("failed", f"{type(exc).__name__}: {exc}")
    if reason:
        return ("skipped", reason)
    return ("processed", "")


def _record(result: BatchResult, pair: SamplePair, outcome: tuple[str, str]) -> None:
    status, detail = outcome
    if status == "processed":
        result.processed.append(pair.source)
    elif status == "skipped":
        logger.info("Skipped %s: %s", pair.source.name, detail)
        result.skipped.append(SampleOutcome(pair.source, "skipped", detail))
    else:
        result.failed.append(SampleOutcome(pair.source, "failed", detail))
