"""Rendering of annotation results to figures.

Replaces ``celltype_image1015.py``.  Which figure to draw follows from what the
annotation stage produced:

* ``obsm["tangram_ct_pred"]`` holds per-spot *proportions*, so each spot is drawn
  as a pie chart whose wedges are the dominant cell types -- the only way to show
  a mixture at a location;
* ``obs["singler_best"]`` holds one *label* per cell, so the sample is drawn as a
  categorical scatter over its spatial coordinates.

The original script mixed this rendering with a bespoke watcher: two threads, a
hand-rolled dedupe dictionary, a signal handler, and a CSV of processed files.
That mode is kept because it is genuinely useful -- figures appear while a
multi-day annotation run is still going -- but it is now opt-in behind
``plot.watch`` and shares the file-stability check with the batch path.
"""

from __future__ import annotations

import time
from itertools import chain
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..config import Config, PlotConfig
from ..io import SamplePair, iter_samples, read_h5ad, relative_output
from ..logging_utils import get_logger
from ._batch import BatchResult, run_batch

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData

__all__ = [
    "estimate_radius",
    "plot_sample",
    "run_plot",
    "top_n_proportions",
    "wait_until_stable",
]

logger = get_logger(__name__)


def _import_pyplot():
    """Import pyplot with a non-interactive backend.

    ``Agg`` is selected before pyplot is first imported: these pipelines run over
    SSH and in batch schedulers where no display exists, and an interactive
    backend fails at import time there.
    """
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    return plt


def wait_until_stable(
    path: Path,
    *,
    checks: int = 3,
    interval: float = 2.0,
    min_size: int = 1024,
    max_wait: float | None = None,
) -> bool:
    """Wait until ``path``'s size stops changing.

    A file being written by another process is readable but truncated, and an
    ``.h5ad`` read from a half-written file fails in ways that look like data
    corruption.  Requiring ``checks`` consecutive identical sizes is a cheap proxy
    for "the writer has finished".

    Args:
        path: File to watch.
        checks: Consecutive equal-size observations required.
        interval: Seconds between observations.
        min_size: Sizes below this are never considered stable, which rejects a
            freshly created empty file.
        max_wait: Give up after this many seconds; defaults to
            ``checks * interval * 8``, so a stalled copy cannot block the queue
            forever.

    Returns:
        ``True`` when the file settled, ``False`` on timeout or disappearance.
    """
    if not path.exists():
        return False
    budget = max_wait if max_wait is not None else checks * interval * 8
    deadline = time.monotonic() + budget

    previous = -1
    stable = 0
    while time.monotonic() < deadline:
        try:
            size = path.stat().st_size
        except FileNotFoundError:
            return False
        if size >= min_size and size == previous:
            stable += 1
            if stable >= checks:
                return True
        else:
            stable = 0
        previous = size
        time.sleep(interval)

    logger.warning("Gave up waiting for %s to stabilise after %.0fs", path.name, budget)
    return False


def top_n_proportions(frame, top_n: int):
    """Zero out all but the ``top_n`` largest values per row, then renormalise.

    A spot's proportion vector has one entry per reference cell type -- often
    thirty or more.  Drawing every entry turns each pie into a ring of slivers, so
    only the dominant types are kept.  Renormalising afterwards keeps each pie a
    full circle, so wedge angles remain comparable between spots.
    """
    import pandas as pd

    values = frame.to_numpy(dtype=np.float64, copy=True)
    n_types = values.shape[1]
    keep = min(top_n, n_types)

    if keep < n_types:
        # argpartition puts the `keep` largest values at the end of each row in
        # O(n) per row, which matters at hundreds of thousands of spots.
        cut = np.partition(values, n_types - keep, axis=1)[:, n_types - keep]
        values[values < cut[:, None]] = 0.0

    totals = values.sum(axis=1, keepdims=True)
    np.divide(values, totals, out=values, where=totals > 0)
    return pd.DataFrame(values, index=frame.index, columns=frame.columns)


def estimate_radius(coordinates: np.ndarray) -> float:
    """Pick a pie radius from nearest-neighbour spacing.

    Spot pitch differs by three orders of magnitude across platforms (55 µm Visium
    spots, 2 µm HD bins, arbitrary units after registration), so a fixed radius is
    wrong everywhere.  Half the median nearest-neighbour distance leaves adjacent
    pies just touching; 0.45 of it leaves a small gap so boundaries stay visible.
    """
    from scipy.spatial import cKDTree

    if len(coordinates) < 2:
        return 1.0
    # k=2 because the nearest neighbour of a point is the point itself.
    distances = cKDTree(coordinates).query(coordinates, k=2)[0][:, 1]
    finite = distances[np.isfinite(distances) & (distances > 0)]
    if finite.size == 0:
        return 1.0
    return float(np.median(finite) * 0.45)


def _spatial_coordinates(adata: AnnData, preferred: str) -> tuple[np.ndarray, str] | None:
    """Find 2-D coordinates in ``obsm``, trying the configured key first."""
    for key in (preferred, "spatial", "X_spatial", "X_umap", "umap"):
        if key in adata.obsm:
            array = np.asarray(adata.obsm[key])
            if array.ndim == 2 and array.shape[1] >= 2:
                return array[:, :2].astype(np.float64, copy=False), key
    return None


def _palette(n_colors: int) -> list:
    """Build a categorical palette of ``n_colors`` distinct colours.

    Chains the three 20-colour tab palettes for up to 60 types, then cycles.  A
    continuous colormap would be wrong here: adjacent cell types are not adjacent
    in any meaningful ordering.
    """
    plt = _import_pyplot()
    maps = [plt.get_cmap("tab20"), plt.get_cmap("tab20b"), plt.get_cmap("tab20c")]
    colors = list(chain.from_iterable([cmap(i) for i in range(cmap.N)] for cmap in maps))
    if n_colors <= len(colors):
        return colors[:n_colors]
    repeats = (n_colors // len(colors)) + 1
    return (colors * repeats)[:n_colors]


def _plot_proportions(adata: AnnData, destination: Path, config: PlotConfig) -> str:
    """Draw a scatter-pie of Tangram proportions."""
    import pandas as pd
    from matplotlib.patches import Wedge

    plt = _import_pyplot()

    raw = adata.obsm[config.proportions_key]
    if isinstance(raw, pd.DataFrame):
        frame = raw.copy()
    else:
        frame = pd.DataFrame(raw, index=adata.obs_names)
    if frame.shape[1] == 0:
        return "proportions matrix has no columns"

    found = _spatial_coordinates(adata, config.spatial_key)
    if found is None:
        return "no 2-D coordinates in obsm"
    coordinates, key_used = found
    if len(coordinates) != frame.shape[0]:
        return f"coordinate count {len(coordinates)} does not match {frame.shape[0]} spots"

    proportions = top_n_proportions(frame, config.top_n_types)
    labels = [str(c) for c in proportions.columns]
    colors = _palette(len(labels))

    figure, axes = plt.subplots(figsize=config.figsize)
    try:
        axes.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            s=1,
            color="lightgrey",
            alpha=0.5,
            zorder=0,
        )

        if len(coordinates) > config.max_pies:
            # Hundreds of thousands of wedges exhaust memory in the PDF backend
            # and would be individually invisible anyway; show the dominant type
            # per spot instead of a pie.
            logger.warning(
                "%d spots exceeds plot.max_pies=%d; drawing dominant cell type only",
                len(coordinates),
                config.max_pies,
            )
            dominant = np.asarray(proportions.to_numpy()).argmax(axis=1)
            axes.scatter(
                coordinates[:, 0],
                coordinates[:, 1],
                s=2,
                c=[colors[i] for i in dominant],
                linewidths=0,
            )
        else:
            radius = config.radius if config.radius is not None else estimate_radius(coordinates)
            values = proportions.to_numpy(dtype=np.float64)
            for (x, y), row in zip(coordinates, values, strict=True):
                total = row.sum()
                if total <= 0:
                    continue
                start = 0.0
                for fraction, color in zip(row / total, colors, strict=False):
                    if fraction <= 0:
                        continue
                    end = start + fraction
                    axes.add_patch(
                        Wedge(
                            (float(x), float(y)),
                            radius,
                            start * 360.0,
                            end * 360.0,
                            facecolor=color,
                            edgecolor="none",
                        )
                    )
                    start = end

        axes.set_aspect("equal")
        if config.invert_y:
            axes.invert_yaxis()
        axes.axis("off")
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c) for c in colors]
        axes.legend(
            handles,
            labels,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=6,
            frameon=False,
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(destination, dpi=config.dpi, bbox_inches="tight")
    finally:
        # Always close: pyplot keeps a global reference, so a leaked figure in a
        # batch of thousands exhausts memory.
        plt.close(figure)

    logger.debug("Rendered scatter-pie using obsm['%s']", key_used)
    return ""


def _plot_labels(adata: AnnData, destination: Path, config: PlotConfig) -> str:
    """Draw a categorical scatter of SingleR labels."""
    scanpy = None
    try:
        from .._deps import require

        scanpy = require("scanpy", feature="plotting cell-type labels")
    except Exception:  # pragma: no cover - handled by falling back to matplotlib
        logger.debug("scanpy unavailable; using the matplotlib fallback")

    found = _spatial_coordinates(adata, config.spatial_key)
    if found is None:
        return "no 2-D coordinates in obsm"
    coordinates, key_used = found

    plt = _import_pyplot()

    if scanpy is not None:
        basis = key_used[2:] if key_used.startswith("X_") else key_used
        scanpy.pl.embedding(
            adata,
            basis=basis,
            color=config.label_key,
            legend_loc="right margin",
            frameon=False,
            title="",
            show=False,
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(destination, dpi=config.dpi, bbox_inches="tight")
        plt.close("all")
        return ""

    categories = adata.obs[config.label_key].astype(str)
    unique = sorted(categories.unique())
    colors = _palette(len(unique))
    lookup = {name: colors[i] for i, name in enumerate(unique)}

    figure, axes = plt.subplots(figsize=config.figsize)
    try:
        axes.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            s=2,
            c=[lookup[v] for v in categories],
            linewidths=0,
        )
        axes.set_aspect("equal")
        if config.invert_y:
            axes.invert_yaxis()
        axes.axis("off")
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=lookup[n]) for n in unique]
        axes.legend(
            handles, unique, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6, frameon=False
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(destination, dpi=config.dpi, bbox_inches="tight")
    finally:
        plt.close(figure)
    return ""


def plot_sample(source: Path, destination: Path, config: PlotConfig) -> str | None:
    """Render one annotated sample.

    Returns:
        ``None`` on success, or a reason string when the sample holds nothing
        plottable.
    """
    adata = read_h5ad(source)

    if config.label_key in adata.obs.columns:
        problem = _plot_labels(adata, destination, config)
    elif config.proportions_key in adata.obsm:
        problem = _plot_proportions(adata, destination, config)
    else:
        return (
            f"no obs['{config.label_key}'] and no obsm['{config.proportions_key}'] "
            f"-- run 'stcompass annotate' first"
        )

    if problem:
        return problem
    logger.info("Wrote %s", destination)
    return None


def run_plot(config: Config, *, exclude: tuple[str, ...] = ()) -> BatchResult:
    """Render every annotated sample under ``paths.annotated`` into ``paths.figures``.

    With ``plot.watch`` enabled the function does not return after the initial
    sweep: it keeps polling for new or modified files and renders them as they
    appear.  Polling is used rather than filesystem events because the atlas lives
    on an NFS mount, where inotify does not fire for writes made on another host.
    """
    src_root = config.paths.require("annotated")
    dst_root = config.paths.require("figures")
    settings = config.plot
    suffix = f".{settings.format.lstrip('.')}"

    def worker(pair: SamplePair) -> str | None:
        if not wait_until_stable(
            pair.source,
            checks=settings.stability_checks,
            interval=settings.stability_interval,
        ):
            return "file did not stabilise (still being written?)"
        return plot_sample(pair.source, pair.destination, settings)

    pairs = iter_samples(
        src_root, dst_root, suffix=suffix, overwrite=config.overwrite, exclude=exclude
    )
    result = run_batch(pairs, worker, n_jobs=config.n_jobs, description="figures")

    if settings.watch:
        _watch_loop(src_root, dst_root, suffix, settings, worker, result, exclude)
    return result


def _watch_loop(
    src_root: Path,
    dst_root: Path,
    suffix: str,
    settings: PlotConfig,
    worker,
    result: BatchResult,
    exclude: tuple[str, ...],
) -> None:
    """Poll for new samples until interrupted.

    Single-threaded on purpose.  The original used a producer thread, a consumer
    thread, a bounded queue that dropped half its contents under load, and a
    debounce dictionary -- machinery whose only job was to avoid reading a file
    twice.  Tracking modification times in one loop achieves that, and rendering
    is I/O-bound on reading the ``.h5ad`` anyway.
    """
    logger.info("Watching %s every %.1fs (Ctrl-C to stop)", src_root, settings.poll_seconds)
    seen: dict[Path, float] = {}
    for path in src_root.rglob("*.h5ad"):
        try:
            seen[path] = path.stat().st_mtime
        except OSError:
            continue

    try:
        while True:
            time.sleep(settings.poll_seconds)
            for source in sorted(src_root.rglob("*.h5ad")):
                if not source.is_file() or source.name in exclude:
                    continue
                try:
                    mtime = source.stat().st_mtime
                except OSError:
                    continue
                if seen.get(source) == mtime:
                    continue
                seen[source] = mtime

                destination = relative_output(source, src_root, dst_root, suffix)
                if destination.exists():
                    continue
                pair = SamplePair(source=source, destination=destination)
                outcome = worker(pair)
                if outcome:
                    logger.info("Skipped %s: %s", source.name, outcome)
                else:
                    result.processed.append(source)
    except KeyboardInterrupt:
        logger.info("Watch stopped: %s", result.summary())
