"""Cell-type annotation against an scRNA-seq reference.

Replaces ``cell_annotation_20251015_human.py`` and
``cell_annotation_20251110_mouse.py``, which were the same pipeline twice: the two
files differed only in a species string, a log file name, and whether the sample
loop ran on one GPU or seven.  Species is a configuration value here, so there is
one implementation to maintain.

The method follows from the platform's resolution rather than from a user choice:

* **spot-based** platforms pool several cells per barcode, so annotation is a
  *deconvolution* problem -- Tangram maps reference cells onto spots and yields a
  proportion vector per spot in ``obsm["tangram_ct_pred"]``;
* **single-cell** platforms segment individual cells, so it is a *classification*
  problem -- SingleR assigns one label per cell in ``obs["singler_best"]``.

Work is grouped by ``(species, tissue)`` because the expensive preparation --
loading the reference and ranking its marker genes -- is shared by every sample in
a group.  The original scripts recomputed markers per tissue too, but re-globbed
the entire atlas tree once per sample to find its file; here the tree is indexed
once up front.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .._deps import require
from ..config import AnnotationConfig, Config
from ..io import SamplePair, load_sample_sheet, read_h5ad, relative_output, write_h5ad
from ..logging_utils import get_logger
from ..platforms import Resolution, canonical_platform, resolution_of
from ._batch import BatchResult, run_batch

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData

__all__ = [
    "AnnotationTask",
    "annotate_with_singler",
    "annotate_with_tangram",
    "index_samples",
    "pick_device",
    "prepare_reference",
    "reference_markers",
    "run_annotation",
]

logger = get_logger(__name__)


@dataclass(frozen=True, slots=True)
class AnnotationTask:
    """One sample plus the reference context needed to annotate it."""

    pair: SamplePair
    platform: str | None
    resolution: Resolution | None
    reference_path: Path
    markers: tuple[str, ...] | None
    """Marker gene set for Tangram; ``None`` for SingleR, which uses all shared genes."""

    device: str
    """Resolved torch device string, e.g. ``cpu`` or ``cuda:2``."""


# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------


def index_samples(root: Path) -> dict[str, list[Path]]:
    """Index every ``.h5ad`` under ``root`` by file stem.

    The original scripts called ``glob(root/**/<sample>.h5ad)`` once per row of the
    sample sheet, which re-walks the whole atlas for every sample -- quadratic in
    practice and the dominant cost for a sheet with thousands of rows.  One walk
    builds the same mapping.

    Returns:
        Stem -> list of matching paths.  A list rather than a single path because
        duplicate sample IDs across platforms do occur and the caller should be
        told rather than silently given the first hit.
    """
    index: dict[str, list[Path]] = {}
    for path in sorted(root.rglob("*.h5ad")):
        if path.is_file():
            index.setdefault(path.stem, []).append(path)
    logger.info("Indexed %d sample file(s) under %s", sum(len(v) for v in index.values()), root)
    return index


def _resolve_sample(index: dict[str, list[Path]], sample: str) -> Path | None:
    matches = index.get(sample)
    if not matches:
        return None
    if len(matches) > 1:
        logger.warning("Sample %s matches %d files; using %s", sample, len(matches), matches[0])
    return matches[0]


def _find_reference(root: Path, species: str, tissue: str) -> Path:
    """Locate the reference ``.h5ad`` for one ``(species, tissue)`` pair.

    Layout is ``<reference_root>/<species>/<tissue>/<anything>.h5ad``.

    Raises:
        FileNotFoundError: The directory or a reference file is missing.
    """
    directory = root / species / tissue
    if not directory.is_dir():
        raise FileNotFoundError(f"no reference directory for {species}/{tissue}: {directory}")
    candidates = sorted(p for p in directory.glob("*.h5ad") if p.is_file())
    if not candidates:
        raise FileNotFoundError(f"no reference .h5ad in {directory}")
    if len(candidates) > 1:
        logger.warning(
            "%d references in %s; using %s", len(candidates), directory, candidates[0].name
        )
    return candidates[0]


# ---------------------------------------------------------------------------
# Reference preparation
# ---------------------------------------------------------------------------


def prepare_reference(adata: AnnData, config: AnnotationConfig) -> AnnData:
    """Drop reference labels with too few cells.

    Singleton labels have no within-group variance, which makes
    ``rank_genes_groups`` fail outright, and a cell type represented by one cell
    cannot support a proportion estimate anyway.

    Raises:
        KeyError: The configured ``label_key`` is not in ``adata.obs``.
        ValueError: No label survives the filter.
    """
    if config.label_key not in adata.obs.columns:
        raise KeyError(
            f"reference has no obs column '{config.label_key}'; "
            f"available: {sorted(adata.obs.columns)[:20]}"
        )

    counts = adata.obs[config.label_key].value_counts()
    keep = counts[counts >= config.min_cells_per_label].index
    if len(keep) == 0:
        raise ValueError(f"no cell type in the reference has >= {config.min_cells_per_label} cells")
    dropped = len(counts) - len(keep)
    if dropped:
        logger.info(
            "Dropping %d reference label(s) with < %d cells", dropped, config.min_cells_per_label
        )
        adata = adata[adata.obs[config.label_key].isin(keep)].copy()
    logger.info("Reference: %d cells x %d genes, %d label(s)", adata.n_obs, adata.n_vars, len(keep))
    return adata


def reference_markers(adata: AnnData, config: AnnotationConfig) -> tuple[str, ...]:
    """Union of the top marker genes of every reference label.

    Tangram matches spatial spots to reference cells on a shared gene set; using
    every gene both slows the mapping and lets housekeeping genes dominate the
    objective, so the gene set is restricted to per-label markers.
    """
    scanpy = require("scanpy", feature="marker gene ranking")
    import pandas as pd

    scanpy.tl.rank_genes_groups(adata, groupby=config.label_key)
    names = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(config.n_marker_genes)
    markers = tuple(sorted(pd.unique(names.to_numpy().ravel())))
    logger.info("Selected %d marker gene(s) from the reference", len(markers))
    return markers


# ---------------------------------------------------------------------------
# Device placement
# ---------------------------------------------------------------------------


def pick_device(requested: str, path: Path | str, n_gpus: int = 1) -> str:
    """Resolve a torch device string for one sample.

    With several GPUs, samples are assigned by a hash of their path rather than
    round-robin, so a resumed run sends each sample to the same device as before
    and two workers never contend for one GPU's memory in a different pattern.

    Falls back to CPU when CUDA is unavailable, so a config written for a GPU box
    still runs on a laptop.

    >>> pick_device("cpu", "/data/s1.h5ad")
    'cpu'
    """
    if requested == "cpu":
        return "cpu"

    try:
        torch = require("torch", feature="GPU device selection")
    except Exception:  # pragma: no cover - torch missing
        logger.warning("torch is unavailable; falling back to CPU")
        return "cpu"

    if not torch.cuda.is_available():
        if requested.startswith("cuda"):
            logger.warning("CUDA requested but unavailable; falling back to CPU")
        return "cpu"

    if requested.startswith("cuda:"):
        return requested

    available = torch.cuda.device_count()
    usable = max(1, min(n_gpus, available))
    if n_gpus > available:
        logger.warning("n_gpus=%d requested but only %d visible", n_gpus, available)
    digest = hashlib.md5(str(path).encode("utf-8"), usedforsecurity=False).hexdigest()
    return f"cuda:{int(digest, 16) % usable}"


# ---------------------------------------------------------------------------
# Annotation methods
# ---------------------------------------------------------------------------


def annotate_with_tangram(
    spatial: AnnData,
    reference: AnnData,
    config: AnnotationConfig,
    *,
    markers: Sequence[str] | None,
    device: str,
) -> None:
    """Deconvolve spot composition with Tangram, writing proportions to ``obsm``.

    ``spatial`` gains ``obsm["tangram_ct_pred"]``: rows are spots, columns are
    reference cell types, values are estimated proportions.

    Both objects are modified by ``tg.pp_adatas`` (it subsets to shared genes), so
    the caller must pass copies it does not need afterwards.
    """
    tangram = require("tangram", feature="Tangram deconvolution")

    gene_list = list(markers) if markers else None
    tangram.pp_adatas(reference, spatial, genes=gene_list)

    n_shared = len(reference.uns.get("training_genes", ()) or ())
    if n_shared == 0:
        raise ValueError(
            "reference and spatial data share no usable genes -- check that both "
            "use the same gene nomenclature (e.g. symbols vs Ensembl IDs)"
        )
    logger.info("Tangram: %d training gene(s) on %s", n_shared, device)

    mapping = tangram.map_cells_to_space(
        reference,
        spatial,
        mode=config.tangram_mode,
        cluster_label=config.label_key,
        density_prior=config.tangram_density_prior,
        device=device,
        num_epochs=config.tangram_epochs,
    )
    tangram.project_cell_annotations(mapping, spatial, annotation=config.label_key)


def annotate_with_singler(
    spatial: AnnData,
    reference: AnnData,
    config: AnnotationConfig,
) -> None:
    """Classify single cells with SingleR, writing labels to ``obs``.

    Adds one ``obs`` column per SingleR output field, prefixed ``singler_``; the
    assigned label lands in ``obs["singler_best"]``.
    """
    singler = require("singler", feature="SingleR annotation")

    # SingleR expects genes x cells, the transpose of the AnnData convention.
    results = singler.annotate_single(
        test_data=spatial.X.T,
        test_features=list(spatial.var_names),
        ref_data=reference.X.T,
        ref_labels=reference.obs[config.label_key].astype(str),
        ref_features=list(reference.var_names),
        num_threads=config.singler_threads,
    )
    frame = results.to_pandas()
    frame.index = spatial.obs_names
    for column in frame.columns:
        spatial.obs[f"singler_{column}"] = frame[column].astype(str)
    logger.info(
        "SingleR: assigned %d label(s) across %d cell(s)",
        spatial.obs["singler_best"].nunique() if "singler_best" in spatial.obs else 0,
        spatial.n_obs,
    )


# ---------------------------------------------------------------------------
# Per-process reference cache
# ---------------------------------------------------------------------------

# Loading a reference costs seconds to minutes.  With n_jobs>1 each worker process
# needs its own copy, so the cache is keyed by path and lives at module level:
# joblib re-imports this module in the worker, and the second task in that worker
# then reuses the reference.  Pickling the AnnData into every task instead -- what
# the original mouse script did -- serialises the whole reference per sample.
_REFERENCE_CACHE: dict[tuple[str, str, int], Any] = {}


def _load_reference_cached(path: Path, config: AnnotationConfig) -> AnnData:
    key = (str(path), config.label_key, config.min_cells_per_label)
    cached = _REFERENCE_CACHE.get(key)
    if cached is None:
        cached = prepare_reference(read_h5ad(path), config)
        _REFERENCE_CACHE.clear()  # one reference per process is enough; bound memory
        _REFERENCE_CACHE[key] = cached
    return cached


def _annotate_one(task: AnnotationTask, config: AnnotationConfig) -> str | None:
    """Annotate a single sample.  Returns a skip reason, or ``None`` on success."""
    spatial = read_h5ad(task.pair.source)
    if spatial.n_obs == 0:
        return "empty sample (0 observations)"

    reference = _load_reference_cached(task.reference_path, config)

    if task.resolution is Resolution.SPOT:
        # Tangram mutates both objects; give it copies so the cached reference and
        # the source file stay intact.
        annotate_with_tangram(
            spatial,
            reference.copy(),
            config,
            markers=task.markers,
            device=task.device,
        )
    else:
        annotate_with_singler(spatial, reference, config)

    write_h5ad(spatial, task.pair.destination)
    logger.info("Wrote %s", task.pair.destination)
    return None


# ---------------------------------------------------------------------------
# Sample sheet driver
# ---------------------------------------------------------------------------


def _iter_groups(frame: Any, config: AnnotationConfig) -> Iterator[tuple[str, str, Any]]:
    """Yield ``(species, tissue, rows)`` groups from the sample sheet."""
    required = [config.species_column, config.tissue_column, config.sample_column]
    missing = [c for c in required if c not in frame.columns]
    if missing:
        raise KeyError(
            f"sample sheet is missing required column(s) {missing}; "
            f"found {sorted(frame.columns)[:20]}"
        )

    if config.category_column and config.category_filter:
        if config.category_column in frame.columns:
            before = len(frame)
            frame = frame[frame[config.category_column] == config.category_filter]
            logger.info(
                "Sample sheet: %d of %d row(s) match %s == %r",
                len(frame),
                before,
                config.category_column,
                config.category_filter,
            )
        else:
            logger.warning(
                "Sample sheet has no '%s' column; skipping the category filter",
                config.category_column,
            )

    if config.species:
        frame = frame[frame[config.species_column].isin(config.species)]
        logger.info("Restricted to species %s: %d row(s)", config.species, len(frame))

    for (species, tissue), rows in frame.groupby(
        [config.species_column, config.tissue_column], sort=True
    ):
        yield str(species), str(tissue), rows


def run_annotation(config: Config) -> BatchResult:
    """Annotate every sample listed in the sample sheet.

    Reads ``paths.metadata`` to learn each sample's species, tissue and platform;
    resolves inputs under ``paths.raw`` and references under ``paths.reference``;
    writes to ``paths.annotated``.

    A ``(species, tissue)`` group whose reference is missing is recorded as a
    skip for every sample in it rather than aborting the run -- an atlas
    legitimately contains tissues with no matching reference.
    """
    src_root = config.paths.require("raw")
    dst_root = config.paths.require("annotated")
    reference_root = config.paths.require("reference")
    sheet_path = config.paths.require("metadata")
    settings = config.annotation

    frame = load_sample_sheet(sheet_path)
    logger.info("Loaded sample sheet with %d row(s) from %s", len(frame), sheet_path)
    index = index_samples(src_root)

    overall = BatchResult()

    for species, tissue, rows in _iter_groups(frame, settings):
        logger.info("=== %s / %s: %d sample(s) ===", species, tissue, len(rows))

        try:
            reference_path = _find_reference(reference_root, species, tissue)
        except FileNotFoundError as exc:
            logger.warning("Skipping %s/%s: %s", species, tissue, exc)
            _record_group_skip(overall, rows, index, settings, f"no reference: {exc}")
            continue

        tasks = _build_tasks(
            rows,
            index,
            settings,
            src_root,
            dst_root,
            reference_path,
            overall,
            overwrite=config.overwrite,
        )
        if not tasks:
            continue

        # Markers are needed only if some sample in this group is spot-based, and
        # ranking them is expensive, so it happens once per group and only on demand.
        if any(task.resolution is Resolution.SPOT for task in tasks):
            try:
                reference = prepare_reference(read_h5ad(reference_path), settings)
                markers = reference_markers(reference, settings)
                del reference
            except Exception as exc:
                logger.error("Could not rank markers for %s/%s: %s", species, tissue, exc)
                _record_group_skip(overall, rows, index, settings, f"marker ranking failed: {exc}")
                continue
            tasks = [
                task if task.resolution is not Resolution.SPOT else _with_markers(task, markers)
                for task in tasks
            ]

        by_source = {task.pair.source: task for task in tasks}

        def worker(pair: SamplePair, _tasks: dict[Path, AnnotationTask] = by_source) -> str | None:
            return _annotate_one(_tasks[pair.source], settings)

        group_result = run_batch(
            [task.pair for task in tasks],
            worker,
            n_jobs=settings.n_jobs,
            description=f"{species}/{tissue}",
        )
        overall.merge(group_result)

    logger.info("Annotation complete: %s", overall.summary())
    return overall


def _with_markers(task: AnnotationTask, markers: tuple[str, ...]) -> AnnotationTask:
    return AnnotationTask(
        pair=task.pair,
        platform=task.platform,
        resolution=task.resolution,
        reference_path=task.reference_path,
        markers=markers,
        device=task.device,
    )


def _build_tasks(
    rows: Any,
    index: dict[str, list[Path]],
    settings: AnnotationConfig,
    src_root: Path,
    dst_root: Path,
    reference_path: Path,
    overall: BatchResult,
    *,
    overwrite: bool = False,
) -> list[AnnotationTask]:
    """Turn sample-sheet rows into tasks, recording unusable rows as skips."""
    from ._batch import SampleOutcome

    tasks: list[AnnotationTask] = []
    for _, row in rows.iterrows():
        sample = str(row[settings.sample_column])
        source = _resolve_sample(index, sample)
        if source is None:
            logger.warning("No file found for sample %s", sample)
            overall.skipped.append(
                SampleOutcome(Path(sample), "skipped", "no matching .h5ad under paths.raw")
            )
            continue

        destination = relative_output(source, src_root, dst_root)
        if destination.exists() and not overwrite:
            logger.debug("Skipping %s (output exists)", sample)
            continue

        raw_platform = str(row.get(settings.platform_column, "") or "")
        platform = canonical_platform(raw_platform)
        resolution = resolution_of(raw_platform)
        if resolution is None:
            logger.warning(
                "Unknown platform %r for sample %s; treating it as single-cell "
                "resolution (SingleR)",
                raw_platform,
                sample,
            )
            resolution = Resolution.SINGLE_CELL

        tasks.append(
            AnnotationTask(
                pair=SamplePair(source=source, destination=destination),
                platform=platform,
                resolution=resolution,
                reference_path=reference_path,
                markers=None,
                device=pick_device(settings.device, source, settings.n_gpus),
            )
        )
    return tasks


def _record_group_skip(
    result: BatchResult,
    rows: Any,
    index: dict[str, list[Path]],
    settings: AnnotationConfig,
    reason: str,
) -> None:
    """Mark every sample in a group as skipped, so counts stay honest."""
    from ._batch import SampleOutcome

    for _, row in rows.iterrows():
        sample = str(row[settings.sample_column])
        source = _resolve_sample(index, sample) or Path(sample)
        result.skipped.append(SampleOutcome(source, "skipped", reason))
