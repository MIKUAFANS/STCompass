"""Quality control, normalisation, embedding and clustering.

Replaces ``QC_ALL0917.py``. The stage decides its thresholds from the platform's
:class:`~stcompass.platforms.Resolution`, because the same cutoff cannot serve
both kinds of assay: a Visium barcode pools several cells and yields thousands of
counts, while a MERFISH cell is measured on a 300-gene panel and may legitimately
carry twenty.

Two behaviours from the original script are preserved deliberately:

* the raw matrix is copied to ``layers["counts"]`` before normalisation, so the
  gene-programs stage can factorise counts after ``X`` has been log-transformed;
* whether ``X`` already holds log-normalised values is *inferred* rather than
  assumed (see :func:`stcompass.matrix.looks_like_counts`), because public
  atlases mix both and carry no reliable flag.

The per-barcode and per-gene filters were commented out in the original script,
so the published atlas kept every barcode.  They are implemented here and
enabled by default; set ``qc.filter_cells``/``qc.filter_genes`` to ``false`` to
reproduce the original output exactly.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .._deps import require
from ..config import Config, QCConfig
from ..io import SamplePair, iter_samples, read_h5ad, write_h5ad
from ..logging_utils import get_logger
from ..matrix import count_nonzero_per_column, count_nonzero_per_row, looks_like_counts, row_sums
from ..platforms import Resolution, infer_platform, resolution_of
from ._batch import BatchResult, run_batch

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData

__all__ = ["QCThresholds", "apply_qc", "qc_sample", "run_qc", "thresholds_for"]

logger = get_logger(__name__)


class QCThresholds:
    """Resolved numeric thresholds for one sample.

    A small holder rather than a dict so the pipeline cannot silently read a
    misspelled key.
    """

    __slots__ = ("min_cells_per_gene", "min_counts", "min_genes", "min_units_after_qc")

    def __init__(
        self,
        min_counts: int,
        min_genes: int,
        min_cells_per_gene: int,
        min_units_after_qc: int,
    ) -> None:
        self.min_counts = min_counts
        self.min_genes = min_genes
        self.min_cells_per_gene = min_cells_per_gene
        self.min_units_after_qc = min_units_after_qc

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"QCThresholds(min_counts={self.min_counts}, min_genes={self.min_genes}, "
            f"min_cells_per_gene={self.min_cells_per_gene}, "
            f"min_units_after_qc={self.min_units_after_qc})"
        )


def thresholds_for(resolution: Resolution | None, config: QCConfig) -> QCThresholds:
    """Pick thresholds for a platform resolution.

    Unknown platforms fall back to the permissive single-cell numbers: discarding
    real data is worse than keeping marginal barcodes, which later stages can
    still filter.

    >>> thresholds_for(Resolution.SPOT, QCConfig()).min_counts
    100
    >>> thresholds_for(None, QCConfig()).min_counts
    20
    """
    if resolution is Resolution.SPOT:
        return QCThresholds(
            min_counts=config.min_counts_spot,
            min_genes=config.min_genes_spot,
            # Prevalence filtering only makes sense on transcriptome-wide assays;
            # on a targeted panel every probe was chosen on purpose.
            min_cells_per_gene=config.min_cells_per_gene,
            min_units_after_qc=config.min_spots_after_qc,
        )
    return QCThresholds(
        min_counts=config.min_counts_single_cell,
        min_genes=config.min_genes_single_cell,
        min_cells_per_gene=0,
        min_units_after_qc=config.min_cells_after_qc,
    )


def _deduplicate(adata: AnnData) -> AnnData:
    """Drop duplicate gene names and make barcodes unique.

    Duplicate ``var_names`` break every downstream ``adata[:, gene]`` lookup, and
    duplicate ``obs_names`` make the object unwritable. The first occurrence of
    each gene is kept, matching the original script.
    """
    _, first_occurrence = np.unique(adata.var_names.to_numpy(), return_index=True)
    if len(first_occurrence) < adata.n_vars:
        logger.info("Dropping %d duplicate gene name(s)", adata.n_vars - len(first_occurrence))
        # np.unique sorts; restore original order so the matrix stays comparable
        # to the source file.
        adata = adata[:, np.sort(first_occurrence)].copy()
    adata.obs_names_make_unique()
    return adata


def _add_qc_metrics(adata: AnnData) -> None:
    """Attach ``n_counts``, ``n_genes`` and, when detectable, ``percent_mt``."""
    scanpy = require("scanpy", feature="quality control")

    adata.obs["n_counts"] = row_sums(adata.X)
    adata.obs["n_genes"] = count_nonzero_per_row(adata.X)

    # Mitochondrial prefix is case-insensitive: human uses MT-, mouse Mt-.
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    n_mt = int(adata.var["mt"].sum())
    if n_mt == 0:
        logger.debug("No mitochondrial genes found; skipping percent_mt")
        return

    scanpy.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
    adata.obs["percent_mt"] = adata.obs["pct_counts_mt"].astype(np.float32)
    logger.debug("percent_mt computed from %d mitochondrial gene(s)", n_mt)


def apply_qc(
    adata: AnnData,
    thresholds: QCThresholds,
    config: QCConfig,
) -> tuple[AnnData, str | None]:
    """Filter barcodes and genes.

    Returns the filtered object and, when the sample should be rejected, a reason
    string. Rejection is a return value rather than an exception because an
    under-sized section is an expected outcome of a heterogeneous atlas, not an
    error.
    """
    if adata.n_obs == 0:
        return adata, "empty sample (0 observations)"

    n_obs_before, n_vars_before = adata.shape

    if config.filter_cells:
        keep = np.ones(adata.n_obs, dtype=bool)
        if thresholds.min_counts > 0:
            keep &= adata.obs["n_counts"].to_numpy() >= thresholds.min_counts
        if thresholds.min_genes > 0:
            keep &= adata.obs["n_genes"].to_numpy() >= thresholds.min_genes
        if not keep.all():
            adata = adata[keep].copy()

    if adata.n_obs == 0:
        return adata, f"no barcodes passed QC (from {n_obs_before})"

    if config.filter_genes and thresholds.min_cells_per_gene > 0:
        per_gene = count_nonzero_per_column(adata.X)
        keep_genes = per_gene >= thresholds.min_cells_per_gene
        if not keep_genes.any():
            return adata, (
                f"no genes passed the min_cells_per_gene={thresholds.min_cells_per_gene} filter"
            )
        if not keep_genes.all():
            adata = adata[:, keep_genes].copy()

    logger.info(
        "QC: %d x %d -> %d x %d",
        n_obs_before,
        n_vars_before,
        adata.n_obs,
        adata.n_vars,
    )

    if adata.n_obs < thresholds.min_units_after_qc:
        return adata, (
            f"only {adata.n_obs} unit(s) left after QC, "
            f"below the minimum of {thresholds.min_units_after_qc}"
        )
    return adata, None


def _normalise_and_embed(adata: AnnData, config: QCConfig) -> None:
    """Normalise, select HVGs, then run PCA, neighbours, clustering and UMAP."""
    scanpy = require("scanpy", feature="quality control")

    if looks_like_counts(adata.X, config.integer_check_sample_size):
        logger.info("Matrix looks like raw counts; normalising and log1p-transforming")
        scanpy.pp.normalize_total(adata, target_sum=config.target_sum)
        scanpy.pp.log1p(adata)
    else:
        logger.info("Matrix appears already normalised; skipping log1p")

    n_top = min(config.n_top_genes, adata.n_vars)
    scanpy.pp.highly_variable_genes(adata, flavor=config.hvg_flavor, n_top_genes=n_top)

    # PCA cannot return more components than min(n_obs, n_vars) - 1.
    n_comps = max(1, min(config.n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    scanpy.tl.pca(adata, n_comps=n_comps)

    n_neighbors = max(2, min(config.n_neighbors, adata.n_obs - 1))
    scanpy.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps)

    cluster = scanpy.tl.leiden if config.cluster_method == "leiden" else scanpy.tl.louvain
    cluster_kwargs = {"resolution": config.resolution, "key_added": config.cluster_method}
    if config.cluster_method == "leiden":
        # Silence the future-default warning and pin the implementation, so
        # cluster labels are reproducible across scanpy versions.
        cluster_kwargs["flavor"] = "igraph"
        cluster_kwargs["n_iterations"] = 2
        cluster_kwargs["directed"] = False
    cluster(adata, **cluster_kwargs)

    if config.compute_umap:
        scanpy.tl.umap(adata)


def qc_sample(
    adata: AnnData,
    *,
    config: QCConfig,
    resolution: Resolution | None,
    sample_name: str = "sample",
) -> tuple[AnnData | None, str | None]:
    """Run the whole QC stage on an in-memory object.

    Separated from the file loop so it can be tested and reused without touching
    the filesystem.

    Returns:
        ``(processed, None)`` on success, or ``(None, reason)`` when the sample is
        rejected.
    """
    logger.info("QC %s (resolution=%s)", sample_name, resolution.value if resolution else "unknown")
    thresholds = thresholds_for(resolution, config)

    adata = _deduplicate(adata)

    # Preserve raw counts before any transformation; the programs stage needs them.
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    _add_qc_metrics(adata)
    adata, reason = apply_qc(adata, thresholds, config)
    if reason is not None:
        return None, reason

    _normalise_and_embed(adata, config)
    return adata, None


def run_qc(config: Config, *, exclude: tuple[str, ...] = ()) -> BatchResult:
    """Run QC over every sample under ``paths.raw``, writing to ``paths.qc``.

    Args:
        config: Full configuration; uses ``paths.raw``, ``paths.qc`` and the ``qc``
            section.
        exclude: File names to skip.

    Returns:
        A :class:`~stcompass.pipelines._batch.BatchResult`.
    """
    src_root = config.paths.require("raw")
    dst_root = config.paths.require("qc")

    def worker(pair: SamplePair) -> str | None:
        platform = infer_platform(pair.source, relative_to=src_root)
        resolution = resolution_of(platform) if platform else None
        if platform is None:
            logger.warning(
                "Could not infer platform for %s; using permissive thresholds",
                pair.source.name,
            )
        adata = read_h5ad(pair.source)
        processed, reason = qc_sample(
            adata,
            config=config.qc,
            resolution=resolution,
            sample_name=pair.source.name,
        )
        if processed is None:
            return reason
        write_h5ad(processed, pair.destination)
        logger.info("Wrote %s", pair.destination)
        return None

    pairs = iter_samples(src_root, dst_root, overwrite=config.overwrite, exclude=exclude)
    return run_batch(pairs, worker, n_jobs=config.n_jobs, description="QC samples")
