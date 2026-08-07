"""Standalone (re-)clustering of QC'd samples.

Replaces ``run_louvain.py``.  The one idea worth keeping from that script is the
size-aware resolution: a 200-spot section clustered at ``resolution=0.2`` collapses
into a single domain, while a 500k-bin Visium HD section clustered at ``1.2``
shatters into hundreds of unusable fragments.  The schedule that maps sample size
to resolution now lives in :class:`~stcompass.config.ClusteringConfig` instead of
an ``if/elif`` chain, so it can be changed without editing code.

The original also caught *every* read error and responded by deleting
``uns/log1p/base`` from the file -- a repair that only applies to one specific
failure.  That workaround now lives in :func:`stcompass.io.read_h5ad`, which
retries once and re-raises anything else.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._deps import require
from ..config import ClusteringConfig, Config
from ..io import SamplePair, iter_samples, read_h5ad, write_h5ad
from ..logging_utils import get_logger
from ._batch import BatchResult, run_batch

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData

__all__ = ["cluster_sample", "run_clustering"]

logger = get_logger(__name__)


def cluster_sample(
    adata: AnnData,
    config: ClusteringConfig,
    *,
    sample_name: str = "sample",
) -> str | None:
    """Cluster ``adata`` in place at a size-appropriate resolution.

    Args:
        adata: QC'd sample; expected to carry a neighbour graph from the QC stage.
        config: Method, resolution schedule and graph parameters.
        sample_name: Used in log messages.

    Returns:
        ``None`` on success, or a reason string when the sample cannot be
        clustered (too few observations to form a graph).
    """
    scanpy = require("scanpy", feature="clustering")

    if adata.n_obs < 3:
        return f"only {adata.n_obs} observation(s); too few to cluster"

    resolution = config.resolution_for(adata.n_obs)
    key_added = config.key_added or config.method

    # Reuse the graph from QC when present: recomputing it is the expensive part,
    # and re-running with a different resolution is the common case.
    if config.recompute_neighbors or "neighbors" not in adata.uns:
        if "X_pca" not in adata.obsm:
            n_comps = max(1, min(config.n_pcs, adata.n_obs - 1, adata.n_vars - 1))
            logger.info("%s: no PCA found; computing %d components", sample_name, n_comps)
            scanpy.tl.pca(adata, n_comps=n_comps)
        n_neighbors = max(2, min(config.n_neighbors, adata.n_obs - 1))
        n_pcs = min(config.n_pcs, adata.obsm["X_pca"].shape[1])
        scanpy.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    logger.info(
        "%s: %s clustering at resolution=%.2f (n_obs=%d) -> obs['%s']",
        sample_name,
        config.method,
        resolution,
        adata.n_obs,
        key_added,
    )

    if config.method == "leiden":
        scanpy.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key_added,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    else:
        scanpy.tl.louvain(adata, resolution=resolution, key_added=key_added)

    n_clusters = adata.obs[key_added].nunique()
    logger.info("%s: found %d cluster(s)", sample_name, n_clusters)
    return None


def run_clustering(config: Config, *, exclude: tuple[str, ...] = ()) -> BatchResult:
    """Cluster every sample under ``paths.qc``, writing to ``paths.clustered``.

    Args:
        config: Uses ``paths.qc``, ``paths.clustered`` and the ``clustering`` section.
        exclude: File names to skip.
    """
    src_root = config.paths.require("qc")
    dst_root = config.paths.require("clustered")

    def worker(pair: SamplePair) -> str | None:
        adata = read_h5ad(pair.source)
        reason = cluster_sample(adata, config.clustering, sample_name=pair.source.name)
        if reason is not None:
            return reason
        write_h5ad(adata, pair.destination)
        logger.info("Wrote %s", pair.destination)
        return None

    pairs = iter_samples(src_root, dst_root, overwrite=config.overwrite, exclude=exclude)
    return run_batch(pairs, worker, n_jobs=config.n_jobs, description="clustering samples")
