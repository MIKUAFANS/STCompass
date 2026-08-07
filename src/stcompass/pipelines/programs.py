"""Gene programs by non-negative matrix factorisation.

Replaces ``spatialVEG_0923.py``.  The factorisation decomposes the expression
matrix ``X`` (barcodes x genes) into ``W`` (barcodes x K) and ``H`` (K x genes):
``W`` says how much of each program a barcode expresses, ``H`` says which genes
make up each program.  Unlike PCA the factors are non-negative and therefore
additive, which is what makes them readable as "programs" rather than axes of
variation.

The implementation is deliberately streaming.  A Visium HD section has millions of
2 µm bins; densifying such a matrix costs ``n_obs * n_vars * 4`` bytes, which is
hundreds of gigabytes.  So the matrix is read from a backed ``.h5ad`` in row
blocks, kept sparse, and fitted with :class:`sklearn.decomposition.MiniBatchNMF`
through repeated ``partial_fit`` calls.

Outputs written back to the sample:

* ``obsm["X_nmf"]`` -- the ``W`` loadings;
* ``var["NNMF_component_<k>"]`` -- the ``H`` weight of each gene in program ``k``,
  aligned to the full gene list with zeros for genes excluded from the fit;
* ``var["NNMF_component_<k>_importance_score"]`` and ``..._rank`` -- per-program
  gene specificity, see :func:`gene_importance`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.sparse as sp

from ..config import Config, ProgramsConfig
from ..io import SamplePair, iter_samples, read_h5ad, write_h5ad
from ..logging_utils import get_logger
from ..matrix import column_variance, nonzero_column_mask, stack_selected_columns
from ._batch import BatchResult, run_batch

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData

__all__ = [
    "factorise",
    "gene_importance",
    "programs_sample",
    "run_programs",
    "select_genes",
]

logger = get_logger(__name__)

COMPONENT_PREFIX = "NNMF_component_"
LOADINGS_KEY = "X_nmf"


def _source_matrix(adata: AnnData, use_counts_layer: bool):
    """Return the matrix to factorise.

    Prefers ``layers["counts"]`` so the factorisation sees raw counts even when
    ``X`` has already been log-normalised by the QC stage.  NMF's multiplicative
    updates assume non-negative, roughly additive input, which counts satisfy.
    """
    if use_counts_layer and "counts" in adata.layers:
        logger.debug("Factorising layers['counts']")
        return adata.layers["counts"]
    logger.debug("Factorising X (no 'counts' layer available)")
    return adata.X


def select_genes(matrix, config: ProgramsConfig) -> np.ndarray:
    """Choose the gene columns to factorise.

    All-zero genes are always dropped: they contribute nothing to the loss and
    would leave a zero row in ``H``.  When :attr:`ProgramsConfig.max_hvg` is set,
    the remaining genes are ranked by variance and the top ones kept, which bounds
    both runtime and memory for transcriptome-wide platforms.

    Returns:
        Sorted indices into the original gene axis.  Sorted so the selection is
        deterministic and so ``H`` columns can be mapped back by position.
    """
    keep = np.flatnonzero(nonzero_column_mask(matrix, config.row_chunk))
    logger.info("Genes: %d total, %d non-empty", matrix.shape[1], keep.size)

    if config.max_hvg is not None and keep.size > config.max_hvg:
        variance = column_variance(matrix, keep, config.row_chunk)
        # argsort ascending then take the tail: the top-variance genes.
        top = np.argsort(variance, kind="stable")[-config.max_hvg :]
        keep = np.sort(keep[top])
        logger.info("Selected %d high-variance gene(s)", keep.size)
    return keep


def factorise(
    matrix: sp.csr_matrix,
    config: ProgramsConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit MiniBatchNMF and return ``(W, H)``.

    ``partial_fit`` is called over row blocks for :attr:`ProgramsConfig.epochs`
    passes.  More than one pass matters: a single pass leaves the factors biased
    towards the rows seen last, because each mini-batch update only partially
    corrects the previous one.

    Args:
        matrix: Barcodes x selected genes, CSR, non-negative.
        config: Factorisation settings.

    Returns:
        ``W`` with shape ``(n_obs, K)`` and ``H`` with shape ``(K, n_selected)``.
    """
    from sklearn.decomposition import MiniBatchNMF

    n_obs, n_genes = matrix.shape
    # NMF cannot extract more components than the smaller matrix dimension.
    n_components = max(1, min(config.n_components, n_obs, n_genes))
    if n_components < config.n_components:
        logger.warning(
            "Reducing n_components from %d to %d for a %d x %d matrix",
            config.n_components,
            n_components,
            n_obs,
            n_genes,
        )

    # 'nndsvda' needs a batch at least as large as n_components to initialise.
    batch_size = max(n_components, min(config.batch_size, n_obs))

    model = MiniBatchNMF(
        n_components=n_components,
        init=config.init,
        batch_size=batch_size,
        random_state=config.random_state,
        max_no_improvement=config.max_no_improvement,
    )

    for epoch in range(config.epochs):
        logger.info("NMF partial_fit epoch %d/%d", epoch + 1, config.epochs)
        for start in range(0, n_obs, batch_size):
            model.partial_fit(matrix[start : start + batch_size])

    dtype = np.float32 if config.save_float32 else np.float64
    loadings = np.zeros((n_obs, n_components), dtype=dtype)
    for start in range(0, n_obs, batch_size):
        stop = min(start + batch_size, n_obs)
        loadings[start:stop] = model.transform(matrix[start:stop]).astype(dtype, copy=False)

    return loadings, model.components_


def gene_importance(weights: np.ndarray) -> np.ndarray:
    """Score how specifically each gene marks each program.

    A gene with a high weight in *every* program is a housekeeping gene and says
    nothing about program identity, so a raw ``H`` weight is a poor marker score.
    This applies the specificity transform from the original script:

    ``score[g, k] = H[g, k] * log(1 + H[g, k] / reference[g, k])``

    where ``reference[g, k]`` is the gene's largest weight across the *other*
    programs -- concretely, the top weight for every program except the top one
    itself, which is compared against the runner-up.  A gene that loads on one
    program only gets a large ratio and scores highly; a flat gene gets a ratio
    near one and scores near zero.

    Args:
        weights: Genes x programs, non-negative.

    Returns:
        Array of the same shape holding the specificity scores.

    >>> import numpy as np
    >>> h = np.array([[10.0, 0.1], [5.0, 5.0]])
    >>> scores = gene_importance(h)
    >>> scores[0, 0] > scores[1, 0]   # specific gene beats the flat one
    True
    """
    weights = np.asarray(weights, dtype=np.float64)
    if weights.ndim != 2:
        raise ValueError(f"weights must be 2-D, got shape {weights.shape}")
    n_genes, n_components = weights.shape
    if n_components == 1:
        # With one program there is nothing to be specific against; the weight
        # itself is the only available ranking.
        return weights.copy()

    # Per row: the largest and second-largest weights.
    order = np.argsort(weights, axis=1, kind="stable")
    top_idx = order[:, -1]
    largest = weights[np.arange(n_genes), top_idx]
    runner_up = weights[np.arange(n_genes), order[:, -2]]

    # Reference for every program is the row max, except for the top program
    # itself, which is compared against the runner-up.
    reference = np.repeat(largest[:, None], n_components, axis=1)
    reference[np.arange(n_genes), top_idx] = runner_up

    # 1e-10 guards the division for genes that are zero everywhere.
    return weights * np.log1p(weights / (reference + 1e-10))


def _write_factors(
    adata: AnnData,
    loadings: np.ndarray,
    components: np.ndarray,
    genes: np.ndarray,
    config: ProgramsConfig,
) -> None:
    """Store ``W``, ``H`` and the importance scores on ``adata``."""
    dtype = np.float32 if config.save_float32 else np.float64
    n_components = components.shape[0]

    adata.obsm[LOADINGS_KEY] = loadings

    # Re-align H to the full gene axis: genes excluded from the fit get zero, so
    # var columns stay comparable across samples with different gene sets.
    full = np.zeros((adata.n_vars, n_components), dtype=dtype)
    full[genes, :] = components.T.astype(dtype, copy=False)

    scores = gene_importance(full)
    for k in range(n_components):
        adata.var[f"{COMPONENT_PREFIX}{k}"] = full[:, k]
        adata.var[f"{COMPONENT_PREFIX}{k}_importance_score"] = scores[:, k].astype(
            dtype, copy=False
        )
        # Rank 1 = most specific gene for this program.
        ranks = np.empty(adata.n_vars, dtype=np.int32)
        ranks[np.argsort(scores[:, k], kind="stable")[::-1]] = np.arange(1, adata.n_vars + 1)
        adata.var[f"{COMPONENT_PREFIX}{k}_rank"] = ranks

    logger.info("Wrote %d gene program(s) to var/obsm", n_components)


def programs_sample(source, config: ProgramsConfig) -> tuple[AnnData | None, str | None]:
    """Factorise one sample, reading it twice by design.

    The first read is *backed* -- the matrix stays on disk while gene selection
    streams over it, so a million-barcode sample costs no more memory than one row
    block.  Only the selected columns are then materialised in memory, and the
    object is re-read unbacked so the factors can be attached and written.

    Returns:
        ``(adata, None)`` on success, or ``(None, reason)`` when the sample cannot
        be factorised.
    """
    backed = read_h5ad(source, backed="r")
    try:
        if backed.n_obs == 0:
            return None, "empty sample (0 observations)"
        if backed.n_obs > config.large_sample_warning:
            logger.warning(
                "%s has %d observations; factorisation will be slow",
                getattr(source, "name", source),
                backed.n_obs,
            )

        matrix = _source_matrix(backed, config.use_counts_layer)
        genes = select_genes(matrix, config)
        if genes.size == 0:
            return None, "no non-empty genes"
        if genes.size < config.n_components:
            return None, (
                f"only {genes.size} usable gene(s), fewer than n_components={config.n_components}"
            )

        selected = stack_selected_columns(matrix, genes, config.row_chunk)
    finally:
        # Release the HDF5 handle before re-opening the file unbacked.
        if backed.isbacked and backed.file is not None:
            backed.file.close()

    if selected.nnz == 0:
        return None, "selected submatrix is all zeros"
    # NMF requires non-negativity; a log-transformed matrix can hold small
    # negatives after scaling, and sklearn's error message does not name the file.
    if selected.data.min() < 0:
        return None, "matrix contains negative values; NMF requires non-negative input"

    logger.info("Factorising a %d x %d matrix", *selected.shape)
    loadings, components = factorise(selected, config)

    adata = read_h5ad(source)
    _write_factors(adata, loadings, components, genes, config)
    return adata, None


def run_programs(config: Config, *, exclude: tuple[str, ...] = ()) -> BatchResult:
    """Factorise every sample under ``paths.qc``, writing to ``paths.programs``.

    Args:
        config: Full configuration; uses ``paths.qc``, ``paths.programs`` and the
            ``programs`` section.
        exclude: File names to skip.
    """
    src_root = config.paths.require("qc")
    dst_root = config.paths.require("programs")
    settings = config.programs

    def worker(pair: SamplePair) -> str | None:
        adata, reason = programs_sample(pair.source, settings)
        if adata is None:
            return reason
        write_h5ad(adata, pair.destination)
        logger.info("Wrote %s", pair.destination)
        return None

    pairs = iter_samples(src_root, dst_root, overwrite=config.overwrite, exclude=exclude)
    return run_batch(pairs, worker, n_jobs=config.n_jobs, description="NMF samples")
