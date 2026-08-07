"""Matrix helpers shared by the pipelines.

Everything here works on ``numpy`` arrays and ``scipy.sparse`` matrices only, with
no ``scanpy``/``anndata`` import, which keeps the numerically interesting parts of
the package testable in a minimal environment.

The block-wise readers exist because spatial samples routinely exceed a million
barcodes.  Densifying such a matrix costs ``n_obs * n_vars * 4`` bytes -- hundreds
of gigabytes for a Visium HD section -- so every routine that touches a whole
matrix streams it in row blocks and keeps the sparse representation.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import numpy as np
import scipy.sparse as sp

__all__ = [
    "as_csr",
    "block_slices",
    "column_variance",
    "count_nonzero_per_column",
    "count_nonzero_per_row",
    "iter_row_blocks",
    "looks_like_counts",
    "nonzero_column_mask",
    "row_sums",
    "stack_selected_columns",
]

# Seeded so that the "are these raw counts?" decision is identical across runs;
# the original scripts drew from the global RNG and could classify the same file
# differently on a re-run.
_DEFAULT_SEED = 0


def as_csr(matrix: Any) -> sp.csr_matrix:
    """Return ``matrix`` as CSR without copying when it already is CSR."""
    if sp.isspmatrix_csr(matrix):
        return matrix
    if sp.issparse(matrix):
        return matrix.tocsr()
    return sp.csr_matrix(np.asarray(matrix))


def block_slices(n_rows: int, block: int) -> Iterator[slice]:
    """Yield ``slice`` objects covering ``n_rows`` in chunks of ``block``.

    >>> [ (s.start, s.stop) for s in block_slices(5, 2) ]
    [(0, 2), (2, 4), (4, 5)]
    """
    if block < 1:
        raise ValueError(f"block must be >= 1, got {block}")
    for start in range(0, n_rows, block):
        yield slice(start, min(start + block, n_rows))


def iter_row_blocks(matrix: Any, block: int) -> Iterator[tuple[slice, sp.csr_matrix]]:
    """Iterate ``(rows, csr_block)`` pairs over ``matrix``.

    Works for in-memory arrays and for the h5py-backed datasets exposed by
    ``AnnData(backed='r')``, which support slicing but not sparse arithmetic.
    """
    n_rows = matrix.shape[0]
    for rows in block_slices(n_rows, block):
        yield rows, as_csr(matrix[rows])


def looks_like_counts(
    matrix: Any,
    sample_size: int = 2000,
    *,
    seed: int = _DEFAULT_SEED,
    atol: float = 1e-6,
) -> bool:
    """Guess whether ``matrix`` holds raw integer counts.

    Public atlases mix raw and already-normalised matrices with no reliable
    metadata flag, so QC has to infer it: applying ``log1p`` twice compresses real
    biological variation, while skipping it leaves counts on a scale that breaks
    PCA.  Integer-valued entries are the practical signal, since normalised or
    log-transformed data is essentially never integral.

    Only ``sample_size`` stored values are inspected, which makes the check
    independent of matrix size.  Sampling is done *with* replacement: drawing
    without replacement forces a permutation of every stored value, which for a
    billion-nonzero matrix costs more memory than the pipeline stage that follows.

    Args:
        matrix: Sparse or dense expression matrix.
        sample_size: Number of stored values to inspect.
        seed: RNG seed; fixed by default so the decision is reproducible.
        atol: Absolute tolerance when comparing against rounded values.

    Returns:
        ``True`` when the sampled values are integral (an empty matrix counts as
        integral -- there is nothing to transform either way).

    >>> import scipy.sparse as sp, numpy as np
    >>> looks_like_counts(sp.csr_matrix(np.array([[0, 3], [2, 0]])))
    True
    >>> looks_like_counts(np.array([[0.0, 1.4], [0.7, 0.0]]))
    False
    """
    rng = np.random.default_rng(seed)

    if sp.issparse(matrix):
        data = matrix.data
        n_stored = data.size
        if n_stored == 0:
            return True
        if n_stored <= sample_size:
            sampled = data
        else:
            sampled = data[rng.integers(0, n_stored, size=sample_size)]
    else:
        flat = np.asarray(matrix).reshape(-1)
        if flat.size == 0:
            return True
        # Prefer non-zero entries: a sparse-but-dense-stored matrix is mostly
        # zeros, and zeros are integral regardless of the underlying scale.
        nonzero = flat[flat != 0]
        pool = nonzero if nonzero.size else flat
        if pool.size <= sample_size:
            sampled = pool
        else:
            sampled = pool[rng.integers(0, pool.size, size=sample_size)]

    sampled = np.asarray(sampled, dtype=np.float64)
    if not np.all(np.isfinite(sampled)):
        return False
    return bool(np.allclose(sampled, np.round(sampled), atol=atol, rtol=0))


def row_sums(matrix: Any) -> np.ndarray:
    """Total per row, as a 1-D float array."""
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64, copy=False)
    return np.asarray(matrix).sum(axis=1).astype(np.float64, copy=False).ravel()


def count_nonzero_per_row(matrix: Any) -> np.ndarray:
    """Number of non-zero entries per row, as a 1-D int array."""
    if sp.issparse(matrix):
        return np.asarray((matrix != 0).sum(axis=1)).ravel().astype(np.int64, copy=False)
    return np.count_nonzero(np.asarray(matrix), axis=1).astype(np.int64, copy=False)


def count_nonzero_per_column(matrix: Any, block: int = 5000) -> np.ndarray:
    """Number of non-zero entries per column, streamed in row blocks."""
    n_cols = matrix.shape[1]
    counts = np.zeros(n_cols, dtype=np.int64)
    for _, chunk in iter_row_blocks(matrix, block):
        counts += np.asarray((chunk != 0).sum(axis=0)).ravel()
    return counts


def nonzero_column_mask(matrix: Any, block: int = 5000) -> np.ndarray:
    """Boolean mask of columns with at least one non-zero entry.

    NMF cannot use an all-zero gene -- it contributes nothing to the loss and
    leaves a zero row in the factor -- so these columns are dropped before
    factorising and re-inserted as zeros afterwards.
    """
    n_cols = matrix.shape[1]
    mask = np.zeros(n_cols, dtype=bool)
    for _, chunk in iter_row_blocks(matrix, block):
        mask |= np.asarray((chunk != 0).sum(axis=0)).ravel() > 0
        if mask.all():
            break  # nothing left to discover
    return mask


def column_variance(
    matrix: Any,
    columns: np.ndarray | None = None,
    block: int = 5000,
) -> np.ndarray:
    """Variance of each (selected) column, computed from streamed row blocks.

    Uses the sum/sum-of-squares identity so the matrix is read once and never
    densified.  That form is numerically weaker than a two-pass algorithm, but it
    is only used to *rank* genes for HVG selection, where a small absolute error
    in a large variance does not change the ordering.
    """
    n_rows = matrix.shape[0]
    if n_rows == 0:
        width = matrix.shape[1] if columns is None else len(columns)
        return np.zeros(width, dtype=np.float64)

    width = matrix.shape[1] if columns is None else len(columns)
    total = np.zeros(width, dtype=np.float64)
    total_sq = np.zeros(width, dtype=np.float64)

    for rows in block_slices(n_rows, block):
        chunk = matrix[rows]
        if columns is not None:
            chunk = chunk[:, columns]
        chunk = as_csr(chunk)
        total += np.asarray(chunk.sum(axis=0)).ravel()
        total_sq += np.asarray(chunk.multiply(chunk).sum(axis=0)).ravel()

    mean = total / n_rows
    variance = total_sq / n_rows - mean**2
    # Cancellation in the identity above can push a near-zero variance slightly
    # negative; clamp so callers never see an impossible value.
    return np.maximum(variance, 0.0)


def stack_selected_columns(
    matrix: Any,
    columns: np.ndarray,
    block: int = 5000,
) -> sp.csr_matrix:
    """Build an in-memory CSR matrix from ``columns`` of ``matrix``.

    Reads row blocks from a (possibly file-backed) matrix and concatenates them.
    Memory scales with the number of stored values in the *selected* columns, not
    with the full matrix.
    """
    blocks: list[sp.csr_matrix] = []
    for rows in block_slices(matrix.shape[0], block):
        chunk = as_csr(matrix[rows][:, columns])
        blocks.append(chunk)
    if not blocks:
        return sp.csr_matrix((0, len(columns)), dtype=np.float32)
    return sp.vstack(blocks, format="csr")
