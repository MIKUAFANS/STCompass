"""Tests for :mod:`stcompass.matrix`.

These cover the numerics that decide pipeline behaviour: whether a matrix is
treated as raw counts (which controls log-transformation), and the streamed
statistics used for gene selection.  Every case is checked against a dense
reference computed with plain numpy, so a regression in the block-wise code shows
up as a disagreement rather than as a plausible-looking wrong number.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from stcompass.matrix import (
    as_csr,
    block_slices,
    column_variance,
    count_nonzero_per_column,
    count_nonzero_per_row,
    iter_row_blocks,
    looks_like_counts,
    nonzero_column_mask,
    row_sums,
    stack_selected_columns,
)


class TestBlockSlices:
    def test_covers_every_row_exactly_once(self) -> None:
        slices = list(block_slices(10, 3))
        covered = [i for s in slices for i in range(s.start, s.stop)]
        assert covered == list(range(10))

    def test_exact_multiple(self) -> None:
        assert [(s.start, s.stop) for s in block_slices(6, 3)] == [(0, 3), (3, 6)]

    def test_block_larger_than_input(self) -> None:
        assert [(s.start, s.stop) for s in block_slices(2, 100)] == [(0, 2)]

    def test_zero_rows_yields_nothing(self) -> None:
        assert list(block_slices(0, 5)) == []

    def test_rejects_nonpositive_block(self) -> None:
        with pytest.raises(ValueError, match="block must be >= 1"):
            list(block_slices(5, 0))


class TestAsCsr:
    def test_csr_passes_through_without_copy(self) -> None:
        original = sp.csr_matrix(np.eye(3))
        assert as_csr(original) is original

    def test_converts_coo(self) -> None:
        converted = as_csr(sp.coo_matrix(np.eye(3)))
        assert sp.isspmatrix_csr(converted)

    def test_converts_dense(self) -> None:
        converted = as_csr(np.array([[1.0, 0.0], [0.0, 2.0]]))
        assert sp.isspmatrix_csr(converted)
        assert converted.sum() == pytest.approx(3.0)


class TestLooksLikeCounts:
    """The decision that controls whether QC applies ``log1p``.

    Getting it wrong is silent and damaging in both directions: transforming
    already-normalised data compresses real variation, while leaving counts
    untransformed breaks PCA.
    """

    def test_integer_sparse_is_counts(self) -> None:
        assert looks_like_counts(sp.csr_matrix(np.array([[0, 3], [2, 0]])))

    def test_integer_dense_is_counts(self) -> None:
        assert looks_like_counts(np.array([[0, 5], [7, 1]]))

    def test_float_valued_is_not_counts(self) -> None:
        assert not looks_like_counts(np.array([[0.0, 1.4], [0.7, 0.0]]))

    def test_log_transformed_is_not_counts(self) -> None:
        counts = np.array([[0, 3, 10], [5, 0, 2]], dtype=np.float64)
        assert not looks_like_counts(np.log1p(counts))

    def test_float_dtype_holding_whole_numbers_is_counts(self) -> None:
        # A float32 matrix of whole numbers is the usual on-disk form for counts.
        assert looks_like_counts(
            sp.csr_matrix(np.array([[0.0, 3.0], [2.0, 0.0]], dtype=np.float32))
        )

    def test_empty_sparse_counts_as_integral(self) -> None:
        # Nothing to transform either way; the permissive answer avoids a crash.
        assert looks_like_counts(sp.csr_matrix((5, 5)))

    def test_empty_dense_counts_as_integral(self) -> None:
        assert looks_like_counts(np.zeros((0, 3)))

    def test_all_zero_dense_counts_as_integral(self) -> None:
        assert looks_like_counts(np.zeros((4, 4)))

    def test_nan_is_not_counts(self) -> None:
        assert not looks_like_counts(np.array([[1.0, np.nan], [2.0, 3.0]]))

    def test_infinity_is_not_counts(self) -> None:
        assert not looks_like_counts(np.array([[1.0, np.inf], [2.0, 3.0]]))

    def test_deterministic_across_calls(self) -> None:
        """Same input, same answer -- the original used the global RNG and could
        classify one file differently on a re-run."""
        rng = np.random.default_rng(1)
        matrix = sp.random(200, 300, density=0.2, format="csr", random_state=1)
        matrix.data = np.round(matrix.data * 10) + 0.5  # deliberately non-integral
        answers = {looks_like_counts(matrix, sample_size=50) for _ in range(20)}
        assert len(answers) == 1
        assert rng is not None  # rng unused beyond documenting intent

    def test_sampling_finds_fractional_values_in_large_matrix(self) -> None:
        matrix = sp.random(500, 500, density=0.3, format="csr", random_state=7)
        matrix.data = matrix.data + 0.25
        assert not looks_like_counts(matrix, sample_size=100)

    def test_dense_prefers_nonzero_entries(self) -> None:
        """A mostly-zero dense matrix must not be judged integral just because
        zeros dominate the sample."""
        matrix = np.zeros((100, 100))
        matrix[::10, ::10] = 1.5
        assert not looks_like_counts(matrix, sample_size=50)


class TestRowStatistics:
    def test_row_sums_matches_numpy(self) -> None:
        dense = np.array([[1.0, 2.0, 0.0], [0.0, 0.0, 5.0]])
        expected = dense.sum(axis=1)
        assert row_sums(dense) == pytest.approx(expected)
        assert row_sums(sp.csr_matrix(dense)) == pytest.approx(expected)

    def test_row_sums_is_one_dimensional(self) -> None:
        assert row_sums(sp.csr_matrix(np.eye(3))).shape == (3,)

    def test_count_nonzero_per_row_matches_numpy(self) -> None:
        dense = np.array([[1.0, 0.0, 3.0], [0.0, 0.0, 0.0]])
        expected = np.count_nonzero(dense, axis=1)
        assert count_nonzero_per_row(dense).tolist() == expected.tolist()
        assert count_nonzero_per_row(sp.csr_matrix(dense)).tolist() == expected.tolist()


class TestColumnStatistics:
    @pytest.fixture
    def dense(self) -> np.ndarray:
        rng = np.random.default_rng(0)
        matrix = rng.poisson(1.5, size=(97, 41)).astype(np.float64)
        matrix[:, 3] = 0.0  # an all-zero gene
        matrix[:, 17] = 0.0
        return matrix

    def test_count_nonzero_per_column_matches_numpy(self, dense: np.ndarray) -> None:
        expected = np.count_nonzero(dense, axis=0)
        for block in (1, 7, 1000):
            got = count_nonzero_per_column(sp.csr_matrix(dense), block=block)
            assert got.tolist() == expected.tolist()

    def test_nonzero_column_mask_identifies_empty_genes(self, dense: np.ndarray) -> None:
        mask = nonzero_column_mask(sp.csr_matrix(dense), block=10)
        assert not mask[3]
        assert not mask[17]
        assert mask.sum() == dense.shape[1] - 2

    def test_nonzero_column_mask_all_true_short_circuits(self) -> None:
        full = np.ones((50, 4))
        assert nonzero_column_mask(sp.csr_matrix(full), block=5).all()

    def test_column_variance_matches_numpy(self, dense: np.ndarray) -> None:
        expected = dense.var(axis=0)
        for block in (1, 13, 500):
            got = column_variance(sp.csr_matrix(dense), block=block)
            assert got == pytest.approx(expected, abs=1e-8)

    def test_column_variance_with_subset(self, dense: np.ndarray) -> None:
        columns = np.array([0, 5, 9, 20])
        expected = dense[:, columns].var(axis=0)
        got = column_variance(sp.csr_matrix(dense), columns=columns, block=11)
        assert got == pytest.approx(expected, abs=1e-8)

    def test_column_variance_never_negative(self) -> None:
        """Cancellation in the sum-of-squares identity can push a near-zero
        variance below zero; callers must never see an impossible value."""
        constant = np.full((200, 5), 1e6)
        assert (column_variance(sp.csr_matrix(constant), block=17) >= 0.0).all()

    def test_column_variance_of_empty_matrix(self) -> None:
        assert column_variance(sp.csr_matrix((0, 4))).tolist() == [0.0, 0.0, 0.0, 0.0]


class TestStackSelectedColumns:
    def test_reproduces_dense_slice(self) -> None:
        rng = np.random.default_rng(3)
        dense = rng.poisson(1.0, size=(63, 20)).astype(np.float64)
        columns = np.array([1, 4, 4, 11])  # duplicates must be preserved in order
        for block in (1, 8, 100):
            got = stack_selected_columns(sp.csr_matrix(dense), columns, block=block)
            assert got.shape == (63, len(columns))
            assert np.allclose(got.toarray(), dense[:, columns])

    def test_empty_input_returns_empty_matrix(self) -> None:
        got = stack_selected_columns(sp.csr_matrix((0, 5)), np.array([0, 2]))
        assert got.shape == (0, 2)


class TestIterRowBlocks:
    def test_blocks_reassemble_into_original(self) -> None:
        rng = np.random.default_rng(5)
        dense = rng.poisson(1.0, size=(50, 6)).astype(np.float64)
        pieces = []
        for rows, chunk in iter_row_blocks(sp.csr_matrix(dense), 7):
            assert sp.isspmatrix_csr(chunk)
            assert chunk.shape[0] == rows.stop - rows.start
            pieces.append(chunk.toarray())
        assert np.allclose(np.vstack(pieces), dense)

    def test_works_on_dense_input(self) -> None:
        dense = np.arange(12, dtype=np.float64).reshape(4, 3)
        pieces = [chunk.toarray() for _, chunk in iter_row_blocks(dense, 2)]
        assert np.allclose(np.vstack(pieces), dense)
