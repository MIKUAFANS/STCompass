"""Tests for the gene-programs stage.

The interesting part is :func:`gene_importance`.  The original script scored genes
in a per-gene Python loop; this package computes the same quantity with array
operations.  A refactor like that is exactly where a silent numerical regression
hides, so the properties of the score are pinned here, and
:func:`test_matches_reference_loop` re-implements the original formula and asserts
agreement on random input.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from stcompass.config import ProgramsConfig
from stcompass.pipelines.programs import factorise, gene_importance, select_genes


def _reference_loop(weights: np.ndarray) -> np.ndarray:
    """The scoring loop from ``spatialVEG_0923.py``, kept as an oracle.

    For each gene, the score of program ``k`` compares its weight against the
    largest weight in *any other* program, so a gene loading on one program alone
    scores high and a ubiquitous gene scores near zero.
    """
    values = np.asarray(weights, dtype=np.float64)
    n_genes, _ = values.shape
    scores = np.zeros_like(values)
    for i in range(n_genes):
        row = values[i, :]
        first = int(np.argmax(row))
        without_max = row.copy()
        without_max[first] = -np.inf
        second = int(np.argmax(without_max))
        competitor = np.full_like(row, row[first])
        competitor[first] = row[second]
        scores[i, :] = row * np.log(1 + row / (competitor + 1e-10))
    return scores


class TestGeneImportance:
    def test_matches_reference_loop(self):
        """Vectorised scoring agrees with the original loop on random input."""
        rng = np.random.default_rng(20260805)
        for _ in range(200):
            n_genes = int(rng.integers(1, 40))
            n_programs = int(rng.integers(2, 8))
            weights = rng.random((n_genes, n_programs)) * rng.choice([1.0, 10.0, 0.01])
            np.testing.assert_allclose(
                gene_importance(weights), _reference_loop(weights), rtol=1e-9, atol=1e-12
            )

    def test_program_specific_gene_outscores_ubiquitous_gene(self):
        """A gene loading on one program must score above a flat gene."""
        weights = np.array([[10.0, 0.01], [5.0, 5.0]])
        scores = gene_importance(weights)
        assert scores[0, 0] > scores[1, 0]

    def test_flat_gene_scores_equally_across_programs(self):
        weights = np.array([[7.0, 7.0, 7.0]])
        scores = gene_importance(weights)
        assert scores[0, 0] == pytest.approx(scores[0, 1])
        assert scores[0, 1] == pytest.approx(scores[0, 2])

    def test_all_zero_gene_scores_zero(self):
        """An unexpressed gene must not produce NaN from the log."""
        scores = gene_importance(np.array([[0.0, 0.0, 0.0]]))
        assert np.all(np.isfinite(scores))
        assert np.allclose(scores, 0.0)

    def test_output_is_finite_for_extreme_ratios(self):
        weights = np.array([[1e6, 1e-12], [1e-12, 1e6]])
        assert np.all(np.isfinite(gene_importance(weights)))

    def test_single_program_is_handled(self):
        """K=1 has no competitor; the score must still be defined and finite."""
        scores = gene_importance(np.array([[3.0], [0.0]]))
        assert scores.shape == (2, 1)
        assert np.all(np.isfinite(scores))

    def test_shape_preserved(self):
        weights = np.random.default_rng(0).random((17, 5))
        assert gene_importance(weights).shape == (17, 5)

    def test_rejects_non_2d_input(self):
        with pytest.raises(ValueError):
            gene_importance(np.array([1.0, 2.0, 3.0]))


class TestSelectGenes:
    def test_drops_all_zero_genes(self):
        matrix = sp.csr_matrix(np.array([[1.0, 0.0, 2.0], [3.0, 0.0, 4.0]]))
        selected = select_genes(matrix, ProgramsConfig(n_components=1))
        assert selected.tolist() == [0, 2]

    def test_hvg_limit_keeps_highest_variance(self):
        # Column 2 varies most, then column 0, then column 1.
        matrix = sp.csr_matrix(
            np.array(
                [
                    [1.0, 5.0, 0.0],
                    [2.0, 5.0, 50.0],
                    [3.0, 5.0, 100.0],
                ]
            )
        )
        selected = select_genes(matrix, ProgramsConfig(n_components=1, max_hvg=2))
        assert selected.tolist() == [0, 2]

    def test_selection_is_sorted(self):
        """Indices must be sorted so H columns map back by position."""
        rng = np.random.default_rng(3)
        matrix = sp.csr_matrix(rng.random((20, 30)))
        selected = select_genes(matrix, ProgramsConfig(n_components=2, max_hvg=7))
        assert selected.tolist() == sorted(selected.tolist())
        assert len(selected) == 7

    def test_hvg_larger_than_gene_count_keeps_all(self):
        matrix = sp.csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]]))
        selected = select_genes(matrix, ProgramsConfig(n_components=1, max_hvg=100))
        assert selected.tolist() == [0, 1]

    def test_returns_empty_when_no_gene_is_expressed(self):
        """An all-zero sample yields no genes rather than raising.

        ``programs_sample`` turns the empty selection into a recorded skip, so one
        dead file does not abort a batch over thousands of samples.
        """
        matrix = sp.csr_matrix((4, 5), dtype=np.float32)
        selected = select_genes(matrix, ProgramsConfig(n_components=1))
        assert selected.size == 0


class TestFactorise:
    def test_returns_nonnegative_factors_of_expected_shape(self):
        rng = np.random.default_rng(11)
        matrix = sp.csr_matrix(rng.random((60, 12)))
        config = ProgramsConfig(n_components=3, batch_size=16, epochs=1, random_state=0)
        loadings, components = factorise(matrix, config)
        assert loadings.shape == (60, 3)
        assert components.shape == (3, 12)
        assert (loadings >= 0).all()
        assert (components >= 0).all()

    def test_is_deterministic_for_a_fixed_seed(self):
        rng = np.random.default_rng(12)
        matrix = sp.csr_matrix(rng.random((40, 10)))
        config = ProgramsConfig(n_components=2, batch_size=8, epochs=1, random_state=7)
        first = factorise(matrix, config)
        second = factorise(matrix, config)
        np.testing.assert_allclose(first[0], second[0])
        np.testing.assert_allclose(first[1], second[1])

    def test_float32_output_when_requested(self):
        rng = np.random.default_rng(13)
        matrix = sp.csr_matrix(rng.random((30, 8)))
        config = ProgramsConfig(n_components=2, batch_size=8, epochs=1, save_float32=True)
        loadings, _ = factorise(matrix, config)
        assert loadings.dtype == np.float32

    def test_components_capped_by_matrix_size(self):
        """Asking for more programs than genes must not crash."""
        matrix = sp.csr_matrix(np.abs(np.random.default_rng(1).random((5, 3))))
        config = ProgramsConfig(n_components=10, batch_size=2, epochs=1)
        loadings, components = factorise(matrix, config)
        assert loadings.shape[1] <= 3
        assert components.shape[1] == 3

    def test_reconstruction_is_better_than_the_mean_baseline(self):
        """Sanity check that the factorisation actually fits something."""
        rng = np.random.default_rng(5)
        # Rank-2 structure plus a little noise.
        basis = np.abs(rng.random((2, 15)))
        loadings_true = np.abs(rng.random((80, 2)))
        dense = loadings_true @ basis + 0.01 * rng.random((80, 15))
        matrix = sp.csr_matrix(dense)

        config = ProgramsConfig(n_components=2, batch_size=20, epochs=3, random_state=0)
        loadings, components = factorise(matrix, config)

        model_error = np.linalg.norm(dense - loadings @ components)
        baseline_error = np.linalg.norm(dense - dense.mean())
        assert model_error < baseline_error
