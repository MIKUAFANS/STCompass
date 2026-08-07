"""Tests for the plotting helpers.

Only the geometry and the proportion arithmetic are covered here; rendering a
figure needs a real ``AnnData`` and is exercised by the ``requires_scanpy`` tests.
The parts tested are the ones that silently produce a *wrong picture* rather than
an error: a mis-normalised pie shows the wrong composition, and a mis-scaled
radius makes every pie overlap or vanish.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from stcompass.pipelines.visualize import (
    estimate_radius,
    top_n_proportions,
    wait_until_stable,
)


class TestTopNProportions:
    def test_keeps_only_the_largest_entries(self):
        frame = pd.DataFrame([[0.1, 0.2, 0.3, 0.25, 0.15]], columns=list("abcde"))
        result = top_n_proportions(frame, 3)
        kept = [c for c in result.columns if result.iloc[0][c] > 0]
        assert kept == ["b", "c", "d"]

    def test_rows_are_renormalised_to_one(self):
        frame = pd.DataFrame([[0.1, 0.2, 0.3, 0.25, 0.15]], columns=list("abcde"))
        result = top_n_proportions(frame, 3)
        assert result.iloc[0].sum() == pytest.approx(1.0)

    def test_renormalisation_preserves_relative_magnitude(self):
        """b:c:d were 0.2:0.3:0.25; after renormalising the ordering must hold."""
        frame = pd.DataFrame([[0.1, 0.2, 0.3, 0.25, 0.15]], columns=list("abcde"))
        row = top_n_proportions(frame, 3).iloc[0]
        assert row["c"] > row["d"] > row["b"]
        assert row["c"] / row["b"] == pytest.approx(0.3 / 0.2)

    def test_all_zero_row_stays_zero_without_dividing_by_zero(self):
        frame = pd.DataFrame([[0.0, 0.0, 0.0]], columns=list("abc"))
        result = top_n_proportions(frame, 2)
        assert result.iloc[0].tolist() == [0.0, 0.0, 0.0]
        assert np.isfinite(result.to_numpy()).all()

    def test_top_n_larger_than_column_count_keeps_everything(self):
        frame = pd.DataFrame([[0.25, 0.75]], columns=list("ab"))
        result = top_n_proportions(frame, 10)
        assert result.iloc[0].tolist() == pytest.approx([0.25, 0.75])

    def test_top_n_of_one_yields_a_single_full_wedge(self):
        frame = pd.DataFrame([[0.3, 0.7, 0.1]], columns=list("abc"))
        row = top_n_proportions(frame, 1).iloc[0]
        assert row.tolist() == pytest.approx([0.0, 1.0, 0.0])

    def test_each_row_is_treated_independently(self):
        """The dominant type differs per row, so the cut must be per-row."""
        frame = pd.DataFrame(
            [[0.7, 0.2, 0.1], [0.0, 0.1, 0.9]],
            columns=list("abc"),
        )
        result = top_n_proportions(frame, 1)
        assert result.iloc[0].tolist() == pytest.approx([1.0, 0.0, 0.0])
        assert result.iloc[1].tolist() == pytest.approx([0.0, 0.0, 1.0])

    def test_ties_are_kept_rather_than_broken_arbitrarily(self):
        """The cut is inclusive: tied values all survive.

        Breaking a tie by position would make the picture depend on column order,
        so a spot split evenly between two types is drawn as two half-wedges even
        when ``top_n`` is 1.
        """
        frame = pd.DataFrame([[0.5, 0.5, 0.0]], columns=list("abc"))
        row = top_n_proportions(frame, 1).iloc[0]
        assert row.tolist() == pytest.approx([0.5, 0.5, 0.0])
        assert row.sum() == pytest.approx(1.0)

    def test_ties_keep_at_least_top_n(self):
        """With equal values the cut is inclusive, so no row loses all its mass."""
        frame = pd.DataFrame([[0.25, 0.25, 0.25, 0.25]], columns=list("abcd"))
        row = top_n_proportions(frame, 2).iloc[0]
        assert row.sum() == pytest.approx(1.0)

    def test_input_frame_is_not_modified(self):
        frame = pd.DataFrame([[0.1, 0.9]], columns=list("ab"))
        before = frame.to_numpy().copy()
        top_n_proportions(frame, 1)
        assert frame.to_numpy().tolist() == before.tolist()

    def test_index_and_columns_are_preserved(self):
        frame = pd.DataFrame([[0.4, 0.6]], index=["spot-1"], columns=["T cell", "B cell"])
        result = top_n_proportions(frame, 2)
        assert list(result.index) == ["spot-1"]
        assert list(result.columns) == ["T cell", "B cell"]


class TestEstimateRadius:
    def test_unit_grid_gives_045(self):
        """Neighbours one unit apart: radius leaves a visible gap."""
        xs, ys = np.meshgrid(np.arange(5.0), np.arange(5.0))
        coordinates = np.column_stack([xs.ravel(), ys.ravel()])
        assert estimate_radius(coordinates) == pytest.approx(0.45)

    def test_scales_linearly_with_spacing(self):
        """A 100x coarser grid must give a 100x larger radius."""
        xs, ys = np.meshgrid(np.arange(5.0), np.arange(5.0))
        fine = np.column_stack([xs.ravel(), ys.ravel()])
        assert estimate_radius(fine * 100.0) == pytest.approx(45.0)

    def test_single_point_falls_back_to_one(self):
        assert estimate_radius(np.array([[0.0, 0.0]])) == 1.0

    def test_empty_input_falls_back_to_one(self):
        assert estimate_radius(np.empty((0, 2))) == 1.0

    def test_all_duplicate_points_fall_back_to_one(self):
        """Zero distances are excluded, so a degenerate layout cannot give radius 0."""
        coordinates = np.zeros((10, 2))
        assert estimate_radius(coordinates) == 1.0

    def test_result_is_always_positive(self):
        rng = np.random.default_rng(0)
        for _ in range(20):
            coordinates = rng.random((30, 2)) * rng.integers(1, 1000)
            assert estimate_radius(coordinates) > 0.0


class TestWaitUntilStable:
    def test_missing_file_returns_false_immediately(self, tmp_path):
        assert wait_until_stable(tmp_path / "absent.h5ad") is False

    def test_settled_file_is_accepted(self, tmp_path):
        path = tmp_path / "sample.h5ad"
        path.write_bytes(b"x" * 4096)
        assert wait_until_stable(path, checks=2, interval=0.01, min_size=1024) is True

    def test_file_below_min_size_is_rejected(self, tmp_path):
        """A freshly created empty file must not be treated as ready."""
        path = tmp_path / "empty.h5ad"
        path.write_bytes(b"")
        assert (
            wait_until_stable(path, checks=2, interval=0.01, min_size=1024, max_wait=0.1) is False
        )

    def test_growing_file_times_out(self, tmp_path):
        """A file still being appended to must not be declared stable."""
        path = tmp_path / "growing.h5ad"
        path.write_bytes(b"x" * 2048)

        original = time.sleep
        state = {"n": 0}

        def grow(seconds):
            state["n"] += 1
            with path.open("ab") as handle:
                handle.write(b"y" * 512)
            original(0)

        time.sleep = grow
        try:
            result = wait_until_stable(path, checks=3, interval=0.01, min_size=1024, max_wait=0.15)
        finally:
            time.sleep = original
        assert result is False
        assert state["n"] > 0

    def test_deleted_midway_returns_false(self, tmp_path):
        path = tmp_path / "vanishing.h5ad"
        path.write_bytes(b"x" * 4096)

        original = time.sleep

        def remove(seconds):
            if path.exists():
                path.unlink()
            original(0)

        time.sleep = remove
        try:
            assert wait_until_stable(path, checks=3, interval=0.01, min_size=1024) is False
        finally:
            time.sleep = original

    def test_default_budget_is_bounded(self, tmp_path, monkeypatch):
        """Without max_wait the helper must still give up rather than hang."""
        path = tmp_path / "stuck.h5ad"
        path.write_bytes(b"x" * 512)  # below min_size, so never stable

        monkeypatch.setattr(time, "sleep", lambda _s: None)
        assert wait_until_stable(path, checks=2, interval=0.001) is False


class TestPathHandling:
    def test_wait_accepts_a_path_object(self, tmp_path):
        path = Path(tmp_path) / "s.h5ad"
        path.write_bytes(b"x" * 2048)
        assert wait_until_stable(path, checks=1, interval=0.01, min_size=10) is True
