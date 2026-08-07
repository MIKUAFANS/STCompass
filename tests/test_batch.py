"""Tests for the shared batch driver.

The behaviour that matters here is failure isolation: a batch over thousands of
samples must not abort because one file is corrupt, and it must still report
honestly that the failure happened.  The original scripts got the first half right
(a bare ``except`` around the loop body) and the second half wrong -- nothing
counted failures, so a run where every sample errored looked identical to a clean
one.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from stcompass.io import SamplePair
from stcompass.pipelines._batch import BatchResult, SampleOutcome, run_batch


def _pair(name: str) -> SamplePair:
    return SamplePair(source=Path(f"/in/{name}.h5ad"), destination=Path(f"/out/{name}.h5ad"))


# ---------------------------------------------------------------------------
# BatchResult bookkeeping
# ---------------------------------------------------------------------------


class TestBatchResult:
    def test_empty_result_is_ok(self):
        assert BatchResult().ok is True

    def test_empty_result_has_zero_total(self):
        assert BatchResult().total == 0

    def test_failure_makes_result_not_ok(self):
        result = BatchResult(failed=[SampleOutcome(Path("a.h5ad"), "failed", "boom")])
        assert result.ok is False

    def test_skips_do_not_make_result_not_ok(self):
        """A skip is an expected outcome, not an error."""
        result = BatchResult(skipped=[SampleOutcome(Path("a.h5ad"), "skipped", "too few cells")])
        assert result.ok is True

    def test_total_counts_all_three_categories(self):
        result = BatchResult(
            processed=[Path("a")],
            skipped=[SampleOutcome(Path("b"), "skipped", "")],
            failed=[SampleOutcome(Path("c"), "failed", "")],
        )
        assert result.total == 3

    def test_summary_mentions_every_count(self):
        result = BatchResult(
            processed=[Path("a"), Path("b")],
            skipped=[SampleOutcome(Path("c"), "skipped", "")],
            failed=[],
        )
        summary = result.summary()
        assert "2 processed" in summary
        assert "1 skipped" in summary
        assert "0 failed" in summary

    def test_merge_folds_in_other_result(self):
        first = BatchResult(processed=[Path("a")])
        second = BatchResult(
            processed=[Path("b")],
            failed=[SampleOutcome(Path("c"), "failed", "boom")],
        )
        first.merge(second)
        assert len(first.processed) == 2
        assert len(first.failed) == 1
        assert first.ok is False


# ---------------------------------------------------------------------------
# run_batch
# ---------------------------------------------------------------------------


class TestRunBatch:
    def test_empty_input_returns_empty_result(self):
        result = run_batch([], lambda pair: None, progress=False)
        assert result.total == 0
        assert result.ok is True

    def test_worker_returning_none_counts_as_processed(self):
        result = run_batch([_pair("a"), _pair("b")], lambda pair: None, progress=False)
        assert len(result.processed) == 2
        assert result.ok is True

    def test_worker_returning_string_counts_as_skipped(self):
        result = run_batch([_pair("a")], lambda pair: "too few cells", progress=False)
        assert len(result.skipped) == 1
        assert result.skipped[0].detail == "too few cells"

    def test_worker_raising_counts_as_failed(self):
        def worker(pair):
            raise ValueError("corrupt file")

        result = run_batch([_pair("a")], worker, progress=False)
        assert len(result.failed) == 1
        assert result.ok is False

    def test_failure_message_includes_exception_type_and_text(self):
        def worker(pair):
            raise KeyError("missing_column")

        result = run_batch([_pair("a")], worker, progress=False)
        assert "KeyError" in result.failed[0].detail
        assert "missing_column" in result.failed[0].detail

    def test_one_failure_does_not_stop_the_batch(self):
        """The whole point of the driver: keep going past a bad sample."""
        seen = []

        def worker(pair):
            seen.append(pair.source.name)
            if pair.source.name == "b.h5ad":
                raise RuntimeError("boom")
            return

        pairs = [_pair("a"), _pair("b"), _pair("c")]
        result = run_batch(pairs, worker, progress=False)

        assert seen == ["a.h5ad", "b.h5ad", "c.h5ad"]
        assert len(result.processed) == 2
        assert len(result.failed) == 1

    def test_mixed_outcomes_are_categorised_correctly(self):
        def worker(pair):
            name = pair.source.name
            if name == "a.h5ad":
                return None
            if name == "b.h5ad":
                return "skip reason"
            raise ValueError("fail reason")

        result = run_batch([_pair("a"), _pair("b"), _pair("c")], worker, progress=False)
        assert [p.name for p in result.processed] == ["a.h5ad"]
        assert len(result.skipped) == 1
        assert len(result.failed) == 1
        assert result.total == 3

    def test_empty_string_from_worker_counts_as_processed(self):
        """Only a non-empty reason marks a skip; '' is falsy and means success."""
        result = run_batch([_pair("a")], lambda pair: "", progress=False)
        assert len(result.processed) == 1
        assert len(result.skipped) == 0

    def test_keyboard_interrupt_propagates(self):
        """Ctrl-C must stop the run, not be recorded as a per-sample failure."""

        def worker(pair):
            raise KeyboardInterrupt

        with pytest.raises(KeyboardInterrupt):
            run_batch([_pair("a")], worker, progress=False)

    def test_source_paths_are_recorded_for_processed_samples(self):
        result = run_batch([_pair("x")], lambda pair: None, progress=False)
        assert result.processed[0] == Path("/in/x.h5ad")

    def test_worker_receives_the_pair(self):
        received = []
        run_batch([_pair("a")], lambda pair: received.append(pair) or None, progress=False)
        assert received[0].destination == Path("/out/a.h5ad")
