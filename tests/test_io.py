"""Tests for sample discovery, path mirroring and sample-sheet loading."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from stcompass.io import (
    SamplePair,
    iter_samples,
    load_sample_sheet,
    relative_output,
    sanitise_frames,
)


class TestRelativeOutput:
    def test_mirrors_the_relative_path(self, tmp_path):
        source = tmp_path / "raw" / "10xVisium" / "Homo sapiens" / "S1.h5ad"
        source.parent.mkdir(parents=True)
        source.touch()
        out = relative_output(source, tmp_path / "raw", tmp_path / "qc")
        assert out == tmp_path / "qc" / "10xVisium" / "Homo sapiens" / "S1.h5ad"

    def test_replaces_the_suffix(self, tmp_path):
        source = tmp_path / "raw" / "x" / "S1.h5ad"
        source.parent.mkdir(parents=True)
        source.touch()
        out = relative_output(source, tmp_path / "raw", tmp_path / "fig", ".pdf")
        assert out == tmp_path / "fig" / "x" / "S1.pdf"

    def test_rejects_a_path_outside_the_root(self, tmp_path):
        (tmp_path / "raw").mkdir()
        outside = tmp_path / "elsewhere" / "S1.h5ad"
        outside.parent.mkdir(parents=True)
        outside.touch()
        with pytest.raises(ValueError):
            relative_output(outside, tmp_path / "raw", tmp_path / "qc")


class TestIterSamples:
    def test_finds_every_h5ad_and_ignores_other_files(self, tmp_atlas, tmp_path):
        pairs = list(iter_samples(tmp_atlas, tmp_path / "qc"))
        assert len(pairs) == 5
        assert all(p.source.suffix == ".h5ad" for p in pairs)

    def test_output_paths_mirror_the_input_tree(self, tmp_atlas, tmp_path):
        pairs = list(iter_samples(tmp_atlas, tmp_path / "qc"))
        by_name = {p.source.name: p for p in pairs}
        assert by_name["S1.h5ad"].destination == (
            tmp_path / "qc" / "10xVisium" / "Homo sapiens" / "S1.h5ad"
        )

    def test_skips_samples_whose_output_exists(self, tmp_atlas, tmp_path):
        destination = tmp_path / "qc" / "10xVisium" / "Homo sapiens" / "S1.h5ad"
        destination.parent.mkdir(parents=True)
        destination.touch()
        names = {p.source.name for p in iter_samples(tmp_atlas, tmp_path / "qc")}
        assert "S1.h5ad" not in names
        assert "S2.h5ad" in names

    def test_overwrite_reprocesses_existing_outputs(self, tmp_atlas, tmp_path):
        destination = tmp_path / "qc" / "10xVisium" / "Homo sapiens" / "S1.h5ad"
        destination.parent.mkdir(parents=True)
        destination.touch()
        names = {p.source.name for p in iter_samples(tmp_atlas, tmp_path / "qc", overwrite=True)}
        assert "S1.h5ad" in names

    def test_exclude_skips_named_files(self, tmp_atlas, tmp_path):
        pairs = iter_samples(tmp_atlas, tmp_path / "qc", exclude=["S1.h5ad", "S4.h5ad"])
        names = {p.source.name for p in pairs}
        assert names == {"S2.h5ad", "S3.h5ad", "S5.h5ad"}

    def test_suffix_applies_to_every_pair(self, tmp_atlas, tmp_path):
        pairs = list(iter_samples(tmp_atlas, tmp_path / "fig", suffix=".pdf"))
        assert all(p.destination.suffix == ".pdf" for p in pairs)

    def test_missing_input_root_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="does not exist"):
            list(iter_samples(tmp_path / "nope", tmp_path / "out"))

    def test_order_is_deterministic(self, tmp_atlas, tmp_path):
        first = [p.source for p in iter_samples(tmp_atlas, tmp_path / "a")]
        second = [p.source for p in iter_samples(tmp_atlas, tmp_path / "b")]
        assert first == second


class TestSamplePair:
    def test_name_strips_the_suffix(self):
        pair = SamplePair(source=Path("/a/S1.h5ad"), destination=Path("/b/S1.h5ad"))
        assert pair.name == "S1"


class TestSanitiseFrames:
    """``_index`` is reserved by AnnData; a data column with that name blocks writing."""

    def test_renames_a_reserved_obs_column(self):
        adata = _FakeAnnData(
            obs=pd.DataFrame({"_index": ["a", "b"], "keep": [1, 2]}),
            var=pd.DataFrame({"gene": ["g1"]}),
        )
        sanitise_frames(adata)
        assert "_index" not in adata.obs.columns
        assert list(adata.obs["cell_index"]) == ["a", "b"]
        assert "keep" in adata.obs.columns

    def test_renames_a_reserved_var_column(self):
        adata = _FakeAnnData(
            obs=pd.DataFrame({"x": [1]}),
            var=pd.DataFrame({"_index": ["g1", "g2"]}),
        )
        sanitise_frames(adata)
        assert list(adata.var["gene_index"]) == ["g1", "g2"]

    def test_avoids_colliding_with_an_existing_column(self):
        adata = _FakeAnnData(
            obs=pd.DataFrame({"_index": ["a"], "cell_index": ["taken"]}),
            var=pd.DataFrame({"g": [1]}),
        )
        sanitise_frames(adata)
        assert "cell_index_2" in adata.obs.columns
        assert list(adata.obs["cell_index"]) == ["taken"]

    def test_leaves_clean_frames_untouched(self):
        adata = _FakeAnnData(
            obs=pd.DataFrame({"a": [1]}),
            var=pd.DataFrame({"b": [2]}),
        )
        sanitise_frames(adata)
        assert list(adata.obs.columns) == ["a"]
        assert list(adata.var.columns) == ["b"]


class _FakeAnnData:
    """Minimal stand-in exposing the ``obs``/``var`` attributes sanitise_frames uses."""

    def __init__(self, obs, var):
        self.obs = obs
        self.var = var


class TestLoadSampleSheet:
    def test_reads_a_csv(self, tmp_path):
        path = tmp_path / "sheet.csv"
        path.write_text("SampleName,Tissue\nS1,Brain\n", encoding="utf-8")
        frame = load_sample_sheet(path)
        assert list(frame["SampleName"]) == ["S1"]

    def test_reads_a_tsv_using_tab_separation(self, tmp_path):
        path = tmp_path / "sheet.tsv"
        path.write_text("SampleName\tTissue\nS1\tBrain\n", encoding="utf-8")
        frame = load_sample_sheet(path)
        assert list(frame.columns) == ["SampleName", "Tissue"]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="sample sheet not found"):
            load_sample_sheet(tmp_path / "absent.csv")

    def test_unsupported_extension_names_the_file(self, tmp_path):
        path = tmp_path / "sheet.parquet"
        path.touch()
        with pytest.raises(ValueError, match="unsupported sample sheet format"):
            load_sample_sheet(path)

    def test_reads_an_xlsx_when_openpyxl_is_available(self, tmp_path):
        pytest.importorskip("openpyxl")
        path = tmp_path / "sheet.xlsx"
        pd.DataFrame({"SampleName": ["S1"], "Tissue": ["Brain"]}).to_excel(path, index=False)
        frame = load_sample_sheet(path)
        assert list(frame["Tissue"]) == ["Brain"]
