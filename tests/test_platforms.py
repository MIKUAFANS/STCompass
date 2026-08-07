"""Tests for the platform registry.

These cover the two behaviours that were wrong in the original scripts: the
substring scan that let ``ST`` match the ``/mnt/cstr/...`` mount point, and the
set iteration whose winner depended on hash randomisation.
"""

from __future__ import annotations

import pytest

from stcompass.platforms import (
    Resolution,
    canonical_platform,
    infer_platform,
    known_platforms,
    resolution_of,
)


class TestCanonicalPlatform:
    @pytest.mark.parametrize(
        ("label", "expected"),
        [
            ("10xVisium", "10xVisium"),
            ("visium", "10xVisium"),
            ("Visium", "10xVisium"),
            ("Visium HD", "10xVisiumHD"),
            ("10xVisiumHD", "10xVisiumHD"),
            ("Stereo Seq", "Stereo-seq"),
            ("stereo-seq", "Stereo-seq"),
            ("StereoSeq", "Stereo-seq"),
            ("Slide seq", "Slide-seq"),
            ("Slide-seq V2", "Slide-seqV2"),
            ("slideseq2", "Slide-seqV2"),
            ("Xenium", "10xXenium"),
            ("MERSCOPE", "MERFISH"),
            ("seqFISH+", "seqFISH+"),
            ("CosMx SMI", "CosMx"),
            ("  merfish  ", "MERFISH"),
        ],
    )
    def test_recognises_spelling_variants(self, label, expected):
        assert canonical_platform(label) == expected

    @pytest.mark.parametrize("label", ["", "Bogus-seq", "scRNA", "unknown"])
    def test_returns_none_for_unknown(self, label):
        assert canonical_platform(label) is None

    def test_separators_are_irrelevant(self):
        variants = ["Well-ST-seq", "Well ST seq", "WellSTseq", "well_st_seq"]
        assert {canonical_platform(v) for v in variants} == {"Well-ST-seq"}


class TestResolutionOf:
    @pytest.mark.parametrize(
        "label", ["10xVisium", "Visium HD", "Stereo Seq", "ST", "Slide-seq", "HDST"]
    )
    def test_spot_platforms(self, label):
        assert resolution_of(label) is Resolution.SPOT

    @pytest.mark.parametrize("label", ["MERFISH", "Xenium", "CosMx", "osmFISH", "STARmap"])
    def test_single_cell_platforms(self, label):
        assert resolution_of(label) is Resolution.SINGLE_CELL

    def test_unknown_returns_default(self):
        assert resolution_of("Bogus") is None
        assert resolution_of("Bogus", Resolution.SPOT) is Resolution.SPOT

    def test_every_registered_platform_has_a_resolution(self):
        for name, resolution in known_platforms().items():
            assert isinstance(resolution, Resolution), name


class TestInferPlatform:
    def test_reads_platform_from_directory(self):
        assert infer_platform("/data/atlas/10xVisium/Homo sapiens/S1.h5ad") == "10xVisium"

    def test_handles_spaced_directory_names(self):
        assert infer_platform("/data/atlas/Stereo Seq/Mus musculus/S2.h5ad") == "Stereo-seq"

    def test_short_label_does_not_match_inside_a_longer_word(self):
        """``ST`` must not match the ``cstr`` mount point.

        This is the original bug: ``"st" in "/mnt/cstr/celldata/..."`` is true, so
        every sample under that mount was classified as spot-based ``ST``.
        """
        assert infer_platform("/mnt/cstr/celldata/liuwenhao/spatial/S3.h5ad") is None

    def test_short_label_matches_a_whole_component(self):
        assert infer_platform("/data/atlas/ST/Homo sapiens/S.h5ad") == "ST"

    def test_longest_match_wins(self):
        """``10xVisiumHD`` must not be truncated to ``10xVisium``."""
        assert infer_platform("/data/atlas/10xVisiumHD/x/S.h5ad") == "10xVisiumHD"
        assert infer_platform("/data/atlas/Slide-seq V2/x/S.h5ad") == "Slide-seqV2"

    def test_deepest_component_wins(self):
        """A nested platform directory overrides one higher up."""
        path = "/data/10xVisium/reprocessed/MERFISH/S.h5ad"
        assert infer_platform(path) == "MERFISH"

    def test_ignores_the_file_stem(self):
        """A sample ID containing a platform name must not leak into the result."""
        assert infer_platform("/data/plain/GSM1234_ST_rep1.h5ad") is None

    def test_relative_to_strips_the_mount_point(self):
        path = "/mnt/cstr/atlas/10xVisium/Homo sapiens/S.h5ad"
        assert infer_platform(path, relative_to="/mnt/cstr/atlas") == "10xVisium"

    def test_relative_to_that_does_not_match_is_ignored(self):
        path = "/data/atlas/MERFISH/x/S.h5ad"
        assert infer_platform(path, relative_to="/somewhere/else") == "MERFISH"

    def test_directory_path_without_suffix(self):
        assert infer_platform("/data/atlas/CosMx/Homo sapiens") == "CosMx"

    def test_is_deterministic(self):
        """Repeated calls must agree; the original iterated an unordered set."""
        path = "/data/atlas/10xVisium/x/S.h5ad"
        assert len({infer_platform(path) for _ in range(50)}) == 1

    def test_unknown_path_returns_none(self):
        assert infer_platform("/data/nothing/here/S.h5ad") is None
