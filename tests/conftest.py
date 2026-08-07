"""Shared fixtures.

The test suite is split in two by dependency weight.  Most tests exercise the
pure-Python layers -- the platform registry, the config loader, the matrix
helpers, the batch driver, the NMF importance scoring -- and run in any
environment that has numpy, scipy and pandas.

Tests that need ``scanpy`` (or anndata's writer) are marked ``requires_scanpy``
and skip themselves when it is absent, so ``pytest`` is green on a laptop and
thorough on a machine with the full stack.
"""

from __future__ import annotations

import importlib.util

import pytest


def _installed(module: str) -> bool:
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):  # pragma: no cover - malformed install
        return False


HAS_SCANPY = _installed("scanpy")
HAS_ANNDATA = _installed("anndata")


def pytest_collection_modifyitems(config, items):
    """Skip ``requires_scanpy`` tests when scanpy is not installed."""
    if HAS_SCANPY:
        return
    skip = pytest.mark.skip(reason="scanpy is not installed")
    for item in items:
        if "requires_scanpy" in item.keywords:
            item.add_marker(skip)


@pytest.fixture
def tmp_atlas(tmp_path):
    """A small input tree shaped like the atlas: ``<platform>/<species>/<sample>.h5ad``.

    Files are empty placeholders; tests that only exercise discovery and path
    mirroring do not need real HDF5 content.
    """
    layout = {
        "10xVisium/Homo sapiens": ["S1.h5ad", "S2.h5ad"],
        "10xVisium/Mus musculus": ["S3.h5ad"],
        "MERFISH/Homo sapiens": ["S4.h5ad"],
        "Stereo Seq/Mus musculus": ["S5.h5ad"],
    }
    root = tmp_path / "raw"
    for relative, names in layout.items():
        directory = root / relative
        directory.mkdir(parents=True, exist_ok=True)
        for name in names:
            (directory / name).write_bytes(b"")
    # A non-.h5ad file, to prove discovery filters by extension.
    (root / "10xVisium" / "README.txt").write_text("not a sample", encoding="utf-8")
    return root
