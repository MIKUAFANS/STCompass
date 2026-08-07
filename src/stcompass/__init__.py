"""STCompass -- reproducible batch pipelines for spatial transcriptomics atlases.

The package turns a directory of ``.h5ad`` samples into an annotated atlas through
five stages, each available as a subcommand of the ``stcompass`` CLI and as a
function in :mod:`stcompass.pipelines`:

===========  ==========================================================
``qc``       Filter barcodes/genes, normalise, embed and cluster.
``cluster``  Re-cluster QC'd samples at a size-aware resolution.
``annotate`` Label cell types against an scRNA-seq reference.
``programs`` Factorise expression into NMF gene programs.
``plot``     Render annotation results to figures.
===========  ==========================================================

Each stage mirrors its input directory tree into its output root and skips samples
that already have an output, so a run over thousands of samples can be interrupted
and resumed.

Heavy dependencies (scanpy, torch, Tangram, SingleR) are imported lazily, so
``import stcompass`` works in a minimal environment and only the stage you invoke
needs its extras installed.

Example:
    >>> from stcompass import load_config
    >>> from stcompass.pipelines import run_qc
    >>> config = load_config("configs/example.yaml")   # doctest: +SKIP
    >>> run_qc(config)                                 # doctest: +SKIP
"""

from __future__ import annotations

from ._deps import MissingDependencyError
from .config import (
    AnnotationConfig,
    ClusteringConfig,
    Config,
    ConfigError,
    PathsConfig,
    PlotConfig,
    ProgramsConfig,
    QCConfig,
    load_config,
)
from .logging_utils import configure as configure_logging
from .platforms import Resolution, canonical_platform, known_platforms, resolution_of

__version__ = "0.1.0"

__all__ = [
    "AnnotationConfig",
    "ClusteringConfig",
    "Config",
    "ConfigError",
    "MissingDependencyError",
    "PathsConfig",
    "PlotConfig",
    "ProgramsConfig",
    "QCConfig",
    "Resolution",
    "__version__",
    "canonical_platform",
    "configure_logging",
    "known_platforms",
    "load_config",
    "resolution_of",
]
