"""Command-line interface.

One entry point, five subcommands -- one per pipeline stage::

    stcompass qc       --config atlas.yaml
    stcompass cluster  --config atlas.yaml
    stcompass annotate --config atlas.yaml
    stcompass programs --config atlas.yaml
    stcompass plot     --config atlas.yaml

Configuration comes from a YAML file; every path and tuning knob can also be
overridden on the command line, with the flag winning over the file.  That split
exists because the two are used differently: the YAML file is the reproducible
artefact you commit next to the results, while the flags are for the one-off
re-run over a subdirectory.

Exit codes are meaningful, so the commands compose in a shell script or a
workflow engine:

====  ==========================================================
``0``  every sample processed or deliberately skipped
``1``  at least one sample failed
``2``  bad usage or configuration; nothing ran
====  ==========================================================
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Callable, Sequence
from typing import Any

from . import __version__
from ._deps import MissingDependencyError
from .config import Config, ConfigError, load_config
from .logging_utils import configure as configure_logging
from .logging_utils import get_logger
from .platforms import known_platforms

__all__ = ["build_parser", "main"]

logger = get_logger(__name__)

EXIT_OK = 0
EXIT_FAILED_SAMPLES = 1
EXIT_USAGE = 2


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def _add_common(parser: argparse.ArgumentParser) -> None:
    """Options accepted by every subcommand."""
    parser.add_argument(
        "-c",
        "--config",
        metavar="FILE",
        help="YAML configuration file. Every setting can also be given as a flag below.",
    )
    parser.add_argument(
        "-j",
        "--n-jobs",
        type=int,
        metavar="N",
        help="Worker processes (default: 1, or the value in the config file).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=None,
        help="Re-process samples whose output already exists. Off by default, so an "
        "interrupted run resumes instead of repeating finished work.",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=None,
        metavar="FILENAME",
        help="Skip this file name; repeatable. Use for samples known to fail.",
    )
    parser.add_argument(
        "--log-file",
        metavar="FILE",
        help="Also write a DEBUG-level log here.",
    )
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print DEBUG-level progress to stderr.",
    )
    verbosity.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Only print warnings and errors.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve the configuration, list the samples that would be processed, "
        "then stop without reading or writing any data.",
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="stcompass",
        description=(
            "Reproducible batch pipelines for spatial transcriptomics atlases. "
            "Each stage mirrors its input directory tree into its output root and "
            "skips samples that already have an output, so runs are resumable."
        ),
        epilog=(
            "Examples:\n"
            "  stcompass qc --config configs/example.yaml\n"
            "  stcompass cluster --qc-dir data/qc --out data/clustered --method leiden\n"
            "  stcompass annotate --config atlas.yaml --species 'Homo sapiens'\n"
            "  stcompass plot --annotated-dir data/annotated --out figures --watch\n"
            "\nSee docs/usage.md for the full pipeline walkthrough."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"stcompass {__version__}")

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    # --- qc ---------------------------------------------------------------
    qc = subparsers.add_parser(
        "qc",
        help="Filter, normalise, embed and cluster raw samples.",
        description=(
            "Quality control. Thresholds follow each sample's platform: spot-based "
            "platforms pool several cells per barcode and tolerate a higher floor, "
            "imaging-based panels measure few genes and need a lower one. Raw counts "
            "are preserved in layers['counts'] for the programs stage."
        ),
    )
    _add_common(qc)
    qc.add_argument("--raw-dir", metavar="DIR", help="Input root (paths.raw).")
    qc.add_argument("--out", metavar="DIR", help="Output root (paths.qc).")
    qc.add_argument("--min-counts-spot", type=int, metavar="N")
    qc.add_argument("--min-genes-spot", type=int, metavar="N")
    qc.add_argument("--min-counts-single-cell", type=int, metavar="N")
    qc.add_argument("--min-genes-single-cell", type=int, metavar="N")
    qc.add_argument(
        "--no-filter-cells",
        dest="filter_cells",
        action="store_false",
        default=None,
        help="Keep every barcode. Reproduces the published atlas, which was built "
        "with the per-barcode thresholds disabled.",
    )
    qc.add_argument(
        "--no-filter-genes",
        dest="filter_genes",
        action="store_false",
        default=None,
        help="Keep every gene regardless of how few barcodes detect it.",
    )
    qc.add_argument("--n-top-genes", type=int, metavar="N", help="Highly variable genes for PCA.")
    qc.add_argument("--resolution", type=float, metavar="R", help="Clustering resolution.")
    qc.add_argument(
        "--cluster-method",
        choices=["leiden", "louvain"],
        help="Community detection algorithm.",
    )
    qc.add_argument(
        "--no-umap",
        dest="compute_umap",
        action="store_false",
        default=None,
        help="Skip UMAP; useful when the embedding dominates runtime.",
    )

    # --- cluster ----------------------------------------------------------
    cluster = subparsers.add_parser(
        "cluster",
        help="Re-cluster QC'd samples at a size-aware resolution.",
        description=(
            "Standalone clustering. Without --resolution, the resolution is chosen "
            "from the number of barcodes: small sections need a coarse graph to "
            "yield interpretable domains, large ones fragment unless it drops."
        ),
    )
    _add_common(cluster)
    cluster.add_argument("--qc-dir", metavar="DIR", help="Input root (paths.qc).")
    cluster.add_argument("--out", metavar="DIR", help="Output root (paths.clustered).")
    cluster.add_argument("--method", choices=["leiden", "louvain"])
    cluster.add_argument(
        "--resolution",
        type=float,
        metavar="R",
        help="Fixed resolution for every sample, overriding the size-aware schedule.",
    )
    cluster.add_argument("--key-added", metavar="NAME", help="Destination obs column.")
    cluster.add_argument(
        "--recompute-neighbors",
        action="store_true",
        default=None,
        help="Rebuild the kNN graph instead of reusing the one stored during QC.",
    )

    # --- annotate ---------------------------------------------------------
    annotate = subparsers.add_parser(
        "annotate",
        help="Label cell types against an scRNA-seq reference.",
        description=(
            "Cell-type annotation. The method follows the platform: spot-based data "
            "is deconvolved with Tangram into per-spot proportions "
            "(obsm['tangram_ct_pred']), single-cell data is classified with SingleR "
            "into labels (obs['singler_best']). Requires a sample sheet describing "
            "each sample's species, tissue and platform."
        ),
    )
    _add_common(annotate)
    annotate.add_argument("--raw-dir", metavar="DIR", help="Input root (paths.raw).")
    annotate.add_argument("--out", metavar="DIR", help="Output root (paths.annotated).")
    annotate.add_argument(
        "--reference-dir",
        metavar="DIR",
        help="Reference root, laid out as <dir>/<species>/<tissue>/*.h5ad.",
    )
    annotate.add_argument(
        "--metadata",
        metavar="FILE",
        help="Sample sheet (.xlsx or .csv) describing species, tissue and platform.",
    )
    annotate.add_argument(
        "--species",
        action="append",
        metavar="NAME",
        help="Restrict to this species; repeatable. Default: every species in the sheet.",
    )
    annotate.add_argument("--label-key", metavar="COLUMN", help="Reference cell-type column.")
    annotate.add_argument("--n-marker-genes", type=int, metavar="N")
    annotate.add_argument("--tangram-epochs", type=int, metavar="N")
    annotate.add_argument(
        "--device",
        metavar="DEV",
        help="'auto', 'cpu', 'cuda' or 'cuda:N'. Falls back to CPU when CUDA is absent.",
    )
    annotate.add_argument(
        "--n-gpus",
        type=int,
        metavar="N",
        help="Spread samples across this many GPUs, assigned by a hash of the path so "
        "a resumed run reuses the same placement.",
    )
    annotate.add_argument("--singler-threads", type=int, metavar="N")

    # --- programs ---------------------------------------------------------
    programs = subparsers.add_parser(
        "programs",
        help="Factorise expression into NMF gene programs.",
        description=(
            "Gene programs via MiniBatchNMF over row blocks read from a backed "
            ".h5ad, so million-cell samples factorise without densifying the matrix. "
            "Writes loadings to obsm['X_nmf'] and per-gene weights, importance "
            "scores and ranks to var."
        ),
    )
    _add_common(programs)
    programs.add_argument("--qc-dir", metavar="DIR", help="Input root (paths.qc).")
    programs.add_argument("--out", metavar="DIR", help="Output root (paths.programs).")
    programs.add_argument("--n-components", type=int, metavar="K", help="Number of programs.")
    programs.add_argument(
        "--max-hvg",
        type=int,
        metavar="N",
        help="Restrict the factorisation to this many high-variance genes.",
    )
    programs.add_argument("--batch-size", type=int, metavar="N")
    programs.add_argument("--epochs", type=int, metavar="N", help="Passes over the matrix.")
    programs.add_argument("--row-chunk", type=int, metavar="N", help="Rows read per block.")
    programs.add_argument("--random-state", type=int, metavar="SEED")
    programs.add_argument(
        "--use-x",
        dest="use_counts_layer",
        action="store_false",
        default=None,
        help="Factorise X instead of preferring layers['counts'].",
    )

    # --- plot -------------------------------------------------------------
    plot = subparsers.add_parser(
        "plot",
        help="Render annotation results to figures.",
        description=(
            "Plot annotated samples. Tangram proportions are drawn as per-spot pie "
            "charts, SingleR labels as a categorical embedding. With --watch the "
            "command keeps running and renders files as an annotation job produces "
            "them."
        ),
    )
    _add_common(plot)
    plot.add_argument("--annotated-dir", metavar="DIR", help="Input root (paths.annotated).")
    plot.add_argument("--out", metavar="DIR", help="Output root (paths.figures).")
    plot.add_argument("--format", metavar="EXT", help="Output extension: pdf, png, svg.")
    plot.add_argument("--dpi", type=int, metavar="N")
    plot.add_argument(
        "--top-n-types",
        type=int,
        metavar="N",
        help="Cell types per pie; the rest are dropped and the remainder renormalised.",
    )
    plot.add_argument(
        "--max-pies",
        type=int,
        metavar="N",
        help="Above this many spots, fall back to a plain scatter of the dominant type.",
    )
    plot.add_argument(
        "--radius",
        type=float,
        metavar="R",
        help="Pie radius in data units. Default: derived from spot spacing.",
    )
    plot.add_argument(
        "--invert-y",
        action="store_true",
        default=None,
        help="Flip the y-axis for image-style coordinates.",
    )
    plot.add_argument(
        "--watch",
        action="store_true",
        default=None,
        help="Keep running and render new files as they appear.",
    )
    plot.add_argument("--poll-seconds", type=float, metavar="S")

    # --- platforms --------------------------------------------------------
    platforms_cmd = subparsers.add_parser(
        "platforms",
        help="List the recognised platforms and their resolution class.",
        description=(
            "Print the platform registry. Useful for checking how a label in your "
            "sample sheet will be interpreted, since the resolution class decides "
            "both the QC thresholds and the annotation method."
        ),
    )
    platforms_cmd.add_argument(
        "--check",
        metavar="LABEL",
        action="append",
        help="Resolve this label instead of listing everything; repeatable.",
    )

    return parser


# ---------------------------------------------------------------------------
# Option -> config mapping
# ---------------------------------------------------------------------------

# Flag name -> dotted config key.  Kept as data so a new option is one line here
# and cannot drift from the parser.
_OVERRIDES: dict[str, dict[str, str]] = {
    "qc": {
        "raw_dir": "paths.raw",
        "out": "paths.qc",
        "min_counts_spot": "qc.min_counts_spot",
        "min_genes_spot": "qc.min_genes_spot",
        "min_counts_single_cell": "qc.min_counts_single_cell",
        "min_genes_single_cell": "qc.min_genes_single_cell",
        "filter_cells": "qc.filter_cells",
        "filter_genes": "qc.filter_genes",
        "n_top_genes": "qc.n_top_genes",
        "resolution": "qc.resolution",
        "cluster_method": "qc.cluster_method",
        "compute_umap": "qc.compute_umap",
    },
    "cluster": {
        "qc_dir": "paths.qc",
        "out": "paths.clustered",
        "method": "clustering.method",
        "resolution": "clustering.resolution",
        "key_added": "clustering.key_added",
        "recompute_neighbors": "clustering.recompute_neighbors",
    },
    "annotate": {
        "raw_dir": "paths.raw",
        "out": "paths.annotated",
        "reference_dir": "paths.reference",
        "metadata": "paths.metadata",
        "species": "annotation.species",
        "label_key": "annotation.label_key",
        "n_marker_genes": "annotation.n_marker_genes",
        "tangram_epochs": "annotation.tangram_epochs",
        "device": "annotation.device",
        "n_gpus": "annotation.n_gpus",
        "singler_threads": "annotation.singler_threads",
    },
    "programs": {
        "qc_dir": "paths.qc",
        "out": "paths.programs",
        "n_components": "programs.n_components",
        "max_hvg": "programs.max_hvg",
        "batch_size": "programs.batch_size",
        "epochs": "programs.epochs",
        "row_chunk": "programs.row_chunk",
        "random_state": "programs.random_state",
        "use_counts_layer": "programs.use_counts_layer",
    },
    "plot": {
        "annotated_dir": "paths.annotated",
        "out": "paths.figures",
        "format": "plot.format",
        "dpi": "plot.dpi",
        "top_n_types": "plot.top_n_types",
        "max_pies": "plot.max_pies",
        "radius": "plot.radius",
        "invert_y": "plot.invert_y",
        "watch": "plot.watch",
        "poll_seconds": "plot.poll_seconds",
    },
}

# Which paths each command reads, so --dry-run can report the sample list without
# importing scanpy.
_STAGE_PATHS: dict[str, tuple[str, str]] = {
    "qc": ("raw", "qc"),
    "cluster": ("qc", "clustered"),
    "annotate": ("raw", "annotated"),
    "programs": ("qc", "programs"),
    "plot": ("annotated", "figures"),
}


def _config_from_args(args: argparse.Namespace) -> Config:
    """Build a :class:`Config` from the config file plus command-line overrides."""
    overrides: dict[str, Any] = {}
    for flag, dotted in _OVERRIDES.get(args.command, {}).items():
        value = getattr(args, flag, None)
        if value is not None:
            overrides[dotted] = value

    # ``annotation.n_jobs`` shadows the global one so a GPU stage can run fewer
    # workers than the CPU stages in the same config file.
    if args.n_jobs is not None:
        overrides["n_jobs"] = args.n_jobs
        if args.command == "annotate":
            overrides["annotation.n_jobs"] = args.n_jobs
    if args.overwrite:
        overrides["overwrite"] = True
    if args.log_file is not None:
        overrides["log_file"] = args.log_file

    return load_config(args.config, **overrides)


def _console_level(args: argparse.Namespace) -> int:
    if getattr(args, "verbose", False):
        return logging.DEBUG
    if getattr(args, "quiet", False):
        return logging.WARNING
    return logging.INFO


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


def _cmd_platforms(args: argparse.Namespace) -> int:
    """Print the platform registry, or resolve specific labels."""
    from .platforms import canonical_platform, resolution_of

    if args.check:
        exit_code = EXIT_OK
        for label in args.check:
            canonical = canonical_platform(label)
            if canonical is None:
                print(f"{label!r}: not recognised")
                exit_code = EXIT_FAILED_SAMPLES
            else:
                print(f"{label!r}: {canonical} ({resolution_of(label).value})")
        return exit_code

    table = known_platforms()
    width = max(len(name) for name in table)
    headings = {
        "spot": "Spot-based (deconvolution)",
        "single_cell": "Single-cell (classification)",
    }
    for resolution in ("spot", "single_cell"):
        names = sorted(n for n, r in table.items() if r.value == resolution)
        print(f"\n{headings[resolution]}:")
        for name in names:
            print(f"  {name.ljust(width)}  {resolution}")
    print()
    return EXIT_OK


def _dry_run(config: Config, command: str, exclude: tuple[str, ...] = ()) -> int:
    """List the samples a stage would process, without touching them.

    ``exclude`` is applied here for the same reason it is applied in the real run:
    a preview that ignored it would list samples the stage will not touch, which
    defeats the point of checking before committing to a long job.
    """
    from .io import iter_samples

    src_name, dst_name = _STAGE_PATHS[command]
    try:
        src_root = config.paths.require(src_name)
        dst_root = config.paths.require(dst_name)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE

    print(f"stage:  {command}")
    print(f"input:  {src_root}")
    print(f"output: {dst_root}")

    if command == "annotate":
        # Annotation is driven by the sample sheet, not by a tree walk, so listing
        # the pairs here would misrepresent what the stage does.
        for name in ("reference", "metadata"):
            try:
                print(f"{name}: {config.paths.require(name)}")
            except ConfigError as exc:
                print(f"error: {exc}", file=sys.stderr)
                return EXIT_USAGE
        print(
            "\n(annotate selects samples from the sample sheet; run without "
            "--dry-run to process them)"
        )
        return EXIT_OK

    suffix = f".{config.plot.format}" if command == "plot" else None
    try:
        pairs = list(
            iter_samples(
                src_root,
                dst_root,
                suffix=suffix,
                overwrite=config.overwrite,
                exclude=exclude,
            )
        )
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE

    print(f"\n{len(pairs)} sample(s) would be processed:")
    for pair in pairs[:50]:
        print(f"  {pair.source}  ->  {pair.destination}")
    if len(pairs) > 50:
        print(f"  ... and {len(pairs) - 50} more")
    return EXIT_OK


def _run_stage(config: Config, command: str, exclude: tuple[str, ...]) -> int:
    """Dispatch to a pipeline and turn its result into an exit code."""
    from . import pipelines

    runners: dict[str, Callable[..., Any]] = {
        "qc": pipelines.run_qc,
        "cluster": pipelines.run_clustering,
        "annotate": pipelines.run_annotation,
        "programs": pipelines.run_programs,
        "plot": pipelines.run_plot,
    }
    runner = runners[command]

    # run_annotation takes its sample list from the sheet, so it has no --exclude.
    if command == "annotate":
        if exclude:
            logger.warning("--exclude is not supported by 'annotate'; ignoring")
        result = runner(config)
    else:
        result = runner(config, exclude=exclude)

    print(f"\n{command}: {result.summary()}")
    if result.failed:
        print(f"\n{len(result.failed)} sample(s) failed:", file=sys.stderr)
        for outcome in result.failed[:20]:
            print(f"  {outcome.sample.name}: {outcome.detail}", file=sys.stderr)
        if len(result.failed) > 20:
            print(f"  ... and {len(result.failed) - 20} more", file=sys.stderr)
        if config.log_file:
            print(f"\nFull log: {config.log_file}", file=sys.stderr)
        return EXIT_FAILED_SAMPLES
    return EXIT_OK


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``stcompass`` command.

    Returns a process exit code rather than calling :func:`sys.exit`, so the
    function is directly testable.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return EXIT_USAGE

    if args.command == "platforms":
        return _cmd_platforms(args)

    try:
        config = _config_from_args(args)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE

    configure_logging(level=_console_level(args), log_file=config.log_file, force=True)

    exclude = tuple(args.exclude or ())

    if args.dry_run:
        return _dry_run(config, args.command, exclude)

    try:
        return _run_stage(config, args.command, exclude)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except MissingDependencyError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
