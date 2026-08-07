"""Typed configuration for the pipelines.

Every parameter that the original scripts hard-coded at module level -- input and
output roots, the sample sheet path, count thresholds, the number of NMF
components, GPU counts -- lives here instead, loaded from YAML.  That is what
makes a run reproducible: the config file is a small artefact you can commit
next to the results, whereas an edited-in-place script is not.

Each pipeline has its own dataclass, and :class:`Config` bundles them so one
YAML file can drive an entire atlas build::

    paths:
      raw: /data/atlas/raw
      qc: /data/atlas/qc
    qc:
      min_counts_spot: 100

Unknown keys are rejected rather than ignored, because a silently-dropped
``min_count`` typo would produce plausible-looking but wrong results.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TypeVar

import yaml

__all__ = [
    "AnnotationConfig",
    "ClusteringConfig",
    "Config",
    "ConfigError",
    "PathsConfig",
    "PlotConfig",
    "ProgramsConfig",
    "QCConfig",
    "load_config",
]


class ConfigError(ValueError):
    """The configuration file is malformed or contains unknown keys."""


T = TypeVar("T")


def _expand(path: str | Path) -> Path:
    """Expand ``~`` and environment variables, then resolve to an absolute path."""
    import os

    return Path(os.path.expandvars(str(path))).expanduser()


@dataclass(slots=True)
class PathsConfig:
    """Filesystem roots for one atlas build.

    All pipelines mirror the input directory tree into their output root, so a
    sample at ``raw/Homo sapiens/Brain/S1.h5ad`` becomes
    ``qc/Homo sapiens/Brain/S1.h5ad``.  Keeping that layout means a stage can be
    re-run over a subtree without bookkeeping.
    """

    raw: Path | None = None
    """Input ``.h5ad`` files, as downloaded/converted from the source atlas."""

    qc: Path | None = None
    """Output of ``stcompass qc``; input to clustering and gene programs."""

    clustered: Path | None = None
    """Output of ``stcompass cluster``."""

    annotated: Path | None = None
    """Output of ``stcompass annotate``."""

    programs: Path | None = None
    """Output of ``stcompass programs`` (NMF gene programs)."""

    figures: Path | None = None
    """Output of ``stcompass plot``."""

    reference: Path | None = None
    """scRNA-seq references, laid out as ``reference/<species>/<tissue>/*.h5ad``."""

    metadata: Path | None = None
    """Sample sheet (``.xlsx`` or ``.csv``) describing species/tissue/platform."""

    def __post_init__(self) -> None:
        for f in dataclasses.fields(self):
            value = getattr(self, f.name)
            if value is not None:
                setattr(self, f.name, _expand(value))

    def require(self, name: str) -> Path:
        """Return the path called ``name`` or raise if it was not configured."""
        value = getattr(self, name, None)
        if value is None:
            raise ConfigError(
                f"paths.{name} is required for this command but was not set. "
                f"Add it to the config file or pass the matching command-line option."
            )
        return value


@dataclass(slots=True)
class QCConfig:
    """Quality control thresholds and the embedding that follows them.

    Defaults reproduce the thresholds from the original ``QC_ALL0917.py``:
    spot-based platforms capture several cells per barcode and so tolerate a
    higher floor, while imaging-based panels measure a few hundred genes and need
    a much lower one.
    """

    min_counts_spot: int = 100
    """Minimum total counts per barcode, spot-based platforms."""

    min_genes_spot: int = 30
    """Minimum detected genes per barcode, spot-based platforms."""

    min_counts_single_cell: int = 20
    """Minimum total counts per cell, imaging-based platforms."""

    min_genes_single_cell: int = 10
    """Minimum detected genes per cell, imaging-based platforms."""

    min_cells_per_gene: int = 5
    """Drop genes seen in fewer barcodes.  Spot-based platforms only: on a
    targeted panel every gene is deliberately chosen, so filtering by prevalence
    would discard real signal."""

    min_spots_after_qc: int = 50
    """Reject a spot-based sample left with fewer barcodes than this."""

    min_cells_after_qc: int = 10
    """Reject an imaging-based sample left with fewer cells than this."""

    filter_cells: bool = True
    """Apply the per-barcode count/gene filters.  The source atlas was built with
    these disabled (the thresholds were commented out) to keep every barcode for
    downstream deconvolution; set ``false`` to reproduce that behaviour."""

    filter_genes: bool = True
    """Apply :attr:`min_cells_per_gene` on spot-based platforms."""

    target_sum: float = 1e4
    """Counts per barcode after library-size normalisation."""

    n_top_genes: int = 3000
    """Highly variable genes to keep before PCA.  Samples with fewer genes than
    this use all of them."""

    hvg_flavor: str = "seurat"
    """``flavor`` passed to :func:`scanpy.pp.highly_variable_genes`."""

    n_pcs: int = 50
    """Principal components for the neighbour graph."""

    n_neighbors: int = 15
    """Neighbours for the kNN graph."""

    cluster_method: str = "leiden"
    """``leiden`` or ``louvain``."""

    resolution: float = 1.0
    """Clustering resolution used during QC."""

    compute_umap: bool = True
    """Compute a UMAP embedding.  Disable for very large samples where the
    embedding dominates runtime."""

    integer_check_sample_size: int = 2000
    """Values sampled from the matrix to decide whether it holds raw counts.
    Sampling keeps the check O(1) in matrix size; see
    :func:`stcompass.matrix.looks_like_counts`."""

    def __post_init__(self) -> None:
        if self.cluster_method not in {"leiden", "louvain"}:
            raise ConfigError(
                f"qc.cluster_method must be 'leiden' or 'louvain', got {self.cluster_method!r}"
            )
        if self.n_top_genes < 1:
            raise ConfigError("qc.n_top_genes must be >= 1")
        if self.target_sum <= 0:
            raise ConfigError("qc.target_sum must be > 0")


@dataclass(slots=True)
class ClusteringConfig:
    """Standalone (re-)clustering of QC'd samples.

    The original ``run_louvain.py`` picked a resolution from the number of
    barcodes on a sliding scale: small sections need a coarse graph to yield
    meaningful domains, large ones fragment unless the resolution drops.  That
    table is expressed here as :attr:`resolution_schedule` so it can be tuned
    without touching code.
    """

    method: str = "louvain"
    """``leiden`` or ``louvain``."""

    resolution: float | None = None
    """Fixed resolution.  When ``None``, :attr:`resolution_schedule` decides."""

    resolution_schedule: list[tuple[int, float]] = field(
        default_factory=lambda: [
            (100, 1.2),
            (500, 0.7),
            (5_000, 0.5),
            (20_000, 0.3),
        ]
    )
    """``(max_cells, resolution)`` pairs, checked in order; the first entry whose
    ``max_cells`` exceeds ``n_obs`` wins."""

    fallback_resolution: float = 0.2
    """Resolution for samples larger than every schedule entry."""

    key_added: str | None = None
    """Destination column in ``adata.obs``.  Defaults to :attr:`method`."""

    recompute_neighbors: bool = False
    """Rebuild the kNN graph instead of reusing the one stored during QC."""

    n_pcs: int = 50
    n_neighbors: int = 15

    def __post_init__(self) -> None:
        if self.method not in {"leiden", "louvain"}:
            raise ConfigError(
                f"clustering.method must be 'leiden' or 'louvain', got {self.method!r}"
            )
        # YAML gives lists-of-lists; normalise to tuples and sort by threshold so
        # ``resolution_for`` can scan in order regardless of how it was written.
        self.resolution_schedule = sorted(
            ((int(a), float(b)) for a, b in self.resolution_schedule),
            key=lambda pair: pair[0],
        )

    def resolution_for(self, n_obs: int) -> float:
        """Return the clustering resolution to use for ``n_obs`` barcodes.

        >>> ClusteringConfig().resolution_for(80)
        1.2
        >>> ClusteringConfig().resolution_for(1_000_000)
        0.2
        """
        if self.resolution is not None:
            return self.resolution
        for max_cells, res in self.resolution_schedule:
            if n_obs < max_cells:
                return res
        return self.fallback_resolution


@dataclass(slots=True)
class AnnotationConfig:
    """Cell-type annotation against an scRNA-seq reference.

    Two methods, chosen by platform resolution rather than by the user: Tangram
    maps reference cells onto spots and yields per-spot *proportions*, which is
    the right model for multi-cell barcodes; SingleR assigns one *label* per
    unit, which is right for segmented cells.
    """

    label_key: str = "cell_ontology_class"
    """Column in the reference ``.obs`` holding cell-type labels."""

    min_cells_per_label: int = 2
    """Drop reference labels with fewer cells.  Singleton labels break
    ``rank_genes_groups`` (no within-group variance) and cannot be validated."""

    n_marker_genes: int = 100
    """Top marker genes per reference label, unioned into the Tangram gene set."""

    tangram_mode: str = "cells"
    """``cells``, ``clusters`` or ``constrained``."""

    tangram_epochs: int = 300
    """Training epochs for the mapping."""

    tangram_density_prior: str | None = "rna_count_based"
    """Spot density prior; ``uniform`` suits platforms with equal-area spots."""

    device: str = "auto"
    """``auto``, ``cpu``, ``cuda`` or an explicit ``cuda:N``."""

    n_gpus: int = 1
    """GPUs to spread samples across when ``device`` resolves to CUDA."""

    n_jobs: int = 1
    """Worker processes.  With multiple GPUs, samples are assigned to devices by
    a hash of their path so re-runs are deterministic."""

    singler_threads: int = 4
    """Threads for SingleR."""

    species_column: str = "OrganismSimple"
    tissue_column: str = "Tissue"
    sample_column: str = "SampleName"
    platform_column: str = "Biotech"
    category_column: str | None = "Biotech Categories"
    """Sample-sheet column names.  Defaults match the atlas sheet shipped with
    the original scripts."""

    category_filter: str | None = "Spatial Transcriptomics"
    """Keep only rows whose :attr:`category_column` equals this value."""

    species: list[str] = field(default_factory=list)
    """Restrict the run to these species; empty means all."""

    def __post_init__(self) -> None:
        if self.tangram_mode not in {"cells", "clusters", "constrained"}:
            raise ConfigError(
                "annotation.tangram_mode must be 'cells', 'clusters' or 'constrained', "
                f"got {self.tangram_mode!r}"
            )
        if self.n_gpus < 1:
            raise ConfigError("annotation.n_gpus must be >= 1")
        if self.n_jobs < 1:
            raise ConfigError("annotation.n_jobs must be >= 1")


@dataclass(slots=True)
class ProgramsConfig:
    """Non-negative matrix factorisation into gene programs.

    Uses :class:`sklearn.decomposition.MiniBatchNMF` over row blocks read from a
    backed ``.h5ad``, which is what lets million-cell samples factorise without
    densifying the matrix.
    """

    n_components: int = 10
    """Number of programs (the ``K`` of the factorisation)."""

    max_hvg: int | None = None
    """Restrict the factorisation to this many high-variance genes; ``None`` uses
    every non-empty gene."""

    batch_size: int = 4096
    """Rows per ``partial_fit`` call."""

    epochs: int = 2
    """Passes over the matrix.  One pass leaves the factors under-fitted."""

    row_chunk: int = 5000
    """Rows read per block from the backed file."""

    use_counts_layer: bool = True
    """Prefer ``layers['counts']`` over ``X``, so the factorisation sees raw
    counts even after ``X`` has been log-normalised."""

    random_state: int = 0
    save_float32: bool = True
    """Store factors as ``float32``; halves file size at irrelevant precision cost."""

    init: str = "nndsvda"
    """NMF initialisation; ``nndsvda`` is the sparse-friendly SVD-based default."""

    max_no_improvement: int = 20
    """Stop after this many batches without improvement."""

    large_sample_warning: int = 1_000_000
    """Log a warning above this many barcodes."""

    def __post_init__(self) -> None:
        if self.n_components < 1:
            raise ConfigError("programs.n_components must be >= 1")
        if self.batch_size < 1:
            raise ConfigError("programs.batch_size must be >= 1")
        if self.epochs < 1:
            raise ConfigError("programs.epochs must be >= 1")
        if self.max_hvg is not None and self.max_hvg < self.n_components:
            raise ConfigError("programs.max_hvg must be >= programs.n_components")


@dataclass(slots=True)
class PlotConfig:
    """Rendering of annotation results.

    Tangram proportions are drawn as pie charts at each spot (each wedge a cell
    type); SingleR labels are drawn as a categorical embedding.
    """

    format: str = "pdf"
    """Output file extension, e.g. ``pdf``, ``png``, ``svg``."""

    dpi: int = 300
    figsize: tuple[float, float] = (8.0, 8.0)

    top_n_types: int = 3
    """Cell types per spot in a scatter-pie.  Showing every type turns each pie
    into noise; the top few carry the signal."""

    proportions_key: str = "tangram_ct_pred"
    """``adata.obsm`` key holding Tangram proportions."""

    label_key: str = "singler_best"
    """``adata.obs`` column holding SingleR labels."""

    spatial_key: str = "spatial"
    """``adata.obsm`` key holding coordinates."""

    radius: float | None = None
    """Pie radius in data units.  ``None`` derives it from nearest-neighbour
    spacing so the same code works for 55 µm Visium spots and 2 µm HD bins."""

    invert_y: bool = False
    """Flip the y-axis, for platforms whose coordinates are image-style."""

    max_pies: int = 20_000
    """Above this many spots, fall back to a plain scatter: rendering hundreds of
    thousands of wedges exhausts memory and produces an unreadable figure."""

    watch: bool = False
    """Keep running and render files as they appear."""

    poll_seconds: float = 5.0
    """Polling interval in watch mode.  Polling runs even when watchdog is
    available, because inotify does not fire for writes on NFS mounts."""

    stability_checks: int = 3
    stability_interval: float = 2.0
    """Consecutive equal-size checks (and the gap between them) before a new file
    is considered fully written."""

    def __post_init__(self) -> None:
        if isinstance(self.figsize, list):
            self.figsize = tuple(float(v) for v in self.figsize)  # type: ignore[assignment]
        if len(self.figsize) != 2:
            raise ConfigError("plot.figsize must have exactly two entries")
        if self.top_n_types < 1:
            raise ConfigError("plot.top_n_types must be >= 1")


@dataclass(slots=True)
class Config:
    """Top-level configuration: paths plus one section per pipeline."""

    paths: PathsConfig = field(default_factory=PathsConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
    programs: ProgramsConfig = field(default_factory=ProgramsConfig)
    plot: PlotConfig = field(default_factory=PlotConfig)

    n_jobs: int = 1
    """Default worker count for pipelines that do not set their own."""

    overwrite: bool = False
    """Re-process samples whose output already exists.  Off by default so an
    interrupted batch resumes instead of redoing days of work."""

    log_file: Path | None = None
    """Optional path for a DEBUG-level log of the run."""

    def __post_init__(self) -> None:
        if self.log_file is not None:
            self.log_file = _expand(self.log_file)
        if self.n_jobs < 1:
            raise ConfigError("n_jobs must be >= 1")

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> Config:
        """Build a :class:`Config` from nested plain data, rejecting stray keys."""
        if not isinstance(data, Mapping):
            raise ConfigError(f"configuration must be a mapping, got {type(data).__name__}")

        sections: dict[str, type] = {
            "paths": PathsConfig,
            "qc": QCConfig,
            "clustering": ClusteringConfig,
            "annotation": AnnotationConfig,
            "programs": ProgramsConfig,
            "plot": PlotConfig,
        }
        kwargs: dict[str, Any] = {}
        scalar_names = {f.name for f in dataclasses.fields(cls) if f.name not in sections}

        unknown = set(data) - set(sections) - scalar_names
        if unknown:
            allowed = sorted(set(sections) | scalar_names)
            raise ConfigError(
                f"unknown top-level configuration key(s): {sorted(unknown)}. Allowed: {allowed}"
            )

        for name, section_cls in sections.items():
            raw = data.get(name)
            if raw is None:
                continue
            if not isinstance(raw, Mapping):
                raise ConfigError(f"'{name}' section must be a mapping")
            kwargs[name] = _build_section(section_cls, raw, name)

        for name in scalar_names:
            if name in data:
                kwargs[name] = data[name]

        return cls(**kwargs)


def _numeric_kind(annotation: str) -> type[int] | type[float] | None:
    """Return ``int``/``float`` when ``annotation`` denotes a plain numeric field.

    ``from __future__ import annotations`` means :func:`dataclasses.fields` hands
    back annotations as strings, so this inspects the text rather than the type.
    Container and union-of-containers annotations return ``None`` -- only scalars
    are coerced.
    """
    base = annotation.replace("| None", "").replace("None |", "").strip()
    if base == "int":
        return int
    if base == "float":
        return float
    return None


def _coerce_scalar(
    value: Any,
    kind: type[int] | type[float],
    section_name: str,
    key: str,
) -> Any:
    """Coerce a YAML scalar to ``kind``, or raise a :class:`ConfigError`.

    YAML 1.1 -- which PyYAML implements -- only reads an exponent as a float when
    the exponent carries an explicit sign, so ``target_sum: 1e4`` arrives as the
    string ``"1e4"``.  Without this, that value reached a dataclass comparison and
    surfaced as ``TypeError: '<=' not supported between 'str' and 'int'``, which
    says nothing about which config key was wrong.
    """
    if value is None:
        # Fields such as ``max_hvg`` and ``radius`` are genuinely optional; the
        # dataclass declares them ``| None`` and validates them itself.
        return value
    if isinstance(value, bool):
        # bool is a subclass of int, so `n_top_genes: true` would otherwise pass
        # every numeric check as 1 and silently produce a one-gene analysis.
        raise ConfigError(f"{section_name}.{key} must be a number, got the boolean {value!r}")
    if isinstance(value, kind):
        return value
    if isinstance(value, (int, float)):
        return kind(value)
    if isinstance(value, str):
        try:
            return kind(float(value.strip()))
        except ValueError:
            raise ConfigError(
                f"{section_name}.{key} must be a number, got {value!r}. "
                f"Note that YAML needs a signed exponent: write 1.0e+4 or 10000."
            ) from None
    raise ConfigError(f"{section_name}.{key} must be a number, got {type(value).__name__}")


def _build_section(section_cls: type[T], raw: Mapping[str, Any], section_name: str) -> T:
    """Instantiate a section dataclass, reporting unknown keys with their section."""
    fields = {f.name: f for f in dataclasses.fields(section_cls)}  # type: ignore[arg-type]
    unknown = set(raw) - set(fields)
    if unknown:
        raise ConfigError(
            f"unknown key(s) in '{section_name}' section: {sorted(unknown)}. "
            f"Allowed: {sorted(fields)}"
        )

    prepared: dict[str, Any] = {}
    for key, value in raw.items():
        kind = _numeric_kind(str(fields[key].type))
        prepared[key] = _coerce_scalar(value, kind, section_name, key) if kind else value

    try:
        return section_cls(**prepared)  # type: ignore[call-arg]
    except ConfigError:
        raise
    except (TypeError, ValueError) as exc:
        # A field-level validation problem; name the section so the user knows
        # where to look.
        raise ConfigError(f"invalid '{section_name}' section: {exc}") from exc


def load_config(path: str | Path | None = None, **overrides: Any) -> Config:
    """Load configuration from ``path``, applying ``overrides`` on top.

    ``overrides`` uses dotted keys so command-line flags can target a nested
    field: ``load_config("atlas.yaml", **{"qc.resolution": 0.5})``.  ``None``
    values are ignored, which lets the CLI pass every flag unconditionally
    without unset ones clobbering the file.
    """
    data: dict[str, Any] = {}
    if path is not None:
        config_path = Path(path).expanduser()
        if not config_path.is_file():
            raise ConfigError(f"configuration file not found: {config_path}")
        with config_path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
        if loaded is None:
            loaded = {}
        if not isinstance(loaded, Mapping):
            raise ConfigError(f"{config_path} must contain a YAML mapping at the top level")
        data = dict(loaded)

    for dotted, value in overrides.items():
        if value is None:
            continue
        target = data
        parts = dotted.split(".")
        for part in parts[:-1]:
            existing = target.get(part)
            if existing is None:
                existing = {}
                target[part] = existing
            elif not isinstance(existing, dict):
                raise ConfigError(f"cannot set {dotted}: '{part}' is not a mapping")
            target = existing
        target[parts[-1]] = value

    return Config.from_mapping(data)
