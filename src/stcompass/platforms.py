"""Registry of spatial transcriptomics platforms.

Every pipeline in this package branches on one question: does a measurement unit
hold one cell or many?  Spot-based platforms (Visium, Slide-seq, Stereo-seq)
capture a mixture of cells per barcode, so they need permissive count filters and
*deconvolution* to estimate cell-type proportions.  Imaging-based platforms
(MERFISH, Xenium, CosMx) segment individual cells but probe only a few hundred
genes, so they need low count thresholds and *classification* against a
reference.

The original scripts each carried their own copy of these platform lists, and the
copies had drifted apart -- ``Stereo Seq`` appeared in one and ``Stereo-seq`` in
another, and the QC script listed platforms the annotation script did not know.
Centralising them here makes the resolution rules explicit and testable.
"""

from __future__ import annotations

from enum import Enum

__all__ = [
    "Resolution",
    "canonical_platform",
    "infer_platform",
    "known_platforms",
    "resolution_of",
]


class Resolution(str, Enum):
    """Whether one measurement unit corresponds to one cell."""

    SPOT = "spot"
    """Multi-cell capture area (Visium, Slide-seq, Stereo-seq, ...)."""

    SINGLE_CELL = "single_cell"
    """Segmented single cells (MERFISH, Xenium, CosMx, ...)."""


# Canonical platform name -> resolution class.
_PLATFORMS: dict[str, Resolution] = {
    # --- spot-based / barcoded arrays -------------------------------------
    "10xVisium": Resolution.SPOT,
    "10xVisiumHD": Resolution.SPOT,
    "ST": Resolution.SPOT,
    "Slide-seq": Resolution.SPOT,
    "Slide-seqV2": Resolution.SPOT,
    "Stereo-seq": Resolution.SPOT,
    "HDST": Resolution.SPOT,
    "Well-ST-seq": Resolution.SPOT,
    "CBSST-seq": Resolution.SPOT,
    "sci-Space": Resolution.SPOT,
    "Pixel-seq": Resolution.SPOT,
    # --- imaging-based / single-cell resolution ---------------------------
    "MERFISH": Resolution.SINGLE_CELL,
    "seqFISH": Resolution.SINGLE_CELL,
    "seqFISH+": Resolution.SINGLE_CELL,
    "osmFISH": Resolution.SINGLE_CELL,
    "EASI-FISH": Resolution.SINGLE_CELL,
    "EEL-FISH": Resolution.SINGLE_CELL,
    "STARmap": Resolution.SINGLE_CELL,
    "10xXenium": Resolution.SINGLE_CELL,
    "CosMx": Resolution.SINGLE_CELL,
    "ExSeq": Resolution.SINGLE_CELL,
}

# Spelling variants seen in the sample sheets and directory names of the source
# atlas.  Keys are normalised (lower-case, separators removed) by ``_normalise``.
_ALIASES: dict[str, str] = {
    "visium": "10xVisium",
    "10xvisium": "10xVisium",
    "visiumhd": "10xVisiumHD",
    "10xvisiumhd": "10xVisiumHD",
    "spatialtranscriptomics": "ST",
    "slideseq": "Slide-seq",
    "slideseqv2": "Slide-seqV2",
    "slideseq2": "Slide-seqV2",
    "stereoseq": "Stereo-seq",
    "cbsstseq": "CBSST-seq",
    "wellstseq": "Well-ST-seq",
    "scispace": "sci-Space",
    "pixelseq": "Pixel-seq",
    "xenium": "10xXenium",
    "10xxenium": "10xXenium",
    "cosmx": "CosMx",
    "cosmxsmi": "CosMx",
    "merfish": "MERFISH",
    "merscope": "MERFISH",
    "seqfish": "seqFISH",
    "seqfishplus": "seqFISH+",
    "osmfish": "osmFISH",
    "easifish": "EASI-FISH",
    "eelfish": "EEL-FISH",
    "starmap": "STARmap",
    "exseq": "ExSeq",
    "hdst": "HDST",
}


def _normalise(name: str) -> str:
    """Strip case, whitespace and separators so ``Stereo Seq`` == ``stereo-seq``."""
    out = []
    for ch in name.strip().lower():
        if ch.isalnum():
            out.append(ch)
        elif ch == "+":
            out.append("plus")
        # '-', '_', ' ', '.' and friends are dropped
    return "".join(out)


_NORMALISED: dict[str, str] = {_normalise(name): name for name in _PLATFORMS}
_NORMALISED.update({_normalise(alias): target for alias, target in _ALIASES.items()})


def known_platforms() -> dict[str, Resolution]:
    """Return a copy of the canonical platform table."""
    return dict(_PLATFORMS)


def canonical_platform(name: str) -> str | None:
    """Map a free-form platform label to its canonical spelling.

    Returns ``None`` when the label is not recognised; callers decide whether an
    unknown platform is fatal or merely falls back to a default.

    >>> canonical_platform("stereo seq")
    'Stereo-seq'
    >>> canonical_platform("Slide-seq V2")
    'Slide-seqV2'
    >>> canonical_platform("Bogus-seq") is None
    True
    """
    if not name:
        return None
    return _NORMALISED.get(_normalise(name))


def resolution_of(name: str, default: Resolution | None = None) -> Resolution | None:
    """Return the :class:`Resolution` for a platform label.

    ``default`` is returned for unrecognised labels, which keeps batch jobs
    running over a heterogeneous atlas instead of aborting on one odd folder.

    >>> resolution_of("Visium")
    <Resolution.SPOT: 'spot'>
    >>> resolution_of("Xenium")
    <Resolution.SINGLE_CELL: 'single_cell'>
    """
    canonical = canonical_platform(name)
    if canonical is None:
        return default
    return _PLATFORMS[canonical]


# Labels short enough to appear inside unrelated words ("ST" sits inside "fastq",
# "cstr", "Stereo").  They are only ever accepted as a whole path component.
_MIN_SUBSTRING_LENGTH = 4

# Longest-first so "10xVisiumHD" is tried before "10xVisium" and "Slide-seqV2"
# before "Slide-seq"; ties broken alphabetically to keep the order deterministic
# across processes.  The original code iterated a set union, whose order varies
# with per-process hash randomisation, so the same file could be classified
# differently on a re-run.
_SUBSTRING_CANDIDATES: list[tuple[str, str]] = sorted(
    ((normalised, target) for normalised, target in _NORMALISED.items()),
    key=lambda pair: (-len(pair[0]), pair[0]),
)


def infer_platform(path: object, *, relative_to: object = None) -> str | None:
    """Infer a platform from the directory names of ``path``.

    Atlases are commonly laid out as ``<root>/<platform>/<species>/<sample>.h5ad``,
    so the platform can be recovered from the path when no sample sheet is
    available.

    Matching is deliberately stricter than the substring scan it replaces, which
    tested ``label.lower() in str(full_path).lower()`` against the *entire*
    absolute path.  That gave two wrong answers in practice: the label ``ST``
    matched the mount point ``/mnt/cstr/...`` and so claimed almost every sample,
    and iteration over a set of labels meant the winner depended on hash
    randomisation.  Here each path component is matched on its own, the deepest
    component wins (it is the most specific), whole-component matches are
    preferred over substrings, and short labels must match a whole component.

    Args:
        path: File or directory path to inspect.
        relative_to: Optional root to strip first, so components above the atlas
            root (``/mnt``, ``/home/user``) cannot contribute a match.

    Returns:
        Canonical platform name, or ``None`` when nothing matches.

    >>> infer_platform("/data/atlas/10xVisium/Homo sapiens/S1.h5ad")
    '10xVisium'
    >>> infer_platform("/data/atlas/Stereo Seq/Mus musculus/S2.h5ad")
    'Stereo-seq'
    >>> infer_platform("/mnt/cstr/celldata/spatial/S3.h5ad") is None
    True
    """
    import contextlib
    from pathlib import Path, PurePath

    candidate = path if isinstance(path, PurePath) else Path(str(path))
    if relative_to is not None:
        root = relative_to if isinstance(relative_to, PurePath) else Path(str(relative_to))
        # Not under root: fall back to inspecting the whole path.
        with contextlib.suppress(ValueError):
            candidate = candidate.relative_to(root)

    # Directory names only: the file's own stem is a sample ID, and sample IDs
    # such as "GSM1234_ST_rep1" would otherwise leak a spurious platform.
    components = list(candidate.parts[:-1]) if candidate.suffix else list(candidate.parts)

    # Deepest first: the directory nearest the file is the most specific.
    for component in reversed(components):
        normalised = _normalise(component)
        if not normalised:
            continue
        exact = _NORMALISED.get(normalised)
        if exact is not None:
            return exact

    for component in reversed(components):
        normalised = _normalise(component)
        if not normalised:
            continue
        for label, target in _SUBSTRING_CANDIDATES:
            if len(label) >= _MIN_SUBSTRING_LENGTH and label in normalised:
                return target

    return None
