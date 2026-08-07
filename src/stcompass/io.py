"""Discovery, reading and writing of ``.h5ad`` samples.

Two concerns live here.

**Tree mirroring.**  Every pipeline stage reads from one root and writes the same
relative path under another, so ``raw/Homo sapiens/Brain/S1.h5ad`` becomes
``qc/Homo sapiens/Brain/S1.h5ad``.  :func:`iter_samples` walks a tree and pairs
each input with its output, skipping work that already exists -- which is what
makes an interrupted batch resumable.

**Hardened round-tripping.**  Files assembled by other tools break AnnData's
readers in two recurring ways, and both were worked around ad hoc in the original
scripts.  :func:`read_h5ad` and :func:`write_h5ad` handle them centrally:

* ``uns/log1p/base`` written as an HDF5 scalar that the reader cannot decode
  (older scanpy versions wrote ``None`` in a form newer ones reject);
* an ``_index`` column inside ``obs``/``var``, which collides with the reserved
  name AnnData uses for the index itself and makes the file unwritable.
"""

from __future__ import annotations

import shutil
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._deps import require
from .logging_utils import get_logger

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData

__all__ = [
    "SamplePair",
    "iter_samples",
    "load_sample_sheet",
    "mirror_path",
    "read_h5ad",
    "relative_output",
    "sanitise_frames",
    "write_h5ad",
]

logger = get_logger(__name__)

H5AD_SUFFIX = ".h5ad"

# Written by AnnData to mark the index column; a data column with this name makes
# the file unwritable ("Variable name '_index' is reserved").
_RESERVED_INDEX = "_index"


@dataclass(frozen=True, slots=True)
class SamplePair:
    """One input file and the output path it maps to."""

    source: Path
    destination: Path

    @property
    def name(self) -> str:
        """File name without the ``.h5ad`` suffix, used for log messages."""
        return self.source.stem


def relative_output(
    source: Path,
    src_root: Path,
    dst_root: Path,
    suffix: str | None = None,
) -> Path:
    """Map ``source`` under ``src_root`` to the same relative path under ``dst_root``.

    ``suffix`` replaces the file extension, which is how the plotting stage turns
    ``sample.h5ad`` into ``sample.pdf`` while keeping the species/tissue folders.

    >>> relative_output(Path("/a/x/y/s.h5ad"), Path("/a"), Path("/b"), ".pdf").as_posix()
    '/b/x/y/s.pdf'
    """
    relative = source.resolve().relative_to(src_root.resolve())
    if suffix is not None:
        relative = relative.with_suffix(suffix)
    return dst_root / relative


# Kept as an alias because "mirror" is the term used in the docs and CLI help.
mirror_path = relative_output


def iter_samples(
    src_root: Path,
    dst_root: Path,
    *,
    suffix: str | None = None,
    overwrite: bool = False,
    exclude: Iterable[str] = (),
) -> Iterator[SamplePair]:
    """Yield input/output pairs for every ``.h5ad`` file under ``src_root``.

    Args:
        src_root: Directory to walk recursively.
        dst_root: Root of the mirrored output tree.
        suffix: Output extension; ``None`` keeps ``.h5ad``.
        overwrite: When ``False``, pairs whose output exists are skipped so a
            re-run resumes rather than repeating completed samples.
        exclude: File names (with extension) to skip -- an escape hatch for
            samples known to fail, so one bad file does not block a batch.

    Yields:
        :class:`SamplePair` in sorted order, which keeps runs comparable.
    """
    if not src_root.exists():
        raise FileNotFoundError(f"input directory does not exist: {src_root}")

    excluded = set(exclude)
    for source in sorted(src_root.rglob(f"*{H5AD_SUFFIX}")):
        if not source.is_file():
            continue
        if source.name in excluded:
            logger.info("Skipping %s (excluded by configuration)", source.name)
            continue
        destination = relative_output(source, src_root, dst_root, suffix)
        if destination.exists() and not overwrite:
            logger.debug("Skipping %s (output exists)", source.name)
            continue
        yield SamplePair(source=source, destination=destination)


def _drop_broken_log1p_base(path: Path) -> bool:
    """Delete an undecodable ``uns/log1p/base`` entry, in place.

    Returns ``True`` when something was removed.  The entry is metadata recording
    the logarithm base used for a previous transformation and is not needed to
    read the matrix, so removing it is lossless for our purposes.
    """
    h5py = require("h5py", feature="repairing a malformed .h5ad file")
    try:
        with h5py.File(path, "r+") as handle:
            group = handle.get("uns/log1p")
            if group is not None and "base" in group:
                del group["base"]
                return True
    except OSError as exc:
        # A read-only mount or a genuinely corrupt file: report and let the
        # caller's original read error surface.
        logger.debug("Could not repair %s: %s", path, exc)
    return False


def read_h5ad(path: str | Path, *, backed: str | None = None, repair: bool = True) -> AnnData:
    """Read an ``.h5ad`` file, repairing the known-bad ``uns/log1p/base`` entry.

    Args:
        path: File to read.
        backed: Pass ``"r"`` to keep the matrix on disk, for samples too large to
            hold in memory.
        repair: Attempt the in-place repair described above and retry once. Set
            ``False`` to treat any read failure as fatal.

    Raises:
        FileNotFoundError: The path does not exist.
    """
    anndata = require("anndata", feature="reading .h5ad files")
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"no such .h5ad file: {file_path}")

    read = getattr(anndata, "read_h5ad", None) or anndata.io.read_h5ad
    try:
        return read(file_path, backed=backed)
    except Exception as exc:
        if not repair:
            raise
        logger.warning("Failed to read %s (%s); attempting repair", file_path.name, exc)
        if not _drop_broken_log1p_base(file_path):
            raise
        return read(file_path, backed=backed)


def sanitise_frames(adata: AnnData) -> AnnData:
    """Rename reserved ``_index`` columns so the object can be written.

    AnnData refuses to write a frame containing a column literally named
    ``_index``; such columns appear when a file is round-tripped through tools
    that reset the index.  Renamed rather than dropped, because the column holds
    the original barcode or gene identifiers.
    """
    for attr, fallback in (("obs", "cell_index"), ("var", "gene_index")):
        frame = getattr(adata, attr, None)
        if frame is None or _RESERVED_INDEX not in frame.columns:
            continue
        target = fallback
        counter = 1
        while target in frame.columns:
            counter += 1
            target = f"{fallback}_{counter}"
        setattr(adata, attr, frame.rename(columns={_RESERVED_INDEX: target}))
        logger.debug("Renamed reserved '%s' column in .%s to '%s'", _RESERVED_INDEX, attr, target)
    return adata


def _strip_log1p_base(adata: AnnData) -> None:
    """Remove ``uns['log1p']['base']`` when it is ``None``.

    scanpy stores ``base=None`` after ``log1p``, and some writer/reader version
    combinations cannot round-trip that value -- the file writes cleanly but fails
    to load.  Dropping it keeps downstream reads working; the remaining
    ``uns['log1p']`` dict still records that the transform was applied.
    """
    log1p = adata.uns.get("log1p")
    if isinstance(log1p, dict) and log1p.get("base", "missing") is None:
        del log1p["base"]


def write_h5ad(
    adata: AnnData,
    path: str | Path,
    *,
    compression: str | None = "gzip",
    atomic: bool = True,
) -> Path:
    """Write ``adata`` to ``path``, creating parents and sanitising frames first.

    Args:
        adata: Object to write.
        path: Destination file.
        compression: Passed to AnnData; ``gzip`` typically halves file size for
            sparse count matrices.
        atomic: Write to a temporary file in the same directory and move it into
            place on success.  This matters for the watch-mode plotter, which
            would otherwise read a partially written file, and it prevents a
            crash from leaving a truncated output that a resumed run would skip.

    Returns:
        The path written.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    sanitise_frames(adata)
    _strip_log1p_base(adata)

    if not atomic:
        adata.write_h5ad(destination, compression=compression)
        return destination

    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        adata.write_h5ad(temporary, compression=compression)
        # os.replace semantics: atomic within a filesystem, overwrites silently.
        shutil.move(str(temporary), str(destination))
    finally:
        if temporary.exists():
            temporary.unlink(missing_ok=True)
    return destination


def load_sample_sheet(path: str | Path, **read_kwargs: Any) -> Any:
    """Read a sample sheet from ``.xlsx``, ``.xls`` or ``.csv``.

    Args:
        path: Sheet describing the samples to process.
        **read_kwargs: Forwarded to the pandas reader.

    Raises:
        FileNotFoundError: The sheet does not exist.
        ValueError: The extension is not a supported tabular format.
    """
    import pandas as pd

    sheet = Path(path)
    if not sheet.is_file():
        raise FileNotFoundError(f"sample sheet not found: {sheet}")

    suffix = sheet.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        require("openpyxl", feature="reading an Excel sample sheet")
        return pd.read_excel(sheet, **read_kwargs)
    if suffix in {".csv", ".tsv", ".txt"}:
        separator = "\t" if suffix == ".tsv" else read_kwargs.pop("sep", ",")
        return pd.read_csv(sheet, sep=separator, **read_kwargs)
    raise ValueError(
        f"unsupported sample sheet format '{suffix}' ({sheet.name}); use .xlsx or .csv"
    )
