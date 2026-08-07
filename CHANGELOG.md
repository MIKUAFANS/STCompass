# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

本项目的所有重要变更都记录在此文件中。

## [0.1.0] - 2026-08-05

First packaged release. The project previously existed as six standalone scripts
(preserved under `legacy/`); this release turns them into an installable package
with a CLI, a configuration layer and a test suite.

首个打包发布版本。此前项目由六个独立脚本组成（已保留在 `legacy/`），本次发布将其
整理为可安装的包，并提供命令行界面、配置层与测试套件。

### Added

- `stcompass` CLI with five stages — `qc`, `cluster`, `annotate`, `programs`,
  `plot` — plus a `platforms` helper for inspecting the platform registry.
- YAML configuration layer (`stcompass.config`) covering every parameter that
  the original scripts hard-coded. Unknown keys and out-of-range values are
  rejected at load time.
- Platform registry (`stcompass.platforms`) mapping ~21 platform labels, with
  spelling aliases, to a spot / single-cell resolution class.
- `--dry-run` on every stage, listing the samples that would be processed.
- Meaningful exit codes: `0` clean, `1` some samples failed, `2` bad usage.
- Structured logging with an optional `DEBUG`-level log file.
- Test suite (271 tests) covering the platform registry, config loader, matrix
  helpers, batch driver, NMF importance scoring, plotting maths and the CLI.
- GitHub Actions CI running lint, format check and tests on Python 3.10–3.12.
- Example configurations in `configs/`, including `reproduce-atlas.yaml` which
  restores the original scripts' behaviour.

### Fixed

Bugs found in the original scripts while porting them:

- **Platform inference matched substrings of the whole path.** `"ST"` matched
  the mount point `/mnt/cstr/...`, so nearly every sample was classified as
  Spatial Transcriptomics and given the wrong QC thresholds. Matching is now
  per-path-component, deepest-first, with whole-component matches preferred.
- **Platform lookup was order-dependent.** Iterating a `set` of labels meant
  `10xVisiumHD` could lose to `10xVisium`, and the winner varied between
  processes because of hash randomisation. Candidates are now sorted
  longest-first, so the result is deterministic.
- **The count-matrix check drew from the global RNG**, so the same file could be
  classified as raw counts on one run and as normalised on the next — changing
  whether `log1p` was applied. The sample is now drawn from a seeded generator.
- **`np.random.choice(..., replace=False)`** on a matrix's stored values forced a
  full permutation, allocating memory proportional to the number of non-zeros.
  Sampling is now with replacement.
- **Sample lookup re-walked the entire atlas tree once per sample sheet row**,
  which is quadratic. The tree is now indexed once.
- **Failures were logged but never counted**, so a run in which every sample
  errored still exited successfully. Outcomes are now aggregated and reflected
  in the exit code.
- **`adata.write_h5ad` wrote in place**, leaving a truncated file if the process
  died — which a resumed run would then skip as "already done". Writes now go
  to a temporary file and are moved into place.
- **Duplicate-gene removal reordered the gene axis**, because `np.unique` sorts.
  Original order is now preserved.
- **YAML `target_sum: 1e4` silently became the string `"1e4"`** (PyYAML follows
  YAML 1.1, where an unsigned exponent is not a float), surfacing later as a
  `TypeError` from inside a dataclass. Numeric fields are now coerced, with an
  error message that names the key.

### Changed

- The per-barcode and per-gene QC filters, which were commented out in the
  original `QC_ALL0917.py`, are implemented and **enabled by default**. Set
  `qc.filter_cells: false` and `qc.filter_genes: false` (or use
  `configs/reproduce-atlas.yaml`) to reproduce the original output.
- Human and mouse annotation are one code path; species is configuration.
- Heavy dependencies (scanpy, torch, Tangram, SingleR, watchdog) are imported
  lazily and declared as optional extras, so `import stcompass` works in a
  minimal environment and each stage only needs its own dependencies.
- `gene_importance` is vectorised rather than looping per gene; verified to
  match the original loop's output on random inputs.

### Notes

The pipeline stages have been verified only to the extent possible without
scanpy, Tangram and SingleR installed: the pure-Python layers are covered by
tests, and the scanpy-dependent paths are exercised by tests marked
`requires_scanpy`, which skip when the stack is absent. They have not been run
against the original atlas data. See the "Status" section of `README.md`.
