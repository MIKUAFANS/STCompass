# STCompass

Reproducible batch pipelines for large-scale spatial transcriptomics atlases:
quality control, clustering, cell-type annotation, NMF gene programs, and figures.

面向大规模空间转录组图谱的可复现批处理流程：质控、聚类、细胞类型注释、NMF 基因程序、可视化。

[![CI](https://github.com/MIKUAFANS/STCompass/actions/workflows/ci.yml/badge.svg)](https://github.com/MIKUAFANS/STCompass/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## What this is / 这是什么

STCompass processes a directory tree of `.h5ad` samples through five stages. Each
stage is a subcommand, mirrors its input directory tree into its output root, and
skips samples that already have an output — so a run over thousands of samples can
be interrupted and resumed without bookkeeping.

STCompass 把一整棵 `.h5ad` 样本目录送过五个阶段。每个阶段都是一个子命令，会把输入目录结构
镜像到输出根目录，并自动跳过已产出的样本 —— 因此处理数千样本的任务可以中断后继续，无需额外记账。

| Stage / 阶段 | Command | Reads | Writes | Replaces |
|---|---|---|---|---|
| Quality control / 质控 | `stcompass qc` | `paths.raw` | `paths.qc` | `QC_ALL0917.py` |
| Clustering / 聚类 | `stcompass cluster` | `paths.qc` | `paths.clustered` | `run_louvain.py` |
| Annotation / 注释 | `stcompass annotate` | `paths.raw` + reference | `paths.annotated` | `cell_annotation_*.py` |
| Gene programs / 基因程序 | `stcompass programs` | `paths.qc` | `paths.programs` | `spatialVEG_0923.py` |
| Figures / 出图 | `stcompass plot` | `paths.annotated` | `paths.figures` | `celltype_image1015.py` |

---

## Install / 安装

```bash
pip install -e .                  # core: scanpy, anndata, sklearn, matplotlib
pip install -e '.[tangram]'       # + Tangram deconvolution (needs CUDA torch)
pip install -e '.[singler]'       # + SingleR classification
pip install -e '.[excel]'         # + .xlsx sample sheets
pip install -e '.[watch]'         # + filesystem watching for `plot --watch`
pip install -e '.[dev]'           # + pytest, ruff
```

Heavy dependencies are imported lazily, so `import stcompass` works without
Tangram, SingleR or torch installed; a stage tells you exactly which extra to
install if you invoke it without its dependency.

重依赖是惰性导入的，因此没有装 Tangram / SingleR / torch 也能 `import stcompass`；
当你调用缺少依赖的阶段时，程序会明确告诉你该装哪个 extra。

Details, including the GPU and R-backend notes: [`docs/installation.md`](docs/installation.md).

---

## Quick start / 快速开始

```bash
# 1. Copy a config and edit the paths
cp configs/example.yaml atlas.yaml

# 2. See what would happen — no files are touched
stcompass qc --config atlas.yaml --dry-run

# 3. Run the stages
stcompass qc       --config atlas.yaml
stcompass annotate --config atlas.yaml
stcompass plot     --config atlas.yaml
```

Every setting can also be given as a flag, which wins over the file:

每个配置项也都能用命令行参数给出，且优先级高于配置文件：

```bash
stcompass cluster --qc-dir data/qc --out data/clustered --method leiden -j 8
```

Check how a platform label in your sample sheet will be interpreted:

检查样本表里的平台名会被如何解析：

```bash
stcompass platforms --check 'Stereo Seq' --check 'Visium HD'
# 'Stereo Seq': Stereo-seq (spot)
# 'Visium HD': 10xVisiumHD (spot)
```

Full walkthrough: [`docs/usage.md`](docs/usage.md).

---

## Why the platform matters / 为什么平台类型是关键

Every stage branches on one question: **does one measurement unit hold one cell,
or many?**

每个阶段都取决于同一个问题：**一个测量单元里是一个细胞，还是很多个？**

- **Spot-based** (Visium, Slide-seq, Stereo-seq): a barcode pools several cells.
  Needs higher count thresholds, and annotation is a *deconvolution* problem —
  Tangram produces per-spot cell-type **proportions**.
  **点阵型**：一个 barcode 混合多个细胞，需要更高的计数阈值；注释是*解卷积*问题 ——
  Tangram 输出每个 spot 的细胞类型**比例**。
- **Imaging-based** (MERFISH, Xenium, CosMx): cells are segmented individually but
  only a few hundred genes are probed. Needs low thresholds, and annotation is a
  *classification* problem — SingleR assigns one **label** per cell.
  **成像型**：细胞被单独分割，但只探测几百个基因，需要较低阈值；注释是*分类*问题 ——
  SingleR 给每个细胞一个**标签**。

The registry in `stcompass/platforms.py` maps 21 platform names (plus spelling
variants such as `Stereo Seq` / `stereo-seq`) onto these two classes, so the
method is chosen from the data rather than by hand.

`stcompass/platforms.py` 中的注册表把 21 个平台名（以及 `Stereo Seq` / `stereo-seq`
这类拼写变体）映射到这两个类别，因此方法由数据决定，而不是手工指定。

---

## What each stage writes / 各阶段写入的字段

| Stage | Key | Contents |
|---|---|---|
| `qc` | `layers["counts"]` | raw counts, preserved before normalisation |
| `qc` | `obs["n_counts"]`, `obs["n_genes"]`, `obs["percent_mt"]` | per-barcode QC metrics |
| `qc` | `obs["leiden"]` / `obs["louvain"]`, `obsm["X_pca"]`, `obsm["X_umap"]` | embedding and clusters |
| `cluster` | `obs["louvain"]` (or `clustering.key_added`) | re-clustered labels |
| `annotate` | `obsm["tangram_ct_pred"]` | spot × cell-type proportions (spot-based) |
| `annotate` | `obs["singler_*"]`, incl. `obs["singler_best"]` | per-cell labels (imaging-based) |
| `programs` | `obsm["X_nmf"]` | barcode × program loadings (*W*) |
| `programs` | `var["NNMF_component_{k}"]` | gene weights (*H*) |
| `programs` | `var["NNMF_component_{k}_importance_score"]`, `..._rank` | program-specificity score and rank |
| `plot` | — | one figure per sample, mirroring the input tree |

Stage-by-stage detail, including the algorithms and their parameters:
[`docs/pipelines.md`](docs/pipelines.md).

---

## Project layout / 目录结构

```
STCompass/
├── src/stcompass/
│   ├── cli.py             # five subcommands, exit codes, --dry-run
│   ├── config.py          # typed YAML config; unknown keys are rejected
│   ├── platforms.py       # platform registry -> spot vs single-cell
│   ├── io.py              # tree mirroring, hardened .h5ad read/write
│   ├── matrix.py          # block-wise sparse helpers (no scanpy import)
│   ├── logging_utils.py   # console + file logging
│   ├── _deps.py           # lazy optional imports with install hints
│   └── pipelines/
│       ├── _batch.py      # per-sample loop, error capture, exit summary
│       ├── qc.py  clustering.py  annotation.py  programs.py  visualize.py
├── tests/                 # 271 tests; scanpy-dependent ones self-skip
├── configs/               # example.yaml, minimal.yaml, reproduce-atlas.yaml
├── docs/                  # installation, usage, configuration, pipelines, migration
├── pyproject.toml         # packaging, ruff, pytest, coverage
└── CHANGELOG.md  CONTRIBUTING.md  LICENSE
```

---

## Exit codes / 退出码

Meaningful, so the commands compose in a shell script or workflow engine:

| Code | Meaning |
|---|---|
| `0` | every sample processed or deliberately skipped |
| `1` | at least one sample failed |
| `2` | bad usage or configuration; nothing ran |

A *skip* (too few cells after QC, no spatial coordinates) is not a failure — it is
an expected outcome for a heterogeneous atlas and is reported separately in the
run summary.

*跳过*（质控后细胞太少、没有空间坐标）不算失败 —— 对异质图谱来说这是预期结果，
会在运行摘要中单独统计。

---

## Development / 开发

```bash
pip install -e '.[dev]'
ruff check src tests
ruff format --check src tests
pytest
```

271 tests cover the platform registry, config loading and rejection, matrix
helpers, the batch driver, NMF importance scoring, plotting maths and the CLI.
Tests needing scanpy are marked `requires_scanpy` and skip themselves when it is
absent, so the suite is green on a laptop and thorough on a full stack.

271 个测试覆盖平台注册表、配置加载与拒绝、矩阵工具、批处理驱动、NMF 重要性打分、
绘图数学和 CLI。需要 scanpy 的测试标记为 `requires_scanpy`，在缺少依赖时自动跳过 ——
因此在笔记本上是全绿，在完整环境里是完整覆盖。

See [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License / 许可

[MIT](LICENSE).

Third-party methods retain their own licences and citation requirements — notably
[Tangram](https://github.com/broadinstitute/Tangram),
[SingleR](https://github.com/SingleR-inc/singler-py) and
[scanpy](https://github.com/scverse/scanpy). If you publish results produced with
this package, cite those methods.

第三方方法保留其各自的许可与引用要求 —— 特别是 Tangram、SingleR 和 scanpy。
若使用本包产出的结果发表论文，请引用这些方法。
