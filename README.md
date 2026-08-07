# STCompassDB

**Reproducible batch pipelines for large-scale spatial transcriptomics atlases** —
quality control, clustering, cell-type annotation, NMF gene programs and figure
generation, driven by a single declarative configuration.

**面向大规模空间转录组图谱的可复现批处理流程** —— 质控、聚类、细胞类型注释、
NMF 基因程序与出图，由单一声明式配置文件驱动。

![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)
![Tests](https://img.shields.io/badge/tests-271-brightgreen.svg)
![Lint](https://img.shields.io/badge/lint-ruff-261230.svg)
![Status](https://img.shields.io/badge/status-beta-orange.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Naming.** `STCompassDB` is the project name. The Python distribution, the
> importable package and the console script remain `stcompass` — `pip install -e .`,
> `import stcompass`, `stcompass qc ...`.
>
> **命名说明.** `STCompassDB` 为项目名称；Python 发行包名、导入名与命令行入口
> 仍为 `stcompass`。

---

## Abstract / 摘要

Spatial transcriptomics atlases aggregate samples from assays whose measurement
units differ by three orders of magnitude in physical size and by two in
gene-panel breadth, so no single parameter set can serve an entire atlas.
STCompassDB makes one property the discriminant from which thresholds, annotation
method and figure type are all derived — *does one measurement unit hold one cell,
or a mixture?* — and resolves that property from a registry of 21 platforms and 27
spelling variants rather than from per-sample manual curation. The five stages are
file-level idempotent: each mirrors its input directory tree into its output root
and skips samples whose output already exists, so an atlas-scale run is
interruptible and resumable with no external work queue or database. Optional
heavy dependencies (torch, Tangram, SingleR) are resolved at call time, and 271
unit tests pin the numerical behaviour of every scoring and thresholding rule so
that a published result cannot change silently between versions.

空间转录组图谱汇集了测量单元物理尺寸相差三个数量级、基因 panel 广度相差两个数量级的
多种技术平台，因此单一参数集无法服务整个图谱。STCompassDB 将**一个测量单元包含单个
细胞还是细胞混合物**作为唯一判别依据，由它推导阈值、注释方法与图形类型；该属性通过
包含 21 个平台与 27 个拼写变体的注册表解析，而非逐样本人工标注。五个阶段在文件级别
幂等：各阶段将输入目录树镜像到输出根目录，并跳过已有输出的样本，因此图谱规模的运行
可中断续跑，无需外部任务队列或数据库。重依赖在调用时惰性解析；271 个单元测试锁定所有
打分与阈值规则的数值行为，确保已发表结果不会在版本间静默改变。

---

## Contents

1. [System overview / 系统总览](#1-system-overview--系统总览)
2. [Design rationale / 设计依据](#2-design-rationale--设计依据)
3. [Installation / 安装](#3-installation--安装)
4. [Quick start / 快速开始](#4-quick-start--快速开始)
5. [Methods / 方法](#5-methods--方法)
6. [Output schema / 输出字段规范](#6-output-schema--输出字段规范)
7. [Architecture / 软件架构](#7-architecture--软件架构)
8. [Scalability / 可扩展性](#8-scalability--可扩展性)
9. [Reproducibility / 可复现性](#9-reproducibility--可复现性)
10. [Testing / 测试](#10-testing--测试)
11. [Limitations / 已知局限](#11-limitations--已知局限)
12. [Citation / 引用](#12-citation--引用)

---

## 1. System overview / 系统总览

The atlas build is a directed acyclic graph over five stages, each a subcommand of
the `stcompass` CLI and an equivalent function in `stcompass.pipelines`:

```
                     ┌──────────────┐
   paths.raw ───────▶│      qc      │──▶ paths.qc ──┬──▶ cluster ──▶ paths.clustered
        │            └──────────────┘               │
        │                                           └──▶ programs ─▶ paths.programs
        │            ┌──────────────┐
        └───────────▶│   annotate   │──▶ paths.annotated ──▶ plot ──▶ paths.figures
                     └──────────────┘
                            ▲
        paths.reference ────┤
        paths.metadata  ────┘
```

Note that `annotate` reads `paths.raw`, not `paths.qc`: deconvolution benefits from
every barcode, including those a QC threshold would discard.

注意 `annotate` 读取 `paths.raw` 而非 `paths.qc`：解卷积受益于全部 barcode，
包括会被质控阈值滤除的那些。

| Stage / 阶段 | Command | Reads | Writes | Parallelism |
|---|---|---|---|---|
| Quality control / 质控 | `stcompass qc` | `paths.raw` | `paths.qc` | `n_jobs` |
| Clustering / 聚类 | `stcompass cluster` | `paths.qc` | `paths.clustered` | `n_jobs` |
| Annotation / 注释 | `stcompass annotate` | `paths.raw` + `reference` + `metadata` | `paths.annotated` | `annotation.n_jobs` × `n_gpus` |
| Gene programs / 基因程序 | `stcompass programs` | `paths.qc` | `paths.programs` | `n_jobs` |
| Figures / 出图 | `stcompass plot` | `paths.annotated` | `paths.figures` | `n_jobs` |

Exit codes are stage-independent, so the commands compose in a shell script or a
workflow engine:

| Code | Meaning |
|---|---|
| `0` | every sample processed or deliberately skipped |
| `1` | at least one sample failed |
| `2` | bad usage or configuration; nothing ran |

A *skip* — too few cells after QC, no spatial coordinates, no reference for a
tissue — is an expected outcome for a heterogeneous atlas, not a failure, and is
counted separately in the run summary.

*跳过*（质控后细胞过少、缺少空间坐标、组织无匹配参考）对异质图谱而言是预期结果，
不计为失败，在运行摘要中单独统计。

---

## 2. Design rationale / 设计依据

Every stage branches on one question: **does one measurement unit hold one cell, or
many?** The answer determines count thresholds, whether prevalence filtering is
meaningful, whether annotation is a regression or a classification problem, and
which figure is honest.

每个阶段都取决于同一个问题：**一个测量单元里是一个细胞，还是很多个？** 答案决定
计数阈值、prevalence 过滤是否有意义、注释是回归还是分类问题，以及何种图形才是诚实的。

| Property | `Resolution.SPOT` (11 platforms) | `Resolution.SINGLE_CELL` (10 platforms) |
|---|---|---|
| Examples | 10xVisium, 10xVisiumHD, Slide-seq(V2), Stereo-seq, HDST, ST, sci-Space, Pixel-seq, Well-ST-seq, CBSST-seq | MERFISH, seqFISH(+), osmFISH, STARmap, 10xXenium, CosMx, ExSeq, EASI-FISH, EEL-FISH |
| One unit contains | a mixture of cells | one segmented cell |
| Gene coverage | transcriptome-wide | targeted panel, 10²–10³ genes |
| `min_counts` / `min_genes` | 100 / 30 | 20 / 10 |
| Prevalence filter | applied (`min_cells_per_gene=5`) | **disabled** — every probe on a targeted panel was chosen deliberately, so filtering by prevalence discards designed signal |
| Minimum units after QC | 50 | 10 |
| Annotation task | **deconvolution** → per-unit cell-type *proportions* (Tangram) | **classification** → one *label* per cell (SingleR) |
| Figure | scatter-pie, one pie per unit | categorical spatial scatter |

Resolution is obtained from the sample sheet's platform column when annotating, and
otherwise inferred from the directory path. Path inference is deliberately stricter
than a substring scan over the absolute path: each path component is matched
independently, the deepest component wins, whole-component matches beat substrings,
and labels shorter than four characters must match a whole component. The naive
scan this replaces matched the label `ST` inside the mount point `/mnt/cstr/…` and
so claimed nearly every sample; it also iterated a set, making the winner depend on
per-process hash randomisation. Candidates are now sorted longest-first with
alphabetical tie-breaking, so `10xVisiumHD` is tried before `10xVisium` and the
result is identical across processes.

平台解析优先取样本表的平台列，其次从目录路径推断。路径推断按路径分量逐一匹配、
最深分量优先、整分量匹配优于子串匹配、长度小于 4 的标签必须整分量匹配 —— 它替代的
朴素全路径子串扫描会让标签 `ST` 命中挂载点 `/mnt/cstr/…` 从而错认几乎所有样本，
且因遍历 set 而受进程 hash 随机化影响。候选按长度降序、字母序稳定排序，
因此 `10xVisiumHD` 先于 `10xVisium` 被尝试，且跨进程结果一致。

Unrecognised platforms do not abort a run. QC falls back to the permissive
single-cell thresholds (discarding real data is worse than keeping marginal
barcodes, which later stages can still filter) and annotation falls back to
SingleR, with a warning in both cases.

未识别的平台不会中断运行：质控回退到较宽松的单细胞阈值（丢弃真实数据比保留边缘
barcode 更糟，后续阶段仍可过滤），注释回退到 SingleR，两者均记录警告。

Check how a label in your sample sheet will be interpreted before committing to a
run:

```bash
stcompass platforms --check 'Stereo Seq' --check 'Visium HD'
# 'Stereo Seq': Stereo-seq (spot)
# 'Visium HD': 10xVisiumHD (spot)
```

---

## 3. Installation / 安装

Requires Python ≥ 3.10.

```bash
pip install -e .                  # core: scanpy, anndata, scikit-learn, matplotlib
pip install -e '.[tangram]'       # + Tangram deconvolution (needs a CUDA build of torch)
pip install -e '.[singler]'       # + SingleR classification
pip install -e '.[excel]'         # + .xlsx sample sheets
pip install -e '.[watch]'         # + filesystem watching for `plot --watch`
pip install -e '.[dev]'           # + pytest, pytest-cov, ruff
pip install -e '.[all]'           # every optional extra
```

Optional dependencies are resolved at call time through `stcompass._deps.require`,
so `import stcompass` succeeds without torch, Tangram or SingleR installed; a stage
invoked without its dependency raises `MissingDependencyError` naming the exact
extra to install rather than failing at import.

可选依赖通过 `stcompass._deps.require` 在调用时解析，因此未安装 torch / Tangram /
SingleR 也能 `import stcompass`；缺少依赖时调用相应阶段会抛出 `MissingDependencyError`
并指明需安装的 extra，而非在导入阶段失败。

Dependencies carry lower bounds rather than exact pins, because this package is
installed alongside a user-managed scientific stack where hard pins reliably produce
unsolvable environments. Reproducible runs should pin via a lock file — see
[`docs/installation.md`](docs/installation.md), which also covers GPU and R-backend
setup.

依赖采用下界而非精确锁定：本包与用户自管的科学计算栈共存，硬锁定极易导致环境无解。
可复现运行请使用 lock 文件锁定版本。

---

## 4. Quick start / 快速开始

```bash
# 1. Copy a configuration template and edit the paths
cp configs/example.yaml atlas.yaml

# 2. Resolve the config and list what would be processed — no file is read or written
stcompass qc --config atlas.yaml --dry-run

# 3. Run the stages
stcompass qc       --config atlas.yaml
stcompass annotate --config atlas.yaml
stcompass plot     --config atlas.yaml
```

Every setting is also available as a flag, which takes precedence over the file. The
split is deliberate: the YAML file is the reproducible artefact you commit next to
the results, while flags are for a one-off re-run over a subtree.

每个配置项也都有对应命令行参数，且优先级高于配置文件。这一划分是有意的：YAML 是与结果
一并提交的可复现产物，命令行参数用于对子目录的一次性重跑。

```bash
stcompass cluster --qc-dir data/qc --out data/clustered --method leiden -j 8
```

Three configuration templates ship with the package:

| File | Purpose |
|---|---|
| [`configs/minimal.yaml`](configs/minimal.yaml) | the smallest config that runs — paths only |
| [`configs/example.yaml`](configs/example.yaml) | every key with its default and an explanation |
| [`configs/reproduce-atlas.yaml`](configs/reproduce-atlas.yaml) | reproduces the original scripts' behaviour (per-barcode filters **off**), for output comparable to the published atlas rather than to the package defaults |

Full walkthrough: [`docs/usage.md`](docs/usage.md).
Key-by-key reference: [`docs/configuration.md`](docs/configuration.md).

---

## 5. Methods / 方法

**Notation.** $X \in \mathbb{R}_{\ge 0}^{n \times p}$ is the expression matrix of one
sample: $n$ measurement units (spots, bins or segmented cells) by $p$ genes. All
defaults below are the shipped values in `stcompass/config.py`.

**记号.** $X \in \mathbb{R}_{\ge 0}^{n \times p}$ 为单样本表达矩阵，$n$ 个测量单元
× $p$ 个基因。下文默认值均为 `stcompass/config.py` 中的实际取值。

### 5.1 `qc` — filtering, normalisation, embedding

1. **Deduplicate.** Duplicate `var_names` break every `adata[:, gene]` lookup
   downstream and duplicate `obs_names` make the object unwritable. The first
   occurrence of each gene is kept and the original column order restored, so the
   matrix stays positionally comparable to the source file.
2. **Preserve counts.** $X$ is copied to `layers["counts"]` *before* any transform,
   which is what allows §5.4 to factorise counts after `X` has been log-normalised.
3. **QC metrics.** `n_counts`, `n_genes`, and `percent_mt` when mitochondrial genes
   are detectable (the `MT-` prefix is matched case-insensitively, so human `MT-` and
   mouse `Mt-` both resolve).
4. **Filter** on the resolution-dependent thresholds of §2. A sample left with fewer
   than `min_units_after_qc` units is *rejected with a reason*, not raised as an error.
5. **Normalise conditionally.** Public atlases mix raw and pre-normalised matrices
   and carry no reliable flag, so the state is inferred rather than assumed:

$$X_{ij} \leftarrow \log\left(1 + \frac{X_{ij}}{\sum_k X_{ik}} \cdot s\right), \quad s = 10^4$$

   applied **only if** `looks_like_counts(X)` holds. Applying `log1p` twice compresses
   real biological variation; skipping it leaves counts on a scale that breaks PCA.
   The test inspects `sample_size=2000` stored values — sampled *with* replacement,
   because sampling without replacement forces a permutation of every stored value,
   which for a matrix with $10^9$ nonzeros costs more memory than the pipeline stage
   that follows — and returns true when they are integral to `atol=1e-6`. The RNG is
   seeded, so the decision is identical on a re-run.
6. **Embed.** `highly_variable_genes` (`flavor="seurat"`, `n_top_genes=3000`) →
   PCA ($n_{\text{comps}} = \min(50,\, n-1,\, p-1)$) → kNN graph
   ($k = \max(2, \min(15, n-1))$) → Leiden at `resolution=1.0` → UMAP.
   Leiden is pinned to `flavor="igraph"`, `n_iterations=2`, `directed=False` so
   cluster labels are reproducible across scanpy versions rather than shifting with
   the library's changing defaults.

### 5.2 `cluster` — size-aware re-clustering

A 200-spot section clustered at `resolution=0.2` collapses into a single domain,
while a 500k-bin Visium HD section clustered at `1.2` shatters into hundreds of
unusable fragments. Resolution is therefore a step function of sample size:

$$
r(n) = \begin{cases}
1.2 & n < 100 \\
0.7 & 100 \le n < 500 \\
0.5 & 500 \le n < 5000 \\
0.3 & 5000 \le n < 20000 \\
0.2 & n \ge 20000
\end{cases}
$$

The schedule lives in `ClusteringConfig.resolution_schedule` as `(max_cells,
resolution)` pairs — YAML-editable, sorted on load, and overridable by a fixed
`clustering.resolution`. The kNN graph computed during QC is reused unless
`recompute_neighbors` is set, since rebuilding it is the expensive part and re-running
at a different resolution is the common case.

分辨率是样本规模的阶梯函数，schedule 以 `(max_cells, resolution)` 对存于配置中，
加载时排序，可被固定的 `clustering.resolution` 覆盖。除显式要求外复用质控阶段的 kNN 图。

### 5.3 `annotate` — cell-type annotation

Driven by a sample sheet (`.xlsx`/`.csv`) supplying species, tissue and platform per
sample. Work is grouped by `(species, tissue)` because the expensive preparation —
loading the reference and ranking its marker genes — is shared by every sample in a
group. References are laid out as `reference/<species>/<tissue>/*.h5ad`.

1. **Index once.** The sample tree is walked a single time into a `stem → paths`
   mapping. The scripts this replaces called `glob(root/**/<sample>.h5ad)` once per
   sheet row, re-walking the whole atlas per sample — $O(\text{rows} \times \text{files})$
   and the dominant cost for a sheet with thousands of rows. One walk is $O(\text{files})$.
2. **Prepare the reference.** Labels with fewer than `min_cells_per_label=2` cells are
   dropped: a singleton label has no within-group variance, which makes
   `rank_genes_groups` fail outright, and one cell cannot support a proportion estimate.
3. **Rank markers** — only when the group contains at least one spot-based sample,
   since ranking is expensive and SingleR does not need it. The gene set is the union
   of the top `n_marker_genes=100` markers per label. Restricting Tangram to markers
   both accelerates the mapping and prevents housekeeping genes from dominating the
   objective.
4. **Deconvolve or classify**, by resolution:
   - **Spot-based → Tangram.** `map_cells_to_space` (`mode="cells"`, 300 epochs,
     `density_prior="rna_count_based"`), then `project_cell_annotations` writes
     spot × cell-type proportions to `obsm["tangram_ct_pred"]`. An empty shared gene
     set raises immediately with a nomenclature hint (gene symbols vs Ensembl IDs),
     which is the usual cause.
   - **Single-cell → SingleR.** `annotate_single` over all shared genes; each output
     field becomes an `obs["singler_*"]` column, with the assigned label in
     `obs["singler_best"]`.
5. **Place on a device.** With several GPUs, samples are assigned by
   $\text{md5}(\text{path}) \bmod n_{\text{gpus}}$ rather than round-robin, so a
   resumed run sends each sample to the same device as before and memory contention
   is reproducible. CUDA-unavailable falls back to CPU with a warning, so a config
   written for a GPU box still runs on a laptop.

References are cached per worker process, keyed by `(path, label_key,
min_cells_per_label)` and bounded to one entry. This matters because joblib re-imports
the module in each worker: pickling the AnnData into every task instead — what the
original scripts did — serialises the entire reference once per sample.

参考数据按 `(path, label_key, min_cells_per_label)` 在每个工作进程内缓存，且仅保留
一条以限制内存。joblib 会在 worker 中重新导入模块，因此这种做法避免了逐样本序列化
整个参考数据集。

### 5.4 `programs` — NMF gene programs

Factorises $X \approx WH$ with $W \in \mathbb{R}_{\ge 0}^{n \times K}$,
$H \in \mathbb{R}_{\ge 0}^{K \times p}$, minimising

$$\tfrac{1}{2}\lVert X - WH \rVert_F^2 \quad \text{s.t.}\ W, H \ge 0$$

Unlike PCA the factors are non-negative and therefore additive, which is what makes
them readable as *programs* rather than as axes of variation.

The implementation is streaming by necessity. A Visium HD section has millions of 2 µm
bins; densifying costs $4np$ bytes — hundreds of gigabytes — so the matrix is read from
a **backed** `.h5ad` in row blocks (`row_chunk=5000`), kept sparse, and fitted with
`MiniBatchNMF.partial_fit` over `epochs=2` passes at `batch_size=4096`. Two passes
matter: a single pass leaves the factors biased towards the rows seen last, because each
mini-batch update only partially corrects the previous one. Initialisation is `nndsvda`,
and `batch_size` is raised to at least $K$ so that initialisation is well-posed;
$K$ is itself clipped to $\min(K, n, p)$ with a warning.

Genes that are all-zero are always dropped — they contribute nothing to the loss and
would leave a zero row in $H$. With `max_hvg` set, the remainder is ranked by variance
and the top genes kept; indices are sorted so the selection is deterministic and $H$
columns map back by position. $H$ is then re-expanded to the full gene axis with zeros
for excluded genes, keeping `var` columns comparable across samples with different gene
sets.

**Program specificity.** A gene with a high weight in *every* program is a housekeeping
gene and says nothing about program identity, so a raw $H$ weight is a poor marker
score. With $\tilde{H} = H^\top$ (genes × programs):

$$
\text{score}_{gk} = \tilde{H}_{gk} \cdot \log\left(1 + \frac{\tilde{H}_{gk}}{\rho_{gk} + \varepsilon}\right),
\qquad
\rho_{gk} = \begin{cases}
\text{2nd-largest}_k \tilde{H}_{gk} & k = \arg\max_k \tilde{H}_{gk} \\
\max_k \tilde{H}_{gk} & \text{otherwise}
\end{cases}
$$

with $\varepsilon = 10^{-10}$ guarding genes that are zero everywhere. A gene loading on
one program only gets a large ratio and scores highly; a flat gene gets a ratio near one
and scores near zero. Ranks are 1-based, rank 1 being the most specific gene for that
program. With $K = 1$ there is nothing to be specific against, so the weight itself is
the ranking.

### 5.5 `plot` — figures

The figure follows from what the annotation stage produced. `obs["singler_best"]` is
checked first; failing that, `obsm["tangram_ct_pred"]`; failing both, the sample is
skipped with a message naming the command to run.

- **Proportions → scatter-pie.** Only the `top_n_types=3` largest proportions per unit
  are kept, then renormalised so every pie remains a full circle and wedge angles stay
  comparable between units. Selection uses `np.partition`, i.e. $O(K)$ per row rather
  than $O(K \log K)$, which matters at $10^5$ spots. Above `max_pies=20000` units the
  renderer falls back to a plain scatter: hundreds of thousands of wedges exhaust memory
  and produce an unreadable figure.
- **Radius.** Spot pitch spans three orders of magnitude (55 µm Visium spots, 2 µm HD
  bins, arbitrary units after registration), so a fixed radius is wrong everywhere.
  The default is $0.45 \times \operatorname{median}(d_{\text{1-NN}})$ via a k-d tree —
  half the spacing would leave adjacent pies exactly touching, so 0.45 leaves a thin gap
  and boundaries stay visible.
- **Labels → categorical scatter**, with a 60-colour palette chained from `tab20`,
  `tab20b` and `tab20c`, cycling beyond that. A continuous colormap would be wrong:
  adjacent cell types are not adjacent in any meaningful ordering.
- **Watch mode** (`plot.watch`) keeps polling and renders figures as annotation output
  appears, which is useful while a multi-day annotation run is still going. Polling is
  used rather than filesystem events because the atlas typically lives on NFS, where
  inotify does not fire for writes made on another host. A new file is read only after
  `stability_checks=3` consecutive equal-size observations `stability_interval=2.0` s
  apart, bounded by a timeout of `checks × interval × 8` so a stalled copy cannot block
  the queue.

`Agg` is selected before pyplot is first imported: these pipelines run over SSH and under
batch schedulers where no display exists, and an interactive backend fails at import time
there.

Stage-by-stage detail: [`docs/pipelines.md`](docs/pipelines.md).

---

## 6. Output schema / 输出字段规范

| Stage | Key | Type | Contents |
|---|---|---|---|
| `qc` | `layers["counts"]` | $n \times p$ sparse | raw counts, preserved before normalisation |
| `qc` | `obs["n_counts"]`, `obs["n_genes"]`, `obs["percent_mt"]` | $n$ | per-unit QC metrics (`percent_mt` only when MT genes exist) |
| `qc` | `obs["leiden"]` \| `obs["louvain"]` | $n$ categorical | clusters at `qc.resolution` |
| `qc` | `obsm["X_pca"]`, `obsm["X_umap"]` | $n \times 50$, $n \times 2$ | embeddings |
| `cluster` | `obs["louvain"]` (or `clustering.key_added`) | $n$ categorical | re-clustered labels at $r(n)$ |
| `annotate` | `obsm["tangram_ct_pred"]` | $n \times C$ | per-spot cell-type proportions (spot-based) |
| `annotate` | `obs["singler_*"]`, incl. `obs["singler_best"]` | $n$ | per-cell labels and scores (imaging-based) |
| `programs` | `obsm["X_nmf"]` | $n \times K$ float32 | program loadings $W$ |
| `programs` | `var["NNMF_component_{k}"]` | $p$ | gene weights $H$, zero-filled for excluded genes |
| `programs` | `var["NNMF_component_{k}_importance_score"]` | $p$ | specificity score, §5.4 |
| `programs` | `var["NNMF_component_{k}_rank"]` | $p$ int32 | 1-based rank by specificity |
| `plot` | — | file | one figure per sample, mirroring the input tree |

`.h5ad` writes are **atomic**: the object is written to a dot-prefixed temporary file in
the destination directory and moved into place on success. This prevents watch-mode from
reading a partially written file, and prevents a crash from leaving a truncated output
that a resumed run would mistake for finished work. Two malformed-file patterns are also
repaired centrally on read/write — an undecodable `uns/log1p/base` scalar, and a reserved
`_index` column in `obs`/`var` (renamed rather than dropped, since it holds the original
barcode or gene identifiers).

`.h5ad` 写入是原子的：先写入目标目录下的临时文件，成功后移动到位。这可避免 watch 模式
读到半写文件，也避免崩溃留下会被续跑误判为已完成的截断输出。

---

## 7. Architecture / 软件架构

```
STCompassDB/
├── src/stcompass/            # 4,141 LOC across 15 modules
│   ├── cli.py            (641)  five subcommands, flag→config overrides, exit codes, --dry-run
│   ├── config.py         (621)  typed YAML config; unknown keys rejected
│   ├── io.py             (286)  tree mirroring, atomic + hardened .h5ad round-trip
│   ├── platforms.py      (233)  platform registry → spot vs single-cell resolution
│   ├── matrix.py         (230)  block-wise sparse helpers (no scanpy/anndata import)
│   ├── logging_utils.py   (79)  console + optional DEBUG file logging
│   ├── _deps.py           (53)  lazy optional imports with install hints
│   └── pipelines/
│       ├── annotation.py (549)  Tangram / SingleR, grouped by (species, tissue)
│       ├── visualize.py  (456)  scatter-pie, categorical scatter, watch loop
│       ├── programs.py   (307)  streaming MiniBatchNMF + specificity scoring
│       ├── qc.py         (299)  filtering, conditional normalisation, embedding
│       ├── _batch.py     (174)  per-sample loop, error capture, exit summary
│       └── clustering.py (116)  size-aware re-clustering
├── tests/                    # 271 tests; scanpy-dependent ones self-skip
├── configs/                  # minimal.yaml, example.yaml, reproduce-atlas.yaml
├── docs/                     # installation, usage, configuration, pipelines
├── pyproject.toml            # packaging, ruff, pytest, coverage
└── CHANGELOG.md  CONTRIBUTING.md  LICENSE
```

Four invariants hold across the layers:

**Configuration is typed and closed.** Each stage has a dataclass; `Config` bundles them.
Unknown keys are *rejected*, not ignored, because a silently dropped `min_count` typo
produces plausible-looking but wrong results. Validation happens in `__post_init__`, and
errors name the section and key. YAML 1.1 — which PyYAML implements — only reads an
exponent as a float when the exponent carries an explicit sign, so `target_sum: 1e4`
arrives as the *string* `"1e4"`; scalars are coerced with an error message that says to
write `1.0e+4` or `10000` rather than surfacing as a `TypeError` deep in a comparison.

**A per-sample problem is a return value, not an exception.** Workers return `None` on
success or a short reason string to record a skip; only genuine failures raise. `run_batch`
converts exceptions into `("failed", message)` tuples rather than propagating them, which is
what keeps the parallel branch simple — joblib would otherwise cancel remaining tasks on the
first exception. The resulting `BatchResult` is what the CLI turns into an exit code, so a
batch where every sample failed cannot report success.

**Optional dependencies are gated at call time.** No module-level import of scanpy, torch,
Tangram or SingleR anywhere in the package. `matrix.py` deliberately imports neither scanpy
nor anndata, which keeps the numerically interesting parts testable in a minimal environment.

**Pipelines do not import each other.** Each depends only on `config`, `io`, `matrix`,
`platforms` and `_batch`, so a stage can be run, tested and reasoned about in isolation.

配置层类型化且封闭（未知键直接拒绝）；单样本问题以返回值而非异常表达；可选依赖在调用时
才解析；各 pipeline 之间互不导入。

Multiprocessing uses joblib with `prefer="processes"`. The loky backend re-imports the
module in a fresh interpreter, which is precisely what isolates a segfaulting native library
to the one sample that triggered it. `n_jobs=1` runs in-process, keeping tracebacks readable
— the right choice for GPU stages that manage their own device placement.

---

## 8. Scalability / 可扩展性

| Concern | Approach |
|---|---|
| Matrices exceeding RAM | backed `.h5ad` reads in row blocks; the sparse representation is never densified. Densifying costs $4np$ bytes — hundreds of GB for a Visium HD section |
| NMF at $10^6$ units | `MiniBatchNMF.partial_fit` over blocks; peak memory is one block, not one matrix. Gene selection streams over the backed file, and only the selected columns are materialised |
| Counts-vs-normalised check | fixed-size sample of stored values, $O(1)$ in matrix size |
| Per-spot top-$K$ selection | `np.partition`, $O(K)$ per row |
| Sample lookup | one tree walk into a dict, replacing a per-row glob over the whole atlas |
| Reference loading | per-process cache bounded to one entry, avoiding per-sample deserialisation |
| Resumability | output existence is the checkpoint; no external state to corrupt or reconcile |
| Interrupted writes | atomic move, so a partial file is never mistaken for a finished one |
| Failure isolation | separate interpreter per worker; one native crash costs one sample |

---

## 9. Reproducibility / 可复现性

Reproducibility is a design constraint, not a documentation exercise. Concretely:

- **The config is the artefact.** Every parameter the original scripts hard-coded at module
  level — input and output roots, the sample sheet path, count thresholds, $K$, GPU counts —
  is a YAML key. An edited-in-place script cannot be committed next to its results; a config
  file can.
- **Deterministic platform resolution.** Candidate labels are sorted longest-first with
  alphabetical tie-breaking. The set-union iteration this replaces was subject to per-process
  hash randomisation, so the same file could be classified differently on a re-run.
- **Seeded decisions.** The counts heuristic draws from a seeded generator
  (`seed=0`), not the global RNG; NMF uses `random_state=0`; Tangram device assignment is a
  hash of the sample path, so a resumed run reproduces the original placement.
- **Pinned algorithm variants.** Leiden is called with `flavor="igraph"`,
  `n_iterations=2`, `directed=False` so labels do not shift when scanpy changes its defaults.
- **Sorted traversal.** `iter_samples` yields pairs in sorted order, so two runs process the
  same samples in the same order.
- **Numerical behaviour is tested.** The specificity score, the top-$N$ renormalisation, the
  resolution schedule and the counts heuristic each have tests pinning their exact output,
  because a silent change in any of them would alter published results without raising an error.
- **Legacy parity.** [`configs/reproduce-atlas.yaml`](configs/reproduce-atlas.yaml) documents
  and restores the original behaviour, including the fact that the per-barcode and per-gene
  filters were commented out when the published atlas was built. Enabling them — the package
  default — produces smaller, cleaner objects that will **not** match the published files.

可复现性是设计约束：配置文件即产物；平台解析确定性排序；随机决策均设种子；算法变体显式
固定以免随上游默认值漂移；遍历顺序排序；关键数值行为均有测试锁定；并提供与原始脚本行为
对齐的配置以复现已发表结果。

---

## 10. Testing / 测试

```bash
pip install -e '.[dev]'
ruff check src tests
ruff format --check src tests
pytest
```

271 tests across 8 modules cover the platform registry and its alias resolution, config
loading and rejection, the block-wise matrix helpers, the batch driver's outcome accounting,
NMF specificity scoring, the plotting mathematics (top-$N$ renormalisation, radius estimation,
file-stability detection) and the CLI's flag-to-config plumbing and exit codes.

`[dev]` deliberately does **not** install scanpy, torch, tangram-sc or singler. Tests needing
scanpy are marked `requires_scanpy` and skip themselves, so the suite is green on a laptop and
thorough on a full stack, and a contributor can work on the config layer, CLI or matrix helpers
without a multi-gigabyte install. Preference is for testing the single-sample functions
(`qc_sample`, `cluster_sample`, `programs_sample`) over the batch runner, since they accept an
in-memory `AnnData` and touch no filesystem.

271 个测试覆盖平台注册表与别名解析、配置加载与拒绝、分块矩阵工具、批处理结果统计、
NMF 特异性打分、绘图数学与 CLI。`[dev]` 有意不安装 scanpy/torch/tangram/singler；
依赖 scanpy 的测试标记为 `requires_scanpy` 并自动跳过。

Conventions, including the comment policy and what a numerical change must ship with:
[`CONTRIBUTING.md`](CONTRIBUTING.md).

---

## 11. Limitations / 已知局限

Stated explicitly, since each affects how results should be interpreted:

- **No CI workflow is committed.** `CONTRIBUTING.md` refers to "what CI runs", but there is no
  `.github/workflows/` in the repository; the three commands in §10 must be run locally.
- **Project URLs in `pyproject.toml` are placeholders** (`github.com/OWNER/STCompass`), as is
  the clone URL in `CONTRIBUTING.md`.
- **No archived release or DOI.** There is nothing citable beyond the repository itself.
- **Clustering is non-spatial.** Both `qc` and `cluster` operate on a kNN graph in PCA space;
  spatial coordinates are used for rendering only. Spatial-domain methods (SpaGCN, BayesSpace,
  hidden-Markov random fields) are not implemented.
- **$K$ is fixed, not selected.** `programs.n_components` defaults to 10 for every sample; no
  stability or reconstruction-error criterion chooses it, and programs are not matched across
  samples, so `NNMF_component_3` in two samples is not the same program.
- **Every sample is processed independently.** There is no cross-sample integration or batch
  correction; clusters and programs are not comparable across samples without a further step.
- **Tangram needs a CUDA build of torch** for realistic runtimes. CPU fallback works but is
  impractical at atlas scale.
- **Annotation quality is bounded by the reference.** A `(species, tissue)` pair with no
  reference is skipped, and a reference whose cell types do not cover the tissue yields
  confident-looking but wrong labels — this is a property of Tangram and SingleR, not of the
  wrapper.
- **`percent_mt` is silently absent** when no `MT-`-prefixed genes exist, which is the normal
  case for targeted panels.

已知局限均直接影响结果解读方式，包括：仓库未提交 CI 配置；`pyproject.toml` 中项目 URL 为占位；
无归档发布与 DOI；聚类不使用空间坐标；$K$ 为固定值且程序不跨样本对齐；样本间无整合与批次校正；
Tangram 实际需要 CUDA；注释质量受参考数据上限约束。

---

## 12. Citation / 引用

There is no archived release yet, so cite the repository and the commit you used:

```bibtex
@software{stcompassdb,
  title    = {{STCompassDB}: Reproducible batch pipelines for large-scale
              spatial transcriptomics atlases},
  author   = {{STCompassDB contributors}},
  version  = {0.1.0},
  year     = {2026},
  url      = {https://github.com/MIKUAFANS/STCompass},
  note     = {Software. Please also cite the underlying methods listed below.}
}
```

STCompassDB orchestrates third-party methods; it does not replace them. If you publish
results produced with this package, cite the methods your run actually invoked — at minimum
scanpy, plus Tangram for spot-based annotation or SingleR for imaging-based annotation.
Each retains its own licence and citation requirements.

本包编排第三方方法而非取代之。若使用本包产出的结果发表论文，请引用运行中实际调用的方法
（至少包括 scanpy，以及点阵型注释所用的 Tangram 或成像型注释所用的 SingleR）。

**Method references / 方法文献**

1. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data
   analysis. *Genome Biology* 19:15 (2018).
2. Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Annotated data.
   *bioRxiv* 2021.12.16.473007 (2021).
3. Biancalani T, Scalia G, Buffoni L, *et al.* Deep learning and alignment of spatially
   resolved single-cell transcriptomes with Tangram. *Nature Methods* 18:1352–1362 (2021).
4. Aran D, Looney AP, Liu L, *et al.* Reference-based analysis of lung single-cell
   sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*
   20:163–172 (2019). [SingleR]
5. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected
   communities. *Scientific Reports* 9:5233 (2019).
6. Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E. Fast unfolding of communities in
   large networks. *Journal of Statistical Mechanics* P10008 (2008). [Louvain]
7. Lee DD, Seung HS. Learning the parts of objects by non-negative matrix factorization.
   *Nature* 401:788–791 (1999).
8. Févotte C, Idier J. Algorithms for nonnegative matrix factorization with the
   β-divergence. *Neural Computation* 23:2421–2456 (2011). [MiniBatchNMF]
9. Boutsidis C, Gallopoulos E. SVD based initialization: a head start for nonnegative
   matrix factorization. *Pattern Recognition* 41:1350–1362 (2008). [`nndsvda`]
10. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection
    for dimension reduction. *arXiv*:1802.03426 (2018).
11. Satija R, Farrell JA, Gennert D, Schier AF, Regev A. Spatial reconstruction of
    single-cell gene expression data. *Nature Biotechnology* 33:495–502 (2015).
    [`hvg_flavor="seurat"`]
12. Pedregosa F, Varoquaux G, Gramfort A, *et al.* Scikit-learn: machine learning in
    Python. *JMLR* 12:2825–2830 (2011).

---

## License / 许可

[MIT](LICENSE).

Third-party methods retain their own licences and citation requirements — notably
[Tangram](https://github.com/broadinstitute/Tangram),
[SingleR](https://github.com/SingleR-inc/singler-py) and
[scanpy](https://github.com/scverse/scanpy).

第三方方法保留其各自的许可与引用要求，特别是 Tangram、SingleR 与 scanpy。

## Documentation / 文档

| Document | Contents |
|---|---|
| [`docs/installation.md`](docs/installation.md) | environment requirements, extras, GPU and R-backend setup, lock files, troubleshooting |
| [`docs/usage.md`](docs/usage.md) | end-to-end walkthrough, directory layout, per-stage invocation, result interpretation, Python API |
| [`docs/configuration.md`](docs/configuration.md) | every configuration key, its default and its rationale; precedence rules; the YAML exponent gotcha |
| [`docs/pipelines.md`](docs/pipelines.md) | algorithms and parameters, stage by stage |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | development setup, conventions, testing policy |
| [`CHANGELOG.md`](CHANGELOG.md) | version history |









