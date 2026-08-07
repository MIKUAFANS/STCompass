# Usage / 使用指南

This page walks through a full atlas build, stage by stage. Every command takes
`--config FILE`; every setting in that file can also be given as a flag, and the
flag wins.

本页按阶段介绍完整的图谱构建流程。所有命令都接受 `--config FILE`；配置文件中的
每一项也都可以用命令行参数覆盖，命令行优先。

---

## 0. Before you start / 开始之前

```bash
stcompass --help              # list the stages
stcompass qc --help           # options for one stage
stcompass platforms           # how platform labels are interpreted
```

Copy a config and edit the paths:

```bash
cp configs/example.yaml atlas.yaml
$EDITOR atlas.yaml
```

Check what a stage *would* do before it does it — this reads no matrices and
needs no scanpy:

```bash
stcompass qc --config atlas.yaml --dry-run
```

`--dry-run` prints the input root, the output root, and the source → destination
pair for each sample that is not already done. `--exclude` is honoured, so the
preview matches the real run.

`--dry-run` 会打印输入/输出根目录，以及每个待处理样本的源路径与目标路径。它同样
遵守 `--exclude`，因此预览结果与实际运行一致。

---

## 1. Directory layout / 目录结构

Every stage mirrors its input tree into its output root. Nothing is flattened,
and nothing is renamed:

```
raw/10xVisium/Homo sapiens/S1.h5ad
  → qc/10xVisium/Homo sapiens/S1.h5ad
  → clustered/10xVisium/Homo sapiens/S1.h5ad
  → programs/10xVisium/Homo sapiens/S1.h5ad
```

Two consequences worth knowing:

- **Platform inference.** When no sample sheet is given, the platform is read
  from the *directory names*. A layout of `<root>/<platform>/<species>/<sample>.h5ad`
  is what the QC stage expects. If a platform cannot be inferred, QC logs a
  warning and falls back to permissive thresholds rather than guessing.
- **Resumability.** A sample whose output already exists is skipped. Interrupt a
  run with Ctrl-C and start it again; it picks up where it stopped. Pass
  `--overwrite` to force reprocessing.

当未提供样本表时，平台类型从**目录名**推断，期望的布局是
`<root>/<platform>/<species>/<sample>.h5ad`。若无法推断，QC 会记录警告并退回到
宽松阈值，而不是猜测。已存在输出的样本会被跳过，因此运行可以中断后续跑；用
`--overwrite` 可强制重算。

---

## 2. Quality control / 质量控制

```bash
stcompass qc --config atlas.yaml
```

What it does, in order:

1. Drop duplicate gene names (keeping the first) and make barcode names unique.
2. Copy the raw matrix to `layers["counts"]` — the programs stage needs counts
   after `X` has been log-transformed.
3. Compute `obs["n_counts"]`, `obs["n_genes"]`, and `obs["percent_mt"]` when
   mitochondrial genes are detectable (`MT-` prefix, case-insensitive).
4. Filter barcodes and genes using platform-dependent thresholds.
5. Normalise to `target_sum` and `log1p` — **only if** the matrix looks like raw
   counts. See [Counts detection](#counts-detection) below.
6. Select highly variable genes, then PCA → neighbours → clustering → UMAP.

Thresholds differ by platform resolution because the same cutoff cannot serve
both assay families:

| | spot-based | single-cell |
|---|---|---|
| min counts | 100 | 20 |
| min genes | 30 | 10 |
| gene prevalence filter | ≥ 5 barcodes | none |
| min units to keep sample | 50 | 10 |

A Visium barcode pools several cells and yields thousands of counts; a MERFISH
cell is measured on a few-hundred-gene panel and may legitimately carry twenty.
The prevalence filter is skipped entirely on targeted panels — every probe there
was chosen deliberately, so filtering by prevalence discards real signal.

阈值按平台分辨率区分：Visium 的一个 barcode 汇集多个细胞、计数上千；而 MERFISH
的单个细胞只测几百个基因，二十个计数可能是真实的。靶向 panel 上完全跳过基因流行度
过滤——那里每个探针都是特意选择的。

### Counts detection

Public atlases mix raw and already-normalised matrices with no reliable metadata
flag. Applying `log1p` twice compresses real biological variation; skipping it
leaves counts on a scale that breaks PCA. So the stage *infers* it: it samples
2000 stored values and checks whether they are integral. Normalised or
log-transformed data is essentially never integral.

The sample is drawn with a fixed seed, so the decision is reproducible — the
original script drew from the global RNG and could classify the same file
differently on a re-run.

公共图谱混杂了原始计数与已归一化的矩阵，且没有可靠的元数据标记。因此该阶段通过
抽样 2000 个存储值、检查是否为整数来*推断*。随机种子固定，故判定可复现——原脚本
使用全局随机数，同一文件重跑可能得到不同结论。

### Reproducing the original atlas

The published atlas was built with the per-barcode and per-gene filters
**disabled** (they were commented out in `QC_ALL0917.py`), keeping every barcode
for downstream deconvolution. To reproduce that:

```bash
stcompass qc --config configs/reproduce-atlas.yaml
# or:
stcompass qc --config atlas.yaml --no-filter-cells --no-filter-genes
```

---

## 3. Clustering / 聚类

```bash
stcompass cluster --config atlas.yaml --method louvain
```

Reads `paths.qc`, writes `paths.clustered`. This stage exists separately from QC
because re-clustering is the most common thing to redo, and it should not require
redoing normalisation and PCA.

Without `--resolution`, the resolution is chosen from the number of barcodes:

| barcodes | resolution |
|---|---|
| < 100 | 1.2 |
| < 500 | 0.7 |
| < 5 000 | 0.5 |
| < 20 000 | 0.3 |
| ≥ 20 000 | 0.2 |

Small sections need a coarse graph to yield interpretable domains; large ones
fragment into hundreds of clusters unless the resolution drops. These are the
thresholds from the original `run_louvain.py`, now in
`clustering.resolution_schedule` where they can be tuned without editing code.

若不指定 `--resolution`，分辨率将根据 barcode 数量自动选择：小切片需要较粗的图才能
得到可解释的区域，大切片若不降低分辨率会碎裂成上百个簇。这些阈值来自原
`run_louvain.py`，现在位于 `clustering.resolution_schedule`，无需改代码即可调整。

By default the neighbour graph computed during QC is reused. Pass
`--recompute-neighbors` to rebuild it.

---

## 4. Cell-type annotation / 细胞类型注释

```bash
stcompass annotate --config atlas.yaml
```

This is the only stage driven by a **sample sheet** rather than a tree walk,
because it needs each sample's species and tissue to pick a reference.

### Requirements

- `paths.metadata` — an `.xlsx` or `.csv` sheet with columns for species, tissue,
  sample name and platform (names configurable; defaults match the original atlas
  sheet).
- `paths.reference` — scRNA-seq references laid out as
  `<reference>/<species>/<tissue>/*.h5ad`, each with a cell-type column
  (`annotation.label_key`, default `cell_ontology_class`).

### Method selection

The method follows the platform's resolution class — it is not a user choice,
because the two cases are different problems:

| platform | problem | method | output |
|---|---|---|---|
| spot-based | one barcode = several cells → **deconvolution** | Tangram | `obsm["tangram_ct_pred"]`, proportions per spot |
| single-cell | one unit = one cell → **classification** | SingleR | `obs["singler_best"]`, one label per cell |

方法由平台分辨率决定，而非用户选择，因为这是两个不同的问题：spot 平台的一个
barcode 包含多个细胞，是**去卷积**问题；单细胞平台一个单位就是一个细胞，是**分类**
问题。

Unknown platform labels fall back to SingleR with a warning. Check how your
labels resolve before a long run:

```bash
stcompass platforms --check 'Stereo Seq' --check 'Visium HD'
```

### Grouping and cost

Work is grouped by `(species, tissue)` so the expensive preparation — loading the
reference and ranking its marker genes — is shared by every sample in the group.
Marker ranking happens only if some sample in the group is spot-based.

A group whose reference is missing is recorded as a *skip* for every sample in
it, not a failure, and the run continues. An atlas legitimately contains tissues
with no matching reference.

按 `(species, tissue)` 分组，使加载参考数据和排序 marker 基因这类昂贵操作在组内共享。
缺少参考数据的分组会把组内所有样本记为**跳过**而非失败，运行继续——图谱中确实存在
没有对应参考的组织。

### GPU placement

```bash
stcompass annotate --config atlas.yaml --device cuda --n-gpus 4 -j 4
```

Samples are assigned to GPUs by a hash of their path, not round-robin. A resumed
run therefore sends each sample to the same device as before, so memory pressure
is reproducible rather than depending on which samples happened to finish first.
Without CUDA, everything falls back to CPU with a warning.

样本通过路径哈希分配到 GPU，而非轮询。因此续跑时每个样本仍会落到同一设备，显存压力
可复现，不依赖于哪些样本先完成。没有 CUDA 时会带警告退回 CPU。

---

## 5. Gene programs / 基因程序

```bash
stcompass programs --config atlas.yaml --n-components 15
```

Factorises expression into non-negative gene programs with `MiniBatchNMF`.

Reads `paths.qc`, writes `paths.programs`. Outputs:

- `obsm["X_nmf"]` — program loadings per barcode, shape `(n_obs, K)`
- `var["NNMF_component_{k}"]` — gene weight in program *k*
- `var["NNMF_component_{k}_importance_score"] `— specificity score (below)
- `var["NNMF_component_{k}_rank"]` — rank 1 = most specific gene for program *k*

The sample is read twice by design. The first read is **backed**: the matrix
stays on disk while gene selection streams over it in row blocks, so a
million-barcode sample costs no more memory than one block. Only the selected
columns are then materialised, and the file is re-read unbacked to attach the
factors.

该阶段有意读取文件两次。第一次是 **backed** 模式：矩阵留在磁盘上，基因选择按行块
流式扫描，因此百万 barcode 的样本占用的内存不超过一个行块。随后只将选中的列载入
内存，再以非 backed 方式重读文件以写回因子。

`epochs` defaults to 2 rather than 1: `partial_fit` biases the factors towards
the rows seen last, because each mini-batch update only partially corrects the
previous one. A single pass leaves them under-fitted.

`epochs` 默认为 2 而非 1：`partial_fit` 会使因子偏向最后见到的行，单次遍历会导致
欠拟合。

### Importance score

The score for gene *i* in program *k* is

```
s[i,k] = w[i,k] * log(1 + w[i,k] / m[i,k])
```

where `m[i,k]` is the largest weight of gene *i* in any *other* program (for the
gene's own top program, it is the runner-up). A gene weighted highly in one
program and near zero elsewhere scores high; a gene weighted equally everywhere
scores low regardless of magnitude. This ranks by *specificity*, not by weight —
which is what makes the top-ranked genes interpretable as program markers.

该分数按**特异性**而非权重排序：在单个程序中权重高、在其他程序中接近零的基因得分高；
在所有程序中权重相同的基因无论量级如何得分都低。这使排名靠前的基因可解释为程序标记
基因。

---

## 6. Figures / 绘图

```bash
stcompass plot --config atlas.yaml
```

Reads `paths.annotated`, writes `paths.figures`, replacing `.h5ad` with
`.{plot.format}`.

The rendering follows what the annotation stage produced:

- **Tangram proportions** → per-spot pie charts. Only the top
  `plot.top_n_types` types are drawn per spot (default 3); drawing thirty types
  turns each pie into a ring of slivers. Wedges are renormalised so each pie is a
  full circle and angles stay comparable between spots.
- **SingleR labels** → a categorical embedding, preferring spatial coordinates
  and falling back to UMAP/PCA.

Pie radius is derived from nearest-neighbour spacing (0.45 × median distance),
because spot pitch differs by three orders of magnitude across platforms — 55 µm
Visium spots versus 2 µm HD bins — so a fixed radius is wrong everywhere.

饼图半径由最近邻间距推导（0.45 × 中位距离），因为不同平台的点间距相差三个数量级
（55 µm 的 Visium 点与 2 µm 的 HD bin），固定半径在任何情况下都不合适。

Above `plot.max_pies` spots (default 20 000) the stage falls back to a scatter
coloured by dominant type: hundreds of thousands of wedges exhaust memory in the
PDF backend and would be individually invisible anyway.

### Watch mode

```bash
stcompass plot --config atlas.yaml --watch
```

Renders files as an annotation job produces them, so figures appear during a long
run instead of after it. Requires `pip install 'stcompass[watch]'` for inotify;
without it, polling alone still works.

Polling runs **even when** watchdog is available, because inotify does not fire
for writes on NFS mounts — a common setup for atlas-scale storage. New files are
size-checked until stable before being read, so a half-written `.h5ad` is never
parsed.

即便 watchdog 可用，轮询也会同时运行，因为 NFS 挂载上的写入不会触发 inotify——而
图谱级存储常用 NFS。新文件会反复检查大小直到稳定后才读取，因此不会解析写入中的
`.h5ad`。

Stop with Ctrl-C.

---

## 7. Interpreting the result / 结果解读

Each stage prints a summary and sets an exit code:

```
qc: 812 processed, 14 skipped, 2 failed (of 828)
```

| code | meaning |
|---|---|
| 0 | every sample processed or deliberately skipped |
| 1 | at least one sample failed |
| 2 | bad usage or configuration; nothing ran |

**Skipped** is not a failure. A sample is skipped when its output already exists,
or when it cannot meaningfully be processed — too few barcodes after QC, no
non-empty genes, no reference for its tissue, no spatial coordinates to plot. The
reason is logged for each one.

**跳过**不是失败。当输出已存在，或样本无法被有意义地处理时（QC 后 barcode 太少、
没有非空基因、组织没有参考数据、没有空间坐标可绘图），样本会被跳过，并记录原因。

**Failed** means an exception. The last 20 failures are printed to stderr with
their reason; the full traceback goes to `log_file` at DEBUG level:

```bash
stcompass qc --config atlas.yaml --log-file logs/qc.log -v
```

Chaining stages in a shell script works as expected, since a failure stops the
chain:

```bash
set -e
stcompass qc       --config atlas.yaml
stcompass cluster  --config atlas.yaml
stcompass annotate --config atlas.yaml
stcompass programs --config atlas.yaml
stcompass plot     --config atlas.yaml
```

---

## Python API / Python 接口

Every stage is also a function, and the per-sample logic is separable from the
file loop — useful in a notebook:

```python
from stcompass import load_config
from stcompass.pipelines import run_qc
from stcompass.pipelines.qc import qc_sample
from stcompass.platforms import Resolution

config = load_config("atlas.yaml")

# whole tree
result = run_qc(config)
print(result.summary())

# one object already in memory
processed, reason = qc_sample(
    adata,
    config=config.qc,
    resolution=Resolution.SPOT,
    sample_name="S1",
)
if processed is None:
    print("rejected:", reason)
```

See [pipelines.md](pipelines.md) for the full API of each stage.

---

## See also

- [configuration.md](configuration.md) — every config key
- [pipelines.md](pipelines.md) — stage internals and API
- [installation.md](installation.md) — extras and dependencies
- [migration.md](migration.md) — mapping from the original scripts
