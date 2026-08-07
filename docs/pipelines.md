# Pipeline reference / 流程详解

What each stage computes, which keys it writes, and why the defaults are what
they are. For flags see [usage.md](usage.md); for config keys see
[configuration.md](configuration.md).

每个阶段计算什么、写入哪些字段、以及默认值的依据。

---

## The one distinction that drives everything / 贯穿全局的核心区分

Every stage branches on a single question: **does one measurement unit hold one
cell, or many?**

| | Spot-based / 点阵式 | Single-cell / 单细胞式 |
|---|---|---|
| Platforms | Visium, Visium HD, ST, Slide-seq(V2), Stereo-seq, HDST, Well-ST-seq, CBSST-seq, sci-Space, Pixel-seq | MERFISH, seqFISH(+), osmFISH, EASI-FISH, EEL-FISH, STARmap, Xenium, CosMx, ExSeq |
| One barcode = | a capture area with several cells | one segmented cell |
| Genes measured | transcriptome-wide (~20 000) | targeted panel (~100–1000) |
| Typical counts/unit | thousands | tens to hundreds |
| QC floor | high (100 counts / 30 genes) | low (20 counts / 10 genes) |
| Gene prevalence filter | yes — a gene in <5 spots is noise | **no** — every probe was chosen deliberately |
| Annotation | **deconvolution** → proportions | **classification** → one label |

`stcompass platforms` prints the registry. `stcompass platforms --check "Stereo Seq"`
resolves a single label.

一个条码代表一个细胞还是多个细胞——这个区分决定了 QC 阈值和注释方法。

---

## `qc` — quality control, normalisation, embedding

Replaces `QC_ALL0917.py`.

### Steps

1. **Deduplicate.** Duplicate `var_names` break every `adata[:, gene]` lookup;
   duplicate `obs_names` make the file unwritable. First occurrence of each gene
   is kept, in original order.
2. **Preserve counts.** `adata.layers["counts"] = adata.X.copy()` *before* any
   transformation, so `programs` can factorise raw counts later.
3. **Metrics.** `n_counts`, `n_genes`, and `percent_mt` when genes matching
   `MT-*` (case-insensitive, so mouse `Mt-` also matches) are present.
4. **Filter.** Per-barcode count/gene floors by resolution; per-gene prevalence
   filter on spot-based platforms only.
5. **Normalise conditionally.** See below.
6. **HVG → PCA → neighbours → cluster → UMAP.**

### The conditional log1p

Public atlases mix raw and already-normalised matrices and carry no reliable
metadata flag. Applying `log1p` twice compresses real biological variation;
skipping it leaves counts on a scale that breaks PCA. So the stage *infers*:

```python
looks_like_counts(matrix, sample_size=2000)  # integral values → raw counts
```

Only 2000 stored values are sampled, making the check O(1) in matrix size. The
RNG is **seeded**, so the same file gets the same verdict on every run — the
original script drew from the global RNG and could classify a borderline file
differently on a re-run.

公开图谱混杂原始与已归一化矩阵且无可靠标记，因此通过抽样判断整数性来推断。抽样已固定随机种子，保证可重现。

### Writes

| Key | Meaning |
|---|---|
| `layers["counts"]` | raw counts before normalisation |
| `obs["n_counts"]`, `obs["n_genes"]` | per-barcode totals |
| `obs["percent_mt"]` | mitochondrial fraction (if detectable) |
| `var["mt"]` | mitochondrial flag |
| `var["highly_variable"]` | HVG selection |
| `obsm["X_pca"]`, `obsm["X_umap"]` | embeddings |
| `obs["leiden"]` / `obs["louvain"]` | cluster labels |

### Behaviour change from the original

The per-barcode and per-gene filters were **commented out** in
`QC_ALL0917.py`, so the published atlas kept every barcode. They are implemented
and enabled here. To reproduce the original output byte-for-byte:

```yaml
qc:
  filter_cells: false
  filter_genes: false
```

`configs/reproduce-atlas.yaml` already sets this.

原脚本中的过滤条件是被注释掉的；本包默认启用。如需复现原输出，设置上述两项为 `false`。

---

## `cluster` — standalone re-clustering

Replaces `run_louvain.py`.

Reuses the neighbour graph stored by `qc` unless `recompute_neighbors: true`.

### Size-aware resolution

A fixed resolution over-fragments large sections and under-splits small ones, so
resolution scales with barcode count — the schedule from the original script:

| n_obs | resolution |
|---|---|
| < 100 | 1.2 |
| < 500 | 0.7 |
| < 5 000 | 0.5 |
| < 20 000 | 0.3 |
| ≥ 20 000 | 0.2 |

Comparisons are strict `<`, so exactly 100 barcodes → 0.7. `--resolution` pins a
fixed value and bypasses the schedule.

固定分辨率会让大样本过度碎片化、小样本欠分割，因此分辨率随条码数变化。

### Writes

`obs["louvain"]` (or `obs["leiden"]`, or `clustering.key_added`).

---

## `annotate` — cell-type annotation

Replaces `cell_annotation_20251015_human.py` and
`cell_annotation_20251110_mouse.py`, which were the same pipeline twice —
differing only in a species string, a log file name, and one-GPU vs seven-GPU
execution. Species is now a config value.

### Reference layout

```
reference/
├── Homo sapiens/
│   └── Brain/ref.h5ad
└── Mus musculus/
    └── Brain/ref.h5ad
```

One `.h5ad` per `<species>/<tissue>`. Labels come from
`annotation.label_key` (default `cell_ontology_class`).

Labels with fewer than `min_cells_per_label` (default 2) cells are dropped:
a singleton has no within-group variance, which makes `rank_genes_groups` fail
outright, and it cannot support a proportion estimate anyway.

单细胞类型标签会被剔除——单个细胞无组内方差，会直接导致标记基因排序失败。

### Method dispatch

**Spot-based → Tangram.** Reference cells are mapped onto spots; the result is a
proportion vector per spot. The gene set is restricted to the union of the top
`n_marker_genes` (default 100) markers per reference label — using every gene
both slows the mapping and lets housekeeping genes dominate the objective.

Markers are ranked **once per `(species, tissue)` group**, not per sample.

**Single-cell → SingleR.** Each cell is classified against the reference,
yielding one label. Note SingleR expects genes × cells, the transpose of the
AnnData convention; the pipeline handles this.

### Grouping and indexing

Work is grouped by `(species, tissue)` because the expensive preparation —
loading the reference, ranking markers — is shared across the group.

The sample tree is indexed **once** into a `stem → paths` dict. The original
scripts called `glob(root/**/<sample>.h5ad)` once per sheet row, re-walking the
entire atlas for every sample; for a sheet with thousands of rows that walk
dominated runtime.

原脚本对每行样本都重新遍历整个目录树；这里改为一次性建立索引。

### GPU assignment

With `n_gpus > 1`, samples are assigned by **hash of their path**, not
round-robin, so a resumed run sends each sample to the same device as before.
Falls back to CPU when CUDA is unavailable, so a GPU config still runs on a
laptop.

### Writes

| Method | Key |
|---|---|
| Tangram | `obsm["tangram_ct_pred"]` — spots × cell types, proportions |
| SingleR | `obs["singler_best"]` (+ one `obs["singler_*"]` per output field) |

### Skips, not failures

A `(species, tissue)` group whose reference is missing is recorded as a **skip
for every sample in it**, and the run continues — an atlas legitimately contains
tissues with no matching reference. Counts stay honest in the summary.

---

## `programs` — NMF gene programs

Replaces `spatialVEG_0923.py`.

### Memory strategy

The sample is read **twice, by design**:

1. **Backed read** (`backed="r"`): the matrix stays on disk while gene selection
   streams over it in row blocks. A million-barcode sample costs no more memory
   than one block.
2. Only the *selected* columns are materialised in memory.
3. **Unbacked re-read** to attach factors and write the output.

Densifying a Visium HD matrix would cost `n_obs × n_vars × 4` bytes — hundreds of
gigabytes. `MiniBatchNMF.partial_fit` over row blocks avoids that entirely.

先以 backed 模式流式扫描选基因，只将所选列载入内存，避免稠密化整个矩阵。

### Why `epochs: 2`

A single pass leaves factors biased towards the rows seen last: each mini-batch
update only partially corrects the previous one. Two passes is the practical
floor.

### Gene importance

Per gene *i* and program *k*:

```
score[i,k] = V[i,k] * log(1 + V[i,k] / (competing[i,k] + 1e-10))
```

where `competing[i,k]` is the largest *other* program's weight for that gene
(for the gene's own top program, it is the runner-up). A gene loaded heavily on
one program and weakly elsewhere scores high; a gene loaded evenly across
programs scores low. That is the intended meaning of "marker of this program".

This is vectorised here; the original looped in Python over every gene. Verified
numerically identical to the original loop across 200 random matrices
(`tests/test_programs.py`).

某基因若集中负载于单一程序则得分高，若均匀分布于各程序则得分低。已验证与原始逐基因循环数值完全一致。

### Writes

| Key | Meaning |
|---|---|
| `obsm["X_nmf"]` | barcode × program loadings (W) |
| `var["NNMF_component_{k}"]` | gene weight in program *k* (H) |
| `var["NNMF_component_{k}_importance_score"]` | specificity score |
| `var["NNMF_component_{k}_rank"]` | rank 1 = most specific gene |

All-zero genes are dropped before factorising (they contribute nothing to the
loss and leave a zero row in H) and re-inserted as zeros afterwards, so the gene
axis still matches the input file.

### Guards

Skipped with a reason, not a crash: empty sample, no non-empty genes, fewer
usable genes than `n_components`, all-zero submatrix, negative values (NMF
requires non-negative input).

---

## `plot` — figures

Replaces `celltype_image1015.py`.

### Dispatch

| Input | Output |
|---|---|
| `obs["singler_best"]` present | categorical embedding (spatial → UMAP → PCA → t-SNE, first available) |
| `obsm["tangram_ct_pred"]` present | per-spot scatter-pie |
| neither | skipped with a reason |

### Scatter-pie details

**Top-N types.** A proportion vector often has 30+ entries; drawing all of them
turns each pie into a ring of slivers. Only the `top_n_types` (default 3)
largest are kept, then **renormalised** so every pie is a full circle and wedge
angles stay comparable between spots. Ties at the cut are all kept.

**Radius.** Spot pitch varies by three orders of magnitude (55 µm Visium spots,
2 µm HD bins, arbitrary post-registration units), so a fixed radius is wrong
everywhere. Default is 0.45 × median nearest-neighbour distance — adjacent pies
nearly touch, with a small gap so boundaries stay visible.

**`max_pies` (default 20 000).** Above this, the stage falls back to a plain
scatter coloured by dominant type: rendering hundreds of thousands of wedges
exhausts memory in the PDF backend and produces an unreadable figure anyway.

点位间距跨三个数量级，故半径由最近邻距离推导而非固定值。

### Watch mode

`--watch` keeps the command running and renders files as an annotation job
produces them.

**Polling always runs**, even when `watchdog` is installed, because inotify does
not fire for writes on NFS mounts — the common case when annotation runs on a
GPU node and plotting on a login node.

New files are only read once their size stops changing (`stability_checks`
consecutive equal sizes). A half-written `.h5ad` is readable but truncated, and
fails in ways that look like data corruption.

轮询始终启用：NFS 挂载上的写入不会触发 inotify 事件。

---

## Cross-cutting behaviour

### Resumability

Every stage mirrors its input tree into its output root and skips samples whose
output already exists. Interrupt and re-run; finished work is not repeated.
`--overwrite` forces reprocessing.

### Atomic writes

Output is written to a temporary file in the destination directory and moved into
place on success. A crash cannot leave a truncated file that a resumed run would
mistake for finished work, and watch mode never reads a partial file.

### Failure isolation

One bad sample does not stop a batch. Failures are caught per sample, logged with
a traceback at `DEBUG`, and counted. The exit code reflects the outcome: `0` all
processed or deliberately skipped, `1` at least one failed, `2` bad usage.

### Hardened h5ad round-tripping

Two recurring corruptions from files assembled by other tools, handled centrally:

- **`uns/log1p/base`** written in a form the reader cannot decode — deleted
  in place (it is metadata, not data) and the read retried.
- **`_index` column** inside `obs`/`var`, colliding with AnnData's reserved
  name — renamed to `cell_index`/`gene_index` rather than dropped, since it
  holds real identifiers.

---

## See also

- [usage.md](usage.md) — commands and flags
- [configuration.md](configuration.md) — every config key
- [migration.md](migration.md) — mapping from the original scripts
