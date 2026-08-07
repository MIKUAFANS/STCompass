# Configuration reference / 配置参考

Every tunable parameter lives in a YAML file. Nothing is hard-coded in the source,
which is the point: the config file is a small artefact you commit next to your
results, whereas an edited-in-place script is not reproducible.

所有可调参数都在 YAML 文件中，源码里没有硬编码路径。配置文件可以和结果一起提交到版本
控制，而"改一改脚本再跑"的做法无法复现。

Unknown keys are **rejected**, not ignored — a silently dropped `min_count` typo
would produce plausible-looking but wrong results across an entire atlas.

未知键会**直接报错**而不是被忽略：一个拼错的 `min_count` 若被静默丢弃，整个图谱都会得到
"看起来合理但其实错误"的结果。

```bash
stcompass qc --config configs/example.yaml   # loads and validates before running
```

Three ready-made files ship in `configs/`:

| File | Purpose |
|------|---------|
| `minimal.yaml` | Smallest thing that runs — just paths. |
| `example.yaml` | Every key, annotated. Copy this and edit. |
| `reproduce-atlas.yaml` | Settings matching the original scripts' behaviour. |

---

## Precedence / 优先级

Command-line flags beat the file; the file beats the defaults.

```
built-in defaults  <  YAML file  <  command-line flag
```

Unset flags do **not** clobber file values, so this works as expected:

```bash
# n_top_genes stays at whatever atlas.yaml says; only n_jobs is overridden
stcompass qc --config atlas.yaml --n-jobs 8
```

---

## `paths` — filesystem roots

Each stage reads one root and mirrors the relative path into another. A sample at
`raw/Homo sapiens/Brain/S1.h5ad` becomes `qc/Homo sapiens/Brain/S1.h5ad`.

每个阶段从一个根目录读取，并把相对路径镜像到另一个根目录下。

| Key | Used by | Meaning |
|-----|---------|---------|
| `raw` | `qc`, `annotate` | Input `.h5ad` files. |
| `qc` | `qc` (out), `cluster`, `programs` | QC output. |
| `clustered` | `cluster` (out) | Clustering output. |
| `annotated` | `annotate` (out), `plot` | Annotation output. |
| `programs` | `programs` (out) | Gene-program output. |
| `figures` | `plot` (out) | Rendered figures. |
| `reference` | `annotate` | scRNA-seq references, laid out `<species>/<tissue>/*.h5ad`. |
| `metadata` | `annotate` | Sample sheet (`.xlsx` or `.csv`). |

`~` and `$VARS` are expanded. A stage that needs a path you did not set fails
immediately with a message naming the key, before any sample is read:

```
error: paths.reference is required for this command but was not set.
```

---

## Top-level keys

| Key | Default | Meaning |
|-----|---------|---------|
| `n_jobs` | `1` | Worker processes for stages that don't override it. |
| `overwrite` | `false` | Re-process samples that already have output. Off by default so an interrupted run **resumes** instead of redoing days of work. |
| `log_file` | unset | Path for a DEBUG-level log of the run. |

---

## `qc` — quality control

Thresholds are chosen **per sample** from the platform's resolution class, because
one cutoff cannot serve both assay types: a Visium barcode pools several cells and
yields thousands of counts, while a MERFISH cell is measured on a ~300-gene panel
and may legitimately carry twenty.

阈值按平台分辨率**逐样本**选择：Visium 一个 barcode 覆盖多个细胞、计数上千；而 MERFISH
单个细胞只测约 300 个基因，计数二十也是正常的。同一套阈值无法同时适用。

| Key | Default | Meaning |
|-----|---------|---------|
| `min_counts_spot` | `100` | Min counts per barcode, spot-based. |
| `min_genes_spot` | `30` | Min genes per barcode, spot-based. |
| `min_counts_single_cell` | `20` | Min counts per cell, imaging-based. |
| `min_genes_single_cell` | `10` | Min genes per cell, imaging-based. |
| `min_cells_per_gene` | `5` | Drop genes seen in fewer barcodes. **Spot-based only** — on a targeted panel every probe was chosen deliberately, so prevalence filtering discards real signal. |
| `min_spots_after_qc` | `50` | Reject a spot-based sample left with fewer barcodes. |
| `min_cells_after_qc` | `10` | Reject an imaging-based sample left with fewer cells. |
| `filter_cells` | `true` | Apply the per-barcode filters at all. |
| `filter_genes` | `true` | Apply `min_cells_per_gene`. |
| `target_sum` | `10000.0` | Counts per barcode after library-size normalisation. |
| `n_top_genes` | `3000` | HVGs before PCA. Samples with fewer genes use all of them. |
| `hvg_flavor` | `seurat` | Passed to `scanpy.pp.highly_variable_genes`. |
| `n_pcs` | `50` | PCs for the neighbour graph. |
| `n_neighbors` | `15` | Neighbours for the kNN graph. |
| `cluster_method` | `leiden` | `leiden` or `louvain`. |
| `resolution` | `1.0` | Clustering resolution during QC. |
| `compute_umap` | `true` | Compute UMAP. Disable for very large samples where it dominates runtime. |
| `integer_check_sample_size` | `2000` | Values sampled to decide whether the matrix holds raw counts. |

### Reproducing the original atlas

`filter_cells` and `filter_genes` were **commented out** in `QC_ALL0917.py`, so the
published atlas kept every barcode. They are implemented and enabled here. To
reproduce the original output:

原脚本中这两个过滤是被注释掉的，因此已发布的图谱保留了所有 barcode。本包实现并默认启用
它们；若要复现原始输出：

```yaml
qc:
  filter_cells: false
  filter_genes: false
```

That is what `configs/reproduce-atlas.yaml` does.

---

## `clustering` — standalone re-clustering

Without an explicit `resolution`, one is chosen from the number of barcodes. Small
sections need a coarse graph to yield interpretable domains; large ones fragment
into noise unless the resolution drops.

不指定 `resolution` 时按 barcode 数自动选择：小切片需要较粗的图才能得到可解释的区域，
大切片若不降低分辨率就会碎裂成噪声。

| Key | Default | Meaning |
|-----|---------|---------|
| `method` | `louvain` | `leiden` or `louvain`. |
| `resolution` | unset | Fixed resolution; overrides the schedule. |
| `resolution_schedule` | see below | `(max_cells, resolution)` pairs, first match wins (strict `<`). |
| `fallback_resolution` | `0.2` | For samples larger than every entry. |
| `key_added` | = `method` | Destination `obs` column. |
| `recompute_neighbors` | `false` | Rebuild the kNN graph instead of reusing QC's. |
| `n_pcs` / `n_neighbors` | `50` / `15` | Used only when recomputing. |

The default schedule is the table hard-coded in the original `run_louvain.py`:

```yaml
clustering:
  resolution_schedule:
    - [100, 1.2]      # n_obs <   100  ->  1.2
    - [500, 0.7]      # n_obs <   500  ->  0.7
    - [5000, 0.5]     # n_obs <  5000  ->  0.5
    - [20000, 0.3]    # n_obs < 20000  ->  0.3
  fallback_resolution: 0.2
```

---

## `annotation` — cell-type labelling

The method follows the platform, not a user choice. Spot-based data is a
*deconvolution* problem (Tangram, proportions); single-cell data is a
*classification* problem (SingleR, labels).

方法由平台决定而非用户选择：spot 数据是*解卷积*问题（Tangram，比例），单细胞数据是
*分类*问题（SingleR，标签）。

| Key | Default | Meaning |
|-----|---------|---------|
| `label_key` | `cell_ontology_class` | Reference `obs` column holding cell types. |
| `min_cells_per_label` | `2` | Drop reference labels with fewer cells. Singletons break `rank_genes_groups` (no within-group variance). |
| `n_marker_genes` | `100` | Top markers per label, unioned into the Tangram gene set. |
| `tangram_mode` | `cells` | `cells`, `clusters` or `constrained`. |
| `tangram_epochs` | `300` | Training epochs. |
| `tangram_density_prior` | `rna_count_based` | Use `uniform` for equal-area spots. |
| `device` | `auto` | `auto`, `cpu`, `cuda`, or explicit `cuda:N`. |
| `n_gpus` | `1` | GPUs to spread samples across. |
| `n_jobs` | `1` | Worker processes; shadows the global `n_jobs`. |
| `singler_threads` | `4` | Threads for SingleR. |
| `species` | `[]` | Restrict to these species; empty means all. |

Sample-sheet column names (defaults match the atlas sheet from the original
scripts):

| Key | Default |
|-----|---------|
| `sample_column` | `SampleName` |
| `species_column` | `OrganismSimple` |
| `tissue_column` | `Tissue` |
| `platform_column` | `Biotech` |
| `category_column` | `Biotech Categories` |
| `category_filter` | `Spatial Transcriptomics` |

Set `category_column: null` to disable the category filter.

### GPU assignment

With `n_gpus > 1`, samples are assigned to devices by a **hash of their path**
rather than round-robin, so a resumed run sends each sample to the same device as
before. Round-robin would reshuffle assignments after a restart and change which
samples contend for memory.

`n_gpus > 1` 时按路径哈希分配 GPU（而非轮询），因此续跑时每个样本仍落到同一张卡上。
轮询会在重启后重新洗牌，改变显存竞争关系。

---

## `programs` — NMF gene programs

`MiniBatchNMF` over row blocks read from a *backed* `.h5ad`, which is what lets
million-cell samples factorise without densifying the matrix.

对 backed `.h5ad` 按行块做 `MiniBatchNMF`，因此百万细胞样本无需将矩阵稠密化即可分解。

| Key | Default | Meaning |
|-----|---------|---------|
| `n_components` | `10` | Number of programs (the `K`). |
| `max_hvg` | unset | Restrict to this many high-variance genes; unset uses every non-empty gene. |
| `batch_size` | `4096` | Rows per `partial_fit` call. |
| `epochs` | `2` | Passes over the matrix. One pass leaves factors biased toward the rows seen last. |
| `row_chunk` | `5000` | Rows read per block from the backed file. |
| `use_counts_layer` | `true` | Prefer `layers['counts']` over `X`, so the factorisation sees raw counts even after `X` was log-normalised. |
| `random_state` | `0` | Seed. |
| `save_float32` | `true` | Store factors as `float32`; halves file size. |
| `init` | `nndsvda` | Sparse-friendly SVD-based initialisation. |
| `max_no_improvement` | `20` | Stop after this many batches without improvement. |
| `large_sample_warning` | `1000000` | Log a warning above this many barcodes. |

---

## `plot` — figure rendering

| Key | Default | Meaning |
|-----|---------|---------|
| `format` | `pdf` | Output extension: `pdf`, `png`, `svg`. |
| `dpi` | `300` | Raster resolution. |
| `figsize` | `[8.0, 8.0]` | Inches. |
| `top_n_types` | `3` | Cell types per spot in a scatter-pie. Showing all thirty turns each pie into a ring of slivers. |
| `proportions_key` | `tangram_ct_pred` | `obsm` key holding Tangram proportions. |
| `label_key` | `singler_best` | `obs` column holding SingleR labels. |
| `spatial_key` | `spatial` | `obsm` key holding coordinates. |
| `radius` | unset | Pie radius in data units. Unset derives it from nearest-neighbour spacing, so the same code works for 55 µm Visium spots and 2 µm HD bins. |
| `invert_y` | `false` | Flip the y-axis for image-style coordinates. |
| `max_pies` | `20000` | Above this many spots, fall back to a plain scatter — hundreds of thousands of wedges exhaust memory and are individually invisible anyway. |
| `watch` | `false` | Keep running and render files as they appear. |
| `poll_seconds` | `5.0` | Polling interval in watch mode. Polling runs even when `watchdog` is installed, because inotify does not fire for writes on NFS mounts. |
| `stability_checks` | `3` | Consecutive equal-size checks before a new file is considered fully written. |
| `stability_interval` | `2.0` | Seconds between those checks. |

---

## A YAML gotcha: exponent notation

PyYAML implements YAML **1.1**, where an exponent without an explicit sign is a
*string*, not a float:

PyYAML 遵循 YAML **1.1**：指数不带符号时被解析为**字符串**而非浮点数：

```yaml
target_sum: 1e4       # -> the string "1e4"
target_sum: 1.0e+4    # -> the float 10000.0   (note the +)
target_sum: 10000.0   # -> the float 10000.0   (clearest)
```

The config loader coerces such strings to numbers rather than failing deep inside a
dataclass, and a genuinely non-numeric value gets an error naming the key:

配置加载器会把这类字符串转换为数字，而不是在 dataclass 内部抛出难以定位的 `TypeError`；
真正非数值的输入会得到指明键名的错误：

```
error: qc.target_sum must be a number, got 'abc'.
       Note that YAML needs a signed exponent: write 1.0e+4 or 10000.
```

Prefer plain decimals (`10000.0`) to sidestep the issue entirely.

---

## Validating a config

`--dry-run` loads and validates the config, resolves paths, and lists the samples a
stage would process — without importing scanpy or touching a single file:

```bash
stcompass qc --config configs/example.yaml --dry-run
```

See [usage.md](usage.md) for the full walkthrough and [pipelines.md](pipelines.md)
for what each stage writes into the `.h5ad`.
