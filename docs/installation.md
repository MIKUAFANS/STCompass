# Installation / 安装

## Requirements / 环境要求

| Item | Requirement |
|------|-------------|
| Python | 3.10 or newer |
| OS | Linux, macOS (Windows untested) |
| Core stack | numpy, scipy, pandas, scanpy, anndata, scikit-learn |
| Optional | torch + tangram-sc (GPU), singler, watchdog, openpyxl |

The core install pulls in `scanpy`, which is the heaviest hard dependency. The
GPU-based annotation path is *not* installed by default, because a CUDA build of
`torch` cannot be resolved on every machine and would make the package
uninstallable on a laptop.

核心安装包含 `scanpy`。基于 GPU 的注释路径不是默认安装的——CUDA 版 `torch`
无法在所有机器上解析，若设为必需依赖会导致本包在笔记本上无法安装。

## Install / 安装

### From a clone / 从克隆仓库安装

```bash
git clone https://github.com/OWNER/STCompass.git
cd STCompass
pip install -e .
```

`-e` (editable) is the right choice while you are still adjusting thresholds:
edits to the source take effect without reinstalling.

### Extras / 可选依赖组

Install only what your data needs:

```bash
pip install -e '.[excel]'      # .xlsx sample sheets
pip install -e '.[tangram]'    # spot deconvolution (needs CUDA torch)
pip install -e '.[singler]'    # single-cell classification
pip install -e '.[watch]'      # `stcompass plot --watch`
pip install -e '.[dev]'        # pytest + ruff, for contributing
pip install -e '.[all]'        # everything optional
```

| Extra | Pulls in | Needed for |
|-------|----------|-----------|
| `excel` | `openpyxl` | reading `.xlsx` sample sheets in `annotate` |
| `tangram` | `tangram-sc`, `torch` | `annotate` on spot-based platforms |
| `singler` | `singler` | `annotate` on imaging-based platforms |
| `watch` | `watchdog` | `plot --watch` (falls back to polling without it) |
| `dev` | `pytest`, `pytest-cov`, `ruff` | running the test suite |

If you skip an extra and then invoke a stage that needs it, the command exits
with code 2 and prints the exact `pip install` line — it does not fail with an
`ImportError` traceback.

若跳过某个 extra 后又调用了需要它的阶段，命令以退出码 2 结束，并打印所需的
`pip install` 命令，而不是抛出 `ImportError` 堆栈。

## GPU setup / GPU 配置

Tangram needs a CUDA-enabled `torch`. Install it from the PyTorch index that
matches your driver *before* installing the extra:

```bash
# check your driver first
nvidia-smi

# example: CUDA 12.1 wheels
pip install torch --index-url https://download.pytorch.org/whl/cu121
pip install -e '.[tangram]'
```

Verify that the GPU is visible to the package:

```bash
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.device_count())"
```

`annotate` falls back to CPU with a warning when CUDA is unavailable, so a config
written for a GPU box still runs elsewhere — just far more slowly.

当 CUDA 不可用时，`annotate` 会带警告回退到 CPU：为 GPU 机器写的配置在别处
依然能跑，只是慢很多。

## Verify the install / 验证安装

```bash
stcompass --version
stcompass platforms                       # lists the recognised platforms
stcompass platforms --check 'Stereo Seq'  # resolve one label
```

Then check a stage's wiring without processing anything:

```bash
stcompass qc --config configs/example.yaml --dry-run
```

`--dry-run` prints the resolved input/output roots and the samples that would be
processed. It does not import scanpy, so it also works as a fast configuration
check on a machine without the full stack.

`--dry-run` 打印解析后的输入/输出根目录与将被处理的样本。它不导入 scanpy，
因此在没有完整依赖栈的机器上也可作为快速配置检查。

## Reproducible environments / 可复现环境

Dependencies are declared with lower bounds, not exact pins. Hard pins across a
stack this large reliably produce unsolvable environments on machines that
already carry a scientific Python install.

For a run you intend to reproduce, capture the resolved versions yourself:

```bash
pip freeze > requirements-lock.txt
```

and commit that file next to the config you used. The config file records the
parameters; the lock file records the code that consumed them.

依赖使用下界而非精确锁定版本——在这种规模的科学栈上硬锁版本极易导致无法求解的
环境。若需复现某次运行，请自行用 `pip freeze` 导出锁文件，并与所用配置一并提交：
配置文件记录参数，锁文件记录消费这些参数的代码。

## Troubleshooting / 故障排查

**`error: ... requires the 'tangram' package`**
The extra is not installed. Run the `pip install` line the message prints.

**`CUDA out of memory` during `annotate`**
Lower `annotation.n_jobs` (fewer samples resident per GPU at once), or set
`annotation.device: cpu`. Sample-to-GPU assignment is by path hash, so re-running
sends each sample to the same device and the failure is reproducible.

**scanpy import warnings about `igraph` / `leidenalg`**
The QC stage pins `flavor="igraph"` for Leiden clustering so labels are stable
across scanpy versions. Install `igraph` and `leidenalg` if scanpy reports them
missing.

**Reading an `.h5ad` fails with an `uns/log1p/base` error**
Handled automatically: the reader removes the undecodable entry in place and
retries once. The entry is metadata about a previous transform, not matrix data,
so removing it is lossless here.

**读取 `.h5ad` 时报 `uns/log1p/base` 错误**
已自动处理：读取器会原地删除该无法解码的条目并重试一次。该条目是关于既往变换的
元数据，并非矩阵数据，因此删除它在此处无损。

## Next / 下一步

- [usage.md](usage.md) — the five stages, end to end
- [configuration.md](configuration.md) — every config key
- [migration.md](migration.md) — moving from the original scripts
