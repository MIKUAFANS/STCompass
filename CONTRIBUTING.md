# Contributing to STCompass

Thanks for your interest. This document covers the development setup, what CI
checks, and the conventions this codebase follows.

感谢关注。本文档说明开发环境搭建、CI 检查项，以及本代码库遵循的约定。

## Development setup / 开发环境

```bash
git clone https://github.com/OWNER/STCompass.git
cd STCompass
python -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'
```

`[dev]` installs `pytest`, `pytest-cov`, `ruff` and `openpyxl`. It deliberately
does **not** install `scanpy`, `torch`, `tangram-sc` or `singler`: the test suite
is designed to pass without them, so a contributor can work on the configuration
layer, the CLI or the matrix helpers without a multi-gigabyte install. Tests that
genuinely need scanpy are marked `requires_scanpy` and skip themselves.

`[dev]` 不安装 `scanpy`/`torch`/`tangram-sc`/`singler`——测试套件设计为在没有它们
的情况下也能通过，因此贡献者无需数 GB 的安装即可开发配置层、CLI 或矩阵工具。
确实需要 scanpy 的测试标记为 `requires_scanpy`，会自动跳过。

To run the full suite including those tests:

```bash
pip install -e '.[dev,all]'   # adds scanpy via the base deps, plus every extra
```

## Before opening a pull request / 提交 PR 前

Run exactly what CI runs:

```bash
ruff check src tests
ruff format --check src tests
pytest
```

All three must pass. If `ruff format --check` complains, run `ruff format src tests`.

三项都必须通过。若 `ruff format --check` 报错，运行 `ruff format src tests`。

## Conventions / 代码约定

**Comments explain why, not what.** The code says what it does. A comment earns
its place by recording a decision that is not recoverable from reading the code:

```python
# Sampling is done *with* replacement: drawing without replacement forces a
# permutation of every stored value, which for a billion-nonzero matrix costs
# more memory than the pipeline stage that follows.
```

**注释解释"为什么"，而非"是什么"。** 代码本身已说明它在做什么；注释的价值在于
记录无法从代码本身还原的决策。

**Heavy dependencies are imported lazily.** Use `stcompass._deps.require`, which
raises a `MissingDependencyError` naming the extra to install:

```python
scanpy = require("scanpy", feature="quality control")
```

Never import `scanpy`, `torch`, `tangram` or `singler` at module level — it would
make `import stcompass` fail for everyone who does not need that stage.

**重依赖必须惰性导入。** 使用 `stcompass._deps.require`，它会抛出指明所需 extra
的 `MissingDependencyError`。切勿在模块顶层导入这些包。

**A per-sample problem is a return value, not an exception.** A section with too
few cells, a sample with no spatial coordinates, a reference with no matching
tissue — these are expected outcomes when processing a heterogeneous atlas, not
errors. Return a short reason string; the batch driver records it as a skip and
keeps going. Reserve exceptions for genuine failures.

**单样本问题应作为返回值，而非异常。** 细胞过少的切片、缺少空间坐标的样本、
没有匹配组织的参考数据——在处理异质图谱时这些都是预期结果。返回简短的原因
字符串，批处理驱动会记为 skip 并继续。异常仅用于真正的失败。

**New configuration goes in a dataclass.** Add the field to the relevant section
in `src/stcompass/config.py` with a docstring explaining what it controls and why
the default is what it is. Validate it in `__post_init__`. Unknown keys are
rejected by design, so a new key needs no other registration — but if it should
be settable from the command line, add it to `_OVERRIDES` in `cli.py` too.

**新增配置项必须放入 dataclass。** 在 `config.py` 相应 section 中添加字段，
用 docstring 说明其作用及默认值的理由，并在 `__post_init__` 中校验。若需支持
命令行覆盖，同时在 `cli.py` 的 `_OVERRIDES` 中登记。

**Numerical changes need a test that would fail without them.** The gene
importance scoring, the top-N proportion renormalisation, the resolution
schedule and the count-matrix heuristic all have tests pinning their exact
behaviour, because a silent change in any of them would alter published results
without any error.

**数值行为的改动需要配套的、缺少改动就会失败的测试。**

## Testing / 测试

Tests live in `tests/`, one module per source module. Prefer testing the
single-sample function (`qc_sample`, `cluster_sample`, `programs_sample`) over
the batch runner, since those take an in-memory `AnnData` and need no
filesystem.

Mark anything needing scanpy:

```python
@pytest.mark.requires_scanpy
def test_something_with_scanpy():
    ...
```

## Reporting bugs / 报告问题

Include the `stcompass` version (`stcompass --version`), the command you ran, the
relevant part of your config file, and the error. For a pipeline failure, run
with `-v` and attach the log — `--log-file run.log` captures DEBUG-level detail
even when the console is quiet.

请附上 `stcompass --version`、执行的命令、配置文件相关片段与错误信息。
流程失败时请加 `-v` 运行并附上日志。

## Licence / 许可

Contributions are accepted under the MIT Licence, the same terms as the project.

贡献内容以与本项目相同的 MIT 许可协议接受。
