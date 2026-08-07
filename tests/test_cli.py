"""Tests for the command-line interface.

These exercise argument parsing, config-override precedence, exit codes and
``--dry-run`` -- everything that happens *before* a pipeline touches scanpy.  The
runners themselves are monkeypatched, because what matters here is that the CLI
translates flags into the right :class:`~stcompass.config.Config` and turns a
:class:`~stcompass.pipelines._batch.BatchResult` into the right exit code.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from stcompass import cli, pipelines
from stcompass._deps import MissingDependencyError
from stcompass.config import ConfigError
from stcompass.pipelines._batch import BatchResult, SampleOutcome

#: CLI command -> the ``stcompass.pipelines`` attribute it dispatches to.
#: ``_run_stage`` looks these up on the package at call time, so patching the
#: attribute is what intercepts a stage without importing scanpy.
RUNNER_ATTRS = {
    "qc": "run_qc",
    "cluster": "run_clustering",
    "annotate": "run_annotation",
    "programs": "run_programs",
    "plot": "run_plot",
}


def _patch_runner(monkeypatch, command: str, function) -> None:
    """Replace the pipeline function that ``command`` dispatches to."""
    import stcompass.pipelines as pipelines

    monkeypatch.setattr(pipelines, RUNNER_ATTRS[command], function)


@pytest.fixture
def stub_runner(monkeypatch):
    """Replace every stage runner with a recorder that returns a clean result.

    Returns a dict that fills in with the config each runner received, so a test
    can assert on the resolved configuration rather than on side effects.
    """
    captured: dict[str, object] = {}

    def make(name):
        # ``annotate`` is called positionally without ``exclude``; the others pass
        # it as a keyword.  One signature covers both.
        def runner(config, exclude=()):
            captured["command"] = name
            captured["config"] = config
            captured["exclude"] = exclude
            return captured.get("result", BatchResult())

        return runner

    for command in RUNNER_ATTRS:
        _patch_runner(monkeypatch, command, make(command))
    return captured


# ---------------------------------------------------------------------------
# Parser structure
# ---------------------------------------------------------------------------


class TestParser:
    def test_no_command_prints_help_and_returns_usage_error(self, capsys):
        assert cli.main([]) == cli.EXIT_USAGE
        assert "COMMAND" in capsys.readouterr().out

    def test_version_flag_exits_zero(self, capsys):
        with pytest.raises(SystemExit) as excinfo:
            cli.main(["--version"])
        assert excinfo.value.code == 0
        assert "stcompass" in capsys.readouterr().out

    @pytest.mark.parametrize(
        "command", ["qc", "cluster", "annotate", "programs", "plot", "platforms"]
    )
    def test_every_command_has_help(self, command, capsys):
        with pytest.raises(SystemExit) as excinfo:
            cli.main([command, "--help"])
        assert excinfo.value.code == 0
        assert capsys.readouterr().out.strip()

    def test_unknown_command_is_a_usage_error(self):
        with pytest.raises(SystemExit) as excinfo:
            cli.main(["frobnicate"])
        assert excinfo.value.code == cli.EXIT_USAGE

    def test_verbose_and_quiet_are_mutually_exclusive(self):
        with pytest.raises(SystemExit):
            cli.main(["qc", "--verbose", "--quiet"])


# ---------------------------------------------------------------------------
# platforms subcommand
# ---------------------------------------------------------------------------


class TestPlatformsCommand:
    def test_lists_both_resolution_classes(self, capsys):
        assert cli.main(["platforms"]) == cli.EXIT_OK
        out = capsys.readouterr().out
        assert "10xVisium" in out
        assert "MERFISH" in out
        assert "Spot-based" in out
        assert "Single-cell" in out

    def test_check_resolves_a_label(self, capsys):
        assert cli.main(["platforms", "--check", "Stereo Seq"]) == cli.EXIT_OK
        assert "Stereo-seq" in capsys.readouterr().out

    def test_check_reports_unknown_label_with_nonzero_exit(self, capsys):
        assert cli.main(["platforms", "--check", "Bogus-seq"]) == cli.EXIT_FAILED_SAMPLES
        assert "not recognised" in capsys.readouterr().out

    def test_check_is_repeatable(self, capsys):
        code = cli.main(["platforms", "--check", "Visium", "--check", "Xenium"])
        out = capsys.readouterr().out
        assert code == cli.EXIT_OK
        assert "10xVisium" in out
        assert "10xXenium" in out


# ---------------------------------------------------------------------------
# Configuration resolution
# ---------------------------------------------------------------------------


class TestConfigResolution:
    def test_flags_alone_are_enough(self, tmp_path, stub_runner):
        raw = tmp_path / "raw"
        raw.mkdir()
        code = cli.main(["qc", "--raw-dir", str(raw), "--out", str(tmp_path / "qc")])
        assert code == cli.EXIT_OK
        config = stub_runner["config"]
        assert config.paths.raw == raw
        assert config.paths.qc == tmp_path / "qc"

    def test_flag_overrides_config_file(self, tmp_path, stub_runner):
        config_file = tmp_path / "atlas.yaml"
        config_file.write_text(
            "paths:\n"
            f"  raw: {tmp_path / 'from_file'}\n"
            f"  qc: {tmp_path / 'qc'}\n"
            "qc:\n"
            "  n_top_genes: 1000\n",
            encoding="utf-8",
        )
        override = tmp_path / "from_flag"
        code = cli.main(
            ["qc", "--config", str(config_file), "--raw-dir", str(override), "--n-top-genes", "500"]
        )
        assert code == cli.EXIT_OK
        config = stub_runner["config"]
        assert config.paths.raw == override
        assert config.qc.n_top_genes == 500

    def test_unset_flags_do_not_clobber_the_file(self, tmp_path, stub_runner):
        config_file = tmp_path / "atlas.yaml"
        config_file.write_text(
            "paths:\n"
            f"  raw: {tmp_path / 'raw'}\n"
            f"  qc: {tmp_path / 'qc'}\n"
            "qc:\n"
            "  n_top_genes: 1234\n"
            "  resolution: 0.42\n",
            encoding="utf-8",
        )
        assert cli.main(["qc", "--config", str(config_file)]) == cli.EXIT_OK
        config = stub_runner["config"]
        assert config.qc.n_top_genes == 1234
        assert config.qc.resolution == pytest.approx(0.42)

    def test_missing_config_file_is_a_usage_error(self, capsys):
        code = cli.main(["qc", "--config", "/nonexistent/atlas.yaml"])
        assert code == cli.EXIT_USAGE
        assert "error:" in capsys.readouterr().err

    def test_unknown_config_key_is_a_usage_error(self, tmp_path, capsys):
        config_file = tmp_path / "bad.yaml"
        config_file.write_text("qc:\n  no_such_setting: 1\n", encoding="utf-8")
        assert cli.main(["qc", "--config", str(config_file)]) == cli.EXIT_USAGE
        assert "no_such_setting" in capsys.readouterr().err

    def test_missing_required_path_is_a_usage_error(self, capsys):
        """paths.raw is required by qc; without it nothing should run."""
        assert cli.main(["qc"]) == cli.EXIT_USAGE
        assert "paths.raw" in capsys.readouterr().err

    def test_negation_flags_disable_defaults(self, tmp_path, stub_runner):
        raw = tmp_path / "raw"
        raw.mkdir()
        code = cli.main(
            [
                "qc",
                "--raw-dir",
                str(raw),
                "--out",
                str(tmp_path / "qc"),
                "--no-filter-cells",
                "--no-filter-genes",
                "--no-umap",
            ]
        )
        assert code == cli.EXIT_OK
        config = stub_runner["config"]
        assert config.qc.filter_cells is False
        assert config.qc.filter_genes is False
        assert config.qc.compute_umap is False

    def test_n_jobs_also_sets_annotation_n_jobs(self, tmp_path, stub_runner):
        """The GPU stage reads its own n_jobs; -j must reach both."""
        raw = tmp_path / "raw"
        raw.mkdir()
        sheet = tmp_path / "samples.csv"
        sheet.write_text("SampleName,Tissue,OrganismSimple,Biotech\n", encoding="utf-8")
        code = cli.main(
            [
                "annotate",
                "--raw-dir",
                str(raw),
                "--out",
                str(tmp_path / "out"),
                "--reference-dir",
                str(tmp_path / "ref"),
                "--metadata",
                str(sheet),
                "-j",
                "4",
            ]
        )
        assert code == cli.EXIT_OK
        config = stub_runner["config"]
        assert config.n_jobs == 4
        assert config.annotation.n_jobs == 4

    def test_species_flag_accumulates(self, tmp_path, stub_runner):
        raw = tmp_path / "raw"
        raw.mkdir()
        sheet = tmp_path / "samples.csv"
        sheet.write_text("SampleName\n", encoding="utf-8")
        code = cli.main(
            [
                "annotate",
                "--raw-dir",
                str(raw),
                "--out",
                str(tmp_path / "out"),
                "--reference-dir",
                str(tmp_path / "ref"),
                "--metadata",
                str(sheet),
                "--species",
                "Homo sapiens",
                "--species",
                "Mus musculus",
            ]
        )
        assert code == cli.EXIT_OK
        assert stub_runner["config"].annotation.species == ["Homo sapiens", "Mus musculus"]

    def test_exclude_is_passed_through_as_a_tuple(self, tmp_path, stub_runner):
        raw = tmp_path / "raw"
        raw.mkdir()
        code = cli.main(
            [
                "qc",
                "--raw-dir",
                str(raw),
                "--out",
                str(tmp_path / "qc"),
                "--exclude",
                "bad1.h5ad",
                "--exclude",
                "bad2.h5ad",
            ]
        )
        assert code == cli.EXIT_OK
        assert stub_runner["exclude"] == ("bad1.h5ad", "bad2.h5ad")

    def test_overwrite_reaches_the_config(self, tmp_path, stub_runner):
        raw = tmp_path / "raw"
        raw.mkdir()
        cli.main(["qc", "--raw-dir", str(raw), "--out", str(tmp_path / "qc"), "--overwrite"])
        assert stub_runner["config"].overwrite is True

    def test_overwrite_defaults_to_false(self, tmp_path, stub_runner):
        raw = tmp_path / "raw"
        raw.mkdir()
        cli.main(["qc", "--raw-dir", str(raw), "--out", str(tmp_path / "qc")])
        assert stub_runner["config"].overwrite is False


# ---------------------------------------------------------------------------
# Stage routing
# ---------------------------------------------------------------------------


class TestStageRouting:
    @pytest.mark.parametrize(
        ("command", "flags", "expected_paths"),
        [
            ("qc", ["--raw-dir", "IN", "--out", "OUT"], ("raw", "qc")),
            ("cluster", ["--qc-dir", "IN", "--out", "OUT"], ("qc", "clustered")),
            ("programs", ["--qc-dir", "IN", "--out", "OUT"], ("qc", "programs")),
            ("plot", ["--annotated-dir", "IN", "--out", "OUT"], ("annotated", "figures")),
        ],
    )
    def test_each_stage_reads_and_writes_its_own_paths(
        self, command, flags, expected_paths, tmp_path, stub_runner
    ):
        source = tmp_path / "in"
        source.mkdir()
        destination = tmp_path / "out"
        resolved = [
            str(source) if f == "IN" else str(destination) if f == "OUT" else f for f in flags
        ]

        assert cli.main([command, *resolved]) == cli.EXIT_OK
        assert stub_runner["command"] == command

        config = stub_runner["config"]
        in_attr, out_attr = expected_paths
        assert getattr(config.paths, in_attr) == source
        assert getattr(config.paths, out_attr) == destination


# ---------------------------------------------------------------------------
# Exit codes
# ---------------------------------------------------------------------------


class TestExitCodes:
    def _run_qc(self, tmp_path):
        raw = tmp_path / "raw"
        raw.mkdir()
        return cli.main(["qc", "--raw-dir", str(raw), "--out", str(tmp_path / "qc")])

    def test_clean_run_returns_zero(self, tmp_path, stub_runner):
        assert self._run_qc(tmp_path) == cli.EXIT_OK

    def test_failed_samples_return_one(self, tmp_path, stub_runner, capsys):
        result = BatchResult()
        result.failed.append(SampleOutcome(Path("bad.h5ad"), "failed", "boom"))
        stub_runner["result"] = result
        assert self._run_qc(tmp_path) == cli.EXIT_FAILED_SAMPLES
        assert "failed" in capsys.readouterr().err

    def test_skipped_samples_still_return_zero(self, tmp_path, stub_runner):
        """A deliberate skip is not a failure."""
        result = BatchResult()
        result.skipped.append(SampleOutcome(Path("small.h5ad"), "skipped", "too few cells"))
        stub_runner["result"] = result
        assert self._run_qc(tmp_path) == cli.EXIT_OK

    def test_summary_is_printed_to_stdout(self, tmp_path, stub_runner, capsys):
        result = BatchResult()
        result.processed.append(Path("a.h5ad"))
        stub_runner["result"] = result
        self._run_qc(tmp_path)
        assert "1 processed" in capsys.readouterr().out

    def test_missing_dependency_is_a_usage_error(self, tmp_path, monkeypatch, capsys):
        def boom(config, exclude=()):
            raise MissingDependencyError("needs scanpy: pip install stcompass")

        monkeypatch.setattr(pipelines, "run_qc", boom)
        assert self._run_qc(tmp_path) == cli.EXIT_USAGE
        assert "scanpy" in capsys.readouterr().err

    def test_missing_input_directory_is_a_usage_error(self, tmp_path, monkeypatch, capsys):
        def boom(config, exclude=()):
            raise FileNotFoundError("input directory does not exist: /nope")

        monkeypatch.setattr(pipelines, "run_qc", boom)
        assert self._run_qc(tmp_path) == cli.EXIT_USAGE
        assert "does not exist" in capsys.readouterr().err

    def test_config_error_from_runner_is_a_usage_error(self, tmp_path, monkeypatch, capsys):
        def boom(config, exclude=()):
            raise ConfigError("paths.reference is required")

        monkeypatch.setattr(pipelines, "run_qc", boom)
        assert self._run_qc(tmp_path) == cli.EXIT_USAGE
        assert "paths.reference" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# --dry-run
# ---------------------------------------------------------------------------


class TestDryRun:
    def test_lists_samples_without_running(self, tmp_atlas, tmp_path, capsys, monkeypatch):
        def explode(config, exclude=()):  # pragma: no cover - must not be called
            raise AssertionError("--dry-run must not invoke the pipeline")

        monkeypatch.setattr(pipelines, "run_qc", explode)
        code = cli.main(
            ["qc", "--raw-dir", str(tmp_atlas), "--out", str(tmp_path / "qc"), "--dry-run"]
        )
        assert code == cli.EXIT_OK
        out = capsys.readouterr().out
        assert "S1.h5ad" in out
        assert "5" in out  # five .h5ad files in the fixture

    def test_reports_an_empty_tree(self, tmp_path, capsys):
        empty = tmp_path / "empty"
        empty.mkdir()
        code = cli.main(["qc", "--raw-dir", str(empty), "--out", str(tmp_path / "qc"), "--dry-run"])
        assert code == cli.EXIT_OK
        assert "0" in capsys.readouterr().out

    def test_missing_input_directory_is_reported(self, tmp_path, capsys):
        code = cli.main(
            ["qc", "--raw-dir", str(tmp_path / "nope"), "--out", str(tmp_path / "qc"), "--dry-run"]
        )
        assert code == cli.EXIT_USAGE
        assert "error:" in capsys.readouterr().err

    def test_dry_run_respects_exclude(self, tmp_atlas, tmp_path, capsys):
        code = cli.main(
            [
                "qc",
                "--raw-dir",
                str(tmp_atlas),
                "--out",
                str(tmp_path / "qc"),
                "--dry-run",
                "--exclude",
                "S1.h5ad",
            ]
        )
        assert code == cli.EXIT_OK
        assert "S1.h5ad" not in capsys.readouterr().out
