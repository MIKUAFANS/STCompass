"""Tests for the YAML configuration layer.

The emphasis is on *rejection*: a typo in a threshold silently ignored would
produce plausible-looking but wrong results across a whole atlas, so unknown keys
and out-of-range values must fail loudly at load time rather than at hour six of a
batch run.
"""

from __future__ import annotations

import pytest

from stcompass.config import (
    AnnotationConfig,
    ClusteringConfig,
    Config,
    ConfigError,
    PathsConfig,
    PlotConfig,
    ProgramsConfig,
    QCConfig,
    load_config,
)

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def test_defaults_load_without_a_file():
    config = load_config(None)
    assert config.n_jobs == 1
    assert config.overwrite is False
    assert config.qc.min_counts_spot == 100


def test_load_from_yaml(tmp_path):
    path = tmp_path / "atlas.yaml"
    path.write_text(
        """
paths:
  raw: /data/raw
  qc: /data/qc
qc:
  min_counts_spot: 250
  cluster_method: louvain
n_jobs: 4
""",
        encoding="utf-8",
    )
    config = load_config(path)
    assert config.paths.raw.as_posix() == "/data/raw"
    assert config.qc.min_counts_spot == 250
    assert config.qc.cluster_method == "louvain"
    assert config.n_jobs == 4
    # Untouched fields keep their defaults.
    assert config.qc.min_genes_spot == 30


def test_empty_yaml_is_all_defaults(tmp_path):
    path = tmp_path / "empty.yaml"
    path.write_text("", encoding="utf-8")
    assert load_config(path).qc.n_top_genes == 3000


def test_missing_file_is_an_error(tmp_path):
    with pytest.raises(ConfigError, match="not found"):
        load_config(tmp_path / "nope.yaml")


def test_non_mapping_yaml_is_an_error(tmp_path):
    path = tmp_path / "list.yaml"
    path.write_text("- a\n- b\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="mapping"):
        load_config(path)


# ---------------------------------------------------------------------------
# Rejection of unknown keys
# ---------------------------------------------------------------------------


def test_unknown_top_level_key_is_rejected():
    with pytest.raises(ConfigError, match="unknown top-level"):
        Config.from_mapping({"qq": 1})


def test_unknown_section_key_is_rejected_and_names_the_section():
    with pytest.raises(ConfigError, match=r"unknown key.*'qc'"):
        Config.from_mapping({"qc": {"min_count": 100}})


def test_unknown_key_error_lists_valid_alternatives():
    with pytest.raises(ConfigError) as info:
        Config.from_mapping({"qc": {"min_count": 100}})
    # The message should help the user find the real name.
    assert "min_counts_spot" in str(info.value)


def test_section_must_be_a_mapping():
    with pytest.raises(ConfigError, match="must be a mapping"):
        Config.from_mapping({"qc": [1, 2]})


def test_top_level_must_be_a_mapping():
    with pytest.raises(ConfigError, match="must be a mapping"):
        Config.from_mapping([1, 2])  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Dotted overrides (how the CLI applies flags)
# ---------------------------------------------------------------------------


def test_dotted_override_targets_a_nested_field():
    config = load_config(None, **{"qc.resolution": 0.75})
    assert config.qc.resolution == 0.75


def test_dotted_override_beats_the_file(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("qc:\n  resolution: 1.0\n", encoding="utf-8")
    config = load_config(path, **{"qc.resolution": 0.3})
    assert config.qc.resolution == 0.3


def test_none_overrides_are_ignored(tmp_path):
    """The CLI passes every flag unconditionally; unset ones must not clobber."""
    path = tmp_path / "c.yaml"
    path.write_text("qc:\n  resolution: 0.9\n", encoding="utf-8")
    config = load_config(path, **{"qc.resolution": None, "qc.n_top_genes": None})
    assert config.qc.resolution == 0.9


def test_dotted_override_creates_missing_sections():
    config = load_config(None, **{"plot.dpi": 150})
    assert config.plot.dpi == 150


def test_dotted_override_through_a_scalar_is_an_error(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("n_jobs: 2\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="not a mapping"):
        load_config(path, **{"n_jobs.deeper": 1})


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------


def test_paths_expand_user_and_env(monkeypatch):
    monkeypatch.setenv("STC_TEST_ROOT", "/mnt/atlas")
    paths = PathsConfig(raw="$STC_TEST_ROOT/raw")
    assert paths.raw.as_posix() == "/mnt/atlas/raw"


def test_require_returns_a_configured_path():
    paths = PathsConfig(raw="/data/raw")
    assert paths.require("raw").as_posix() == "/data/raw"


def test_require_explains_a_missing_path():
    with pytest.raises(ConfigError, match=r"paths\.raw is required"):
        PathsConfig().require("raw")


def test_unset_paths_stay_none():
    assert PathsConfig().qc is None


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_qc_rejects_an_unknown_cluster_method():
    with pytest.raises(ConfigError, match="cluster_method"):
        QCConfig(cluster_method="kmeans")


@pytest.mark.parametrize("field,value", [("n_top_genes", 0), ("target_sum", 0)])
def test_qc_rejects_out_of_range_values(field, value):
    with pytest.raises(ConfigError):
        QCConfig(**{field: value})


def test_clustering_rejects_an_unknown_method():
    with pytest.raises(ConfigError, match=r"clustering\.method"):
        ClusteringConfig(method="kmeans")


def test_annotation_rejects_an_unknown_tangram_mode():
    with pytest.raises(ConfigError, match="tangram_mode"):
        AnnotationConfig(tangram_mode="magic")


@pytest.mark.parametrize("field", ["n_gpus", "n_jobs"])
def test_annotation_rejects_non_positive_worker_counts(field):
    with pytest.raises(ConfigError):
        AnnotationConfig(**{field: 0})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n_components": 0},
        {"batch_size": 0},
        {"epochs": 0},
        {"n_components": 10, "max_hvg": 5},  # fewer genes than programs
    ],
)
def test_programs_rejects_impossible_settings(kwargs):
    with pytest.raises(ConfigError):
        ProgramsConfig(**kwargs)


def test_programs_allows_max_hvg_none():
    assert ProgramsConfig(max_hvg=None).max_hvg is None


def test_plot_rejects_a_bad_figsize():
    with pytest.raises(ConfigError, match="figsize"):
        PlotConfig(figsize=[8.0])  # type: ignore[arg-type]


def test_plot_rejects_zero_types_per_pie():
    with pytest.raises(ConfigError, match="top_n_types"):
        PlotConfig(top_n_types=0)


def test_config_rejects_zero_jobs():
    with pytest.raises(ConfigError, match="n_jobs"):
        Config(n_jobs=0)


# ---------------------------------------------------------------------------
# YAML type coercion
# ---------------------------------------------------------------------------


def test_figsize_from_yaml_list_becomes_a_tuple(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("plot:\n  figsize: [10, 6]\n", encoding="utf-8")
    figsize = load_config(path).plot.figsize
    assert figsize == (10.0, 6.0)
    assert isinstance(figsize, tuple)


def test_resolution_schedule_from_yaml_is_normalised(tmp_path):
    """YAML gives lists-of-lists, in arbitrary order; both must be handled."""
    path = tmp_path / "c.yaml"
    path.write_text(
        "clustering:\n  resolution_schedule:\n    - [5000, 0.5]\n    - [100, 1.5]\n",
        encoding="utf-8",
    )
    schedule = load_config(path).clustering.resolution_schedule
    assert schedule == [(100, 1.5), (5000, 0.5)]


# ---------------------------------------------------------------------------
# Size-aware resolution schedule
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n_obs,expected",
    [
        (10, 1.2),
        (99, 1.2),
        (100, 0.7),  # boundary: schedule uses a strict <
        (499, 0.7),
        (500, 0.5),
        (4_999, 0.5),
        (5_000, 0.3),
        (19_999, 0.3),
        (20_000, 0.2),  # past every entry -> fallback
        (5_000_000, 0.2),
    ],
)
def test_resolution_schedule_matches_the_original_script(n_obs, expected):
    """These are the thresholds hard-coded in the original ``run_louvain.py``."""
    assert ClusteringConfig().resolution_for(n_obs) == expected


def test_explicit_resolution_overrides_the_schedule():
    config = ClusteringConfig(resolution=0.8)
    assert config.resolution_for(10) == 0.8
    assert config.resolution_for(10_000_000) == 0.8


def test_custom_schedule_is_honoured():
    config = ClusteringConfig(resolution_schedule=[(50, 2.0)], fallback_resolution=0.1)
    assert config.resolution_for(10) == 2.0
    assert config.resolution_for(100) == 0.1


# ---------------------------------------------------------------------------
# YAML 1.1 numeric coercion
# ---------------------------------------------------------------------------
# PyYAML implements YAML 1.1, where an exponent without an explicit sign is a
# *string*: `target_sum: 1e4` loads as "1e4".  Before coercion that reached a
# dataclass comparison and surfaced as
# `TypeError: '<=' not supported between 'str' and 'int'`, naming no key.


@pytest.mark.parametrize("written", ["1e4", "1.0e4", " 1e4 ", "10000", "10000.0"])
def test_unsigned_exponent_strings_are_coerced_to_float(written):
    config = Config.from_mapping({"qc": {"target_sum": written}})
    assert config.qc.target_sum == 10000.0
    assert isinstance(config.qc.target_sum, float)


def test_numeric_strings_are_coerced_to_int_for_int_fields():
    config = Config.from_mapping({"qc": {"n_top_genes": "2000"}})
    assert config.qc.n_top_genes == 2000
    assert isinstance(config.qc.n_top_genes, int)


@pytest.mark.parametrize(
    "section,key,value",
    [
        ("qc", "target_sum", "abc"),
        ("qc", "n_top_genes", "many"),
        ("programs", "n_components", "ten"),
    ],
)
def test_non_numeric_strings_name_the_offending_key(section, key, value):
    with pytest.raises(ConfigError, match=rf"{section}\.{key} must be a number"):
        Config.from_mapping({section: {key: value}})


def test_error_message_explains_the_yaml_exponent_pitfall():
    with pytest.raises(ConfigError, match=r"1\.0e\+4"):
        Config.from_mapping({"qc": {"target_sum": "abc"}})


@pytest.mark.parametrize(
    "section,key",
    [("qc", "n_top_genes"), ("qc", "target_sum"), ("programs", "n_components")],
)
def test_booleans_are_rejected_on_numeric_fields(section, key):
    """`bool` is an int subclass; accepting `true` as 1 would hide a real typo."""
    with pytest.raises(ConfigError, match="boolean"):
        Config.from_mapping({section: {key: True}})


@pytest.mark.parametrize(
    "section,key,value",
    [
        ("qc", "filter_cells", False),
        ("qc", "compute_umap", True),
        ("plot", "invert_y", True),
    ],
)
def test_booleans_still_work_on_boolean_fields(section, key, value):
    config = Config.from_mapping({section: {key: value}})
    assert getattr(getattr(config, section), key) is value


@pytest.mark.parametrize("section,key", [("programs", "max_hvg"), ("plot", "radius")])
def test_none_is_preserved_on_optional_numeric_fields(section, key):
    config = Config.from_mapping({section: {key: None}})
    assert getattr(getattr(config, section), key) is None


def test_coercion_does_not_disturb_list_or_tuple_fields():
    config = Config.from_mapping(
        {
            "annotation": {"species": ["Homo sapiens"]},
            "clustering": {"resolution_schedule": [[100, 1.2]]},
            "plot": {"figsize": [6, 6]},
        }
    )
    assert config.annotation.species == ["Homo sapiens"]
    assert config.clustering.resolution_schedule == [(100, 1.2)]
    assert config.plot.figsize == (6.0, 6.0)


def test_range_validation_still_runs_after_coercion():
    """Coercion must not bypass __post_init__ checks."""
    with pytest.raises(ConfigError, match=r"qc\.target_sum must be > 0"):
        Config.from_mapping({"qc": {"target_sum": "0"}})
