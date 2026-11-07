"""Sanity tests for configuration and metadata wiring (fast; no large data loaded)."""

import pandas as pd

from opengo import config


def test_project_root_points_at_package():
    assert (config.PROJECT_ROOT / "FILE_MANIFEST.tsv").exists()
    assert (config.PROJECT_ROOT / "pyproject.toml").exists()


def test_ecm_panel_matches_supplied_list():
    supplied = pd.read_csv(config.ECM_PROGRAM_PATH, sep="\t")["gene_symbol"].tolist()
    assert sorted(supplied) == sorted(config.ECM_GENES)
    assert len(config.ECM_GENES) == 10


def test_metadata_paths_exist():
    for p in [
        config.SAMPLE_MANIFEST_PATH,
        config.MULTIMODAL_CROSSWALK_PATH,
        config.LABEL_MAPPING_PATH,
        config.BULK_SAMPLES_PATH,
        config.BULK_COLUMN_MANIFEST_PATH,
        config.ANCHOR_ANNOTATIONS_PATH,
        config.COHORT_B_ANNOTATIONS_PATH,
        config.COHORT_C_SAMPLES_PATH,
    ]:
        assert p.exists(), p


def test_sample_manifest_has_73_specimens():
    m = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t")
    assert m.groupby("cohort").size().to_dict() == {"A": 34, "B": 26, "C": 13}
    assert m["sample_id"].is_unique
