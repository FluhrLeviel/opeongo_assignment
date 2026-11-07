"""Stage 1 reconciliation checks. Skipped until `opengo load all` and `opengo genes` have run."""

import pandas as pd
import pytest

from opengo import config

SPEC = config.TABLES_DIR / "stage1_specimens.tsv"
GENES = config.TABLES_DIR / "gene_correspondence.tsv"
EXPECTED_CELLS = {"A": 141_950, "B": 157_881, "C": 29_499}  # DATA_DICTIONARY.md


@pytest.fixture(scope="module")
def spec():
    if not SPEC.exists():
        pytest.skip("run `opengo load all` first")
    return pd.read_csv(SPEC, sep="\t")


def test_all_73_specimens_loaded(spec):
    assert spec.groupby("cohort").size().to_dict() == {"A": 34, "B": 26, "C": 13}


def test_cells_per_cohort_match_documentation(spec):
    assert spec.groupby("cohort")["cells_loaded"].sum().to_dict() == EXPECTED_CELLS


def test_cells_per_specimen_match_manifest(spec):
    bad = spec[spec["cells_loaded"] != spec["manifest_cells"]]
    assert bad.empty, bad[["sample_id", "cells_loaded", "manifest_cells"]]


def test_every_author_annotation_found(spec):
    ab = spec[spec["cohort"].isin(["A", "B"])]
    assert (ab["annotated_matched"] == ab["annotated_expected"]).all()
    assert (ab["loaded_not_annotated"] == 0).all()


def test_counts_are_nonnegative_integers(spec):
    assert spec["integer_counts"].all()
    assert (spec["min_count"] >= 0).all()


def test_ecm_panel_in_all_cohorts():
    if not GENES.exists():
        pytest.skip("run `opengo genes` first")
    ecm = pd.read_csv(config.TABLES_DIR / "ecm_panel_presence.tsv", sep="\t")
    assert len(ecm) == 10
    assert ecm["shared_ABC"].all()


def test_trm_override_applied_and_traceable():
    path = config.TABLES_DIR / "stage1_population_overrides.tsv"
    if not path.exists():
        pytest.skip("run `opengo genes` first")
    ov = pd.read_csv(path, sep="\t")
    assert len(ov) == 1
    row = ov.iloc[0]
    assert (row["cohort"], row["population_suggested"], row["population"]) == (
        "A",
        "Lymphoid",
        "Myeloid",
    )
    assert row["cells"] == 510
