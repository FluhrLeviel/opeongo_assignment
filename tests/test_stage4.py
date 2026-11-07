"""Stage 4 checks. Skipped until `opengo overview` has run."""

import pandas as pd
import pytest

from opengo import config

MED = config.TABLES_DIR / "stage4_ecm_median_logcpm.tsv"


@pytest.fixture(scope="module")
def med():
    if not MED.exists():
        pytest.skip("run `opengo overview` first")
    return pd.read_csv(MED, sep="\t").set_index(["cohort", "population"])


def test_stroma_is_top_ecm_population_in_both_cohorts(med):
    panel = med[config.ECM_GENES].mean(axis=1)
    for cohort in config.ANALYSIS_COHORTS:
        p = panel.xs(cohort).drop(["ALL", "Unresolved"], errors="ignore")
        assert p.idxmax() in {"Stromal_fibro_myo", "Stromal (fibro+mural)"}


def test_focus_compartment_coverage():
    cov = pd.read_csv(config.TABLES_DIR / "stage4_focus_coverage.tsv", sep="\t")
    s = cov[cov["population"] == "Stromal (fibro+mural)"].set_index("cohort")
    assert s.loc["A", "paired_donors_ge10"] >= 3 and s.loc["B", "paired_donors_ge10"] >= 9


def test_summary_choices_agree():
    for c in config.ANALYSIS_COHORTS:
        s = pd.read_csv(
            config.TABLES_DIR / f"stage4_summary_agreement_{c}.tsv", sep="\t", index_col=0
        )
        assert s.loc["mean_z", "pc1"] > 0.9
