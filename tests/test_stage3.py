"""Stage 3 checks. Skipped until `opengo integrate` has run."""

import pandas as pd
import pytest

from opengo import config

M = config.TABLES_DIR / "stage3_integration_metrics.tsv"


@pytest.fixture(scope="module")
def metrics():
    if not M.exists():
        pytest.skip("run `opengo integrate` first")
    return pd.read_csv(M, sep="\t").set_index("embedding")


def test_all_embeddings_evaluated(metrics):
    assert set(metrics.index) == {"A_baseline", "AB_baseline", "AB_harmony", "AB_harmony_spec"}


def test_harmony_actually_corrects(metrics):
    # Guards against the harmonypy 2.1.0 regression (1 iteration, no correction).
    assert metrics.loc["AB_harmony", "batch_ASW"] > metrics.loc["AB_baseline", "batch_ASW"] + 0.1


def test_harmony_does_not_erase_condition_structure():
    c = pd.read_csv(config.TABLES_DIR / "stage3_condition_structure.tsv", sep="\t")
    p = c.pivot_table(index=["cohort", "population"], columns="embedding", values="condition_ASW")
    stromal = p.xs("Stromal (fibro+mural)", level="population")
    assert (stromal["AB_harmony"] > 0.5 * stromal["AB_baseline"]).all()
