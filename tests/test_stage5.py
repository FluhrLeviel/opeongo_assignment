"""Stage 5 checks. Skipped until `opengo pseudobulk` has run."""

import pandas as pd
import pytest

from opengo import config
from opengo.stage5 import sign_flip_p

E = config.TABLES_DIR / "stage5_effects.tsv"


def test_sign_flip_exact_minimum():
    # n positive differences -> smallest attainable two-sided p is 2 / 2^n
    assert sign_flip_p([1, 2, 3]) == pytest.approx(2 / 8)
    assert sign_flip_p([1.0] * 9) == pytest.approx(2 / 512)


@pytest.fixture(scope="module")
def eff():
    if not E.exists():
        pytest.skip("run `opengo pseudobulk` first")
    return pd.read_csv(E, sep="\t")


def test_donors_are_the_unit(eff):
    p = eff[eff["analysis"].str.startswith("primary") & (eff["feature"] == "panel_z")]
    n = p.set_index("cohort")["n_donors"]
    assert n["A"] == 3 and n["B"] == 9  # paired donors, repeated regions averaged


def test_every_gene_reported_with_bh(eff):
    p = eff[eff["analysis"].str.startswith("primary") & eff["cohort"].isin(["A", "B"])]
    genes = p[p["feature"].isin(config.ECM_GENES)]
    assert len(genes) == 20 and genes["q_BH"].notna().all()
