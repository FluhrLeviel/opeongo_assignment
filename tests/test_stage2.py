"""Stage 2 checks. Skipped until `opengo qc` and `opengo explore` have run."""

import anndata as ad
import pandas as pd
import pytest

from opengo import config

QC = config.TABLES_DIR / "stage2_qc_by_specimen.tsv"
PB = config.CHECKPOINTS_DIR / "pseudobulk_raw.h5ad"


def test_qc_covers_every_cell():
    if not QC.exists():
        pytest.skip("run `opengo qc` first")
    q = pd.read_csv(QC, sep="\t")
    assert len(q) == 73
    assert q["cells"].sum() == 141_950 + 157_881 + 29_499


def test_pseudobulk_has_whole_specimen_rows():
    if not PB.exists():
        pytest.skip("run `opengo explore` first")
    pb = ad.read_h5ad(PB, backed="r")
    obs = pb.obs.copy()
    pb.file.close()
    whole = obs[obs["population"] == "ALL"]
    assert whole["sample_id"].nunique() == 73
    assert set(config.ECM_GENES) <= set(pb.var["symbol"])


def test_marker_labelling_validated_before_use():
    acc = config.TABLES_DIR / "stage2_marker_label_accuracy.tsv"
    if not acc.exists():
        pytest.skip("run `opengo explore` first")
    a = pd.read_csv(acc, sep="\t").set_index("cohort")["accuracy"]
    assert a["B"] > 0.9 and a["A"] > 0.75
