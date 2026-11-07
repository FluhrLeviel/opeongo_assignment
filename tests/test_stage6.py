"""Stage 6 checks. Skipped until `opengo multimodal` has run."""

import pandas as pd
import pytest

from opengo import config

V = config.TABLES_DIR / "stage6_crosswalk_validation.tsv"


@pytest.fixture(scope="module")
def v():
    if not V.exists():
        pytest.skip("run `opengo multimodal` first")
    return pd.read_csv(V, sep="\t")


def test_crosswalk_counts_match_documentation(v):
    assert len(v) == 34
    assert v["linked"].sum() == 29
    p = v[v["primary"]]
    assert len(p) == 27 and (p["condition"] == "Adjacent").sum() == 4


def test_every_link_resolves_to_a_real_bulk_column(v):
    linked = v[v["linked"]]
    assert linked["in_RNA_file"].all() and linked["in_protein_file"].all()
    assert linked["protein_aliquot_matches_anchor"].all()


def test_no_pair_invented_from_donor_identity():
    md = pd.read_csv(config.TABLES_DIR / "stage6_matched_donors.tsv", sep="\t")
    assert set(md.loc[md["primary_pair"], "donor_id"]) == {"C3L-00079", "C3L-00088", "C3N-01200"}
    assert set(md.loc[~md["primary_pair"], "donor_id"]) == {"C3N-00242"}
