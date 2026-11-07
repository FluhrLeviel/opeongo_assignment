"""Stage 1: explicit gene correspondence across cohorts.

A and C carry Ensembl gene IDs (different annotation releases); B carries symbols only (an older
annotation). Correspondence rules, in order:

1. A and C are joined on Ensembl ID.
2. A B symbol is assigned an Ensembl ID if it equals exactly one A symbol, else exactly one C
   symbol. Symbols matching several IDs are `ambiguous`; no match is `unmatched`.
3. If two B symbols land on the same ID, both are flagged `collision` and left unassigned.

Genes absent from a cohort are *not measured* there; they are never filled with zeros.
"""

from __future__ import annotations

import anndata as ad
import pandas as pd

from opengo import config


def cohort_var(cohort: str) -> pd.DataFrame:
    a = ad.read_h5ad(config.cohort_raw_path(cohort), backed="r")
    var = a.var.copy()
    a.file.close()
    for col in var.columns:  # h5ad stores strings as categoricals
        if isinstance(var[col].dtype, pd.CategoricalDtype):
            var[col] = var[col].astype(str)
    if "gene_id" not in var:
        var["gene_id"] = None
    return var


def _unique_symbol_map(var: pd.DataFrame) -> dict[str, str]:
    counts = var["symbol"].value_counts()
    uniq = var[var["symbol"].map(counts) == 1]
    return dict(zip(uniq["symbol"], uniq["gene_id"], strict=True))


def build_gene_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    va, vb, vc = (cohort_var(c) for c in config.COHORTS)
    a = va[["gene_id", "symbol"]].rename(columns={"symbol": "symbol_A"})
    c = vc[["gene_id", "symbol"]].rename(columns={"symbol": "symbol_C"})
    table = a.merge(c, on="gene_id", how="outer")
    table["in_A"] = table["symbol_A"].notna()
    table["in_C"] = table["symbol_C"].notna()

    map_a, map_c = _unique_symbol_map(va), _unique_symbol_map(vc)
    dup_a = set(va["symbol"][va["symbol"].duplicated(keep=False)])
    dup_c = set(vc["symbol"][vc["symbol"].duplicated(keep=False)])
    b = pd.DataFrame({"symbol_B": vb["symbol"].astype(str).values})
    ids, how = [], []
    for s in b["symbol_B"]:
        if s in map_a:
            ids.append(map_a[s])
            how.append("A_symbol")
        elif s in map_c:
            ids.append(map_c[s])
            how.append("C_symbol")
        elif s in dup_a or s in dup_c:
            ids.append(None)
            how.append("ambiguous")
        else:
            ids.append(None)
            how.append("unmatched")
    b["gene_id"], b["B_match"] = ids, how
    collided = b["gene_id"].notna() & b["gene_id"].duplicated(keep=False)
    b.loc[collided, "B_match"] = "collision"
    b.loc[collided, "gene_id"] = None

    mapped = b[b["gene_id"].notna()]
    table = table.merge(mapped[["gene_id", "symbol_B", "B_match"]], on="gene_id", how="outer")
    unmapped_b = b[b["gene_id"].isna()][["symbol_B", "B_match"]]
    table = pd.concat([table, unmapped_b], ignore_index=True)
    table["in_B"] = table["symbol_B"].notna()
    for col in ["in_A", "in_C"]:
        table[col] = table[col].fillna(False).astype(bool)
    table["symbol"] = table["symbol_A"].fillna(table["symbol_C"]).fillna(table["symbol_B"])
    table["n_cohorts"] = table[["in_A", "in_B", "in_C"]].sum(axis=1)
    table["shared_ABC"] = table["n_cohorts"].eq(3) & table["gene_id"].notna()

    summary = pd.DataFrame(
        {
            "metric": [
                "genes_A",
                "genes_B",
                "genes_C",
                "A_and_C_by_id",
                "B_mapped_via_A_symbol",
                "B_mapped_via_C_symbol",
                "B_ambiguous",
                "B_collision",
                "B_unmatched",
                "shared_A_B_C",
                "symbol_differs_A_vs_C_same_id",
            ],
            "value": [
                len(va),
                len(vb),
                len(vc),
                int((table["in_A"] & table["in_C"]).sum()),
                int((b["B_match"] == "A_symbol").sum()),
                int((b["B_match"] == "C_symbol").sum()),
                int((b["B_match"] == "ambiguous").sum()),
                int((b["B_match"] == "collision").sum()),
                int((b["B_match"] == "unmatched").sum()),
                int(table["shared_ABC"].sum()),
                int(
                    (table["in_A"] & table["in_C"] & (table["symbol_A"] != table["symbol_C"])).sum()
                ),
            ],
        }
    )
    cols = ["gene_id", "symbol", "symbol_A", "symbol_B", "symbol_C", "in_A", "in_B", "in_C"]
    cols += ["B_match", "n_cohorts", "shared_ABC"]
    return table[cols], summary


def ecm_check(table: pd.DataFrame) -> pd.DataFrame:
    """Each panel gene, located by the supplied Ensembl ID, in every cohort."""
    panel = pd.read_csv(config.ECM_PROGRAM_PATH, sep="\t")
    out = panel[["gene_symbol", "ensembl_gene_id"]].merge(
        table, left_on="ensembl_gene_id", right_on="gene_id", how="left"
    )
    return out[
        ["gene_symbol", "ensembl_gene_id", "symbol_A", "symbol_B", "symbol_C", "B_match"]
        + ["in_A", "in_B", "in_C", "shared_ABC"]
    ]
