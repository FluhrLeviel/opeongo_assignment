"""Stage 6: linked modalities (cohort A snRNA <-> CPTAC bulk RNA and protein).

1. Crosswalk validation: every A snRNA specimen's bulk link is resolved through bulk_samples.tsv
   to a real column of the condition-specific bulk file; protein aliquot / parent specimen are
   checked against the author sample table. Matched subset = primary (same parent specimen) links.
2. Within-assay tumor - adjacent comparison, paired by donor, for each ECM gene and the panel:
     snRNA stroma (linked specimens), snRNA whole specimen (bridge to bulk), bulk RNA, protein;
   on the matched donors; C3N-00242 (flagged tumor protein link) as sensitivity.
   Context (NOT the matched subset): all CPTAC donors with both tumor and normal columns.
3. Composition: bulk RNA stromal marker genes (non-panel) as a proxy for stromal content; does the
   bulk ECM difference track it? snRNA stromal fraction of the matched specimens.

Bulk values are used as supplied (log2, already normalised); no re-normalisation, no imputation:
missing protein values stay missing and the panel is averaged over available genes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from opengo import config
from opengo.stage4 import COMBINED_STROMAL, STAGE4_DIR
from opengo.stage5 import describe

GENES = config.ECM_GENES
BULK_FILES = {
    ("RNA", "Tumor"): "CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    ("RNA", "Adjacent"): "CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt",
    (
        "protein",
        "Tumor",
    ): "CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    (
        "protein",
        "Adjacent",
    ): "CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt",
}
COMPOSITION_MARKERS = {
    "mural": ["RGS5", "NOTCH3", "PDGFRB", "MYH11"],
    "fibroblast": ["PDGFRA", "C7", "FBLN1"],
    "endothelial": ["VWF", "CDH5", "PECAM1"],
}


# ----------------------------------------------------------------------------- bulk io


def ecm_ids() -> dict[str, str]:
    p = pd.read_csv(config.ECM_PROGRAM_PATH, sep="\t")
    return dict(zip(p["ensembl_gene_id"], p["gene_symbol"], strict=True))


def marker_ids() -> dict[str, str]:
    g = pd.read_csv(config.TABLES_DIR / "gene_correspondence.tsv", sep="\t", low_memory=False)
    want = {s for v in COMPOSITION_MARKERS.values() for s in v}
    g = g[g["symbol_A"].isin(want)].drop_duplicates("symbol_A")
    return dict(zip(g["gene_id"], g["symbol_A"], strict=True))


def load_bulk(modality: str, condition: str, genes: dict[str, str]) -> pd.DataFrame:
    """Rows = gene symbols (requested genes only), columns = donor IDs. Version suffix stripped."""
    d = pd.read_csv(config.BULK_DIR / BULK_FILES[(modality, condition)], sep="\t", index_col=0)
    d.index = d.index.str.split(".").str[0]
    d = d[d.index.isin(genes)]
    d = d[~d.index.duplicated()]
    d.index = d.index.map(genes)
    return d


def bulk_columns(modality: str, condition: str) -> list[str]:
    path = config.BULK_DIR / BULK_FILES[(modality, condition)]
    with path.open() as fh:
        return fh.readline().rstrip("\n").split("\t")[1:]


# ----------------------------------------------------------------------------- crosswalk


def validate_crosswalk() -> pd.DataFrame:
    cw = pd.read_csv(config.MULTIMODAL_CROSSWALK_PATH, sep="\t", dtype=str)
    bs = pd.read_csv(config.BULK_SAMPLES_PATH, sep="\t", dtype=str)
    anchor = pd.read_csv(config.ANCHOR_SAMPLES_PATH, sep="\t", dtype=str)
    cols = {(m, c): set(bulk_columns(m, c)) for m, c in BULK_FILES}
    v = cw.merge(
        bs[["bulk_sample_id", "source_column", "source_condition"]], on="bulk_sample_id", how="left"
    )
    v["linked"] = v["bulk_sample_id"].notna() & (v["bulk_sample_id"] != "")
    v["primary"] = v["primary_specimen_link"].str.lower().eq("true")
    for mod in ("RNA", "protein"):
        v[f"in_{mod}_file"] = [
            bool(r.linked and r.source_column in cols[(mod, r.condition)]) for r in v.itertuples()
        ]
    a = anchor.set_index("Aliquot.snRNA")
    v["anchor_bulk_aliquot"] = v["sample_id"].map(a["Aliquot.bulk"])
    v["anchor_parent"] = v["sample_id"].map(a["Sample"])
    v["protein_aliquot_matches_anchor"] = v["protein_aliquot"] == v["anchor_bulk_aliquot"]
    v["snRNA_parent_matches_protein_parent"] = (
        v["snRNA_parent_specimen"] == v["protein_parent_specimen"]
    )
    v["condition_matches_bulk_file"] = (
        v["source_condition"].replace({"Normal": "Adjacent"}) == v["condition"]
    ) | ~v["linked"]
    return v


def matched_donors(v: pd.DataFrame, include_flagged: bool = False) -> list[str]:
    ok = v["linked"] & v["in_RNA_file"] & v["in_protein_file"]
    ok &= v["primary"] | (
        include_flagged & (v["matching_status"] == "multiple_protein_specimen_ids")
    )
    s = v[ok]
    both = s.groupby("donor_id")["condition"].nunique()
    return sorted(both[both == 2].index)


# ----------------------------------------------------------------------------- comparisons


def _paired(tumor: pd.DataFrame, adj: pd.DataFrame, donors: list[str]) -> pd.DataFrame:
    """donor x gene differences (tumor - adjacent), only donors present in both."""
    d = [x for x in donors if x in tumor.columns and x in adj.columns]
    return (tumor[d] - adj[d]).T


def _panel(diff: pd.DataFrame, sd: pd.Series) -> pd.DataFrame:
    z = diff[GENES] / sd[GENES]
    out = diff.copy()
    out["panel_z"] = z.mean(axis=1, skipna=True)
    out["n_genes_available"] = z.notna().sum(axis=1)
    return out


def snrna_specimen_levels(v: pd.DataFrame, population: str) -> pd.DataFrame:
    """log2-CPM of linked snRNA specimens: donor x condition x gene (A only)."""
    long = pd.read_csv(STAGE4_DIR / "ecm_pseudobulk_long.tsv.gz", sep="\t")
    linked = v[
        v["linked"] & (v["primary"] | (v["matching_status"] == "multiple_protein_specimen_ids"))
    ]
    d = long[(long["population"] == population) & long["sample_id"].isin(linked["sample_id"])]
    d = d[d["n_cells"] >= config.PSEUDOBULK_MIN_CELLS]
    d = d.merge(linked[["sample_id", "donor_id"]], on="sample_id")
    sd_all = (
        long[
            (long["cohort"] == "A")
            & (long["population"] == population)
            & (long["n_cells"] >= config.PSEUDOBULK_MIN_CELLS)
        ]
        .groupby("gene")["logcpm"]
        .std()
    )
    wide = d.pivot_table(index="gene", columns=["condition", "donor_id"], values="logcpm")
    return wide, sd_all


def matched_comparison(v: pd.DataFrame, donors: list[str]) -> dict[str, pd.DataFrame]:
    out = {}
    for label, pop in [("snRNA stroma", COMBINED_STROMAL), ("snRNA whole specimen", "ALL")]:
        wide, sd = snrna_specimen_levels(v, pop)
        t = wide["Tumor"] if "Tumor" in wide.columns.get_level_values(0) else pd.DataFrame()
        a = wide["Adjacent"] if "Adjacent" in wide.columns.get_level_values(0) else pd.DataFrame()
        out[label] = _panel(_paired(t.reindex(GENES), a.reindex(GENES), donors), sd)
    ids = ecm_ids()
    for mod in ("RNA", "protein"):
        t, a = load_bulk(mod, "Tumor", ids), load_bulk(mod, "Adjacent", ids)
        sd = pd.concat([t, a], axis=1).std(axis=1)
        out[f"bulk {mod}"] = _panel(_paired(t.reindex(GENES), a.reindex(GENES), donors), sd)
    return out


def full_cohort() -> dict[str, pd.DataFrame]:
    """Context: every CPTAC donor with both tumor and normal columns (not the matched subset)."""
    ids = ecm_ids()
    out = {}
    for mod in ("RNA", "protein"):
        t, a = load_bulk(mod, "Tumor", ids), load_bulk(mod, "Adjacent", ids)
        donors = sorted(set(t.columns) & set(a.columns))
        sd = pd.concat([t, a], axis=1).std(axis=1)
        out[f"bulk {mod} (all paired donors)"] = _panel(
            _paired(t.reindex(GENES), a.reindex(GENES), donors), sd
        )
    return out


def summarise(diffs: dict[str, pd.DataFrame], scope: str) -> pd.DataFrame:
    rows = []
    for assay, d in diffs.items():
        for f in GENES + ["panel_z"]:
            x = d[f].dropna() if f in d else pd.Series(dtype=float)
            r = {"scope": scope, "assay": assay, "feature": f, **describe(x)}
            if len(x) > 2:
                r["p_wilcoxon"] = float(stats.wilcoxon(x).pvalue) if (x != 0).any() else np.nan
            rows.append(r)
    t = pd.DataFrame(rows)
    t["q_BH"] = np.nan
    for assay in t["assay"].unique():
        m = (t["assay"] == assay) & t["feature"].isin(GENES) & t["p_signflip"].notna()
        if m.sum() > 1:
            t.loc[m, "q_BH"] = multipletests(t.loc[m, "p_signflip"], method="fdr_bh")[1]
    return t


# ----------------------------------------------------------------------------- composition


def composition_proxy() -> dict[str, pd.DataFrame]:
    """Bulk RNA (all paired donors): does the ECM difference track stromal marker differences?"""
    ids = ecm_ids()
    mids = marker_ids()
    t = pd.concat([load_bulk("RNA", "Tumor", ids), load_bulk("RNA", "Tumor", mids)])
    a = pd.concat([load_bulk("RNA", "Adjacent", ids), load_bulk("RNA", "Adjacent", mids)])
    donors = sorted(set(t.columns) & set(a.columns))
    sd = pd.concat([t, a], axis=1).std(axis=1)
    diff = (t[donors] - a[donors]).T
    z = diff / sd
    scores = pd.DataFrame(
        {
            "ECM_panel_z": z[GENES].mean(axis=1),
            "ECM_noDCN_LUM_z": z[[g for g in GENES if g not in ("DCN", "LUM")]].mean(axis=1),
            **{
                f"{k}_markers_z": z[[g for g in v if g in z]].mean(axis=1)
                for k, v in COMPOSITION_MARKERS.items()
            },
        }
    )
    corr = []
    for k in COMPOSITION_MARKERS:
        rho, p = stats.spearmanr(scores["ECM_panel_z"], scores[f"{k}_markers_z"])
        rho8, p8 = stats.spearmanr(scores["ECM_noDCN_LUM_z"], scores[f"{k}_markers_z"])
        corr.append(
            {
                "marker_set": k,
                "spearman_rho": rho,
                "p": p,
                "spearman_rho_noDCN_LUM": rho8,
                "p_noDCN_LUM": p8,
                "n_donors": len(scores),
                "mean_marker_change_z": scores[f"{k}_markers_z"].mean(),
                "frac_donors_marker_down": float((scores[f"{k}_markers_z"] < 0).mean()),
            }
        )
    # ECM change left after removing the part explained by stromal (mural + fibroblast) markers
    x = scores[["mural_markers_z", "fibroblast_markers_z"]].to_numpy()
    xx = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(xx, scores["ECM_panel_z"].to_numpy(), rcond=None)
    resid_intercept = beta[0]
    adj = pd.DataFrame(
        [
            {
                "ECM_change_unadjusted": scores["ECM_panel_z"].mean(),
                "ECM_noDCN_LUM_change": scores["ECM_noDCN_LUM_z"].mean(),
                "ECM_change_at_zero_stromal_marker_change": resid_intercept,
                "note": "intercept extrapolates to zero marker change (rare in data)",
                "beta_mural": beta[1],
                "beta_fibroblast": beta[2],
            }
        ]
    )
    return {"scores": scores.reset_index(names="donor_id"), "corr": pd.DataFrame(corr), "adj": adj}


def snrna_stromal_fraction(v: pd.DataFrame) -> pd.DataFrame:
    long = pd.read_csv(STAGE4_DIR / "ecm_pseudobulk_long.tsv.gz", sep="\t")
    b = long[(long["gene"] == GENES[0]) & (long["cohort"] == "A")]
    n = b.pivot_table(index="sample_id", columns="population", values="n_cells", fill_value=0)
    frac = (n[COMBINED_STROMAL] / n["ALL"]).rename("stromal_fraction").reset_index()
    linked = v[v["linked"]][["sample_id", "donor_id", "condition", "matching_status"]]
    return linked.merge(frac, on="sample_id", how="left")


# ----------------------------------------------------------------------------- agreement


def agreement_table(effects: pd.DataFrame) -> pd.DataFrame:
    """Tumor - adjacent mean per gene/panel across all evidence lines (direction comparison)."""
    s5 = pd.read_csv(config.TABLES_DIR / "stage5_effects.tsv", sep="\t")
    s5 = s5[s5["analysis"].str.startswith("primary") & s5["cohort"].isin(["A", "B"])]
    rows = []
    for _, r in s5[s5["feature"].isin(GENES + ["panel_z"])].iterrows():
        rows.append(
            {
                "evidence": f"snRNA/scRNA stroma, cohort {r['cohort']} (Stage 5)",
                "feature": r["feature"],
                "mean": r["mean"],
                "n_donors": r["n_donors"],
                "n_up": r["n_up"],
                "q_BH": r.get("q_BH"),
            }
        )
    for _, r in effects.iterrows():
        rows.append(
            {
                "evidence": f"{r['assay']} | {r['scope']}",
                "feature": r["feature"],
                "mean": r["mean"],
                "n_donors": r["n_donors"],
                "n_up": r["n_up"],
                "q_BH": r.get("q_BH"),
            }
        )
    t = pd.DataFrame(rows)
    t["direction"] = np.sign(t["mean"]).map({1.0: "up", -1.0: "down", 0.0: "none"})
    return t


# ----------------------------------------------------------------------------- record


def specimen_inclusion(v: pd.DataFrame, primary: list[str], flagged: list[str]) -> pd.DataFrame:
    """Every cohort A snRNA specimen: its bulk link and role in the matched comparison."""
    r = v[["sample_id", "donor_id", "condition", "matching_status", "bulk_sample_id"]].copy()
    in_primary = v["primary"] & v["donor_id"].isin(primary)
    in_sens = (v["primary"] | (v["matching_status"] == "multiple_protein_specimen_ids")) & v[
        "donor_id"
    ].isin(flagged)
    r["matched_primary"] = in_primary.to_numpy()
    r["matched_sensitivity"] = in_sens.to_numpy()
    r["reason"] = np.select(
        [
            ~v["linked"],
            v["matching_status"].eq("different_parent_specimen"),
            v["matching_status"].eq("multiple_protein_specimen_ids"),
            ~v["donor_id"].isin(flagged),
            ~v["donor_id"].isin(primary),
        ],
        [
            "no author bulk link (no pair invented)",
            "protein from a different parent specimen: excluded",
            "protein from two parent specimens: sensitivity only",
            "linked, but donor has no linked adjacent specimen",
            "sensitivity only (this donor's tumor link is flagged)",
        ],
        default="matched pair (primary link)",
    )
    return r


# ----------------------------------------------------------------------------- driver


def run() -> dict[str, pd.DataFrame]:
    v = validate_crosswalk()
    primary = matched_donors(v)
    with_flagged = matched_donors(v, include_flagged=True)
    matched = matched_comparison(v, primary)
    matched_s = matched_comparison(v, with_flagged)
    full = full_cohort()
    comp = composition_proxy()

    per_donor = []
    for assay, d in {**matched, **{f"{k} [+C3N-00242]": x for k, x in matched_s.items()}}.items():
        x = d.reset_index(names="donor_id")
        x.insert(0, "assay", assay)
        per_donor.append(x)

    tables = {
        "stage6_crosswalk_validation": v,
        "stage6_matched_donors": pd.DataFrame(
            {"donor_id": with_flagged, "primary_pair": [d in primary for d in with_flagged]}
        ),
        "stage6_effects": pd.concat(
            [
                summarise(matched, f"matched subset ({len(primary)} donors, primary links)"),
                summarise(matched_s, f"matched + C3N-00242 ({len(with_flagged)} donors)"),
                summarise(full, "full bulk cohort (context, not matched)"),
            ],
            ignore_index=True,
        ),
        "stage6_matched_per_donor": pd.concat(per_donor, ignore_index=True),
        "stage6_agreement": None,
        "stage6_composition_scores": comp["scores"],
        "stage6_composition_corr": comp["corr"],
        "stage6_composition_adjusted": comp["adj"],
        "stage6_snrna_stromal_fraction": snrna_stromal_fraction(v),
        "stage6_protein_missingness": pd.DataFrame(
            [
                {
                    "condition": c,
                    "gene": g,
                    "missing_frac": float(
                        load_bulk("protein", c, ecm_ids()).reindex(GENES).loc[g].isna().mean()
                    ),
                }
                for c in ("Tumor", "Adjacent")
                for g in GENES
            ]
        ),
    }
    tables["stage6_agreement"] = agreement_table(tables["stage6_effects"])
    tables["stage6_specimen_inclusion"] = specimen_inclusion(v, primary, with_flagged)
    for name, df in tables.items():
        df.to_csv(config.TABLES_DIR / f"{name}.tsv", sep="\t", index=False)
    return tables
