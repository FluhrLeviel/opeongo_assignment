"""Stage 4: ECM panel overview across populations and choice of focus population.

Values:
- Cell level (overview only): fraction of cells with >= 1 UMI and mean log1p(CP10k) per gene,
  library size from all genes of the cohort. Computed per specimen and aggregated.
- Donor level: pseudobulk raw counts per specimen x population (QC-passing cells), log2-CPM with
  library size from all genes of the cohort (`lib_size_all_genes`), pseudocount 1 CPM.

Panel summary (justified in `panel_coherence`): for each pseudobulk sample, the mean over the 10
genes of the per-gene z-score of log2-CPM, z computed within cohort and population. Equal gene
weights; interpretable as an average standardized shift; robust to one highly expressed gene
dominating (unlike log of summed CPM). Individual genes are always shown next to it.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from opengo import config
from opengo.explore import PSEUDOBULK_PATH, qc_cells

STAGE4_DIR = config.CHECKPOINTS_DIR / "stage4"
# Stromal split by author label (A: Fibroblasts / Myofibroblasts; B: Fibroblast / Pericytes)
FINE_STROMAL = {"Stromal_fibro_myo", "Pericyte"}
COMBINED_STROMAL = "Stromal (fibro+mural)"


# ----------------------------------------------------------------------------- cell level


def cell_level_stats(progress=None) -> pd.DataFrame:
    """Per specimen x group x gene: n cells, n expressing, sum of log1p(CP10k)."""
    rows = []
    for cohort in config.ANALYSIS_COHORTS:
        cells = qc_cells(cohort)
        backed = ad.read_h5ad(config.cohort_raw_path(cohort), backed="r")
        sym = backed.var["symbol"].astype(str).to_numpy()
        cols = [int(np.where(sym == g)[0][0]) for g in config.ECM_GENES]
        obs = backed.obs
        for sid, idx in obs.groupby("sample_id", observed=True).indices.items():
            idx = np.sort(idx)
            sub = backed[idx].to_memory()
            ok = cells.loc[sub.obs_names, "qc_pass"].to_numpy(bool)
            x = sp.csr_matrix(sub.X)[ok]
            lib = np.asarray(x.sum(axis=1)).ravel()
            ecm = x[:, cols].toarray()
            lognorm = np.log1p(ecm / lib[:, None] * 1e4)
            o = sub.obs[ok]
            pop = o["population"].astype(str).to_numpy()
            fine = np.where(
                np.isin(pop, list(FINE_STROMAL)), o["author_broad"].astype(str).to_numpy(), pop
            )
            for grp in np.unique(fine):
                m = fine == grp
                for j, g in enumerate(config.ECM_GENES):
                    rows.append(
                        {
                            "cohort": cohort,
                            "sample_id": sid,
                            "condition": str(o["condition"].iloc[0]),
                            "population": pop[m][0],
                            "group": grp,
                            "gene": g,
                            "n": int(m.sum()),
                            "n_pos": int((ecm[m, j] > 0).sum()),
                            "sum_lognorm": float(lognorm[m, j].sum()),
                        }
                    )
            if progress:
                progress(cohort, sid)
        backed.file.close()
    df = pd.DataFrame(rows)
    STAGE4_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(STAGE4_DIR / "ecm_cell_stats.tsv.gz", sep="\t", index=False, compression="gzip")
    return df


def dotplot_table(stats: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    g = stats.groupby(by + ["gene"], observed=True)[["n", "n_pos", "sum_lognorm"]].sum()
    g["frac"] = g["n_pos"] / g["n"]
    g["mean"] = g["sum_lognorm"] / g["n"]
    return g.reset_index()


# ----------------------------------------------------------------------------- pseudobulk


def ecm_pseudobulk(min_cells: int = config.PSEUDOBULK_MIN_CELLS) -> pd.DataFrame:
    """Long table: specimen x population x gene log2-CPM (+ combined stromal compartment)."""
    pb = ad.read_h5ad(PSEUDOBULK_PATH)
    pb = pb[pb.obs["cohort"].isin(config.ANALYSIS_COHORTS).to_numpy()].copy()
    sym = pb.var["symbol"].astype(str)
    cols = [int(np.where(sym.to_numpy() == g)[0][0]) for g in config.ECM_GENES]
    x = pb.X.toarray() if sp.issparse(pb.X) else np.asarray(pb.X)
    counts = pd.DataFrame(x[:, cols], index=pb.obs_names, columns=config.ECM_GENES)
    obs = pb.obs.copy()
    for col in obs.columns:  # h5ad categoricals -> plain values
        if isinstance(obs[col].dtype, pd.CategoricalDtype):
            obs[col] = obs[col].astype(str)
    obs.index = obs.index.astype(str)
    counts.index = obs.index

    # combined stromal compartment = sum of fibroblast/myofibroblast + pericyte pseudobulks
    st = obs[obs["population"].isin(FINE_STROMAL)]
    comb_c = counts.loc[st.index].groupby(st["sample_id"]).sum()
    comb_o = (
        st.groupby("sample_id")
        .agg(
            cohort=("cohort", "first"),
            donor_key=("donor_key", "first"),
            condition=("condition", "first"),
            assay=("assay", "first"),
            n_cells=("n_cells", "sum"),
            lib_size_all_genes=("lib_size_all_genes", "sum"),
        )
        .assign(population=COMBINED_STROMAL)
    )
    comb_o["sample_id"] = comb_o.index
    comb_o.index = comb_o.index + "|" + COMBINED_STROMAL
    comb_c.index = comb_o.index
    counts = pd.concat([counts, comb_c])
    obs = pd.concat([obs, comb_o[obs.columns]])

    logcpm = np.log2(counts.div(obs["lib_size_all_genes"], axis=0) * 1e6 + 1)
    long = logcpm.join(obs).melt(
        id_vars=list(obs.columns), value_vars=config.ECM_GENES, var_name="gene", value_name="logcpm"
    )
    long["enough_cells"] = long["n_cells"] >= min_cells
    STAGE4_DIR.mkdir(parents=True, exist_ok=True)
    long.to_csv(STAGE4_DIR / "ecm_pseudobulk_long.tsv.gz", sep="\t", index=False)
    return long


def panel_score(long: pd.DataFrame) -> pd.DataFrame:
    """Mean of per-gene z-scores (z within cohort x population) per pseudobulk sample."""
    d = long[long["enough_cells"]].copy()
    grp = d.groupby(["cohort", "population", "gene"])["logcpm"]
    d["z"] = (d["logcpm"] - grp.transform("mean")) / grp.transform("std").replace(0, np.nan)
    keys = ["cohort", "population", "sample_id", "donor_key", "condition", "n_cells"]
    return d.groupby(keys, observed=True)["z"].mean().rename("panel_z").reset_index()


def panel_coherence(long: pd.DataFrame, population: str = COMBINED_STROMAL) -> dict:
    """Do the 10 genes co-vary across specimens? Gene-gene Spearman + PC1 loadings, per cohort."""
    out = {}
    for cohort in config.ANALYSIS_COHORTS:
        d = long[
            (long["cohort"] == cohort) & (long["population"] == population) & long["enough_cells"]
        ]
        wide = d.pivot_table(index="sample_id", columns="gene", values="logcpm")[config.ECM_GENES]
        corr = wide.corr(method="spearman")
        z = (wide - wide.mean()) / wide.std()
        _, s, vt = np.linalg.svd(z.fillna(0).to_numpy(), full_matrices=False)
        load = pd.Series(vt[0] * np.sign(vt[0].sum()), index=config.ECM_GENES)
        summaries = pd.DataFrame(
            {
                "mean_z": z.mean(axis=1),
                "log_sum_cpm": np.log2((2**wide - 1).sum(axis=1) + 1),
                "pc1": z.fillna(0).to_numpy() @ load.to_numpy(),
            }
        )
        out[cohort] = {
            "corr": corr,
            "pc1_loadings": load,
            "pc1_var": float(s[0] ** 2 / (s**2).sum()),
            "n": len(wide),
            "summary_corr": summaries.corr(method="spearman"),
        }
    return out


# ----------------------------------------------------------------------------- focus choice


def coverage_table(long: pd.DataFrame, thresholds=(10, 20, 50)) -> pd.DataFrame:
    base = long[long["gene"] == config.ECM_GENES[0]]
    rows = []
    for (cohort, pop), d in base.groupby(["cohort", "population"]):
        per_donor = d.groupby(["donor_key", "condition"])["n_cells"].sum().unstack(fill_value=0)
        row = {"cohort": cohort, "population": pop}
        for t in thresholds:
            ok = (per_donor.get("Tumor", 0) >= t) & (per_donor.get("Adjacent", 0) >= t)
            row[f"paired_donors_ge{t}"] = int(ok.sum())
        rows.append(row)
    return pd.DataFrame(rows)


def expression_table(long: pd.DataFrame) -> pd.DataFrame:
    """Median log2-CPM per population x gene (specimens with enough cells), per cohort."""
    d = long[long["enough_cells"]]
    return (
        d.groupby(["cohort", "population", "gene"])["logcpm"]
        .median()
        .unstack("gene")[config.ECM_GENES]
    )


def run(progress=None) -> dict[str, pd.DataFrame]:
    stats = cell_level_stats(progress=progress)
    long = ecm_pseudobulk()
    score = panel_score(long)
    coh = panel_coherence(long)
    tables = {
        "stage4_ecm_dotplot_population": dotplot_table(stats, ["cohort", "population"]),
        "stage4_ecm_dotplot_stromal_fine": dotplot_table(
            stats[stats["population"].isin(FINE_STROMAL)], ["cohort", "group"]
        ),
        "stage4_ecm_median_logcpm": expression_table(long).reset_index(),
        "stage4_panel_score": score,
        "stage4_focus_coverage": coverage_table(long),
    }
    coh_rows = []
    for cohort, c in coh.items():
        c["corr"].to_csv(config.TABLES_DIR / f"stage4_ecm_gene_corr_{cohort}.tsv", sep="\t")
        c["summary_corr"].to_csv(
            config.TABLES_DIR / f"stage4_summary_agreement_{cohort}.tsv", sep="\t"
        )
        off = c["corr"].to_numpy()[np.triu_indices(10, 1)]
        for g, v in c["pc1_loadings"].items():
            gene_corr = c["corr"].loc[g].drop(g)
            coh_rows.append(
                {
                    "cohort": cohort,
                    "gene": g,
                    "pc1_loading": v,
                    "mean_corr_with_others": gene_corr.mean(),
                    "pc1_var_explained": c["pc1_var"],
                    "median_pairwise_corr": float(np.median(off)),
                    "n_specimens": c["n"],
                }
            )
    tables["stage4_panel_coherence"] = pd.DataFrame(coh_rows)
    for name, df in tables.items():
        df.to_csv(config.TABLES_DIR / f"{name}.tsv", sep="\t", index=False)
    return tables
