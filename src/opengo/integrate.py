"""Stage 3: unintegrated baseline vs Harmony on cohorts A and B.

Embeddings are built on a stratified subsample (up to `cells_per_specimen` random QC-passing cells
per specimen plus *all* QC-passing endothelial, fibroblast/myofibroblast and pericyte cells, so rare
stromal populations are fully represented). Embeddings are for visualisation, evaluation and label
transfer only; no test is run on them (testing uses pseudobulk counts, Stages 4-5).

Embeddings compared:
- A_baseline      : cohort A alone, PCA (the anchor cohort's own structure)
- AB_baseline     : A + B, PCA, no correction
- AB_harmony      : A + B, Harmony on PCA, batch = cohort (primary)
- AB_harmony_spec : A + B, Harmony on PCA, batch = specimen (sensitivity; specimen is confounded
                    with condition, so this is expected to erode tumor/adjacent structure)

Metrics (scIB-style; Luecken et al. 2022):
- batch mixing   : iLISI over cohort (rescaled 0-1), batch ASW within each population
- identity       : cLISI over population (rescaled 0-1), cell-type ASW, ARI/NMI of Leiden vs labels,
                   kNN label purity per population
- cross-cohort   : kNN label transfer A -> B and B -> A (accuracy, balanced accuracy)
- biology        : tumor/adjacent silhouette within each (cohort, population) — does integration
                   erase the condition structure we later test?
"""

from __future__ import annotations

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import (
    adjusted_rand_score,
    balanced_accuracy_score,
    normalized_mutual_info_score,
    silhouette_samples,
    silhouette_score,
)
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors

from opengo import config
from opengo.explore import _cohort_gene_index, _CsrKnn, qc_cells, shared_genes

P = config.INTEGRATION
STAGE3_DIR = config.CHECKPOINTS_DIR / "stage3"
EVAL_MAP = {"Stromal_fibro_myo": "Stromal (fibro+mural)", "Pericyte": "Stromal (fibro+mural)"}


# ----------------------------------------------------------------------------- data


def load_subsample() -> ad.AnnData:
    rng = np.random.default_rng(config.RANDOM_SEED)
    shared = shared_genes()
    parts = []
    for cohort in config.ANALYSIS_COHORTS:
        cells = qc_cells(cohort)
        backed = ad.read_h5ad(config.cohort_raw_path(cohort), backed="r")
        obs = backed.obs
        ok = cells.reindex(backed.obs_names)["qc_pass"].to_numpy(bool)
        pop = obs["population"].astype(str).to_numpy()
        rare = np.isin(pop, P["keep_all"])
        keep = []
        for _, idx in obs.groupby("sample_id", observed=True).indices.items():
            idx = np.asarray(idx)
            keep.extend(idx[ok[idx] & rare[idx]].tolist())
            common = idx[ok[idx] & ~rare[idx]]
            if len(common) > P["cells_per_specimen"]:
                common = rng.choice(common, P["cells_per_specimen"], replace=False)
            keep.extend(common.tolist())
        keep = np.sort(np.array(keep))
        sub = backed[keep].to_memory()
        gidx = _cohort_gene_index(cohort, backed.var, shared)
        backed.file.close()
        cols = ["cohort", "sample_id", "donor_key", "condition", "assay", "population"]
        cols += ["author_broad", "author_detailed"]
        part = ad.AnnData(
            X=sp.csr_matrix(sub.X)[:, gidx],
            obs=sub.obs[cols].astype(str),
            var=pd.DataFrame(index=shared["gene_id"].astype(str).to_numpy()),
        )
        parts.append(part)
    adata = ad.concat(parts, join="inner")
    adata.var["symbol"] = shared.set_index("gene_id").loc[adata.var_names, "symbol"].to_numpy()
    adata.obs["pop_eval"] = adata.obs["population"].replace(EVAL_MAP)
    adata.obs.loc[adata.obs["population"] == "Unresolved", "pop_eval"] = np.nan
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def _pca(adata: ad.AnnData, batch_key: str | None) -> np.ndarray:
    a = adata.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=P["n_hvg"], batch_key=batch_key, flavor="seurat")
    sc.pp.pca(a, n_comps=P["n_pcs"], mask_var="highly_variable", random_state=config.RANDOM_SEED)
    return a.obsm["X_pca"]


def _graph_umap_leiden(adata: ad.AnnData, rep: str, tag: str) -> None:
    sc.pp.neighbors(
        adata,
        use_rep=rep,
        n_neighbors=P["n_neighbors"],
        transformer=_CsrKnn(P["n_neighbors"]),
        key_added=tag,
        random_state=config.RANDOM_SEED,
    )
    sc.tl.leiden(
        adata,
        neighbors_key=tag,
        key_added=f"leiden_{tag}",
        resolution=P["leiden_resolution"],
        flavor="igraph",
        n_iterations=2,
        random_state=config.RANDOM_SEED,
    )
    sc.tl.umap(adata, neighbors_key=tag, random_state=config.RANDOM_SEED)
    adata.obsm[f"X_umap_{tag}"] = adata.obsm.pop("X_umap")


def build_embeddings(adata: ad.AnnData, progress=print) -> dict[str, ad.AnnData]:
    out = {}
    a_only = adata[adata.obs["cohort"] == "A"].copy()
    a_only.obsm["X_A_baseline"] = _pca(a_only, batch_key=None)
    progress("A_baseline: PCA done")
    _graph_umap_leiden(a_only, "X_A_baseline", "A_baseline")
    out["A"] = a_only

    pcs = _pca(adata, batch_key="cohort")
    adata.obsm["X_AB_baseline"] = pcs
    progress("AB PCA done")
    for key, tag in zip(P["harmony_keys"], ["AB_harmony", "AB_harmony_spec"], strict=True):
        ho = harmonypy.run_harmony(
            pcs, adata.obs, [key], random_state=config.RANDOM_SEED, verbose=False
        )
        z = np.asarray(ho.Z_corr)
        adata.obsm[f"X_{tag}"] = z.T if z.shape[0] != adata.n_obs else z
        progress(f"{tag}: Harmony done")
    for tag in ["AB_baseline", "AB_harmony", "AB_harmony_spec"]:
        _graph_umap_leiden(adata, f"X_{tag}", tag)
        progress(f"{tag}: graph/UMAP/Leiden done")
    out["AB"] = adata
    return out


# ----------------------------------------------------------------------------- metrics


def _rescale_lisi(lisi: np.ndarray, n: int, inverse: bool) -> float:
    m = float(np.median(lisi))
    return (n - m) / (n - 1) if inverse else (m - 1) / (n - 1)


def embedding_metrics(adata: ad.AnnData, tag: str) -> dict:
    rng = np.random.default_rng(config.RANDOM_SEED)
    emb = adata.obsm[f"X_{tag}"]
    obs = adata.obs
    lab = obs["pop_eval"]
    has = lab.notna().to_numpy()
    idx = np.where(has)[0]
    if len(idx) > P["metric_subsample"]:
        idx = np.sort(rng.choice(idx, P["metric_subsample"], replace=False))
    e, o = emb[idx], obs.iloc[idx]
    n_pop = o["pop_eval"].nunique()
    res = {"embedding": tag, "cells": adata.n_obs}

    keys = ["pop_eval"] + (["cohort"] if obs["cohort"].nunique() > 1 else [])
    lisi = harmonypy.compute_lisi(e, o[keys].reset_index(drop=True), keys)
    res["cLISI_scaled"] = _rescale_lisi(lisi[:, 0], n_pop, inverse=True)
    if "cohort" in keys:
        res["iLISI_scaled"] = _rescale_lisi(lisi[:, 1], 2, inverse=False)
        basw = []
        for _pop, g in o.groupby("pop_eval", observed=True):
            if g["cohort"].nunique() == 2 and g["cohort"].value_counts().min() >= 20:
                s = silhouette_samples(e[o.index.get_indexer(g.index)], g["cohort"])
                basw.append(float(np.mean(1 - np.abs(s))))
        res["batch_ASW"] = float(np.mean(basw)) if basw else np.nan
    res["celltype_ASW"] = (silhouette_score(e, o["pop_eval"]) + 1) / 2
    lab_all = lab[has]
    leiden = obs.loc[has, f"leiden_{tag}"]
    res["ARI"] = adjusted_rand_score(lab_all, leiden)
    res["NMI"] = normalized_mutual_info_score(lab_all, leiden)

    # condition structure within (cohort, population): tumor/adjacent silhouette
    cond = []
    for (coh, pop), g in obs[has].groupby(["cohort", "pop_eval"], observed=True):
        vc = g["condition"].value_counts()
        if len(vc) == 2 and vc.min() >= 50:
            gi = g.index
            if len(gi) > 4000:
                gi = pd.Index(rng.choice(gi, 4000, replace=False))
            s = silhouette_score(emb[obs.index.get_indexer(gi)], obs.loc[gi, "condition"])
            cond.append({"embedding": tag, "cohort": coh, "population": pop, "condition_ASW": s})
    return {"summary": res, "condition": pd.DataFrame(cond)}


def knn_purity(adata: ad.AnnData, tag: str, k: int = 15) -> pd.DataFrame:
    has = adata.obs["pop_eval"].notna().to_numpy()
    e = adata.obsm[f"X_{tag}"][has]
    lab = adata.obs.loc[has, "pop_eval"].to_numpy()
    coh = adata.obs.loc[has, "cohort"].to_numpy()
    nn = NearestNeighbors(n_neighbors=k + 1).fit(e)
    _, ind = nn.kneighbors(e)
    same = (lab[ind[:, 1:]] == lab[:, None]).mean(axis=1)
    return (
        pd.DataFrame({"cohort": coh, "population": lab, "purity": same})
        .groupby(["cohort", "population"])["purity"]
        .mean()
        .rename(tag)
        .reset_index()
    )


def label_transfer(adata: ad.AnnData, tag: str, k: int = 15) -> tuple[pd.DataFrame, list]:
    rows, confs = [], []
    for src, dst in [("A", "B"), ("B", "A")]:
        ms = (adata.obs["cohort"] == src) & adata.obs["pop_eval"].notna()
        md = (adata.obs["cohort"] == dst) & adata.obs["pop_eval"].notna()
        clf = KNeighborsClassifier(n_neighbors=k).fit(
            adata.obsm[f"X_{tag}"][ms.to_numpy()], adata.obs.loc[ms, "pop_eval"]
        )
        truth = adata.obs.loc[md, "pop_eval"]
        pred = clf.predict(adata.obsm[f"X_{tag}"][md.to_numpy()])
        rows.append(
            {
                "embedding": tag,
                "direction": f"{src}->{dst}",
                "accuracy": float((pred == truth.to_numpy()).mean()),
                "balanced_accuracy": balanced_accuracy_score(truth, pred),
            }
        )
        conf = pd.crosstab(truth, pd.Series(pred, index=truth.index, name="predicted"))
        conf.attrs.update({"embedding": tag, "direction": f"{src}->{dst}"})
        confs.append(conf)
    return pd.DataFrame(rows), confs


# ----------------------------------------------------------------------------- driver


def run(progress=print) -> dict[str, pd.DataFrame]:
    STAGE3_DIR.mkdir(parents=True, exist_ok=True)
    adata = load_subsample()
    progress(f"subsample: {adata.n_obs:,} cells x {adata.n_vars:,} shared genes")
    emb = build_embeddings(adata, progress=progress)

    summaries, conds, purity, transfer = [], [], [], []
    for key, tags in [
        ("A", ["A_baseline"]),
        ("AB", ["AB_baseline", "AB_harmony", "AB_harmony_spec"]),
    ]:
        a = emb[key]
        for tag in tags:
            m = embedding_metrics(a, tag)
            summaries.append(m["summary"])
            conds.append(m["condition"])
            purity.append(knn_purity(a, tag).set_index(["cohort", "population"]))
            if key == "AB":
                t, confs = label_transfer(a, tag)
                transfer.append(t)
                for conf in confs:
                    d = conf.attrs["direction"].replace("->", "to")
                    conf.to_csv(config.TABLES_DIR / f"stage3_transfer_{tag}_{d}.tsv", sep="\t")
            progress(f"{tag}: metrics done")

    # save obs + coordinates for figures (no expression; small)
    for key, a in emb.items():
        df = a.obs.copy()
        for name, arr in a.obsm.items():
            if name.startswith("X_umap_"):
                df[[f"{name[7:]}_umap1", f"{name[7:]}_umap2"]] = arr
        df.to_csv(STAGE3_DIR / f"embedding_obs_{key}.tsv.gz", sep="\t", compression="gzip")

    out = {
        "stage3_integration_metrics": pd.DataFrame(summaries),
        "stage3_condition_structure": pd.concat(conds, ignore_index=True),
        "stage3_knn_purity": pd.concat(purity, axis=1).reset_index(),
        "stage3_label_transfer": pd.concat(transfer, ignore_index=True),
        "stage3_subsample_composition": adata.obs.groupby(
            ["cohort", "population", "condition"], observed=True
        )
        .size()
        .rename("cells")
        .reset_index(),
    }
    for name, df in out.items():
        df.to_csv(config.TABLES_DIR / f"{name}.tsv", sep="\t", index=False)
    return out
