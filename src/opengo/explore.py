"""Stage 2 exploratory structure: pseudobulk, marker-based labels, stromal identity checks.

1. Pseudobulk raw counts per specimen x population (QC-passing cells), over the genes shared by
   all cohorts. Saved as a checkpoint and reused in Stages 4-6. "ALL" = whole specimen.
2. Specimen-level pseudobulk PCA: how much variation follows cohort vs condition vs donor?
3. Marker-based coarse labelling (cohort C has no labels). The same procedure is first run on
   subsamples of A and B and scored against the author labels, so its accuracy is known before
   it is trusted on C.
4. Stromal sub-analysis of A and B: are A "Myofibroblasts" partly mural cells? Where do B's
   fibroblast/pericyte/endothelial label conflicts ("Unresolved") fall?

ECM panel genes are never used as markers (circularity).
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import statsmodels.formula.api as smf

from opengo import config
from opengo.qc import qc_obs_path

MARKERS: dict[str, list[str]] = {
    "Malignant_epithelial": ["CA9", "NDUFA4L2", "ANGPTL4", "EGLN3"],
    "Nonmalignant_epithelial": ["LRP2", "CUBN", "UMOD", "SLC34A1", "AQP2", "SLC12A1"],
    "Endothelial": ["VWF", "CDH5", "RAMP2", "PLVAP", "KDR"],
    "Fibroblast": ["PDGFRA", "C7", "FBLN1", "MFAP5"],
    "Mural": ["RGS5", "NOTCH3", "MYH11", "PDGFRB"],
    "Myeloid": ["CD68", "CD14", "LYZ", "CSF1R", "C1QB"],
    "Lymphoid": ["CD3E", "CD2", "NKG7", "MS4A1", "IL7R"],
}
STROMAL_MARKERS = {k: MARKERS[k] for k in ("Fibroblast", "Mural", "Endothelial")}
PSEUDOBULK_PATH = config.CHECKPOINTS_DIR / "pseudobulk_raw.h5ad"


# ----------------------------------------------------------------------------- helpers


class _CsrKnn:
    """pynndescent kNN transformer returning scipy `csr_matrix`.

    Works around scanpy 1.11.5 calling `.getnnz` on the `csr_array` that pynndescent now returns
    for large inputs (> ~8k cells). Same algorithm as scanpy's default for large n.
    """

    def __init__(self, n_neighbors: int = 15, random_state: int = config.RANDOM_SEED):
        self.n_neighbors = n_neighbors
        self.random_state = random_state

    def fit_transform(self, x, y=None):
        from pynndescent import PyNNDescentTransformer

        t = PyNNDescentTransformer(n_neighbors=self.n_neighbors, random_state=self.random_state)
        return sp.csr_matrix(t.fit_transform(x))

    def get_params(self, deep: bool = True) -> dict:
        return {"n_neighbors": self.n_neighbors, "random_state": self.random_state}


def shared_genes() -> pd.DataFrame:
    g = pd.read_csv(config.TABLES_DIR / "gene_correspondence.tsv", sep="\t", low_memory=False)
    return g[g["shared_ABC"]].reset_index(drop=True)


def _cohort_gene_index(cohort: str, var: pd.DataFrame, shared: pd.DataFrame) -> np.ndarray:
    key = shared["symbol_B"] if cohort == "B" else shared["gene_id"]
    idx = pd.Index(var.index.astype(str)).get_indexer(key.astype(str))
    if (idx < 0).any():
        raise ValueError(f"{cohort}: {(idx < 0).sum()} shared genes not found")
    return idx


def qc_cells(cohort: str) -> pd.DataFrame:
    return pd.read_csv(qc_obs_path(cohort), sep="\t", index_col=0, low_memory=False)


def _symbols_unique(adata: ad.AnnData) -> None:
    adata.var_names = adata.var["symbol"].astype(str).to_numpy()
    adata.var_names_make_unique()


# ----------------------------------------------------------------------------- pseudobulk


def build_pseudobulk(progress=None) -> ad.AnnData:
    shared = shared_genes()
    rows, mats = [], []
    for cohort in config.COHORTS:
        cells = qc_cells(cohort)
        backed = ad.read_h5ad(config.cohort_raw_path(cohort), backed="r")
        gidx = _cohort_gene_index(cohort, backed.var, shared)
        obs = backed.obs
        for sid, idx in obs.groupby("sample_id", observed=True).indices.items():
            idx = np.sort(idx)
            sub = backed[idx].to_memory()
            ok = cells.loc[sub.obs_names, "qc_pass"].to_numpy(bool)
            pops = sub.obs["population"].astype(str).to_numpy()
            x = sp.csr_matrix(sub.X)
            lib_all = np.asarray(x.sum(axis=1)).ravel()
            xs = x[:, gidx]
            meta = sub.obs.iloc[0]
            for pop in ["ALL", *sorted(set(pops))]:
                m = ok if pop == "ALL" else ok & (pops == pop)
                if m.sum() == 0:
                    continue
                mats.append(sp.csr_matrix(xs[m].sum(axis=0)))
                rows.append(
                    {
                        "cohort": cohort,
                        "sample_id": sid,
                        "donor_key": meta["donor_key"],
                        "condition": meta["condition"],
                        "assay": meta["assay"],
                        "population": pop,
                        "n_cells": int(m.sum()),
                        "lib_size_all_genes": float(lib_all[m].sum()),
                    }
                )
            if progress:
                progress(cohort, sid)
        backed.file.close()
    obs = pd.DataFrame(rows)
    obs.index = (obs["sample_id"] + "|" + obs["population"]).to_numpy()
    var = shared.set_index("gene_id")[["symbol", "symbol_A", "symbol_B", "symbol_C"]]
    pb = ad.AnnData(X=sp.vstack(mats).tocsr().astype(np.float64), obs=obs, var=var)
    pb.write_h5ad(PSEUDOBULK_PATH)
    return pb


def logcpm(pb: ad.AnnData) -> np.ndarray:
    x = pb.X.toarray() if sp.issparse(pb.X) else np.asarray(pb.X)
    return np.log2(x / x.sum(axis=1, keepdims=True) * 1e6 + 1)


def pseudobulk_pca(pb: ad.AnnData, n_hvg: int = 2000, n_pcs: int = 10) -> dict:
    """PCA on log-CPM of the most variable genes + R² of cohort/condition/donor per PC."""
    y = logcpm(pb)
    keep = y.mean(axis=0) > 1  # expressed genes
    y = y[:, keep]
    top = np.argsort(y.var(axis=0))[::-1][:n_hvg]
    z = y[:, top] - y[:, top].mean(axis=0)
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    pcs = u[:, :n_pcs] * s[:n_pcs]
    var_exp = (s**2 / (s**2).sum())[:n_pcs]
    scores = pd.DataFrame(pcs, columns=[f"PC{i + 1}" for i in range(n_pcs)], index=pb.obs_names)
    scores = scores.join(pb.obs[["cohort", "condition", "donor_key", "sample_id", "n_cells"]])
    r2 = []
    for i in range(n_pcs):
        pc = f"PC{i + 1}"
        for term in ["cohort", "condition", "donor_key"]:
            fit = smf.ols(f"{pc} ~ C({term})", data=scores).fit()
            r2.append({"PC": pc, "var_explained": var_exp[i], "term": term, "R2": fit.rsquared})
    return {"scores": scores, "r2": pd.DataFrame(r2)}


# ----------------------------------------------------------------------------- marker labels


def _subsample(cohort: str, n: int, mask_fn=None) -> ad.AnnData:
    rng = np.random.default_rng(config.RANDOM_SEED)
    cells = qc_cells(cohort)
    backed = ad.read_h5ad(config.cohort_raw_path(cohort), backed="r")
    ok = cells.reindex(backed.obs_names)["qc_pass"].to_numpy(bool)
    if mask_fn is not None:
        ok &= mask_fn(backed.obs)
    idx = np.where(ok)[0]
    if n and len(idx) > n:
        idx = np.sort(rng.choice(idx, n, replace=False))
    sub = backed[idx].to_memory()
    backed.file.close()
    return sub


def marker_label(adata: ad.AnnData, markers: dict[str, list[str]], resolution: float = 1.0):
    """Cluster, score marker sets, assign each cluster to its top-scoring set.

    Returns (adata with obs['marker_label', 'leiden', 'label_margin'], cluster table).
    """
    adata = adata.copy()
    _symbols_unique(adata)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    batch = "sample_id" if adata.obs["sample_id"].nunique() > 1 else None
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=batch, flavor="seurat")
    sc.pp.pca(adata, n_comps=30, mask_var="highly_variable", random_state=config.RANDOM_SEED)
    sc.pp.neighbors(adata, n_neighbors=15, transformer=_CsrKnn(15), random_state=config.RANDOM_SEED)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        random_state=config.RANDOM_SEED,
    )
    sc.tl.umap(adata, random_state=config.RANDOM_SEED)
    used = {}
    for name, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        used[name] = present
        x = adata[:, present].X
        x = x.toarray() if sp.issparse(x) else np.asarray(x)
        z = (x - x.mean(axis=0)) / (x.std(axis=0) + 1e-9)
        adata.obs[f"score_{name}"] = z.mean(axis=1)
    score_cols = [f"score_{k}" for k in markers]
    cl = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    top2 = np.sort(cl.to_numpy(), axis=1)[:, ::-1][:, :2]
    cl["label"] = cl[score_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    cl["margin"] = top2[:, 0] - top2[:, 1]
    cl["n_cells"] = adata.obs["leiden"].value_counts().reindex(cl.index)
    adata.obs["marker_label"] = adata.obs["leiden"].map(cl["label"]).astype(str)
    adata.obs["label_margin"] = adata.obs["leiden"].map(cl["margin"]).astype(float)
    cl.attrs["markers_used"] = used
    return adata, cl


def compare_to_author(adata: ad.AnnData) -> pd.DataFrame:
    """Confusion table: author population vs marker label (Fibroblast+Mural -> stromal)."""
    lab = adata.obs["marker_label"].replace({"Fibroblast": "Stromal", "Mural": "Stromal"})
    auth = (
        adata.obs["population"]
        .astype(str)
        .replace({"Stromal_fibro_myo": "Stromal", "Pericyte": "Stromal"})
    )
    return pd.crosstab(auth, lab)


def label_accuracy(conf: pd.DataFrame) -> float:
    conf = conf.drop(index=[i for i in conf.index if i in ("Unresolved",)], errors="ignore")
    hits = sum(conf.loc[i, i] for i in conf.index if i in conf.columns)
    return float(hits / conf.to_numpy().sum())
