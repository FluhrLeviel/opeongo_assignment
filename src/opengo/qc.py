"""Stage 2: per-cell QC metrics and per-specimen, assay-aware outlier flags.

Nothing is removed here. Every cell gets metrics and flags; the inclusion rule is decided at
Check-in 2 from the figures and then applied downstream from `qc_pass`.

Flag rules (per specimen, so thresholds adapt to each specimen's depth):
- low_counts / low_genes: log1p value below median - QC_NMADS_LOW * MAD, or below absolute floors.
- high_counts: log1p total counts above median + QC_NMADS_HIGH * MAD (possible multiplets).
- high_mt: % mitochondrial above median + QC_NMADS_MT * MAD **and** above an assay floor.
  Nuclei (snRNA) should carry ~0 % mitochondrial reads, so a high value there signals ambient RNA /
  incomplete cytoplasm stripping; whole cells (scRNA) from kidney can be biologically high
  (proximal tubule), hence a permissive floor.
- doublet: Scrublet call per specimen (B also carries the authors' doublet score for comparison).
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from opengo import config


def qc_obs_path(cohort: str):
    return config.CHECKPOINTS_DIR / f"qc_cells_{cohort}.tsv.gz"


def _mad(x: np.ndarray) -> float:
    return float(np.median(np.abs(x - np.median(x)))) * 1.4826


def _gene_masks(symbols: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    s = symbols.astype(str).str.upper()
    mt = s.str.startswith("MT-").to_numpy()
    ribo = s.str.match(r"^RP[SL]\d").to_numpy()
    return mt, ribo


def specimen_metrics(adata: ad.AnnData, assay: str) -> pd.DataFrame:
    """Metrics + flags for one specimen (adata holds raw counts, cells x genes)."""
    x = sp.csr_matrix(adata.X)
    mt, ribo = _gene_masks(adata.var["symbol"])
    total = np.asarray(x.sum(axis=1)).ravel()
    n_genes = np.diff(x.indptr)
    pct_mt = 100 * np.asarray(x[:, mt].sum(axis=1)).ravel() / np.maximum(total, 1)
    pct_ribo = 100 * np.asarray(x[:, ribo].sum(axis=1)).ravel() / np.maximum(total, 1)
    df = pd.DataFrame(
        {
            "total_counts": total,
            "n_genes": n_genes,
            "pct_mt": pct_mt,
            "pct_ribo": pct_ribo,
        },
        index=adata.obs_names,
    )
    lt, lg = np.log1p(total), np.log1p(n_genes)
    p = config.QC_PARAMS
    df["low_counts"] = (lt < np.median(lt) - p["nmads_low"] * _mad(lt)) | (total < p["min_counts"])
    df["low_genes"] = (lg < np.median(lg) - p["nmads_low"] * _mad(lg)) | (n_genes < p["min_genes"])
    df["high_counts"] = lt > np.median(lt) + p["nmads_high"] * _mad(lt)
    floor = p["mt_floor_nuclei"] if assay.startswith("snRNA") else p["mt_floor_cells"]
    df["high_mt"] = (pct_mt > np.median(pct_mt) + p["nmads_mt"] * _mad(pct_mt)) & (pct_mt > floor)

    # Scrublet on raw counts of this specimen only (doublets form within a capture).
    tmp = ad.AnnData(X=x.copy(), obs=pd.DataFrame(index=adata.obs_names))
    try:
        sc.pp.scrublet(tmp, random_state=config.RANDOM_SEED)
        df["doublet_score"] = tmp.obs["doublet_score"].to_numpy()
        df["doublet"] = tmp.obs["predicted_doublet"].astype(bool).to_numpy()
    except Exception as err:  # very small specimens can fail; record, don't crash
        df["doublet_score"] = np.nan
        df["doublet"] = False
        df["doublet_note"] = str(err)[:80]
    flags = ["low_counts", "low_genes", "high_counts", "high_mt", "doublet"]
    df["qc_flags"] = df[flags].apply(lambda r: ",".join(f for f in flags if r[f]), axis=1)
    df["qc_pass"] = ~df[flags].any(axis=1)
    return df


def run_qc(cohort: str, progress=None) -> pd.DataFrame:
    """Compute metrics for every cell of a cohort, one specimen at a time (memory-safe)."""
    backed = ad.read_h5ad(config.cohort_raw_path(cohort), backed="r")
    obs = backed.obs.copy()
    var = backed.var.copy()
    out = []
    for sid, idx in obs.groupby("sample_id", observed=True).indices.items():
        sub = backed[np.sort(idx)].to_memory()
        sub.var = var
        assay = str(obs["assay"].iloc[idx[0]])
        m = specimen_metrics(sub, assay)
        out.append(m)
        if progress:
            progress(sid, m)
    backed.file.close()
    metrics = pd.concat(out)
    keep = [
        "cohort",
        "sample_id",
        "donor_key",
        "condition",
        "assay",
        "population",
        "author_broad",
    ]
    extra = [c for c in obs.columns if c.startswith("author_") and c not in keep]
    cols = [c for c in keep + extra if c in obs.columns]
    res = obs[cols].join(metrics)
    res.to_csv(qc_obs_path(cohort), sep="\t", compression="gzip")
    return res


def summarise(cells: pd.DataFrame) -> dict[str, pd.DataFrame]:
    flags = ["low_counts", "low_genes", "high_counts", "high_mt", "doublet"]
    by_spec = cells.groupby(["cohort", "sample_id", "condition"], observed=True).agg(
        cells=("qc_pass", "size"),
        pass_frac=("qc_pass", "mean"),
        median_counts=("total_counts", "median"),
        median_genes=("n_genes", "median"),
        median_pct_mt=("pct_mt", "median"),
        **{f"{f}_frac": (f, "mean") for f in flags},
    )
    by_pop = cells.groupby(["cohort", "population"], observed=True).agg(
        cells=("qc_pass", "size"),
        pass_frac=("qc_pass", "mean"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_mt", "median"),
        **{f"{f}_frac": (f, "mean") for f in flags},
    )
    return {
        "stage2_qc_by_specimen": by_spec.reset_index(),
        "stage2_qc_by_population": by_pop.reset_index(),
    }
