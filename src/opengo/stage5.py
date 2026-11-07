"""Stage 5: tumor vs adjacent in the stromal (fibroblast + mural) compartment.

Tested values: pseudobulk raw counts summed per specimen over QC-passing cells of the
compartment -> log2-CPM (library size = all genes of the cohort, pseudocount 1 CPM). Never the
integration embedding.

Design:
- Unit = donor. Repeated tumor regions of one donor are averaged (log2-CPM) so each donor counts
  once; adjacent tissue has one specimen per donor.
- Primary analysis = paired: per donor, tumor - adjacent (log2 fold change, log2FC) for each gene
  and for the panel summary. Panel summary per donor = mean over genes of the per-gene difference
  divided by that gene's SD across the cohort's specimens (the paired form of the mean-z score;
  standardized units). The mean log2FC over genes is reported too (interpretable units).
- Per cohort: mean difference, 95 % t-interval, exact sign-flip p (all 2^n sign patterns),
  direction count; Benjamini-Hochberg across the 10 genes.
- Combined: inverse-variance meta-analysis of cohort estimates (fixed and DerSimonian-Laird
  random effects, I^2), and a pooled sign-flip over all paired donors (standardized scale).
- Robustness (main uncertainty = fibroblast:mural mix of the compartment): 8-gene summary without
  the fibroblast-specific DCN/LUM; pericyte-only (B); fibroblast share tumor vs adjacent (B).
- Sensitivities: min cells 20/50; one tumor region per donor instead of the mean; mixed model on all
  specimens (donor random intercept, includes unpaired donors); leave-one-donor-out.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

from opengo import config
from opengo.stage4 import COMBINED_STROMAL, STAGE4_DIR

FOCUS = COMBINED_STROMAL
GENES = config.ECM_GENES
NO_FIB = [g for g in GENES if g not in ("DCN", "LUM")]


def load_long() -> pd.DataFrame:
    return pd.read_csv(STAGE4_DIR / "ecm_pseudobulk_long.tsv.gz", sep="\t")


# ----------------------------------------------------------------------------- donor table


def donor_pairs(
    long: pd.DataFrame,
    population: str = FOCUS,
    min_cells: int = config.PSEUDOBULK_MIN_CELLS,
    region: str = "mean",
) -> dict[str, pd.DataFrame]:
    """Per cohort: paired donors x genes, tumor - adjacent log2FC; plus per-gene SDs."""
    d = long[(long["population"] == population) & (long["n_cells"] >= min_cells)]
    d = d[d["cohort"].isin(config.ANALYSIS_COHORTS)]
    out = {}
    for cohort, dc in d.groupby("cohort"):
        sd = dc.groupby("gene")["logcpm"].std()
        if region == "first":  # one tumor region per donor (lowest specimen id)
            dc = dc.sort_values("sample_id")
            first = dc.groupby(["donor_key", "condition"])["sample_id"].transform("first")
            dc = dc[dc["sample_id"] == first]
        donor = dc.groupby(["donor_key", "condition", "gene"])["logcpm"].mean().unstack("gene")
        donor = donor[GENES]
        conds = donor.index.get_level_values("condition")
        t = donor[conds == "Tumor"].droplevel("condition")
        a = donor[conds == "Adjacent"].droplevel("condition")
        both = t.index.intersection(a.index)
        diff = t.loc[both] - a.loc[both]
        diff.attrs["sd"] = sd[GENES]
        out[cohort] = diff
    return out


def add_summaries(diff: pd.DataFrame) -> pd.DataFrame:
    sd = diff.attrs["sd"]
    z = diff[GENES] / sd
    out = diff.copy()
    out["panel_z"] = z.mean(axis=1)
    out["panel_z_noDCN_LUM"] = z[NO_FIB].mean(axis=1)
    out["panel_log2FC"] = diff[GENES].mean(axis=1)
    out.attrs = {}  # drop metadata so pandas concat works
    return out


# ----------------------------------------------------------------------------- statistics


def sign_flip_p(x: np.ndarray) -> float:
    """Exact two-sided sign-flip p for the mean of paired differences."""
    x = np.asarray(x, float)
    n = len(x)
    if n == 0:
        return np.nan
    obs = abs(x.mean())
    if n <= 20:
        signs = np.array(list(itertools.product([1, -1], repeat=n)))
    else:
        signs = np.random.default_rng(config.RANDOM_SEED).choice([1, -1], size=(100_000, n))
    null = np.abs((signs * x).mean(axis=1))
    return float((null >= obs - 1e-12).mean())


def describe(x: pd.Series) -> dict:
    x = x.dropna().to_numpy(float)
    n = len(x)
    m = x.mean() if n else np.nan
    sd = x.std(ddof=1) if n > 1 else np.nan
    se = sd / np.sqrt(n) if n > 1 else np.nan
    tcrit = stats.t.ppf(0.975, n - 1) if n > 1 else np.nan
    return {
        "n_donors": n,
        "mean": m,
        "sd": sd,
        "se": se,
        "ci_low": m - tcrit * se if n > 1 else np.nan,
        "ci_high": m + tcrit * se if n > 1 else np.nan,
        "n_up": int((x > 0).sum()),
        "p_signflip": sign_flip_p(x),
        "min_possible_p": 2 / 2**n if n else np.nan,
    }


def meta(est: pd.DataFrame) -> dict:
    """Inverse-variance fixed + DerSimonian-Laird random effects over cohort estimates."""
    e = est.dropna(subset=["mean", "se"])
    e = e[e["se"] > 0]
    if len(e) < 2:
        return {}
    w = 1 / e["se"] ** 2
    fe = float((w * e["mean"]).sum() / w.sum())
    se_fe = float(np.sqrt(1 / w.sum()))
    q = float((w * (e["mean"] - fe) ** 2).sum())
    k = len(e)
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - (w**2).sum() / w.sum()))
    wr = 1 / (e["se"] ** 2 + tau2)
    re = float((wr * e["mean"]).sum() / wr.sum())
    se_re = float(np.sqrt(1 / wr.sum()))
    return {
        "fixed": fe,
        "fixed_low": fe - 1.96 * se_fe,
        "fixed_high": fe + 1.96 * se_fe,
        "random": re,
        "random_low": re - 1.96 * se_re,
        "random_high": re + 1.96 * se_re,
        "Q": q,
        "Q_p": float(1 - stats.chi2.cdf(q, k - 1)),
        "I2": max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0,
        "tau2": tau2,
    }


def effects_table(pairs: dict[str, pd.DataFrame], label: str) -> pd.DataFrame:
    rows = []
    cols = GENES + ["panel_z", "panel_z_noDCN_LUM", "panel_log2FC"]
    for cohort, diff in pairs.items():
        diff = add_summaries(diff)
        for col in cols:
            rows.append(
                {"analysis": label, "cohort": cohort, "feature": col, **describe(diff[col])}
            )
    t = pd.DataFrame(rows)
    # BH across the 10 genes within cohort
    t["q_BH"] = np.nan
    for cohort in t["cohort"].unique():
        m = (t["cohort"] == cohort) & t["feature"].isin(GENES) & t["p_signflip"].notna()
        if m.any():
            t.loc[m, "q_BH"] = multipletests(t.loc[m, "p_signflip"], method="fdr_bh")[1]
    # combined
    comb = []
    for col in cols:
        e = t[t["feature"] == col].set_index("cohort")
        mm = meta(e)
        pooled = (
            pd.concat([add_summaries(d)[col] for d in pairs.values()])
            if col.startswith("panel_z")
            else None
        )
        row = {"analysis": label, "cohort": "A+B", "feature": col, **mm}
        if pooled is not None:
            row.update(
                {
                    "n_donors": len(pooled),
                    "n_up": int((pooled > 0).sum()),
                    "p_signflip": sign_flip_p(pooled.to_numpy()),
                }
            )
        comb.append(row)
    return pd.concat([t, pd.DataFrame(comb)], ignore_index=True)


def mixed_model(long: pd.DataFrame, population: str = FOCUS) -> pd.DataFrame:
    """Specimen-level panel z ~ condition + (1 | donor), all donors incl. unpaired, per cohort."""
    d = long[(long["population"] == population) & (long["n_cells"] >= config.PSEUDOBULK_MIN_CELLS)]
    rows = []
    for cohort, dc in d.groupby("cohort"):
        grp = dc.groupby("gene")["logcpm"]
        dc = dc.assign(z=(dc["logcpm"] - grp.transform("mean")) / grp.transform("std"))
        spec = dc.groupby(["sample_id", "donor_key", "condition"])["z"].mean().reset_index()
        spec["tumor"] = (spec["condition"] == "Tumor").astype(int)
        fit = smf.mixedlm("z ~ tumor", spec, groups=spec["donor_key"]).fit(reml=True)
        ci = fit.conf_int().loc["tumor"]
        rows.append(
            {
                "cohort": cohort,
                "specimens": len(spec),
                "donors": spec["donor_key"].nunique(),
                "estimate_z": fit.params["tumor"],
                "ci_low": ci[0],
                "ci_high": ci[1],
                "p_wald": fit.pvalues["tumor"],
            }
        )
    return pd.DataFrame(rows)


def leave_one_out(pairs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for cohort, diff in pairs.items():
        d = add_summaries(diff)["panel_z"]
        for donor in d.index:
            rows.append({"cohort": cohort, "left_out": donor, "panel_z": d.drop(donor).mean()})
    return pd.DataFrame(rows)


def fibroblast_share(long: pd.DataFrame) -> pd.DataFrame:
    """B only: fibroblasts / (fibroblasts + pericytes) per donor and condition."""
    b = long[(long["cohort"] == "B") & (long["gene"] == GENES[0])]
    n = b.pivot_table(
        index=["donor_key", "condition", "sample_id"],
        columns="population",
        values="n_cells",
        fill_value=0,
    )
    n = n.groupby(["donor_key", "condition"])[["Stromal_fibro_myo", "Pericyte"]].sum()
    n["fib_share"] = n["Stromal_fibro_myo"] / (n["Stromal_fibro_myo"] + n["Pericyte"])
    return n.reset_index()


# ----------------------------------------------------------------------------- driver


def run() -> dict[str, pd.DataFrame]:
    long = load_long()
    primary = donor_pairs(long)
    tables = {}
    eff = [effects_table(primary, "primary: stroma (fibro+mural), >=10 cells, regions averaged")]
    eff.append(effects_table(donor_pairs(long, min_cells=20), "sensitivity: >=20 cells"))
    eff.append(effects_table(donor_pairs(long, min_cells=50), "sensitivity: >=50 cells"))
    eff.append(effects_table(donor_pairs(long, region="first"), "sensitivity: one tumor region"))
    eff.append(
        effects_table(donor_pairs(long, population="Pericyte"), "robustness: pericytes only")
    )
    eff.append(effects_table(donor_pairs(long, population="Endothelial"), "contrast: endothelial"))
    tables["stage5_effects"] = pd.concat(eff, ignore_index=True)

    per_donor = []
    for cohort, diff in primary.items():
        d = add_summaries(diff).reset_index()
        d.insert(0, "cohort", cohort)
        per_donor.append(d)
    tables["stage5_per_donor_log2FC"] = pd.concat(per_donor, ignore_index=True)

    levels = long[(long["population"] == FOCUS) & (long["n_cells"] >= config.PSEUDOBULK_MIN_CELLS)]
    tables["stage5_donor_levels"] = (
        levels.groupby(["cohort", "donor_key", "condition", "gene"])["logcpm"].mean().reset_index()
    )
    tables["stage5_mixed_model"] = mixed_model(long)
    tables["stage5_leave_one_out"] = leave_one_out(primary)
    fs = fibroblast_share(long)
    tables["stage5_B_fibroblast_share"] = fs
    w = fs.pivot(index="donor_key", columns="condition", values="fib_share").dropna()
    tables["stage5_B_fibroblast_share_test"] = pd.DataFrame(
        [{"feature": "fibroblast share (T - Adj)", **describe(w["Tumor"] - w["Adjacent"])}]
    )
    for name, df in tables.items():
        df.to_csv(config.TABLES_DIR / f"{name}.tsv", sep="\t", index=False)
    return tables
