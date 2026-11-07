"""Stage 2 figures -> outputs/step_2_output/. Run with `uv run opengo figures 2`."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from opengo import config  # noqa: E402
from opengo.explore import MARKERS, STROMAL_MARKERS  # noqa: E402
from opengo.qc import qc_obs_path  # noqa: E402
from opengo.stage2 import STAGE2_DIR  # noqa: E402

OUT_DIR = config.OUTPUTS_DIR / "step_2_output"
COND = {"Tumor": "#c0392b", "Adjacent": "#2471a3"}
COHORT_COLORS = {"A": "#1f77b4", "B": "#ff7f0e", "C": "#2ca02c"}
FLAGS = ["low_counts", "low_genes", "high_counts", "high_mt", "doublet"]
T = config.TABLES_DIR


def _save(fig, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / name
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return p


def _cells() -> pd.DataFrame:
    cols = ["cohort", "sample_id", "condition", "population", "total_counts", "n_genes"]
    cols += ["pct_mt", "doublet_score", "qc_pass", *FLAGS]
    parts = []
    for c in config.COHORTS:
        d = pd.read_csv(qc_obs_path(c), sep="\t", index_col=0, low_memory=False)
        parts.append(d[[x for x in cols if x in d.columns]])
    return pd.concat(parts)


# --------------------------------------------------------------------- 1. QC distributions


def fig_qc_distributions(cells: pd.DataFrame) -> Path:
    metrics = [
        ("total_counts", "counts per cell (log10)", True),
        ("n_genes", "genes per cell (log10)", True),
        ("pct_mt", "% mitochondrial", False),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(20, 11), gridspec_kw={"width_ratios": [34, 26, 13]})
    for col, cohort in enumerate(config.COHORTS):
        d = cells[cells["cohort"] == cohort]
        order = (
            d.groupby("sample_id")
            .agg(cond=("condition", "first"), med=("total_counts", "median"))
            .sort_values(["cond", "med"])
            .index
        )
        for row, (m, label, logy) in enumerate(metrics):
            ax = axes[row, col]
            vals = np.log10(d[m].clip(lower=1)) if logy else d[m]
            data = [vals[d["sample_id"] == s].to_numpy() for s in order]
            parts = ax.violinplot(data, showextrema=False, widths=0.85)
            conds = d.groupby("sample_id")["condition"].first().reindex(order)
            for body, c in zip(parts["bodies"], conds, strict=True):
                body.set_facecolor(COND[c])
                body.set_alpha(0.7)
            ax.scatter(
                np.arange(1, len(order) + 1), [np.median(x) for x in data], s=6, c="k", zorder=3
            )
            ax.set_xticks([])
            if col == 0:
                ax.set_ylabel(label)
            if row == 0:
                assay = "snRNA" if cohort == "A" else "scRNA"
                ax.set_title(f"Cohort {cohort} ({assay}), {len(order)} specimens")
            if m == "pct_mt":
                ax.set_ylim(0, min(100, np.percentile(d[m], 99.5) * 1.3 + 1))
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.7) for c in COND.values()]
    fig.legend(handles, COND.keys(), loc="upper right")
    fig.suptitle(
        "Per-specimen QC distributions (sorted by condition, then median depth; dot = median)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "01_qc_distributions.png")


# --------------------------------------------------------------------- 2. flag rates


def fig_flag_rates(cells: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(20, 10), gridspec_kw={"height_ratios": [1.1, 1]})
    ax = axes[0]
    spec = cells.groupby(["cohort", "condition", "sample_id"])[FLAGS + ["qc_pass"]].mean()
    spec = spec.sort_index()
    bottom = np.zeros(len(spec))
    colors = sns.color_palette("Set2", len(FLAGS))
    for f, c in zip(FLAGS, colors, strict=True):
        ax.bar(np.arange(len(spec)), spec[f], bottom=bottom, color=c, label=f, width=0.85)
        bottom += spec[f].to_numpy()
    ax.plot(np.arange(len(spec)), 1 - spec["qc_pass"], "k_", ms=8, label="any flag (union)")
    ax.set_xticks(
        np.arange(len(spec)),
        [f"{c}:{s}" for c, _, s in spec.index],
        rotation=90,
        fontsize=6,
    )
    for i, (_, cond, _) in enumerate(spec.index):
        ax.get_xticklabels()[i].set_color(COND[cond])
    ax.set_ylabel("fraction of cells flagged (stacked)")
    ax.set_title("Flag rates per specimen (red label = tumor, blue = adjacent)")
    ax.legend(ncol=6, fontsize=8, loc="upper left")

    ax = axes[1]
    pop = cells.groupby(["cohort", "population"])[FLAGS + ["qc_pass"]].mean()
    n = cells.groupby(["cohort", "population"]).size()
    pop["any flag"] = 1 - pop.pop("qc_pass")
    pop.index = [f"{c}: {p} (n={n[(c, p)]:,})" for c, p in pop.index]
    sns.heatmap(
        pop.T * 100,
        annot=True,
        fmt=".0f",
        cmap="Reds",
        vmin=0,
        vmax=40,
        ax=ax,
        cbar_kws={"label": "% cells flagged"},
        annot_kws={"fontsize": 7},
    )
    ax.set_title(
        "Flag rates per population — are flags technical, or concentrated in one cell type?"
    )
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    return _save(fig, "02_flag_rates.png")


# --------------------------------------------------------------------- 3. B author QC check


def fig_b_author_check() -> Path:
    b = pd.read_csv(qc_obs_path("B"), sep="\t", index_col=0, low_memory=False)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.2))
    pairs = [
        ("author_total_umi", "total_counts", "TotalUMI (author)", "counts (recomputed)", True),
        ("author_n_genes", "n_genes", "gene.num (author)", "genes (recomputed)", True),
        ("author_mt_ratio", "pct_mt", "MT.ratio ×100 (author)", "% mt (recomputed)", False),
        (
            "author_doublet_score",
            "doublet_score",
            "doublet score (author)",
            "Scrublet (ours)",
            False,
        ),
    ]
    rng = np.random.default_rng(0)
    idx = rng.choice(len(b), min(20000, len(b)), replace=False)
    s = b.iloc[idx]
    for ax, (x, y, xl, yl, log) in zip(axes, pairs, strict=True):
        xv = s[x] * (100 if x == "author_mt_ratio" else 1)
        ax.scatter(xv, s[y], s=2, alpha=0.3)
        rho = b[x].corr(b[y], method="spearman")
        if log:
            ax.set_xscale("log")
            ax.set_yscale("log")
        lim = [min(xv.min(), s[y].min()), max(xv.max(), s[y].max())]
        if x != "author_doublet_score":
            ax.plot(lim, lim, "k--", lw=0.7)
        ax.set(xlabel=xl, ylabel=yl, title=f"Spearman ρ = {rho:.3f}")
    fig.suptitle("Cohort B: author QC fields vs our recomputation from the deposited counts")
    fig.tight_layout()
    return _save(fig, "03_B_author_qc_vs_recomputed.png")


# --------------------------------------------------------------------- 4. pseudobulk PCA


def fig_pseudobulk_pca() -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    for row, label in enumerate(["ALL", "Endothelial"]):
        sc_ = pd.read_csv(STAGE2_DIR / f"pca_scores_{label}.tsv", sep="\t", index_col=0)
        r2 = pd.read_csv(T / f"stage2_pseudobulk_pca_r2_{label}.tsv", sep="\t")
        ve = r2.drop_duplicates("PC").set_index("PC")["var_explained"]
        for col, (pcx, pcy) in enumerate([("PC1", "PC2"), ("PC3", "PC4")]):
            ax = axes[row, col]
            for (coh, cond), d in sc_.groupby(["cohort", "condition"]):
                ax.scatter(
                    d[pcx],
                    d[pcy],
                    c=COHORT_COLORS[coh],
                    marker="o" if cond == "Tumor" else "^",
                    s=45,
                    edgecolor="k" if cond == "Adjacent" else "none",
                    label=f"{coh} {cond}",
                    alpha=0.85,
                )
            ax.set(
                xlabel=f"{pcx} ({ve[pcx]:.0%})",
                ylabel=f"{pcy} ({ve[pcy]:.0%})",
                title=f"{'Whole-specimen' if label == 'ALL' else 'Endothelial'} pseudobulk",
            )
            if col == 0:
                ax.legend(fontsize=7)
        ax = axes[row, 2]
        piv = r2.pivot(index="PC", columns="term", values="R2").reindex(ve.index)[:6]
        piv.plot.bar(ax=ax, width=0.8, color=["#7f7f7f", "#c0392b", "#8e44ad"])
        ax.set(ylabel="R² of PC explained by term", ylim=(0, 1), title="What drives each PC?")
        ax.set_xticklabels([f"{p}\n({ve[p]:.0%})" for p in piv.index], rotation=0)
    fig.suptitle(
        "Exploratory pseudobulk PCA (log-CPM, 2,000 most variable shared genes; "
        "circle = tumor, triangle = adjacent)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "04_pseudobulk_pca.png")


# --------------------------------------------------------------------- 5. marker validation


def fig_marker_validation() -> Path:
    acc = pd.read_csv(T / "stage2_marker_label_accuracy.tsv", sep="\t").set_index("cohort")
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    for ax, c in zip(axes, ("A", "B"), strict=True):
        conf = pd.read_csv(T / f"stage2_marker_vs_author_{c}.tsv", sep="\t", index_col=0)
        frac = conf.div(conf.sum(axis=1), axis=0)
        sns.heatmap(frac, annot=conf, fmt="d", cmap="Blues", vmin=0, vmax=1, ax=ax, cbar=False)
        ax.set(
            xlabel="marker-based label (our procedure)",
            ylabel="author population",
            title=f"Cohort {c}: agreement {acc.loc[c, 'accuracy']:.1%} "
            f"(n={acc.loc[c, 'cells']:,}, Unresolved excluded)",
        )
        plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
    fig.suptitle(
        "Validating the marker-based labelling on labelled cohorts before using it on C "
        "(colour = row fraction)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "05_marker_label_validation.png")


# --------------------------------------------------------------------- 6. cohort C labels


def fig_c_labels() -> Path:
    obs = pd.read_csv(STAGE2_DIR / "C_marker_labels.tsv", sep="\t", index_col=0)
    dot = pd.read_csv(STAGE2_DIR / "C_marker_dot.tsv", sep="\t")
    fig = plt.figure(figsize=(20, 6.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.6])
    ax = fig.add_subplot(gs[0])
    pal = dict(zip(sorted(obs["marker_label"].unique()), sns.color_palette("tab10"), strict=False))
    for lab, d in obs.groupby("marker_label"):
        ax.scatter(d["umap_1"], d["umap_2"], s=1, color=pal[lab], label=f"{lab} ({len(d):,})")
    ax.legend(markerscale=8, fontsize=7, loc="best")
    ax.set(title="Cohort C UMAP: marker-based labels", xticks=[], yticks=[])
    ax = fig.add_subplot(gs[1])
    for cond, d in obs.groupby("condition"):
        ax.scatter(d["umap_1"], d["umap_2"], s=1, color=COND[cond], label=cond, alpha=0.5)
    ax.legend(markerscale=8)
    ax.set(title="Cohort C UMAP: condition", xticks=[], yticks=[])
    ax = fig.add_subplot(gs[2])
    genes = [g for gs_ in MARKERS.values() for g in gs_ if g in set(dot["gene"])]
    groups = [g for g in MARKERS if g in set(dot["group"])]
    d = dot[dot["group"].isin(groups)]
    xi = {g: i for i, g in enumerate(genes)}
    yi = {g: i for i, g in enumerate(groups)}
    sca = ax.scatter(
        d["gene"].map(xi),
        d["group"].map(yi),
        s=d["frac"] * 220,
        c=d["mean"],
        cmap="viridis",
    )
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=8)
    ax.set_yticks(
        range(len(groups)), [f"{g} (n={int(d[d.group == g].n.iloc[0]):,})" for g in groups]
    )
    ax.invert_yaxis()
    plt.colorbar(sca, ax=ax, label="mean log-norm expr")
    ax.set_title("Marker dot plot (size = fraction expressing)")
    fig.tight_layout()
    return _save(fig, "06_cohortC_marker_labels.png")


# --------------------------------------------------------------------- 7. stromal identity


def fig_stromal() -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(19, 12), gridspec_kw={"width_ratios": [1.3, 1]})
    genes = [g for gs_ in STROMAL_MARKERS.values() for g in gs_]
    for row, c in enumerate(("A", "B")):
        dot = pd.read_csv(STAGE2_DIR / f"stromal_{c}_marker_dot.tsv", sep="\t")
        piv = dot.pivot(index="group", columns="gene", values="frac")[
            [g for g in genes if g in set(dot["gene"])]
        ]
        n = dot.drop_duplicates("group").set_index("group")["n"]
        piv.index = [f"{g} (n={n[g]:,})" for g in piv.index]
        sns.heatmap(
            piv * 100,
            annot=True,
            fmt=".0f",
            cmap="Greens",
            vmin=0,
            vmax=80,
            ax=axes[row, 0],
            cbar_kws={"label": "% cells expressing"},
            annot_kws={"fontsize": 7},
        )
        axes[row, 0].set(title=f"Cohort {c}: stromal markers by author label", xlabel="", ylabel="")
        ct = pd.read_csv(T / f"stage2_stromal_{c}_author_vs_marker.tsv", sep="\t", index_col=0)
        frac = ct.div(ct.sum(axis=1), axis=0)
        sns.heatmap(
            frac,
            annot=ct,
            fmt="d",
            cmap="Purples",
            vmin=0,
            vmax=1,
            ax=axes[row, 1],
            cbar=False,
            annot_kws={"fontsize": 7},
        )
        axes[row, 1].set(
            title=f"Cohort {c}: author label vs marker-based sub-cluster label",
            xlabel="marker label",
            ylabel="",
        )
        plt.setp(axes[row, 1].get_yticklabels(), fontsize=7)
    fig.suptitle(
        "Stromal identity check with non-panel markers (fibroblast: PDGFRA, C7, FBLN1, MFAP5; "
        "mural: RGS5, NOTCH3, MYH11, PDGFRB; endothelial: VWF, CDH5, RAMP2, PLVAP, KDR)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "07_stromal_identity.png")


# --------------------------------------------------------------------- 8. cohort C integrity


def _panel_detection(cohort: str, n: int = 20000) -> pd.Series:
    """Fraction of cells with >=1 UMI for each ECM panel gene (random subsample)."""
    import anndata as ad

    rng = np.random.default_rng(config.RANDOM_SEED)
    a = ad.read_h5ad(config.cohort_raw_path(cohort), backed="r")
    idx = np.arange(a.n_obs)
    if a.n_obs > n:
        idx = np.sort(rng.choice(a.n_obs, n, replace=False))
    sym = a.var["symbol"].astype(str).to_numpy()
    cols = [int(np.where(sym == g)[0][0]) for g in config.ECM_GENES]
    x = a[idx].to_memory().X[:, cols]
    a.file.close()
    return pd.Series(np.asarray((x > 0).mean(axis=0)).ravel(), index=config.ECM_GENES)


def fig_c_integrity(cells: pd.DataFrame) -> Path:
    cells = cells.assign(nonmt=cells["total_counts"] * (1 - cells["pct_mt"] / 100))
    fig, axes = plt.subplots(1, 3, figsize=(20, 5.2), gridspec_kw={"width_ratios": [1, 1.2, 1]})
    ax = axes[0]
    bins = np.logspace(0.5, 4.5, 80)
    for c in config.COHORTS:
        d = cells.loc[cells["cohort"] == c, "nonmt"]
        ax.hist(
            d,
            bins=bins,
            histtype="step",
            lw=2,
            density=True,
            color=COHORT_COLORS[c],
            label=f"{c}: max {d.max():,.0f}",
        )
    ax.set_xscale("log")
    ax.set(
        xlabel="non-mitochondrial counts per cell (log)",
        ylabel="density",
        title="Non-mitochondrial counts: C has a hard ceiling (~400)",
    )
    ax.legend()
    ax = axes[1]
    d = cells[cells["cohort"] == "C"]
    order = d.groupby("sample_id")["pct_mt"].median().sort_values().index
    sns.boxplot(
        data=d,
        x="sample_id",
        y="pct_mt",
        order=order,
        hue="condition",
        palette=COND,
        ax=ax,
        fliersize=0,
        dodge=False,
    )
    ax.set(xlabel="", ylabel="% mitochondrial", title="Cohort C: % mitochondrial per specimen")
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
    ax = axes[2]
    det = pd.DataFrame({c: _panel_detection(c) for c in config.COHORTS})
    det.plot.barh(ax=ax, color=[COHORT_COLORS[c] for c in det.columns], width=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="fraction of all cells with ≥1 UMI", title="ECM panel gene detection by cohort")
    fig.suptitle(
        "Cohort C integrity: mitochondrial-dominated, capped non-mitochondrial signal",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "08_cohortC_integrity.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="notebook")
    cells = _cells()
    return [
        fig_qc_distributions(cells),
        fig_flag_rates(cells),
        fig_b_author_check(),
        fig_pseudobulk_pca(),
        fig_marker_validation(),
        fig_c_labels(),
        fig_stromal(),
        fig_c_integrity(cells),
    ]
