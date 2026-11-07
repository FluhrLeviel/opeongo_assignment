"""Report figures (one per assignment section) -> report/figures/. `uv run opengo figures report`.

Built only from tables and small checkpoints written by Stages 1-6, so the report can be
regenerated without recomputing anything heavy.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from opengo import config  # noqa: E402

OUT_DIR = config.PROJECT_ROOT / "report" / "figures"
T = config.TABLES_DIR
GENES = config.ECM_GENES
COH = {"A": "#1f77b4", "B": "#ff7f0e", "C": "#2ca02c"}
COND = {"Tumor": "#c0392b", "Adjacent": "#2471a3"}
POP_COLORS = {
    "Malignant_epithelial": "#7b241c",
    "Nonmalignant_epithelial": "#e59866",
    "Endothelial": "#1f618d",
    "Stromal_fibro_myo": "#196f3d",
    "Pericyte": "#82e0aa",
    "Myeloid": "#7d3c98",
    "Lymphoid": "#d4ac0d",
    "Unresolved": "#bfc9ca",
}
SHORT = {
    "Malignant_epithelial": "Malignant epi.",
    "Nonmalignant_epithelial": "Non-malig. epi.",
    "Stromal_fibro_myo": "Fibro/myofibro",
    "Stromal (fibro+mural)": "Stroma (fibro+mural)",
}


def _label(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def _save(fig, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / name
    fig.savefig(p, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return p


# ----------------------------------------------------------------------------- figure 1


def figure1() -> Path:
    fig = plt.figure(figsize=(11, 6.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.45, wspace=0.6)

    ax = fig.add_subplot(gs[0, 0])
    bins = np.logspace(0.5, 4.6, 70)
    for c in config.COHORTS:
        q = pd.read_csv(
            config.CHECKPOINTS_DIR / f"qc_cells_{c}.tsv.gz",
            sep="\t",
            usecols=["total_counts", "pct_mt"],
        )
        nonmt = q["total_counts"] * (1 - q["pct_mt"] / 100)
        ax.hist(
            nonmt,
            bins=bins,
            histtype="step",
            lw=1.6,
            density=True,
            color=COH[c],
            label=f"{c} (max {nonmt.max():,.0f})",
        )
    ax.set_xscale("log")
    ax.set(xlabel="non-mito counts / cell", ylabel="density", yticks=[])
    ax.set_title("QC: cohort C is capped", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    _label(ax, "a")

    ab = pd.read_csv(
        config.CHECKPOINTS_DIR / "stage3" / "embedding_obs_AB.tsv.gz", sep="\t", index_col=0
    )
    rng = np.random.default_rng(0)
    ab = ab.iloc[rng.permutation(len(ab))]
    panels = [
        ("AB_baseline", "cohort", COH, "A+B unintegrated (cohort)", "b"),
        ("AB_harmony", "cohort", COH, "Harmony, batch = cohort", "c"),
        ("AB_harmony", "population", POP_COLORS, "Harmony (population)", "d"),
    ]
    for i, (tag, col, pal, title, lab) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i + 1])
        ax.scatter(
            ab[f"{tag}_umap1"],
            ab[f"{tag}_umap2"],
            c=ab[col].map(pal).fillna("#ddd"),
            s=0.15,
            linewidths=0,
            rasterized=True,
        )
        ax.set(xticks=[], yticks=[])
        ax.set_title(title, fontsize=9)
        _label(ax, lab)
    h = [
        plt.Line2D([], [], marker="o", ls="", color=v, label=SHORT.get(k, k), ms=4)
        for k, v in POP_COLORS.items()
    ]
    ax.legend(handles=h, fontsize=6, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)

    m = pd.read_csv(T / "stage3_integration_metrics.tsv", sep="\t").set_index("embedding")
    tr = pd.read_csv(T / "stage3_label_transfer.tsv", sep="\t")
    cond = pd.read_csv(T / "stage3_condition_structure.tsv", sep="\t")
    names = {
        "AB_baseline": "unintegrated",
        "AB_harmony": "Harmony (cohort)",
        "AB_harmony_spec": "Harmony (specimen)",
    }
    rows = {}
    for e in names:
        st = cond[(cond["embedding"] == e) & (cond["population"] == "Stromal (fibro+mural)")]
        rows[names[e]] = {
            "cohort mixing\n(batch ASW)": m.loc[e, "batch_ASW"],
            "cell-type\nseparation (ASW)": m.loc[e, "celltype_ASW"],
            "label transfer\nA↔B (bal. acc.)": tr[tr["embedding"] == e]["balanced_accuracy"].mean(),
            "tumor/adjacent\nstructure, stroma": st["condition_ASW"].mean(),
        }
    d = pd.DataFrame(rows).T
    ax = fig.add_subplot(gs[1, :3])
    d.T.plot.bar(ax=ax, rot=0, color=["#95a5a6", "#27ae60", "#c0392b"], width=0.75)
    for cont in ax.containers:
        ax.bar_label(cont, fmt="%.2f", fontsize=7, padding=1)
    ax.set(ylim=(0, 1.15), ylabel="score")
    ax.set_title(
        "Integration evaluation (A-alone tumor/adjacent structure in stroma = 0.26)", fontsize=9
    )
    ax.legend(fontsize=7, frameon=False, ncol=3, loc="upper left")
    plt.setp(ax.get_xticklabels(), fontsize=8)
    _label(ax, "e")

    ax = fig.add_subplot(gs[1, 3])
    cov = pd.read_csv(T / "stage4_focus_coverage.tsv", sep="\t")
    cov = cov[
        cov["population"].isin(
            ["Stromal (fibro+mural)", "Stromal_fibro_myo", "Pericyte", "Endothelial"]
        )
    ]
    w = cov.pivot_table(
        index="population", columns="cohort", values="paired_donors_ge10", fill_value=0
    )
    w = w.reindex(["Stromal (fibro+mural)", "Stromal_fibro_myo", "Pericyte", "Endothelial"])
    w.index = ["stroma\n(fibro+mural)", "fibro/\nmyofibro", "pericyte", "endothelial"]
    w.plot.barh(ax=ax, stacked=True, color=[COH["A"], COH["B"]], width=0.7)
    ax.invert_yaxis()
    ax.set(xlabel="paired donors (≥10 cells)")
    ax.set_title("What combining adds", fontsize=9)
    ax.legend(fontsize=7, frameon=False, title="cohort", title_fontsize=7, loc="lower right")
    plt.setp(ax.get_yticklabels(), fontsize=7)
    _label(ax, "f")
    return _save(fig, "figure1_qc_integration.png")


# ----------------------------------------------------------------------------- figure 2


def figure2() -> Path:
    fig = plt.figure(figsize=(11, 6.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 0.8, 0.8], hspace=0.5, wspace=0.42)

    ax = fig.add_subplot(gs[:, 0])
    dp = pd.read_csv(T / "stage4_ecm_dotplot_population.tsv", sep="\t")
    pops = [
        "Stromal_fibro_myo",
        "Pericyte",
        "Endothelial",
        "Nonmalignant_epithelial",
        "Malignant_epithelial",
        "Myeloid",
        "Lymphoid",
    ]
    rows = [
        (c, p) for c in ("A", "B") for p in pops if ((dp.cohort == c) & (dp.population == p)).any()
    ]
    xi = {g: i for i, g in enumerate(GENES)}
    for y, (c, p) in enumerate(rows):
        d = dp[(dp["cohort"] == c) & (dp["population"] == p)]
        sca = ax.scatter(
            d["gene"].map(xi),
            [y] * len(d),
            s=8 + d["frac"] * 230,
            c=d["mean"],
            cmap="magma_r",
            vmin=0,
            vmax=2.1,
            edgecolor="k",
            linewidth=0.3,
        )
    ax.set_yticks(range(len(rows)), [f"{c}: {SHORT.get(p, p)}" for c, p in rows], fontsize=8)
    ax.set_xticks(range(len(GENES)), GENES, rotation=60, fontsize=8)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.axhline(len([r for r in rows if r[0] == "A"]) - 0.5, color="grey", lw=0.6)
    cb = fig.colorbar(sca, ax=ax, shrink=0.3, pad=0.03, location="right", aspect=12)
    cb.set_label("mean log1p(CP10k)\n(size = fraction expressing)", fontsize=7)
    ax.set_title("ECM panel by population (QC-passing cells)", fontsize=9)
    _label(ax, "a")

    lv = pd.read_csv(T / "stage5_donor_levels.tsv", sep="\t")
    g = lv.groupby(["cohort", "gene"])["logcpm"]
    lv["z"] = (lv["logcpm"] - g.transform("mean")) / g.transform("std")
    sc = lv.groupby(["cohort", "donor_key", "condition"])["z"].mean().unstack("condition")
    eff = pd.read_csv(T / "stage5_effects.tsv", sep="\t")
    prim = eff[eff["analysis"].str.startswith("primary")]
    for i, c in enumerate(["B", "A"]):
        ax = fig.add_subplot(gs[0, 1 + i])
        s = sc.xs(c).dropna()
        for _, r in s.iterrows():
            ax.plot(
                [0, 1],
                [r["Adjacent"], r["Tumor"]],
                "-o",
                ms=3,
                lw=1.2,
                color="#c0392b" if r["Tumor"] > r["Adjacent"] else "#7f8c8d",
            )
        st = prim[(prim["cohort"] == c) & (prim["feature"] == "panel_z")].iloc[0]
        ax.set(xticks=[0, 1], xticklabels=["Adj.", "Tumor"], xlim=(-0.25, 1.25))
        ax.set_title(
            f"Stroma, cohort {c}: {int(st['n_up'])}/{int(st['n_donors'])} up\n"
            f"Δ {st['mean']:+.2f} [{st['ci_low']:+.2f}, {st['ci_high']:+.2f}]",
            fontsize=8,
        )
        if i == 0:
            ax.set_ylabel("panel score (mean z)", fontsize=8)
        _label(ax, "b" if c == "B" else "c")

    ax = fig.add_subplot(gs[1, 1])
    for c, off in [("A", -0.17), ("B", 0.17)]:
        d = (
            prim[(prim["cohort"] == c) & prim["feature"].isin(GENES)]
            .set_index("feature")
            .reindex(GENES)
        )
        yy = np.arange(len(GENES)) + off
        lo = np.clip(d["ci_low"], -6, 6)
        hi = np.clip(d["ci_high"], -6, 6)
        ax.errorbar(
            d["mean"].clip(-6, 6),
            yy,
            xerr=[d["mean"].clip(-6, 6) - lo, hi - d["mean"].clip(-6, 6)],
            fmt="o",
            ms=3,
            color=COH[c],
            capsize=1.5,
            lw=0.8,
            label=c,
        )
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(range(len(GENES)), GENES, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlim=(-6, 6))
    ax.set_xlabel("log2FC tumor − adjacent (95 % CI)", fontsize=7)
    ax.set_title("Per gene (stroma)", fontsize=8)
    ax.legend(
        fontsize=7,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=2,
        title="cohort",
        title_fontsize=7,
    )
    _label(ax, "d")

    ax = fig.add_subplot(gs[1, 2])
    b = eff[(eff["cohort"] == "B") & (eff["feature"].isin(["panel_z", "panel_z_noDCN_LUM"]))]
    items = [
        ("primary", "panel_z", "primary"),
        ("primary", "panel_z_noDCN_LUM", "no DCN/LUM"),
        ("robustness: pericytes only", "panel_z", "pericytes only"),
        ("sensitivity: >=50 cells", "panel_z", "≥50 cells"),
        ("sensitivity: one tumor region", "panel_z", "one region"),
        ("contrast: endothelial", "panel_z", "endothelial"),
    ]
    for y, (an, f, lab) in enumerate(items):
        r = b[b["analysis"].str.startswith(an) & (b["feature"] == f)].iloc[0]
        ax.errorbar(
            r["mean"],
            y,
            xerr=[[r["mean"] - r["ci_low"]], [r["ci_high"] - r["mean"]]],
            fmt="o",
            ms=3,
            color="#7f8c8d" if lab == "endothelial" else COH["B"],
            capsize=2,
        )
    mm = pd.read_csv(T / "stage5_mixed_model.tsv", sep="\t").set_index("cohort").loc["B"]
    y = len(items)
    ax.errorbar(
        mm["estimate_z"],
        y,
        xerr=[[mm["estimate_z"] - mm["ci_low"]], [mm["ci_high"] - mm["estimate_z"]]],
        fmt="s",
        ms=3,
        color=COH["B"],
        capsize=2,
    )
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(range(y + 1), [i[2] for i in items] + ["mixed model"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("panel Δ, cohort B (95 % CI)", fontsize=7)
    ax.set_title("Robustness (composition) +\nsensitivities", fontsize=8)
    _label(ax, "e")
    return _save(fig, "figure2_ecm_program.png")


# ----------------------------------------------------------------------------- figure 3


def figure3() -> Path:
    fig = plt.figure(figsize=(11, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.7, 0.75, 0.75], wspace=0.45)
    a = pd.read_csv(T / "stage6_agreement.tsv", sep="\t")
    cols = {
        "snRNA/scRNA stroma, cohort B (Stage 5)": "B scRNA\nstroma (9)",
        "snRNA/scRNA stroma, cohort A (Stage 5)": "A snRNA\nstroma (3)",
    }
    for ev in a["evidence"].unique():
        if "matched subset" in ev:
            short = ev.split(" — ")[0].replace("snRNA whole specimen", "snRNA whole")
            cols[ev] = f"{short}\n(matched 3)"
        elif ev.startswith("bulk RNA (all"):
            cols[ev] = "bulk RNA\n(72 pairs)"
        elif ev.startswith("bulk protein (all"):
            cols[ev] = "protein\n(80 pairs)"
    a = a[a["evidence"].isin(cols)].assign(evidence=lambda d: d["evidence"].map(cols))
    piv = a.pivot_table(index="feature", columns="evidence", values="mean").reindex(
        GENES + ["panel_z"]
    )[list(dict.fromkeys(cols.values()))]
    ax = fig.add_subplot(gs[0])
    sns.heatmap(
        piv,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        annot=True,
        fmt="+.1f",
        annot_kws={"fontsize": 6.5},
        ax=ax,
        cbar_kws={"shrink": 0.7, "label": "T − Adj"},
    )
    ax.axhline(len(GENES), color="k", lw=1.2)
    for x in (2, 6):
        ax.axvline(x, color="k", lw=1.2)
    ax.set(xlabel="", ylabel="")
    ax.set_title("Tumor − adjacent within each assay (log2; panel in SD units)", fontsize=9)
    ax.set_xticklabels(
        [t.get_text().replace("\n", " ") for t in ax.get_xticklabels()],
        fontsize=7,
        rotation=40,
        ha="right",
    )
    plt.setp(ax.get_yticklabels(), fontsize=8)
    _label(ax, "a")

    ax = fig.add_subplot(gs[1])
    s = pd.read_csv(T / "stage6_composition_scores.tsv", sep="\t")
    corr = pd.read_csv(T / "stage6_composition_corr.tsv", sep="\t").set_index("marker_set")
    ax.scatter(s["fibroblast_markers_z"], s["ECM_noDCN_LUM_z"], s=10, alpha=0.7, color="#34495e")
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("fibroblast marker change", fontsize=8)
    ax.set_ylabel("ECM change, 8 genes (no DCN/LUM)", fontsize=8)
    ax.set_title(
        f"Bulk RNA composition (72 pairs)\nρ = "
        f"{corr.loc['fibroblast', 'spearman_rho_noDCN_LUM']:.2f}; fibroblast markers "
        f"down in {corr.loc['fibroblast', 'frac_donors_marker_down']:.0%}",
        fontsize=8,
    )
    _label(ax, "b")

    ax = fig.add_subplot(gs[2])
    pdn = pd.read_csv(T / "stage6_matched_per_donor.tsv", sep="\t")
    pdn = pdn[~pdn["assay"].str.contains(r"\[")]
    assays = ["snRNA stroma", "snRNA whole specimen", "bulk RNA", "bulk protein"]
    for donor, d in pdn.groupby("donor_id"):
        d = d.set_index("assay").reindex(assays)
        ax.plot(range(4), d["panel_z"], "-o", ms=3, label=donor)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(4), ["sn\nstroma", "sn\nwhole", "bulk\nRNA", "protein"], fontsize=7)
    ax.set_ylabel("panel Δ (T − Adj)", fontsize=8)
    ax.set_title("Matched donors (n = 3)", fontsize=8)
    ax.legend(fontsize=6, frameon=False)
    _label(ax, "c")
    return _save(fig, "figure3_multimodal.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="paper")
    return [figure1(), figure2(), figure3()]
