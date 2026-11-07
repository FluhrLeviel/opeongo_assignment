"""Report figures (one per assignment section) -> <report dir>/figures/. `opengo figures report`.

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

OUT_DIR = config.REPORT_DIR / "figures"
T = config.TABLES_DIR
GENES = config.ECM_GENES
COH = {"A": "#1f77b4", "B": "#ff7f0e", "C": "#2ca02c"}
UP, DOWN = "#c0392b", "#7f8c8d"
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


def _label(ax, s, x=-0.12):
    ax.text(x, 1.06, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def _save(fig, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / name
    fig.savefig(p, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return p


def _ci(ax, est, lo, hi, y, color, marker="o"):
    err = [[est - lo], [hi - est]] if pd.notna(lo) else None
    ax.errorbar(est, y, xerr=err, fmt=marker, ms=3, color=color, capsize=2)


# ----------------------------------------------------------------------------- figure 1


def figure1() -> Path:
    fig = plt.figure(figsize=(11, 6.4))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05], hspace=0.6, wspace=0.55)

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
            label=f"{c}: max {nonmt.max():,.0f}",
        )
    ax.set_xscale("log")
    ax.set(xlabel="non-mito counts / cell", ylabel="density", yticks=[])
    ax.set_title("QC: cohort C is capped", fontsize=9)
    ax.legend(
        fontsize=6.5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.32),
        ncol=3,
        handlelength=1,
        columnspacing=0.8,
    )
    _label(ax, "a")

    ab = pd.read_csv(
        config.CHECKPOINTS_DIR / "stage3" / "embedding_obs_AB.tsv.gz", sep="\t", index_col=0
    )
    ab = ab.iloc[np.random.default_rng(0).permutation(len(ab))]
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
        _label(ax, lab, -0.05)
        if i == 0:
            h = [plt.Line2D([], [], marker="o", ls="", color=COH[k], label=k, ms=4) for k in "AB"]
            ax.legend(handles=h, fontsize=7, frameon=False, loc="lower left")
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
            "cell-type sep.\n(ASW)": m.loc[e, "celltype_ASW"],
            "label transfer\nA↔B": tr[tr["embedding"] == e]["balanced_accuracy"].mean(),
            "tumor/adjacent\nin stroma": st["condition_ASW"].mean(),
        }
    d = pd.DataFrame(rows).T
    ax = fig.add_subplot(gs[1, :2])
    d.T.plot.bar(ax=ax, rot=0, color=["#95a5a6", "#27ae60", "#c0392b"], width=0.78)
    for cont in ax.containers:
        ax.bar_label(cont, fmt="%.2f", fontsize=6.5, padding=1)
    ax.set(ylim=(0, 1.12), ylabel="score")
    ax.set_title("Integration evaluation", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, ncol=3, loc="upper left")
    plt.setp(ax.get_xticklabels(), fontsize=7)
    _label(ax, "e", -0.08)

    ax = fig.add_subplot(gs[1, 2])
    qp = pd.read_csv(T / "stage2_qc_by_population.tsv", sep="\t")
    qp["removed"] = 100 * (1 - qp["pass_frac"])
    order = [
        "Stromal_fibro_myo",
        "Pericyte",
        "Endothelial",
        "Nonmalignant_epithelial",
        "Malignant_epithelial",
        "Myeloid",
        "Lymphoid",
        "Unresolved",
        "Unlabelled",
    ]
    hm = qp.pivot_table(index="population", columns="cohort", values="removed")
    hm = hm.reindex([o for o in order if o in hm.index])
    hm.index = [SHORT.get(i, i) for i in hm.index]
    sns.heatmap(
        hm,
        annot=True,
        fmt=".0f",
        cmap="Reds",
        vmin=0,
        vmax=30,
        cbar=False,
        ax=ax,
        annot_kws={"fontsize": 7},
        linewidths=0.4,
        linecolor="white",
    )
    ax.set(xlabel="cohort", ylabel="")
    ax.set_title("% cells removed by QC", fontsize=9)
    plt.setp(ax.get_yticklabels(), fontsize=7)
    _label(ax, "f", -0.75)

    ax = fig.add_subplot(gs[1, 3])
    cov = pd.read_csv(T / "stage4_focus_coverage.tsv", sep="\t")
    keep = ["Stromal (fibro+mural)", "Stromal_fibro_myo", "Pericyte", "Endothelial"]
    w = (
        cov[cov["population"].isin(keep)]
        .pivot_table(
            index="population", columns="cohort", values="paired_donors_ge10", fill_value=0
        )
        .reindex(keep)
    )
    w.index = ["stroma\n(fibro+mural)", "fibro/\nmyofibro", "pericyte", "endothelial"]
    w.plot.barh(ax=ax, stacked=True, color=[COH["A"], COH["B"]], width=0.7, legend=False)
    ax.invert_yaxis()
    ax.set(xlabel="paired donors (≥10 cells)")
    ax.set_title("What combining adds", fontsize=9)
    h = [plt.Rectangle((0, 0), 1, 1, color=COH[k], label=f"cohort {k}") for k in "AB"]
    ax.legend(
        handles=h,
        fontsize=6.5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.32),
        ncol=2,
    )
    plt.setp(ax.get_yticklabels(), fontsize=7)
    _label(ax, "g", -0.55)
    return _save(fig, "figure1_qc_integration.png")


# ----------------------------------------------------------------------------- figure 2


def figure2() -> Path:
    fig = plt.figure(figsize=(11, 7.0))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.45, 0.8, 0.8, 0.8], hspace=0.6, wspace=0.55)

    ax = fig.add_subplot(gs[:, 0])
    dp = pd.read_csv(T / "stage4_ecm_dotplot_population.tsv", sep="\t")
    fine = pd.read_csv(T / "stage4_ecm_dotplot_stromal_fine.tsv", sep="\t")
    blocks = []
    other = [
        "Endothelial",
        "Nonmalignant_epithelial",
        "Malignant_epithelial",
        "Myeloid",
        "Lymphoid",
    ]
    for c, stromal in [
        ("A", ["Fibroblasts", "Myofibroblasts"]),
        ("B", ["Fibroblast", "Pericytes"]),
    ]:
        for grp in stromal:
            blocks.append((f"{c}: {grp}", fine[(fine["cohort"] == c) & (fine["group"] == grp)]))
        for p in other:
            blocks.append(
                (f"{c}: {SHORT.get(p, p)}", dp[(dp["cohort"] == c) & (dp["population"] == p)])
            )
    xi = {g: i for i, g in enumerate(GENES)}
    for y, (_, d) in enumerate(blocks):
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
    ax.set_yticks(range(len(blocks)), [b[0] for b in blocks], fontsize=7.5)
    ax.set_xticks(range(len(GENES)), GENES, rotation=60, fontsize=7.5)
    ax.set_ylim(len(blocks) - 0.5, -0.6)
    ax.axhline(len(other) + 1.5, color="grey", lw=0.6)
    cb = fig.colorbar(sca, ax=ax, shrink=0.25, pad=0.02, location="right", aspect=10)
    cb.set_label("mean log1p(CP10k)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    for f in (0.1, 0.5, 0.9):
        ax.scatter([], [], s=8 + f * 230, c="grey", label=f"{f:.0%}")
    ax.legend(
        title="cells expressing",
        title_fontsize=6.5,
        fontsize=6.5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=3,
    )
    ax.set_title("ECM panel by population (QC-passing cells)", fontsize=9)
    _label(ax, "a", -0.4)

    lv = pd.read_csv(T / "stage5_donor_levels.tsv", sep="\t")
    g = lv.groupby(["cohort", "gene"])["logcpm"]
    lv["z"] = (lv["logcpm"] - g.transform("mean")) / g.transform("std")
    sc = lv.groupby(["cohort", "donor_key", "condition"])["z"].mean().unstack("condition")
    eff = pd.read_csv(T / "stage5_effects.tsv", sep="\t")
    prim = eff[eff["analysis"].str.startswith("primary")]
    for i, c in enumerate(["B", "A"]):
        ax = fig.add_subplot(gs[0, 1 + i])
        for _, r in sc.xs(c).dropna().iterrows():
            ax.plot(
                [0, 1],
                [r["Adjacent"], r["Tumor"]],
                "-o",
                ms=3,
                lw=1.2,
                color=UP if r["Tumor"] > r["Adjacent"] else DOWN,
            )
        st = prim[(prim["cohort"] == c) & (prim["feature"] == "panel_z")].iloc[0]
        ax.set(xticks=[0, 1], xticklabels=["Adj.", "Tumor"], xlim=(-0.25, 1.25))
        ax.set_ylabel("panel score (mean z)", fontsize=7)
        ax.set_title(
            f"Stroma, cohort {c}: {int(st['n_up'])}/{int(st['n_donors'])} up\n"
            f"Δ {st['mean']:+.2f} [{st['ci_low']:+.2f}, {st['ci_high']:+.2f}]",
            fontsize=8,
        )
        _label(ax, "b" if c == "B" else "c")
    ax.legend(
        handles=[
            plt.Line2D([], [], color=UP, label="tumor > adjacent"),
            plt.Line2D([], [], color=DOWN, label="tumor < adjacent"),
        ],
        fontsize=6.5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.17),
        ncol=1,
    )

    ax = fig.add_subplot(gs[0, 3])
    fs = pd.read_csv(T / "stage5_B_fibroblast_share.tsv", sep="\t")
    w = fs.pivot(index="donor_key", columns="condition", values="fib_share").dropna()
    for _, r in w.iterrows():
        ax.plot(
            [0, 1],
            [r["Adjacent"], r["Tumor"]],
            "-o",
            ms=3,
            lw=1.1,
            color=UP if r["Tumor"] > r["Adjacent"] else DOWN,
        )
    t = pd.read_csv(T / "stage5_B_fibroblast_share_test.tsv", sep="\t").iloc[0]
    ax.set(xticks=[0, 1], xticklabels=["Adj.", "Tumor"], xlim=(-0.25, 1.25))
    ax.set_ylabel("fibroblast share of compartment", fontsize=7)
    ax.set_title(
        f"B: fibroblast share, {int(t['n_up'])}/{int(t['n_donors'])} up\n"
        f"Δ {t['mean']:+.2f} [{t['ci_low']:+.2f}, {t['ci_high']:+.2f}]",
        fontsize=8,
    )
    _label(ax, "f")

    ax = fig.add_subplot(gs[1, 1])
    for c, off in [("A", -0.17), ("B", 0.17)]:
        d = prim[(prim["cohort"] == c) & prim["feature"].isin(GENES)].set_index("feature")
        d = d.reindex(GENES)
        yy = np.arange(len(GENES)) + off
        mean = d["mean"].clip(-6, 6)
        lo, hi = np.clip(d["ci_low"], -6, 6), np.clip(d["ci_high"], -6, 6)
        ax.errorbar(
            mean,
            yy,
            xerr=[mean - lo, hi - mean],
            fmt="o",
            ms=3,
            color=COH[c],
            capsize=1.5,
            lw=0.8,
            label=f"cohort {c}",
        )
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(range(len(GENES)), GENES, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlim=(-6, 6))
    ax.set_xlabel("log2FC tumor − adjacent (95 % CI)", fontsize=7)
    ax.set_title("Per gene (stroma)", fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    _label(ax, "d")

    ax = fig.add_subplot(gs[1, 2:])
    b = eff[eff["cohort"] == "B"]
    labels = []

    def row(analysis, feature, label, color=COH["B"]):
        r = b[b["analysis"].str.startswith(analysis) & (b["feature"] == feature)].iloc[0]
        _ci(ax, r["mean"], r["ci_low"], r["ci_high"], len(labels), color)
        labels.append(label)

    row("primary", "panel_z", "primary (fibro + mural)")
    row("primary", "panel_z_noDCN_LUM", "without DCN/LUM")
    row("robustness: pericytes only", "panel_z", "pericytes only")
    row("sensitivity: >=50 cells", "panel_z", "≥50 cells / specimen")
    row("sensitivity: one tumor region", "panel_z", "one tumor region / donor")
    mm = pd.read_csv(T / "stage5_mixed_model.tsv", sep="\t").set_index("cohort").loc["B"]
    _ci(ax, mm["estimate_z"], mm["ci_low"], mm["ci_high"], len(labels), COH["B"], "s")
    labels.append("mixed model, all specimens")
    lo = pd.read_csv(T / "stage5_leave_one_out.tsv", sep="\t")
    lo = lo[lo["cohort"] == "B"]["panel_z"]
    ax.plot([lo.min(), lo.max()], [len(labels)] * 2, "-", color=COH["B"], lw=4, alpha=0.6)
    labels.append("leave-one-donor-out range")
    row("contrast: endothelial", "panel_z", "contrast: endothelial cells", DOWN)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.yaxis.tick_right()
    ax.invert_yaxis()
    ax.set_xlabel("panel Δ tumor − adjacent, cohort B (95 % CI)", fontsize=7)
    ax.set_title("Robustness (composition) and sensitivities", fontsize=8)
    _label(ax, "e", -0.05)
    return _save(fig, "figure2_ecm_program.png")


# ----------------------------------------------------------------------------- figure 3


def figure3() -> Path:
    fig = plt.figure(figsize=(11, 4.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.75, 0.62, 0.62, 0.7], wspace=0.55)
    a = pd.read_csv(T / "stage6_agreement.tsv", sep="\t")
    cols = {"snRNA/scRNA stroma, cohort B (Stage 5)": "B scRNA stroma (9)"}
    for ev in a["evidence"].unique():
        if ev.startswith("bulk RNA (all"):
            cols[ev] = "bulk RNA (72)"
        elif ev.startswith("bulk protein (all"):
            cols[ev] = "protein (80)"
        elif "matched subset" in ev:
            short = ev.split(" — ")[0].replace("snRNA whole specimen", "snRNA whole")
            cols[ev] = short.replace("bulk protein", "protein")
    order = [
        "B scRNA stroma (9)",
        "snRNA stroma",
        "snRNA whole",
        "bulk RNA",
        "protein",
        "bulk RNA (72)",
        "protein (80)",
    ]
    a = a[a["evidence"].isin(cols)].assign(evidence=lambda d: d["evidence"].map(cols))
    feats = GENES + ["panel_z"]
    piv = a.pivot_table(index="feature", columns="evidence", values="mean").reindex(feats)[order]
    qv = a.pivot_table(index="feature", columns="evidence", values="q_BH").reindex(feats)
    qv = qv.reindex(columns=order)
    annot = piv.map(lambda v: f"{v:+.1f}") + qv.map(lambda q: "•" if q < 0.05 else "")
    ax = fig.add_subplot(gs[0])
    sns.heatmap(
        piv,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        annot=annot.to_numpy(),
        fmt="",
        annot_kws={"fontsize": 6.5},
        ax=ax,
        cbar_kws={"shrink": 0.6, "label": "tumor − adjacent", "pad": 0.02},
    )
    ax.axhline(len(GENES), color="k", lw=1.2)
    for x in (1, 5):
        ax.axvline(x, color="k", lw=1.2)
    for x0, x1, lab in [
        (0, 1, "single-\ncell"),
        (1, 5, "matched cohort A\ndonors (n = 3)"),
        (5, 7, "whole CPTAC\n(context)"),
    ]:
        ax.text((x0 + x1) / 2, -0.35, lab, ha="center", va="bottom", fontsize=7)
    ax.set(xlabel="", ylabel="")
    ax.set_xticklabels(order, fontsize=7, rotation=35, ha="right")
    plt.setp(ax.get_yticklabels(), fontsize=7.5)
    _label(ax, "a", -0.12)

    s = pd.read_csv(T / "stage6_composition_scores.tsv", sep="\t")
    corr = pd.read_csv(T / "stage6_composition_corr.tsv", sep="\t").set_index("marker_set")
    for i, (k, lab) in enumerate([("fibroblast", "b"), ("mural", "c")]):
        ax = fig.add_subplot(gs[1 + i])
        ax.scatter(s[f"{k}_markers_z"], s["ECM_noDCN_LUM_z"], s=9, alpha=0.7, color="#34495e")
        ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xlabel(f"{k} marker change", fontsize=7.5)
        if i == 0:
            ax.set_ylabel("ECM change, 8 genes", fontsize=7.5)
        ax.set_title(
            f"bulk RNA (72): {k}\nρ = {corr.loc[k, 'spearman_rho_noDCN_LUM']:.2f}; "
            f"marker Δ {corr.loc[k, 'mean_marker_change_z']:+.2f}",
            fontsize=7.5,
        )
        ax.tick_params(labelsize=7)
        _label(ax, lab, -0.25)

    ax = fig.add_subplot(gs[3])
    pdn = pd.read_csv(T / "stage6_matched_per_donor.tsv", sep="\t")
    pdn = pdn[~pdn["assay"].str.contains(r"\[")]
    assays = ["snRNA stroma", "snRNA whole specimen", "bulk RNA", "bulk protein"]
    donors = sorted(pdn["donor_id"].unique())
    for j, (donor, color) in enumerate(
        zip(donors, ["#1f77b4", "#d62728", "#2ca02c"], strict=False)
    ):
        d = pdn[pdn["donor_id"] == donor].set_index("assay").reindex(assays)
        ax.scatter(np.arange(4) + (j - 1) * 0.18, d["panel_z"], s=18, color=color, label=donor)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(4), ["sn\nstroma", "sn\nwhole", "bulk\nRNA", "protein"], fontsize=7)
    ax.set_ylabel("panel Δ (T − Adj)", fontsize=7.5)
    ax.set_title("Matched donors (n = 3)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    _label(ax, "d", -0.25)
    return _save(fig, "figure3_multimodal.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="paper")
    return [figure1(), figure2(), figure3()]
