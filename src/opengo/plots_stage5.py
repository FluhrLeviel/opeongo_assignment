"""Stage 5 figures -> outputs/step_5_output/. Run with `uv run opengo figures 5`."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from opengo import config  # noqa: E402

OUT_DIR = config.OUTPUTS_DIR / "step_5_output"
T = config.TABLES_DIR
GENES = config.ECM_GENES
COHORT = {"A": "#1f77b4", "B": "#ff7f0e", "A+B": "#2c3e50"}


def _save(fig, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / name
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return p


def _effects() -> pd.DataFrame:
    return pd.read_csv(T / "stage5_effects.tsv", sep="\t")


def fig_paired() -> Path:
    lv = pd.read_csv(T / "stage5_donor_levels.tsv", sep="\t")
    grp = lv.groupby(["cohort", "gene"])["logcpm"]
    lv["z"] = (lv["logcpm"] - grp.transform("mean")) / grp.transform("std")
    score = lv.groupby(["cohort", "donor_key", "condition"])["z"].mean().unstack("condition")
    e = _effects()
    prim = e[e["analysis"].str.startswith("primary") & (e["feature"] == "panel_z")]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), gridspec_kw={"width_ratios": [1, 1, 1.3]})
    for ax, cohort in zip(axes[:2], config.ANALYSIS_COHORTS, strict=True):
        s = score.xs(cohort).dropna()
        for donor, r in s.iterrows():
            up = r["Tumor"] > r["Adjacent"]
            ax.plot(
                [0, 1],
                [r["Adjacent"], r["Tumor"]],
                "-o",
                color="#c0392b" if up else "#7f8c8d",
                lw=1.5,
                ms=5,
                alpha=0.85,
            )
            ax.annotate(
                donor.split(":")[1],
                (1, r["Tumor"]),
                xytext=(5, 0),
                textcoords="offset points",
                fontsize=7,
                va="center",
            )
        st = prim[prim["cohort"] == cohort].iloc[0]
        ax.set(
            xticks=[0, 1],
            xticklabels=["Adjacent", "Tumor"],
            xlim=(-0.3, 1.5),
            ylabel="ECM panel score (mean z, donor level)",
            title=f"Cohort {cohort}: {int(st['n_up'])}/{int(st['n_donors'])} donors up\n"
            f"Δ = {st['mean']:+.2f} [{st['ci_low']:+.2f}, {st['ci_high']:+.2f}], "
            f"sign-flip p = {st['p_signflip']:.3f}",
        )
    ax = axes[2]
    rows = e[(e["feature"] == "panel_z") & e["cohort"].isin(["A", "B"])]
    rows = rows[rows["analysis"].str.startswith(("primary", "robustness", "contrast"))]
    labels, y = [], 0
    for _, r in rows.iterrows():
        ax.errorbar(
            r["mean"],
            y,
            xerr=[[r["mean"] - r["ci_low"]], [r["ci_high"] - r["mean"]]]
            if pd.notna(r["ci_low"])
            else None,
            fmt="o",
            color=COHORT[r["cohort"]],
            capsize=3,
        )
        labels.append(
            f"{r['cohort']}: {r['analysis'].split(':')[1].split(',')[0].strip()} "
            f"(n={int(r['n_donors'])})"
        )
        y += 1
    meta = e[
        (e["cohort"] == "A+B")
        & (e["feature"] == "panel_z")
        & e["analysis"].str.startswith("primary")
    ].iloc[0]
    ax.errorbar(
        meta["random"],
        y,
        xerr=[[meta["random"] - meta["random_low"]], [meta["random_high"] - meta["random"]]],
        fmt="D",
        color=COHORT["A+B"],
        capsize=3,
    )
    labels.append(f"A+B random-effects (I² = {meta['I2']:.0%})")
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.invert_yaxis()
    ax.set(xlabel="tumor − adjacent, panel score (95 % CI)", title="Estimates by cohort")
    fig.suptitle(
        "Stromal (fibroblast + mural) ECM panel: paired tumor vs adjacent, donors as replicates",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "01_paired_panel_score.png")


def fig_gene_forest() -> Path:
    e = _effects()
    p = e[e["analysis"].str.startswith("primary")]
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    off = {"A": -0.18, "B": 0.18}
    for cohort in config.ANALYSIS_COHORTS:
        d = (
            p[(p["cohort"] == cohort) & p["feature"].isin(GENES)]
            .set_index("feature")
            .reindex(GENES)
        )
        yy = np.arange(len(GENES)) + off[cohort]
        ax.errorbar(
            d["mean"],
            yy,
            xerr=[d["mean"] - d["ci_low"], d["ci_high"] - d["mean"]],
            fmt="o",
            color=COHORT[cohort],
            capsize=2,
            label=f"cohort {cohort}",
        )
        for yv, (_, r) in zip(yy, d.iterrows(), strict=True):
            if pd.notna(r["q_BH"]) and r["q_BH"] < 0.05:
                ax.text(
                    r["ci_high"] + 0.2,
                    yv,
                    f"q={r['q_BH']:.2f}",
                    fontsize=7,
                    va="center",
                    color=COHORT[cohort],
                )
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(range(len(GENES)), GENES)
    ax.invert_yaxis()
    ax.set_xlim(-8, 8)
    ax.set(
        xlabel="log2 fold change, tumor − adjacent (95 % CI; A intervals clipped)",
        title="Per gene (BH q < 0.05 annotated)",
    )
    ax.legend()
    ax = axes[1]
    pd_ = pd.read_csv(T / "stage5_per_donor_log2FC.tsv", sep="\t")
    pd_["donor"] = pd_["cohort"] + ": " + pd_["donor_key"].str.split(":").str[1]
    sns.heatmap(
        pd_.set_index("donor")[GENES].T,
        cmap="RdBu_r",
        center=0,
        vmin=-5,
        vmax=5,
        annot=True,
        fmt=".1f",
        annot_kws={"fontsize": 6},
        ax=ax,
        cbar_kws={"label": "log2FC tumor − adjacent"},
    )
    ax.set(title="Per donor × gene (donor variability)", xlabel="", ylabel="")
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
    fig.suptitle(
        "Individual genes: fibrillar collagens, FN1, SPARC rise in B; DCN/LUM do not",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "02_gene_effects.png")


def fig_robustness() -> Path:
    e = _effects()
    s = e[(e["feature"] == "panel_z") & e["cohort"].isin(["A", "B"])].copy()
    mm = pd.read_csv(T / "stage5_mixed_model.tsv", sep="\t")
    loo = pd.read_csv(T / "stage5_leave_one_out.tsv", sep="\t")
    fs = pd.read_csv(T / "stage5_B_fibroblast_share.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), gridspec_kw={"width_ratios": [1.5, 0.8, 0.8]})
    ax = axes[0]
    labels, ypos, y = [], [], 0.0

    def row(est, lo, hi, label, color, marker="o"):
        nonlocal y
        err = [[est - lo], [hi - est]] if pd.notna(lo) else None
        ax.errorbar(est, y, xerr=err, fmt=marker, color=color, capsize=3)
        labels.append(label)
        ypos.append(y)
        y += 1

    for cohort in ["B", "A"]:
        for _, r in s[(s["cohort"] == cohort) & s["mean"].notna()].iterrows():
            row(
                r["mean"],
                r["ci_low"],
                r["ci_high"],
                f"{cohort}: {r['analysis'].split(',')[0]} (n={int(r['n_donors'])})",
                COHORT[cohort],
            )
        m = mm[mm["cohort"] == cohort].iloc[0]
        row(
            m["estimate_z"],
            m["ci_low"],
            m["ci_high"],
            f"{cohort}: mixed model, all {int(m['donors'])} donors / {int(m['specimens'])} "
            "specimens",
            COHORT[cohort],
            marker="s",
        )
        lo = loo[loo["cohort"] == cohort]["panel_z"]
        ax.plot([lo.min(), lo.max()], [y, y], "-", color=COHORT[cohort], lw=4, alpha=0.5)
        labels.append(f"{cohort}: leave-one-donor-out range")
        ypos.append(y)
        y += 1.5
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(ypos, labels, fontsize=8)
    ax.invert_yaxis()
    ax.set(
        xlabel="tumor − adjacent, panel score (95 % CI)",
        title="Robustness and sensitivity of the main estimate",
    )
    ax = axes[1]
    w = fs.pivot(index="donor_key", columns="condition", values="fib_share").dropna()
    for _, r in w.iterrows():
        ax.plot([0, 1], [r["Adjacent"], r["Tumor"]], "-o", color="#7f8c8d", ms=4)
    ax.set(
        xticks=[0, 1],
        xticklabels=["Adjacent", "Tumor"],
        xlim=(-0.3, 1.3),
        ylabel="fibroblast share of compartment",
        title="B: fibroblast:mural mix\n(does composition explain it?)",
    )
    ax = axes[2]
    pc = e[
        (e["feature"] == "panel_z")
        & e["analysis"].str.startswith(("robustness", "primary"))
        & (e["cohort"] == "B")
    ]
    ax.bar(
        ["fibro + mural", "pericytes only", "8 genes\n(no DCN/LUM)"],
        [
            pc[pc["analysis"].str.startswith("primary")]["mean"].iloc[0],
            pc[pc["analysis"].str.startswith("robustness")]["mean"].iloc[0],
            e[
                (e["feature"] == "panel_z_noDCN_LUM")
                & (e["cohort"] == "B")
                & e["analysis"].str.startswith("primary")
            ]["mean"].iloc[0],
        ],
        color=["#ff7f0e", "#f5b041", "#fad7a0"],
    )
    ax.set(ylabel="tumor − adjacent (panel score)", title="B: effect without fibroblast\ncontent")
    fig.suptitle(
        "Robustness check targeting the main uncertainty (compartment composition) + sensitivities",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "03_robustness.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="notebook")
    return [fig_paired(), fig_gene_forest(), fig_robustness()]
