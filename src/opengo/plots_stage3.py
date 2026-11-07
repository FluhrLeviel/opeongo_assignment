"""Stage 3 figures -> outputs/step_3_output/. Run with `uv run opengo figures 3`."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from opengo import config  # noqa: E402
from opengo.integrate import STAGE3_DIR  # noqa: E402

OUT_DIR = config.OUTPUTS_DIR / "step_3_output"
T = config.TABLES_DIR
COND = {"Tumor": "#c0392b", "Adjacent": "#2471a3"}
COHORT = {"A": "#1f77b4", "B": "#ff7f0e"}
POP_ORDER = [
    "Malignant_epithelial",
    "Nonmalignant_epithelial",
    "Endothelial",
    "Stromal_fibro_myo",
    "Pericyte",
    "Myeloid",
    "Lymphoid",
    "Unresolved",
]
POP_COLORS = dict(
    zip(
        POP_ORDER,
        ["#7b241c", "#e59866", "#1f618d", "#196f3d", "#82e0aa", "#7d3c98", "#d4ac0d", "#bfc9ca"],
        strict=True,
    )
)
TITLES = {
    "A_baseline": "A alone, unintegrated",
    "AB_baseline": "A + B, unintegrated",
    "AB_harmony": "A + B, Harmony (batch = cohort)",
    "AB_harmony_spec": "A + B, Harmony (batch = specimen)",
}


def _save(fig, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / name
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return p


def _obs(key: str) -> pd.DataFrame:
    return pd.read_csv(STAGE3_DIR / f"embedding_obs_{key}.tsv.gz", sep="\t", index_col=0)


def _scatter(ax, df, tag, col, palette, s=0.4):
    rng = np.random.default_rng(0)
    df = df.iloc[rng.permutation(len(df))]
    ax.scatter(
        df[f"{tag}_umap1"],
        df[f"{tag}_umap2"],
        c=df[col].map(palette).fillna("#dddddd"),
        s=s,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xticks([])
    ax.set_yticks([])


def fig_umaps() -> Path:
    a, ab = _obs("A"), _obs("AB")
    rows = [("A_baseline", a), ("AB_baseline", ab), ("AB_harmony", ab), ("AB_harmony_spec", ab)]
    cols = [("cohort", COHORT, "cohort"), ("population", POP_COLORS, "population")]
    cols += [("condition", COND, "condition")]
    fig, axes = plt.subplots(4, 3, figsize=(16, 21))
    for r, (tag, df) in enumerate(rows):
        for c, (col, pal, label) in enumerate(cols):
            ax = axes[r, c]
            _scatter(ax, df, tag, col, pal)
            if r == 0:
                ax.set_title(f"colour = {label}", fontsize=13)
            if c == 0:
                ax.set_ylabel(TITLES[tag], fontsize=12)
    for c, (_, pal, _) in enumerate(cols):
        handles = [plt.Line2D([], [], marker="o", ls="", color=v, label=k) for k, v in pal.items()]
        axes[3, c].legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.03),
            ncol=2,
            fontsize=9,
            markerscale=1.5,
        )
    fig.suptitle(
        "Stage 3 embeddings (subsample: ≤1,000 cells/specimen + all stromal/endothelial)",
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout()
    return _save(fig, "01_umap_baseline_vs_harmony.png")


def fig_metrics() -> Path:
    m = pd.read_csv(T / "stage3_integration_metrics.tsv", sep="\t").set_index("embedding")
    tr = pd.read_csv(T / "stage3_label_transfer.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), gridspec_kw={"width_ratios": [1, 1.4, 1]})
    batch = m.loc[[i for i in m.index if i.startswith("AB")], ["iLISI_scaled", "batch_ASW"]]
    batch.rename(index=TITLES).plot.barh(ax=axes[0], color=["#7f8c8d", "#bdc3c7"])
    axes[0].set(xlim=(0, 1), title="Batch mixing (higher = cohorts more mixed)", xlabel="score")
    axes[0].invert_yaxis()
    bio = m[["cLISI_scaled", "celltype_ASW", "ARI", "NMI"]].rename(index=TITLES)
    bio.plot.barh(ax=axes[1], color=sns.color_palette("Greens", 4))
    axes[1].set(xlim=(0, 1), title="Cell identity conservation (higher = better)", xlabel="score")
    axes[1].invert_yaxis()
    piv = tr.pivot(index="embedding", columns="direction", values="balanced_accuracy")
    piv.rename(index=TITLES).plot.barh(ax=axes[2], color=["#2e86c1", "#e67e22"])
    axes[2].set(
        xlim=(0, 1),
        title="Cross-cohort kNN label transfer\n(balanced accuracy)",
        xlabel="balanced accuracy",
    )
    axes[2].invert_yaxis()
    for ax in axes:
        for cont in ax.containers:
            ax.bar_label(cont, fmt="%.2f", fontsize=7, padding=2)
    fig.suptitle("Integration evaluation (scIB-style metrics)", fontweight="bold")
    fig.tight_layout()
    return _save(fig, "02_integration_metrics.png")


def fig_condition_and_purity() -> Path:
    cond = pd.read_csv(T / "stage3_condition_structure.tsv", sep="\t")
    pur = pd.read_csv(T / "stage3_knn_purity.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(19, 6.5), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    cond["group"] = cond["cohort"] + ": " + cond["population"]
    order = [t for t in TITLES if t in set(cond["embedding"])]
    piv = cond.pivot(index="group", columns="embedding", values="condition_ASW")[order]
    piv.rename(columns=TITLES).plot.barh(
        ax=ax, width=0.8, color=["#95a5a6", "#34495e", "#27ae60", "#c0392b"]
    )
    ax.axvline(0, color="k", lw=0.6)
    ax.set(
        xlabel="tumor-vs-adjacent silhouette within (cohort, population)",
        title="Is the tumor/adjacent structure kept?\n(the signal Stage 5 will test)",
    )
    ax.legend(fontsize=8)
    ax = axes[1]
    pur["group"] = pur["cohort"] + ": " + pur["population"]
    cols = [t for t in TITLES if t in pur.columns]
    sns.heatmap(
        pur.set_index("group")[cols].rename(columns=TITLES),
        annot=True,
        fmt=".2f",
        cmap="Greens",
        vmin=0,
        vmax=1,
        ax=ax,
        cbar_kws={"label": "kNN label purity"},
    )
    ax.set(
        title="Are populations kept distinct? (fraction of 15 neighbours, same population)",
        ylabel="",
    )
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    fig.tight_layout()
    return _save(fig, "03_condition_structure_and_purity.png")


def fig_transfer() -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    for r, tag in enumerate(["AB_baseline", "AB_harmony"]):
        for c, d in enumerate(["AtoB", "BtoA"]):
            conf = pd.read_csv(T / f"stage3_transfer_{tag}_{d}.tsv", sep="\t", index_col=0)
            frac = conf.div(conf.sum(axis=1), axis=0)
            sns.heatmap(
                frac,
                annot=conf,
                fmt="d",
                cmap="Blues",
                vmin=0,
                vmax=1,
                cbar=False,
                ax=axes[r, c],
                annot_kws={"fontsize": 7},
            )
            axes[r, c].set(
                title=f"{TITLES[tag]}: labels {d.replace('to', ' → ')}",
                ylabel="true (target cohort author label)",
                xlabel="predicted",
            )
            plt.setp(axes[r, c].get_xticklabels(), rotation=35, ha="right")
    fig.suptitle(
        "Cross-cohort label transfer (stromal = fibroblast/myofibroblast + pericyte)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "04_label_transfer.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="notebook")
    return [fig_umaps(), fig_metrics(), fig_condition_and_purity(), fig_transfer()]
