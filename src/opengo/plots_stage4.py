"""Stage 4 figures -> outputs/step_4_output/. Run with `uv run opengo figures 4`."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from opengo import config  # noqa: E402
from opengo.stage4 import COMBINED_STROMAL  # noqa: E402

OUT_DIR = config.OUTPUTS_DIR / "step_4_output"
T = config.TABLES_DIR
GENES = config.ECM_GENES
POP_ORDER = [
    COMBINED_STROMAL,
    "Stromal_fibro_myo",
    "Pericyte",
    "Endothelial",
    "Nonmalignant_epithelial",
    "Malignant_epithelial",
    "Myeloid",
    "Lymphoid",
    "Unresolved",
]


def _save(fig, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / name
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return p


def _dot(ax, d: pd.DataFrame, row_key: str, rows: list[str], vmax: float):
    xi = {g: i for i, g in enumerate(GENES)}
    yi = {r: i for i, r in enumerate(rows)}
    d = d[d[row_key].isin(rows)]
    sca = ax.scatter(
        d["gene"].map(xi),
        d[row_key].map(yi),
        s=20 + d["frac"] * 600,
        c=d["mean"],
        cmap="magma_r",
        vmin=0,
        vmax=vmax,
        edgecolor="k",
        linewidth=0.3,
    )
    ax.set_xticks(range(len(GENES)), GENES, rotation=45, ha="right")
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xlim(-0.6, len(GENES) - 0.4)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    return sca


def fig_dotplot() -> Path:
    pop = pd.read_csv(T / "stage4_ecm_dotplot_population.tsv", sep="\t")
    fine = pd.read_csv(T / "stage4_ecm_dotplot_stromal_fine.tsv", sep="\t")
    vmax = max(pop["mean"].max(), fine["mean"].max())
    fig, axes = plt.subplots(
        1, 3, figsize=(24, 6.5), gridspec_kw={"width_ratios": [1, 1, 0.8], "wspace": 0.75}
    )
    for ax, cohort in zip(axes[:2], config.ANALYSIS_COHORTS, strict=True):
        d = pop[pop["cohort"] == cohort]
        n = d.drop_duplicates("population").set_index("population")["n"]
        rows = [p for p in POP_ORDER if p in set(d["population"])]
        d = d.assign(population=d["population"])
        sca = _dot(ax, d, "population", rows, vmax)
        ax.set_yticklabels([f"{r} (n={n[r]:,})" for r in rows])
        ax.set_title(f"Cohort {cohort} ({'snRNA' if cohort == 'A' else 'scRNA'}): all populations")
    ax = axes[2]
    fine = fine.assign(row=fine["cohort"] + ": " + fine["group"])
    n = fine.drop_duplicates("row").set_index("row")["n"]
    rows = sorted(fine["row"].unique())
    sca = _dot(ax, fine, "row", rows, vmax)
    ax.set_yticklabels([f"{r} (n={n[r]:,})" for r in rows])
    ax.set_title("Stromal compartment by author label")
    cb = fig.colorbar(sca, ax=axes, shrink=0.6, pad=0.01)
    cb.set_label("mean log1p(CP10k)")
    for f in (0.1, 0.5, 0.9):
        axes[2].scatter([], [], s=20 + f * 600, c="grey", label=f"{f:.0%} cells")
    axes[2].legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.25),
        ncol=3,
        frameon=False,
        title="fraction expressing",
    )
    fig.suptitle(
        "ECM panel across major populations (QC-passing cells; dot size = fraction "
        "expressing, colour = mean expression)",
        fontweight="bold",
    )
    return _save(fig, "01_ecm_dotplot.png")


def fig_donor_level() -> Path:
    long = pd.read_csv(config.CHECKPOINTS_DIR / "stage4" / "ecm_pseudobulk_long.tsv.gz", sep="\t")
    long = long[long["enough_cells"] & (long["population"] != "ALL")]
    fig, axes = plt.subplots(2, 5, figsize=(22, 9), sharey=True)
    pops = [p for p in POP_ORDER if p in set(long["population"]) and p != "Unresolved"]
    pal = {"A": "#1f77b4", "B": "#ff7f0e"}
    for ax, g in zip(axes.ravel(), GENES, strict=True):
        d = long[long["gene"] == g]
        sns.boxplot(
            data=d,
            y="population",
            x="logcpm",
            hue="cohort",
            order=pops,
            ax=ax,
            palette=pal,
            fliersize=0,
            linewidth=0.7,
        )
        sns.stripplot(
            data=d,
            y="population",
            x="logcpm",
            hue="cohort",
            order=pops,
            ax=ax,
            palette=pal,
            dodge=True,
            size=2,
            alpha=0.6,
            legend=False,
        )
        ax.set(title=g, xlabel="log2 CPM (pseudobulk)", ylabel="")
        if ax is not axes[0, 0]:
            ax.get_legend().remove()
    fig.suptitle(
        "Donor-level ECM expression: each point = one specimen's pseudobulk "
        f"(≥{config.PSEUDOBULK_MIN_CELLS} cells)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "02_ecm_donor_level.png")


def fig_coherence() -> Path:
    coh = pd.read_csv(T / "stage4_panel_coherence.tsv", sep="\t")
    fig, axes = plt.subplots(
        1, 4, figsize=(24, 5.8), gridspec_kw={"width_ratios": [1, 1, 0.7, 0.7]}
    )
    for ax, cohort in zip(axes[:2], config.ANALYSIS_COHORTS, strict=True):
        c = pd.read_csv(T / f"stage4_ecm_gene_corr_{cohort}.tsv", sep="\t", index_col=0)
        n = int(coh.loc[coh["cohort"] == cohort, "n_specimens"].iloc[0])
        sns.heatmap(
            c,
            annot=True,
            fmt=".2f",
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            ax=ax,
            cbar=False,
            annot_kws={"fontsize": 7},
            square=True,
        )
        ax.set_title(f"Cohort {cohort}: gene–gene Spearman across {n} stromal pseudobulks")
    ax = axes[2]
    piv = coh.pivot(index="gene", columns="cohort", values="pc1_loading").reindex(GENES)
    piv.plot.barh(ax=ax, color=["#1f77b4", "#ff7f0e"])
    ax.invert_yaxis()
    ve = coh.drop_duplicates("cohort").set_index("cohort")["pc1_var_explained"]
    ax.set(
        xlabel="PC1 loading",
        title="Shared axis (PC1)\n" + ", ".join(f"{k}: {v:.0%} var" for k, v in ve.items()),
    )
    ax.axvline(1 / np.sqrt(10), ls="--", c="grey", lw=0.8)
    ax = axes[3]
    rows = []
    for cohort in config.ANALYSIS_COHORTS:
        s = pd.read_csv(T / f"stage4_summary_agreement_{cohort}.tsv", sep="\t", index_col=0)
        rows.append(
            {
                "cohort": cohort,
                "mean_z vs PC1": s.loc["mean_z", "pc1"],
                "mean_z vs log Σ CPM": s.loc["mean_z", "log_sum_cpm"],
            }
        )
    pd.DataFrame(rows).set_index("cohort").plot.bar(ax=ax, rot=0, color=["#27ae60", "#95a5a6"])
    ax.set(ylim=(0, 1), ylabel="Spearman ρ", title="Do summary choices agree?")
    for cont in ax.containers:
        ax.bar_label(cont, fmt="%.2f", fontsize=8)
    fig.suptitle(
        f"Panel coherence within the {COMBINED_STROMAL} compartment (log2-CPM pseudobulk)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "03_panel_coherence.png")


def fig_focus_choice() -> Path:
    cov = pd.read_csv(T / "stage4_focus_coverage.tsv", sep="\t")
    med = pd.read_csv(T / "stage4_ecm_median_logcpm.tsv", sep="\t")
    med["panel_mean_logcpm"] = med[GENES].mean(axis=1)
    m = cov.merge(med[["cohort", "population", "panel_mean_logcpm"]], on=["cohort", "population"])
    m = m[~m["population"].isin(["ALL", "Unresolved", "Malignant_epithelial"])]
    fig, axes = plt.subplots(1, 2, figsize=(17, 5.5))
    ax = axes[0]
    w = m.pivot_table(
        index="population", columns="cohort", values="paired_donors_ge10", fill_value=0
    )
    w = w.reindex([p for p in POP_ORDER if p in w.index])
    w.plot.barh(ax=ax, stacked=True, color=["#1f77b4", "#ff7f0e"])
    ax.invert_yaxis()
    ax.set(
        xlabel="paired donors with ≥10 cells in tumor AND adjacent",
        title="Coverage for a paired tumor-vs-adjacent comparison",
    )
    ax = axes[1]
    agg = m.groupby("population").agg(
        donors=("paired_donors_ge10", "sum"), expr=("panel_mean_logcpm", "mean")
    )
    ax.scatter(agg["donors"], agg["expr"], s=120, c="#2c3e50")
    for p, r in agg.iterrows():
        ax.annotate(p, (r["donors"], r["expr"]), xytext=(6, 4), textcoords="offset points")
    ax.set(
        xlabel="paired donors (A + B, ≥10 cells)",
        ylabel="mean panel log2-CPM (median over specimens)",
        title="Choice: high ECM expression AND adequate paired coverage",
    )
    fig.suptitle(
        "Focus population decision (malignant cells excluded: no adjacent counterpart)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "04_focus_population_choice.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="notebook")
    return [fig_dotplot(), fig_donor_level(), fig_coherence(), fig_focus_choice()]
