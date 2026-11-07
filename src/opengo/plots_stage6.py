"""Stage 6 figures -> outputs/step_6_output/. Run with `uv run opengo figures 6`."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from opengo import config  # noqa: E402

OUT_DIR = config.OUTPUTS_DIR / "step_6_output"
T = config.TABLES_DIR
GENES = config.ECM_GENES


def _save(fig, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUT_DIR / name
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_crosswalk() -> Path:
    v = pd.read_csv(T / "stage6_crosswalk_validation.tsv", sep="\t")
    checks = {
        "linked": "bulk link given",
        "primary": "primary (same parent)",
        "in_RNA_file": "column in bulk RNA file",
        "in_protein_file": "column in protein file",
        "protein_aliquot_matches_anchor": "protein aliquot = author table",
        "snRNA_parent_matches_protein_parent": "snRNA parent = protein parent",
    }
    v = v.sort_values(["condition", "donor_id", "sample_id"])
    m = v[list(checks)].astype(float)
    m.index = v["donor_id"] + " " + v["condition"].str[:3] + " (" + v["sample_id"] + ")"
    fig, axes = plt.subplots(1, 2, figsize=(16, 11), gridspec_kw={"width_ratios": [1.3, 1]})
    sns.heatmap(
        m.rename(columns=checks),
        cmap=["#e74c3c", "#27ae60"],
        cbar=False,
        ax=axes[0],
        linewidths=0.5,
        linecolor="white",
    )
    axes[0].set(title="All 34 cohort A snRNA specimens: link checks (green = pass)", ylabel="")
    plt.setp(axes[0].get_xticklabels(), rotation=30, ha="right")
    plt.setp(axes[0].get_yticklabels(), fontsize=7)
    for i, r in enumerate(v.itertuples()):
        axes[0].text(len(checks) + 0.1, i + 0.5, r.matching_status, fontsize=6, va="center")
    ax = axes[1]
    ax.axis("off")
    md = pd.read_csv(T / "stage6_matched_donors.tsv", sep="\t")
    counts = v["matching_status"].value_counts()
    txt = ["Matching status (34 specimens):"]
    txt += [f"  {k}: {n}" for k, n in counts.items()]
    txt += ["", "Donors with linked tumor AND adjacent:"]
    txt += [
        f"  {r.donor_id}: {'primary pair' if r.primary_pair else 'sensitivity only'}"
        for r in md.itertuples()
    ]
    txt += [
        "",
        "Notes:",
        "  • bulk RNA is a donor/condition column;",
        "    its exact RNA aliquot is not reconciled",
        "  • C3N-00242 tumor: protein from two parent",
        "    specimens (-01, -02) → sensitivity only",
        "  • C3L-00908 tumor: protein from a different",
        "    parent (-03) → excluded (no adjacent anyway)",
        "  • unlinked snRNA specimens: no bulk pair",
        "    invented from donor identity",
    ]
    ax.text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=10)
    fig.suptitle("Multimodal crosswalk validation", fontweight="bold")
    fig.tight_layout()
    return _save(fig, "01_crosswalk_validation.png")


def fig_agreement() -> Path:
    a = pd.read_csv(T / "stage6_agreement.tsv", sep="\t")
    keep = {
        "snRNA/scRNA stroma, cohort B (Stage 5)": "B scRNA stroma (9 pairs)",
        "snRNA/scRNA stroma, cohort A (Stage 5)": "A snRNA stroma (3 pairs)",
    }
    for ev in a["evidence"].unique():
        if ev.startswith("bulk RNA (all"):
            keep[ev] = "bulk RNA, all CPTAC pairs (72)"
        elif ev.startswith("bulk protein (all"):
            keep[ev] = "protein, all CPTAC pairs (80)"
        elif "matched subset" in ev:
            keep[ev] = ev.split(" — ")[0] + " (matched, 3)"
    a = a[a["evidence"].isin(keep)].assign(evidence=lambda d: d["evidence"].map(keep))
    order = list(dict.fromkeys(keep.values()))
    piv = a.pivot_table(index="feature", columns="evidence", values="mean").reindex(
        GENES + ["panel_z"]
    )[order]
    fig, ax = plt.subplots(figsize=(14, 7.5))
    sns.heatmap(
        piv,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        annot=True,
        fmt="+.2f",
        annot_kws={"fontsize": 8},
        ax=ax,
        cbar_kws={"label": "tumor − adjacent (log2; panel_z in SD units)"},
    )
    ax.axhline(len(GENES), color="k", lw=1.5)
    ax.axvline(2, color="k", lw=1.5)
    ax.axvline(4, color="k", lw=1.5)
    ax.set(
        xlabel="",
        ylabel="",
        title="Direction of the tumor-vs-adjacent difference across "
        "assays (each estimated within its own assay)",
    )
    plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
    fig.tight_layout()
    return _save(fig, "02_cross_assay_agreement.png")


def fig_matched() -> Path:
    pdn = pd.read_csv(T / "stage6_matched_per_donor.tsv", sep="\t")
    pdn = pdn[~pdn["assay"].str.contains(r"\[")]
    assays = ["snRNA stroma", "snRNA whole specimen", "bulk RNA", "bulk protein"]
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.5), gridspec_kw={"width_ratios": [1, 1.4]})
    ax = axes[0]
    pal = dict(zip(sorted(pdn["donor_id"].unique()), sns.color_palette("Set1"), strict=False))
    for donor, d in pdn.groupby("donor_id"):
        d = d.set_index("assay").reindex(assays)
        ax.plot(range(len(assays)), d["panel_z"], "-o", color=pal[donor], label=donor)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(range(len(assays)), assays, rotation=20, ha="right")
    ax.set(ylabel="panel score, tumor − adjacent", title="Matched donors: same donor, four assays")
    ax.legend()
    ax = axes[1]
    sf = pd.read_csv(T / "stage6_snrna_stromal_fraction.tsv", sep="\t")
    sf = sf[sf["donor_id"].isin(pdn["donor_id"].unique()) | (sf["donor_id"] == "C3N-00242")]
    sf = sf[sf["matching_status"].isin(["same_parent_specimen", "multiple_protein_specimen_ids"])]
    w = sf.pivot_table(index="donor_id", columns="condition", values="stromal_fraction")
    w.plot.bar(ax=ax, color={"Adjacent": "#2471a3", "Tumor": "#c0392b"}, rot=0)
    ax.set(
        ylabel="stromal (fibro + mural) fraction of nuclei",
        title="snRNA composition of the linked specimens\n(different aliquot from bulk)",
    )
    fig.suptitle(
        "Matched subset (3 primary pairs; C3N-00242 shown in composition only)", fontweight="bold"
    )
    fig.tight_layout()
    return _save(fig, "03_matched_subset.png")


def fig_composition() -> Path:
    sc = pd.read_csv(T / "stage6_composition_scores.tsv", sep="\t")
    corr = pd.read_csv(T / "stage6_composition_corr.tsv", sep="\t").set_index("marker_set")
    adj = pd.read_csv(T / "stage6_composition_adjusted.tsv", sep="\t").iloc[0]
    fig, axes = plt.subplots(1, 4, figsize=(22, 5))
    for ax, k in zip(axes[:3], ["fibroblast", "mural", "endothelial"], strict=True):
        ax.scatter(sc[f"{k}_markers_z"], sc["ECM_panel_z"], s=18, alpha=0.7)
        ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="k", lw=0.5)
        ax.set(
            xlabel=f"{k} marker change (z)",
            ylabel="ECM panel change (z)",
            title=f"{k}: ρ = {corr.loc[k, 'spearman_rho']:.2f} "
            f"(8-gene panel {corr.loc[k, 'spearman_rho_noDCN_LUM']:.2f})\n"
            f"mean marker change {corr.loc[k, 'mean_marker_change_z']:+.2f}",
        )
    ax = axes[3]
    ax.bar(
        ["10-gene panel", "8 genes\n(no DCN/LUM)", "10-gene, extrapolated\nto no marker change"],
        [
            adj["ECM_change_unadjusted"],
            adj["ECM_noDCN_LUM_change"],
            adj["ECM_change_at_zero_stromal_marker_change"],
        ],
        color=["#95a5a6", "#7f8c8d", "#2c3e50"],
    )
    ax.set(
        ylabel="ECM panel change (z)",
        title="Mean ECM change\n(dark bar = regression extrapolation; interpret cautiously)",
    )
    fig.suptitle(
        "Composition proxy in bulk RNA (72 CPTAC paired donors; markers are non-panel "
        "genes: fibroblast PDGFRA/C7/FBLN1, mural RGS5/NOTCH3/PDGFRB/MYH11, "
        "endothelial VWF/CDH5/PECAM1)",
        fontweight="bold",
        fontsize=11,
    )
    fig.tight_layout()
    return _save(fig, "04_composition_proxy.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="notebook")
    return [fig_crosswalk(), fig_agreement(), fig_matched(), fig_composition()]
