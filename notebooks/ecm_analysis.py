# ---
# jupyter:
#   jupytext:
#     formats: py:percent,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # An ECM gene program in ccRCC stroma: analysis walkthrough
#
# **Leviel Fluhr** · Opeongo computational biology home assignment
#
# ## How this notebook fits the repository
#
# The analysis is computed by the `opengo` Python package (`src/opengo/`) through one command:
#
# ```bash
# uv sync
# uv run opengo all   # ≈ 25 min from the raw inputs, peak ≈ 2.6 GB RAM
# ```
#
# That command writes every result table to `outputs/tables/` and every figure to
# `outputs/step_N_output/` and `report/figures/`.
#
# **This notebook is the reproducible walkthrough of those results.** It reads only the
# saved tables and figures, so it runs in seconds from a clean session. It does not recompute
# anything, and it shows the evidence behind each decision in the same order as the report.
# The computation is kept in the package rather than in notebook cells so that each stage can be
# re-run on its own, tested (`uv run pytest`), and rebuilt end to end. A clean rebuild was verified:
# 63 of 65 tables are identical (see Section 7).
#
# Supporting records:
# - `outputs/tables/decisions_log.tsv`: every decision, with its evidence and the alternatives
#   considered.
# - `README.md`: the per-stage reasoning.
#
# **Question.** Where is the 10-gene ECM panel (COL1A1, COL1A2, COL3A1, COL5A1, COL6A1, COL6A2, FN1,
# SPARC, DCN, LUM) expressed? Does it differ between tumor and adjacent tissue in one population?
# Does combining cohorts change that, and do bulk RNA and protein agree?

# %%
import io
from pathlib import Path

import pandas as pd
from IPython.display import Image, display
from PIL import Image as PILImage

from opengo import config

T = config.TABLES_DIR
OUT = config.OUTPUTS_DIR
pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 20)
pd.set_option("display.precision", 3)
pd.set_option("display.max_colwidth", 80)


def table(name: str) -> pd.DataFrame:
    """Read one result table written by the pipeline."""
    return pd.read_csv(T / f"{name}.tsv", sep="\t")


def show(path: Path, width: int = 950) -> None:
    """Embed a downscaled copy of a figure (full-resolution files are in the repository)."""
    img = PILImage.open(path)
    img.thumbnail((width * 1.5, 10_000))
    buf = io.BytesIO()
    img.convert("RGB").save(buf, format="JPEG", quality=85)
    display(Image(data=buf.getvalue(), width=width))


# %% [markdown]
# ## 0 · Inputs
#
# Before anything else, every supplied file (70 files, ≈ 4.8 GB) is checked against the package's
# `FILE_MANIFEST.tsv` by size and SHA-256. All files are present and unchanged.

# %%
table("input_verification")["status"].value_counts().rename("files").to_frame()

# %% [markdown]
# ## 1 · Loading and reconciliation
#
# - **Cohort A** (single-nucleus, GSE227898) is deposited as raw droplet matrices, so it was
#   restricted to the nuclei the authors annotated.
# - **Cohorts B and C** (single-cell) were loaded in full.
# - **Check:** every count must equal the package documentation (141,950 / 157,881 / 29,499 cells).

# %%
spec = table("stage1_specimens")
spec.groupby("cohort").agg(
    specimens=("sample_id", "size"),
    cells=("cells_loaded", "sum"),
    annotations_matched=("annotated_matched", "sum"),
    all_integer_counts=("integer_counts", "all"),
)

# %% [markdown]
# Genes were matched by Ensembl ID in A and C, and by unique symbol in B, which uses an older gene
# annotation. Ambiguous symbols were left out rather than guessed. A gene missing from a cohort is
# treated as unmeasured, never as zero.

# %%
table("gene_correspondence_summary").set_index("metric")

# %% [markdown]
# **A label correction, decided by marker genes.** The supplied label guide maps cohort A's "TRM"
# nuclei to lymphoid, while the authors group them with macrophages.
# - Lineage markers settle it: TRM nuclei express the myeloid marker CSF1R as often as macrophages
#   do, and almost never express T-cell markers.
# - They were therefore remapped to myeloid, and the original label is kept in the data.

# %%
trm = table("stage1_TRM_marker_check").pivot_table(
    index="marker", columns="author_label", values="frac_pos"
)
(100 * trm[["TRM", "Macrophages", "CD8 CTL", "CD4 T-cells"]]).round(0).loc[
    ["CSF1R", "CD68", "CD163", "CD3D", "CD3E", "CD2", "CD8A"]
].rename_axis(columns="% nuclei expressing")

# %% [markdown]
# ## 2 · Quality control
#
# **How cells were flagged.** Flags are set per specimen and per assay, not with global cut-offs. A
# cell is flagged for any of:
# - low counts or genes (median − 5 MAD);
# - high counts (median + 5 MAD);
# - a high mitochondrial fraction: median + 3 MAD **and** above 5 % in nuclei (where it signals
#   ambient RNA) or above 25 % in whole kidney cells;
# - a Scrublet doublet call.
#
# **Cohorts A and B were already author-filtered,** so few cells were removed. Removal was uneven
# across cell types, which is one reason population proportions are never compared between groups.

# %%
qp = table("stage2_qc_by_population")
(
    qp.assign(removed_pct=100 * (1 - qp["pass_frac"]))
    .pivot_table(index="population", columns="cohort", values="removed_pct")
    .round(1)
    .fillna("")
    .rename_axis(columns="% removed")
)

# %% [markdown]
# **Cohort C was excluded.** Its non-mitochondrial counts are capped at about 400 per cell (A and B
# reach 10⁴–10⁵). Most of its reads are mitochondrial, and the ECM genes are barely detected.

# %%
show(OUT / "step_2_output" / "08_cohortC_integrity.png")

# %% [markdown]
# In whole-specimen pseudobulk, the first principal component is almost entirely cohort
# (single-nucleus vs single-cell). Tumor-vs-adjacent effects are therefore always estimated
# **within** a cohort, never by pooling raw expression across cohorts.

# %%
r2 = table("stage2_pseudobulk_pca_r2_ALL")
r2.pivot_table(index="PC", columns="term", values="R2").join(
    r2.drop_duplicates("PC").set_index("PC")["var_explained"]
).loc[["PC1", "PC2", "PC3", "PC4"]].round(2)

# %% [markdown]
# ## 3 · Combining cohorts: unintegrated baseline vs Harmony
#
# Four embeddings were compared on a stratified subsample: A alone, A + B unintegrated, and A + B
# with Harmony using cohort or specimen as the batch.
#
# **Harmony with batch = cohort was chosen.**
# - It improves cohort mixing within each cell type, cell-type separation, and cross-cohort label
#   transfer.
# - It keeps the tumor/adjacent structure in stroma.
#
# **Batch = specimen was rejected.** Specimen is confounded with condition, so it erases that
# structure.
#
# The embedding is used only for pictures and label transfer, never for testing.

# %%
m = table("stage3_integration_metrics").set_index("embedding")
tr = table("stage3_label_transfer").groupby("embedding")["balanced_accuracy"].mean()
cs = table("stage3_condition_structure")
cs = cs[cs["population"] == "Stromal (fibro+mural)"].groupby("embedding")["condition_ASW"].mean()
pd.DataFrame(
    {
        "cohort mixing (batch ASW)": m["batch_ASW"],
        "cell-type separation (ASW)": m["celltype_ASW"],
        "label transfer A↔B": tr,
        "tumor/adjacent in stroma": cs,
    }
).round(2)

# %%
show(config.REPORT_DIR / "figures" / "figure1_qc_integration.png")

# %% [markdown]
# ## 4 · Where the ECM panel is expressed; summary score; focus population
#
# **Where it's expressed.** The stromal compartment has the highest expression of every panel gene
# in both cohorts. Endothelial cells express mainly SPARC and FN1.
#
# **Summary score.** The panel is summarised as the mean of per-gene z-scored log2-CPM. It agrees
# almost perfectly with the alternatives (first principal component; log summed CPM). DCN and LUM
# are fibroblast-specific and are always reported separately.
#
# **Focus population.** The **combined fibroblast + mural compartment** was chosen. It expresses the
# panel most strongly and has the best paired-donor coverage.

# %%
med = table("stage4_ecm_median_logcpm").set_index(["cohort", "population"])
med = med.drop(index=["ALL", "Unresolved", "Stromal_fibro_myo", "Pericyte"], level="population")
med.round(1).rename_axis(columns="median log2-CPM")

# %%
cov = table("stage4_focus_coverage")
cov[
    cov["population"].isin(
        ["Stromal (fibro+mural)", "Stromal_fibro_myo", "Pericyte", "Endothelial"]
    )
].pivot_table(
    index="population", columns="cohort", values="paired_donors_ge10", fill_value=0
).astype(int).rename_axis(columns="paired donors (≥10 cells)")

# %% [markdown]
# ## 5 · Tumor vs adjacent in the fibroblast + mural compartment
#
# **What is tested.** Pseudobulk raw counts per specimen, converted to log2-CPM.
#
# **Design.**
# - Donors are the replicates; repeated tumor regions are averaged per donor.
# - Each donor's tumor is compared with its own adjacent sample.
# - Each cohort gets an exact sign-flip p-value. With 3 pairs the smallest possible p is 0.25.
# - Cohorts are combined by fixed- and random-effects meta-analysis.
#
# **Result.**
# - **Cohort B:** the panel is higher in tumor in all 9 donors.
# - **Cohort A:** with 3 pairs it is too imprecise to confirm or refute B; its interval contains B's
#   estimate.

# %%
e = table("stage5_effects")
prim = e[e["analysis"].str.startswith("primary") & (e["feature"] == "panel_z")].set_index("cohort")
rows = [
    (
        "Cohort B (single-cell)",
        prim.loc["B", "n_donors"],
        prim.loc["B", "mean"],
        prim.loc["B", "ci_low"],
        prim.loc["B", "ci_high"],
        prim.loc["B", "n_up"],
        prim.loc["B", "p_signflip"],
    ),
    (
        "Cohort A (single-nucleus)",
        prim.loc["A", "n_donors"],
        prim.loc["A", "mean"],
        prim.loc["A", "ci_low"],
        prim.loc["A", "ci_high"],
        prim.loc["A", "n_up"],
        prim.loc["A", "p_signflip"],
    ),
    (
        "A + B, fixed effect",
        prim.loc["A+B", "n_donors"],
        prim.loc["A+B", "fixed"],
        prim.loc["A+B", "fixed_low"],
        prim.loc["A+B", "fixed_high"],
        prim.loc["A+B", "n_up"],
        prim.loc["A+B", "p_signflip"],
    ),
    (
        "A + B, random effects",
        prim.loc["A+B", "n_donors"],
        prim.loc["A+B", "random"],
        prim.loc["A+B", "random_low"],
        prim.loc["A+B", "random_high"],
        None,
        None,
    ),
]
pd.DataFrame(
    rows,
    columns=[
        "",
        "paired donors",
        "difference (SD)",
        "95% CI low",
        "95% CI high",
        "donors up",
        "sign-flip p",
    ],
).set_index("").round(3).astype({"paired donors": "Int64", "donors up": "Int64"})

# %% [markdown]
# **Per gene, cohort B.** Fibrillar collagens, FN1 and SPARC drive the increase. DCN and LUM do not
# follow; they are fibroblast-specific in a compartment that is mostly pericytes.

# %%
e[
    e["analysis"].str.startswith("primary")
    & (e["cohort"] == "B")
    & e["feature"].isin(config.ECM_GENES)
][["feature", "n_donors", "mean", "ci_low", "ci_high", "n_up", "q_BH"]].rename(
    columns={"mean": "log2FC"}
).astype({"n_donors": int, "n_up": int}).round(3).set_index("feature")

# %% [markdown]
# **Robustness check (main uncertainty: compartment composition).** Could the effect come from more
# fibroblasts in tumor stroma? Three checks say no:
# - The effect holds in pericytes alone.
# - It holds without the fibroblast genes (DCN, LUM).
# - The fibroblast share does not rise.
#
# It is also stable to every analysis choice tried, and endothelial cells show no change.

# %%
checks = [
    ("primary (fibroblast + mural)", "primary", "panel_z"),
    ("without DCN/LUM", "primary", "panel_z_noDCN_LUM"),
    ("pericytes only", "robustness: pericytes only", "panel_z"),
    ("≥ 20 cells per specimen", "sensitivity: >=20 cells", "panel_z"),
    ("≥ 50 cells per specimen", "sensitivity: >=50 cells", "panel_z"),
    ("one tumor region per donor", "sensitivity: one tumor region", "panel_z"),
    ("contrast: endothelial cells", "contrast: endothelial", "panel_z"),
]
b = e[e["cohort"] == "B"]
rob = pd.DataFrame(
    [
        b[b["analysis"].str.startswith(a) & (b["feature"] == f)]
        .iloc[0][["n_donors", "mean", "ci_low", "ci_high", "n_up"]]
        .rename(label)
        for label, a, f in checks
    ]
)
mm = table("stage5_mixed_model").set_index("cohort").loc["B"]
rob.loc["mixed model, all specimens"] = [
    mm["donors"],
    mm["estimate_z"],
    mm["ci_low"],
    mm["ci_high"],
    None,
]
fs = table("stage5_B_fibroblast_share_test").iloc[0]
rob.loc["fibroblast share (tumor − adjacent)"] = [
    fs["n_donors"],
    fs["mean"],
    fs["ci_low"],
    fs["ci_high"],
    fs["n_up"],
]
rob.rename(columns={"mean": "difference"}).round(3).astype(
    {"n_donors": "Int64", "n_up": "Int64"}
).rename_axis("cohort B check")

# %%
show(config.REPORT_DIR / "figures" / "figure2_ecm_program.png")

# %% [markdown]
# ## 6 · Linked bulk RNA and protein
#
# **Crosswalk.** Every cohort A specimen was traced through the crosswalk to real bulk columns. No
# pair was invented from donor identity.
# - Only **3 donors** have clean tumor and adjacent links. C3N-00242 is a sensitivity case.
# - The matched subset is therefore too small to decide.
#
# **Context from the whole CPTAC cohort.** It comes from the same study as cohort A (72 RNA and
# 80 protein pairs) and shows the same gene pattern as cohort B.
#
# **Composition.** Bulk differences depend strongly on composition. Tumors lose fibroblast markers
# but gain mural and endothelial markers, so the bulk rise may partly reflect expansion of those
# tissues.

# %%
table("stage6_specimen_inclusion")["reason"].value_counts().rename("specimens").to_frame()

# %%
e6 = table("stage6_effects")
e6 = e6[(e6["feature"] == "panel_z") & ~e6["scope"].str.contains("C3N-00242")]
e6.assign(scope=e6["scope"].str.split(" \\(").str[0]).set_index(["scope", "assay"])[
    ["n_donors", "mean", "ci_low", "ci_high", "n_up", "p_signflip"]
].rename(columns={"mean": "difference (SD)"}).round(3).assign(
    # n >= 72 uses 100,000 random sign patterns, so 0 means p < 1e-5
    p_signflip=lambda d: d["p_signflip"].map(lambda p: "< 1e-5" if p == 0 else f"{p:.3f}")
)

# %%
table("stage6_composition_corr")[
    [
        "marker_set",
        "mean_marker_change_z",
        "frac_donors_marker_down",
        "spearman_rho",
        "spearman_rho_noDCN_LUM",
    ]
].round(2).set_index("marker_set")

# %%
show(config.REPORT_DIR / "figures" / "figure3_multimodal.png")

# %% [markdown]
# ## 7 · Reproducibility check
#
# `opengo all` was re-run from the raw inputs into a separate folder and every table was compared
# with the committed results:
# - 63 of 65 tables are identical;
# - one differs at the 1e-9 level (floating point);
# - in one, only the Leiden-based ARI/NMI on the Harmony embeddings varies.
#
# No number quoted in the report changed.

# %%
table("rebuild_comparison")["status"].str.split(" ").str[0].value_counts().rename(
    "tables"
).to_frame()
