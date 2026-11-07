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
# # ECM gene program in ccRCC — analysis notebook
#
# This notebook walks through the whole analysis in order and shows the numbers and figures the
# report is built from. Heavy computation lives in the `opengo` package and its CLI
# (`uv run opengo all` regenerates everything from the supplied inputs, ≈ 45 min, peak ≈ 3 GB);
# this notebook reads the resulting tables and figures, so it runs in seconds from a clean session.
#
# Every decision, with its evidence and the alternatives considered, is in
# `outputs/tables/decisions_log.tsv`; the reasoning per stage is in `README.md`.

# %%
import io
from pathlib import Path

import pandas as pd
from PIL import Image as PILImage
from IPython.display import Image, Markdown, display

from opengo import config

T = config.TABLES_DIR
pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 20)


def show(path: Path, width: int = 1000) -> None:
    """Embed a downscaled copy (keeps the notebook small); full-size files are in outputs/."""
    img = PILImage.open(path)
    img.thumbnail((width * 1.5, 10_000))
    buf = io.BytesIO()
    img.convert("RGB").save(buf, format="JPEG", quality=82)
    display(Image(data=buf.getvalue(), width=width))


def table(name: str) -> pd.DataFrame:
    return pd.read_csv(T / f"{name}.tsv", sep="\t")


# %% [markdown]
# ## 0 · Inputs and environment
# All 70 supplied files are checked against `FILE_MANIFEST.tsv` (size + SHA-256).

# %%
table("input_verification")["status"].value_counts()

# %% [markdown]
# ## 1 · Loading and reconciliation
# Cohort A raw droplet matrices are restricted to the author-annotated nuclei; B and C are loaded
# in full. Cells per cohort must equal the documentation (141,950 / 157,881 / 29,499).

# %%
spec = table("stage1_specimens")
spec.groupby("cohort").agg(
    specimens=("sample_id", "size"),
    cells=("cells_loaded", "sum"),
    annotations_matched=("annotated_matched", "sum"),
    integer_counts=("integer_counts", "all"),
)

# %%
table("gene_correspondence_summary")

# %% [markdown]
# Label check: the suggested guide maps A "TRM" to lymphoid; lineage markers say myeloid.

# %%
trm = table("stage1_TRM_marker_check")
trm.pivot_table(index="marker", columns="author_label", values="frac_pos").round(2)[
    ["TRM", "Macrophages", "CD8 CTL", "CD4 T-cells"]
]

# %%
show(config.OUTPUTS_DIR / "step_1_output" / "07_TRM_marker_check.png", 700)

# %% [markdown]
# ## 2 · QC and exploratory structure
# Per-specimen, assay-aware flags. Flagged cells are removed in A and B; cohort C is excluded
# (non-mitochondrial counts capped at ≤ 409 per cell).

# %%
q = table("stage2_qc_by_specimen")
q.groupby("cohort").apply(
    lambda d: pd.Series({"cells": d["cells"].sum(),
                         "removed_frac": 1 - (d["cells"] * d["pass_frac"]).sum() / d["cells"].sum(),
                         "median_pct_mt": d["median_pct_mt"].median()}),
    include_groups=False,
).round(3)

# %%
show(config.OUTPUTS_DIR / "step_2_output" / "08_cohortC_integrity.png")

# %%
table("stage2_pseudobulk_pca_r2_ALL").query("term == 'cohort'").head(3).round(2)

# %% [markdown]
# ## 3 · Combining cohorts: baseline vs Harmony

# %%
table("stage3_integration_metrics").round(3)

# %%
table("stage3_label_transfer").round(3)

# %%
show(config.OUTPUTS_DIR / "step_3_output" / "01_umap_baseline_vs_harmony.png", 800)

# %% [markdown]
# ## 4 · ECM program across populations; summary and focus population

# %%
med = table("stage4_ecm_median_logcpm").set_index(["cohort", "population"])
med.drop(index="ALL", level="population").round(1)

# %%
table("stage4_focus_coverage").query("population != 'ALL'")

# %%
show(config.OUTPUTS_DIR / "step_4_output" / "03_panel_coherence.png")

# %% [markdown]
# ## 5 · Tumor vs adjacent in the fibroblast + mural compartment
# Pseudobulk log2-CPM per specimen; donors are replicates (tumor regions averaged; paired within
# donor); exact sign-flip p; random-effects meta-analysis across cohorts.

# %%
e = table("stage5_effects")
cols = ["analysis", "cohort", "feature", "n_donors", "mean", "ci_low", "ci_high", "n_up",
        "p_signflip", "q_BH", "random", "random_low", "random_high", "I2"]
e[e["feature"].isin(["panel_z", "panel_z_noDCN_LUM"])][cols].round(3)

# %%
e[e["analysis"].str.startswith("primary") & e["cohort"].eq("B")][
    ["feature", "mean", "ci_low", "ci_high", "n_up", "q_BH"]
].round(3)

# %%
show(config.OUTPUTS_DIR / "step_5_output" / "03_robustness.png")

# %% [markdown]
# ## 6 · Linked modalities

# %%
v = table("stage6_crosswalk_validation")
v[["linked", "primary", "in_RNA_file", "in_protein_file", "protein_aliquot_matches_anchor",
   "snRNA_parent_matches_protein_parent"]].sum()

# %%
table("stage6_matched_donors")

# %%
e6 = table("stage6_effects")
e6[e6["feature"] == "panel_z"][["scope", "assay", "n_donors", "mean", "ci_low", "ci_high", "n_up",
                                "p_signflip"]].round(3)

# %%
table("stage6_composition_corr").round(3)

# %%
show(config.OUTPUTS_DIR / "step_6_output" / "02_cross_assay_agreement.png")

# %% [markdown]
# ## Report figures
# The three report figures (one per assignment section) are built by `opengo figures report`.

# %%
for f in sorted((config.PROJECT_ROOT / "report" / "figures").glob("figure*.png")):
    display(Markdown(f"**{f.name}**"))
    show(f)
