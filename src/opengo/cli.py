"""Command-line entry point: `uv run opengo <command>`.

Each stage of PLAN.md gets one command. Commands read inputs/checkpoints and write
outputs under outputs/, so the pipeline can be re-run from a clean session.
"""

from __future__ import annotations

import sys
import time

import click
from rich.console import Console
from rich.table import Table

from opengo import config

console = Console()


@click.group()
def main() -> None:
    """ECM program in ccRCC: reproducible pipeline (see PLAN.md for the stages)."""


@main.command()
def verify() -> None:
    """Stage 0: check every input file against FILE_MANIFEST.tsv (size + SHA-256)."""
    from opengo.verify import verify_inputs

    t0 = time.perf_counter()
    with console.status("Hashing inputs (≈4.8 GB)…"):
        result = verify_inputs()
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out = config.TABLES_DIR / "input_verification.tsv"
    result.to_csv(out, sep="\t", index=False)

    counts = result["status"].value_counts()
    table = Table(title="Input verification")
    table.add_column("status")
    table.add_column("files", justify="right")
    for status, n in counts.items():
        table.add_row(str(status), str(n))
    console.print(table)
    bad = result[result["status"] != "ok"]
    if not bad.empty:
        console.print(bad.to_string(index=False), style="red")
    console.print(
        f"Wrote {out.relative_to(config.PROJECT_ROOT)} in {time.perf_counter() - t0:.0f}s"
    )
    sys.exit(0 if bad.empty else 1)


def _not_yet(stage: str) -> None:
    console.print(f"[yellow]Not implemented yet — built at {stage} after check-in.[/yellow]")
    sys.exit(2)


@main.command()
@click.argument("cohort", type=click.Choice(["A", "B", "C", "all"]))
def load(cohort: str) -> None:
    """Stage 1: load native files into raw-count checkpoints (per specimen + per cohort)."""
    import pandas as pd

    from opengo.io import load_cohort

    cohorts = config.COHORTS if cohort == "all" else (cohort,)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    for c in cohorts:
        t0 = time.perf_counter()
        console.rule(f"Cohort {c}")

        def progress(rec, t0=t0):
            console.print(
                f"  {rec.sample_id:<18} cells={rec.cells_loaded:>6} "
                f"annotated={rec.annotated_matched}/{rec.annotated_expected} "
                f"genes={rec.n_genes} int={rec.integer_counts} "
                f"[{time.perf_counter() - t0:.0f}s]"
            )

        table = load_cohort(c, progress=progress)
        out = config.TABLES_DIR / f"stage1_specimens_{c}.tsv"
        table.to_csv(out, sep="\t", index=False)
        console.print(
            f"Cohort {c}: {len(table)} specimens, {table['cells_loaded'].sum():,} cells "
            f"-> {config.cohort_raw_path(c).relative_to(config.PROJECT_ROOT)} "
            f"({time.perf_counter() - t0:.0f}s)"
        )
    parts = [
        pd.read_csv(p, sep="\t")
        for c in config.COHORTS
        if (p := config.TABLES_DIR / f"stage1_specimens_{c}.tsv").exists()
    ]
    pd.concat(parts).to_csv(config.TABLES_DIR / "stage1_specimens.tsv", sep="\t", index=False)


@main.command()
def genes() -> None:
    """Stage 1 wrap-up: gene correspondence, ECM panel check, design/coverage tables."""
    from opengo.genes import build_gene_table, ecm_check
    from opengo.io import stage1_summaries

    summaries = stage1_summaries()
    console.print(summaries["stage1_population_overrides"].to_string(index=False))

    table, summary = build_gene_table()
    table.to_csv(config.TABLES_DIR / "gene_correspondence.tsv", sep="\t", index=False)
    summary.to_csv(config.TABLES_DIR / "gene_correspondence_summary.tsv", sep="\t", index=False)
    ecm = ecm_check(table)
    ecm.to_csv(config.TABLES_DIR / "ecm_panel_presence.tsv", sep="\t", index=False)
    console.print(summary.to_string(index=False))
    console.print(ecm.to_string(index=False))
    if not ecm["shared_ABC"].fillna(False).all():
        console.print("[red]ECM panel gene missing from a cohort[/red]")
        sys.exit(1)


@main.command()
@click.argument("stage", type=click.Choice(["1", "2", "3", "4", "5", "6"]))
def figures(stage: str) -> None:
    """Render a stage's figures into outputs/step_<stage>_output/."""
    import importlib

    module = importlib.import_module(f"opengo.plots_stage{stage}")
    for path in module.make_all():
        console.print(f"  {path.relative_to(config.PROJECT_ROOT)}")


@main.command()
def qc() -> None:
    """Stage 2: per-cell QC metrics and per-specimen flags (nothing removed)."""
    from opengo import stage2

    def progress(sid, m):
        console.print(
            f"  {sid:<18} n={len(m):>6} pass={m['qc_pass'].mean():.1%} "
            f"med_counts={m['total_counts'].median():.0f} med_mt={m['pct_mt'].median():.2f}%"
        )

    tables = stage2.run_qc(progress=progress)
    for name, df in tables.items():
        if name != "stage2_qc_by_specimen":
            console.rule(name)
            console.print(df.to_string(index=False))


@main.command()
def explore() -> None:
    """Stage 2: pseudobulk, exploratory PCA, marker labels (C), stromal identity checks."""
    from opengo import stage2

    out = stage2.run_explore(progress=lambda c, s: None)
    for name, df in out.items():
        console.rule(name)
        console.print(df.to_string(index=False))


@main.command()
def integrate() -> None:
    """Stage 3: A alone, A+B unintegrated, A+B Harmony (cohort / specimen); evaluation metrics."""
    from opengo import integrate as integ

    t0 = time.perf_counter()
    out = integ.run(
        progress=lambda msg: console.print(f"  [{time.perf_counter() - t0:.0f}s] {msg}")
    )
    for name in ["stage3_integration_metrics", "stage3_label_transfer"]:
        console.rule(name)
        console.print(out[name].round(3).to_string(index=False))


@main.command()
def overview() -> None:
    """Stage 4: ECM panel across populations, panel coherence, focus-population coverage."""
    from opengo import stage4

    t0 = time.perf_counter()
    out = stage4.run()
    for name in ["stage4_ecm_median_logcpm", "stage4_focus_coverage", "stage4_panel_coherence"]:
        console.rule(name)
        console.print(out[name].round(2).to_string(index=False))
    console.print(f"done in {time.perf_counter() - t0:.0f}s")


@main.command()
def pseudobulk() -> None:
    """Stage 5: tumor vs adjacent (paired, donor-level) in the stromal compartment."""
    from opengo import stage5

    out = stage5.run()
    e = out["stage5_effects"]
    show = e[e["analysis"].str.startswith("primary") & e["feature"].str.startswith("panel")]
    cols = ["cohort", "feature", "n_donors", "mean", "ci_low", "ci_high", "n_up", "p_signflip"]
    cols += ["random", "random_low", "random_high", "I2"]
    console.print(show[[c for c in cols if c in show]].round(3).to_string(index=False))


@main.command()
def multimodal() -> None:
    """Stage 6: crosswalk validation, snRNA vs bulk RNA vs protein, composition proxy."""
    from opengo import stage6

    out = stage6.run()
    v = out["stage6_crosswalk_validation"]
    checks = [
        "linked",
        "primary",
        "in_RNA_file",
        "in_protein_file",
        "protein_aliquot_matches_anchor",
        "snRNA_parent_matches_protein_parent",
        "condition_matches_bulk_file",
    ]
    console.print(v[checks].sum().to_string())
    console.print(out["stage6_matched_donors"].to_string(index=False))
    e = out["stage6_effects"]
    show = e[e["feature"] == "panel_z"][
        ["scope", "assay", "n_donors", "mean", "ci_low", "ci_high", "n_up", "p_signflip"]
    ]
    console.print(show.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
