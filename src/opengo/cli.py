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
@click.argument("cohort", type=click.Choice(["A", "B", "C"]))
def load(cohort: str) -> None:
    """Stage 1: load one cohort's native files into a raw-count checkpoint."""
    _not_yet("Stage 1")


@main.command()
def qc() -> None:
    """Stage 2: per-specimen QC metrics and exploratory structure."""
    _not_yet("Stage 2")


@main.command()
def integrate() -> None:
    """Stage 3: unintegrated baseline vs Harmony, with evaluation metrics."""
    _not_yet("Stage 3")


@main.command()
def pseudobulk() -> None:
    """Stages 4-5: ECM overview and tumor vs adjacent pseudobulk testing."""
    _not_yet("Stage 4-5")


@main.command()
def multimodal() -> None:
    """Stage 6: crosswalk validation and bulk RNA / protein comparison."""
    _not_yet("Stage 6")


if __name__ == "__main__":
    main()
