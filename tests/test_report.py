"""Deliverable checks: 4-page report, report figures, executed notebook."""

import json

import pytest
from pypdf import PdfReader

from opengo import config

R = config.PROJECT_ROOT / "report"


def test_report_within_four_pages():
    pdf = R / "report.pdf"
    if not pdf.exists():
        pytest.skip("run `opengo report` first")
    assert 1 <= len(PdfReader(pdf).pages) <= 4


def test_report_figures_referenced_exist():
    html = (R / "report.html").read_text()
    for name in ["figure1_qc_integration", "figure2_ecm_program", "figure3_multimodal"]:
        assert f"figures/{name}.png" in html
        assert (R / "figures" / f"{name}.png").exists()


def test_notebook_executed_without_errors():
    nb_path = config.PROJECT_ROOT / "notebooks" / "ecm_analysis.ipynb"
    if not nb_path.exists():
        pytest.skip("execute the notebook first")
    nb = json.loads(nb_path.read_text())
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    assert code and all(c.get("execution_count") for c in code)
    assert not [o for c in code for o in c.get("outputs", []) if o.get("output_type") == "error"]
