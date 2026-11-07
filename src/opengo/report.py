"""Stage 7: render the HTML report to PDF with headless Chrome/Chromium and count its pages.

The report source is report/report.html (figures in report/figures/, built by
`opengo figures report`). Set CHROME=/path/to/chrome if it is not found automatically.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from pypdf import PdfReader

from opengo import config

REPORT_DIR = config.PROJECT_ROOT / "report"
CANDIDATES = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "google-chrome",
    "chromium",
    "chromium-browser",
]


def find_chrome() -> str:
    for c in [os.environ.get("CHROME", "")] + CANDIDATES:
        if c and (Path(c).exists() or shutil.which(c)):
            return c
    raise FileNotFoundError("Chrome/Chromium not found; set CHROME=/path/to/chrome")


def render_pdf(
    html: Path = REPORT_DIR / "report.html", pdf: Path = REPORT_DIR / "report.pdf"
) -> tuple[Path, int]:
    subprocess.run(
        [
            find_chrome(),
            "--headless=new",
            "--disable-gpu",
            "--no-pdf-header-footer",
            f"--print-to-pdf={pdf}",
            html.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
    )
    return pdf, len(PdfReader(pdf).pages)
