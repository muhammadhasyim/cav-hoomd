"""CLI tests for dipole FDR flags on 05_advanced_run.py."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "examples" / "05_advanced_run.py"


def test_help_lists_dipole_fdr_flags() -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    help_text = result.stdout
    assert "--enable-dipole-fdr" in help_text
    assert "--dipole-fdr-output-period-ps" in help_text
    assert "--dipole-fdr-max-correlation-time-ps" in help_text
