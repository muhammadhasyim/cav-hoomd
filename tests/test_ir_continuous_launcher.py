"""Tests for continuous IR production launcher."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = REPO_ROOT / "aging_weak_lambda_ir" / "run_local_ir_continuous.sh"


def test_ir_continuous_dry_run_queue_and_defaults(tmp_path: Path) -> None:
    base = tmp_path / "aging_weak_lambda_ir"
    base.mkdir()
    ic_dir = tmp_path / "aging_weak_lambda_ic_library"
    ic_dir.mkdir()
    ic = ic_dir / "init-500-from-ic-library.gsd"
    ic.write_bytes(b"")
    env = {
        **dict(**{k: v for k, v in __import__("os").environ.items()}),
        "REPO": str(tmp_path),
        "PY": sys.executable,
        "DRY_RUN": "1",
        "N_REPLICAS": "2",
        "RUNTIME_PS": "400.0",
        "CYCLES": "3",
    }
    result = subprocess.run(
        ["bash", str(LAUNCHER)],
        check=True,
        capture_output=True,
        text=True,
        env=env,
        cwd=str(REPO_ROOT),
    )
    queue = base / "logs" / "ir_continuous" / "ir_queue.tsv"
    assert queue.is_file()
    lines = queue.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 5 * 2 * 3
    assert "400.0ps" in result.stdout
