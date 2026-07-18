"""Bootstrap sys.path for scripts under reproduction/."""
from __future__ import annotations

import sys
from pathlib import Path

REPRO_ROOT = Path(__file__).resolve().parents[1]
SHARED_DIR = REPRO_ROOT / "shared"


def setup_repro_path(*extra_subdirs: str) -> Path:
    """Insert reproduction/shared (and optional subdirs) onto sys.path."""
    for path in (SHARED_DIR, *(REPRO_ROOT / sub for sub in extra_subdirs)):
        s = str(path)
        if s not in sys.path:
            sys.path.insert(0, s)
    return REPRO_ROOT


def setup_profile(args, default: str = "paper"):
    """Load a dataset profile after ensuring shared modules are importable."""
    setup_repro_path()
    from dataset_profile import setup_profile as _setup_profile

    return _setup_profile(args, default=default)
