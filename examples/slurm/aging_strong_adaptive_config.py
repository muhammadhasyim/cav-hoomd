"""Shared constants for the strong-coupling adaptive aging campaign."""

from __future__ import annotations

from pathlib import Path

# Duplicate λ=0.03 at higher coupling plus three stronger values.
STRONG_ADAPTIVE_LAMBDAS = (0.03, 0.06, 0.09, 0.12)

DEFAULT_OUTPUT_BASE = Path(
    "/scratch/mh7373/projects/cav-hoomd/aging_strong_lambda_adaptive"
)

# Match 05_advanced_run.py / FINITE_SIZE_SCALING_SETUP.md / submit_aging_adaptive.sh.
ADAPTIVE_ERROR_TOLERANCE = 5.0
ADAPTIVE_INITIAL_FRACTION = 1e-5
ADAPTIVE_TIME_CONSTANT_PS = 50.0

# Match packed weak-lambda campaign SLURM settings (submit_aging_campaign.sh).
ADAPTIVE_PACKED_WALLTIME = "02:00:00"
PACKED_CONCURRENT_TASKS = 24
PACKED_CPUS = 8
PACKED_MEMORY = "32G"
PACKED_REPLICAS_PER_TASK = 1
