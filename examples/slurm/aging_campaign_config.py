"""Shared constants for the weak-coupling aging SLURM campaign."""

from __future__ import annotations

# Primary production length used for the first campaign wave (ps).
CAMPAIGN_PRIMARY_RUNTIME_PS = 1600.0

# Target runtime for new submissions and extension goal (ps).
CAMPAIGN_TARGET_RUNTIME_PS = 2000.0

# Extension segment length when continuing from a checkpoint (ps).
CAMPAIGN_EXTENSION_RUNTIME_PS = (
    CAMPAIGN_TARGET_RUNTIME_PS - CAMPAIGN_PRIMARY_RUNTIME_PS
)

CAMPAIGN_SWITCH_TIME_PS = 200.0
CAMPAIGN_FKT_REFERENCE_INTERVAL_PS = 200.0
CAMPAIGN_FKT_MAX_REFS_PRIMARY = 8
CAMPAIGN_FKT_MAX_REFS_TARGET = 10
CAMPAIGN_HDF5_PERIOD_PS = 1.0
CAMPAIGN_FKT_KMAG = 1.0

# Packed-task walltime for 2000 ps fixed-dt production (HH:MM:SS).
CAMPAIGN_PACKED_WALLTIME = "02:00:00"

# Extension jobs are shorter.
CAMPAIGN_EXTENSION_WALLTIME = "00:45:00"

CHECKPOINT_GSD_TEMPLATE = "checkpoint_replica_{replica}.gsd"
FKT_STATE_TEMPLATE = "{job_name}-{replica}_fkt_state.npz"
