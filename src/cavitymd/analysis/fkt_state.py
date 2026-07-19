"""Serialize and deserialize F(k,t) reference state for runtime extensions."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def serialize_fkt_references(references: list[dict[str, Any]]) -> dict[str, np.ndarray]:
    """Convert in-memory F(k,t) reference frames to an NPZ payload."""
    count = len(references)
    payload: dict[str, np.ndarray] = {
        "reference_count": np.array([count], dtype=np.int64),
        "reference_times_ps": np.zeros(count, dtype=np.float64),
        "reference_timesteps": np.zeros(count, dtype=np.int64),
    }
    if count == 0:
        return payload

    rhok_real = np.stack([np.asarray(ref["rhok_real"], dtype=np.float64) for ref in references])
    rhok_imag = np.stack([np.asarray(ref["rhok_imag"], dtype=np.float64) for ref in references])
    payload["rhok_real"] = rhok_real
    payload["rhok_imag"] = rhok_imag
    for index, ref in enumerate(references):
        payload["reference_times_ps"][index] = float(ref["time_ps"])
        payload["reference_timesteps"][index] = int(ref["timestep"])
    if "wavevectors" in references[0]:
        payload["wavevectors"] = np.asarray(references[0]["wavevectors"], dtype=np.float64)
    return payload


def deserialize_fkt_references(payload: dict[str, np.ndarray]) -> list[dict[str, Any]]:
    """Rebuild F(k,t) reference frames from an NPZ payload."""
    count = int(payload["reference_count"][0])
    if count == 0:
        return []

    wavevectors = payload.get("wavevectors")
    references: list[dict[str, Any]] = []
    for index in range(count):
        ref: dict[str, Any] = {
            "timestep": int(payload["reference_timesteps"][index]),
            "time_ps": float(payload["reference_times_ps"][index]),
            "rhok_real": np.asarray(payload["rhok_real"][index], dtype=np.float64),
            "rhok_imag": np.asarray(payload["rhok_imag"][index], dtype=np.float64),
        }
        if wavevectors is not None:
            ref["wavevectors"] = np.asarray(wavevectors, dtype=np.float64)
        references.append(ref)
    return references


def save_fkt_state(
    path: Path,
    *,
    references: list[dict[str, Any]],
    last_reference_time: float | None,
    last_output_time: float | None,
    kmag: float,
    reference_interval_ps: float,
) -> None:
    """Persist serialized F(k,t) tracker state for a runtime extension."""
    payload = serialize_fkt_references(references)
    payload["kmag"] = np.array([kmag], dtype=np.float64)
    payload["reference_interval_ps"] = np.array([reference_interval_ps], dtype=np.float64)
    payload["last_reference_time"] = np.array(
        [-1.0 if last_reference_time is None else float(last_reference_time)],
        dtype=np.float64,
    )
    payload["last_output_time"] = np.array(
        [-1.0 if last_output_time is None else float(last_output_time)],
        dtype=np.float64,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **payload)


def load_fkt_state(path: Path) -> dict[str, Any]:
    """Load serialized F(k,t) tracker state from disk."""
    with np.load(path, allow_pickle=False) as archive:
        payload = {key: archive[key] for key in archive.files}
    references = deserialize_fkt_references(payload)
    last_reference = float(payload["last_reference_time"][0])
    last_output = float(payload["last_output_time"][0])
    return {
        "references": references,
        "last_reference_time": None if last_reference < 0.0 else last_reference,
        "last_output_time": None if last_output < 0.0 else last_output,
        "kmag": float(payload["kmag"][0]),
        "reference_interval_ps": float(payload["reference_interval_ps"][0]),
    }
