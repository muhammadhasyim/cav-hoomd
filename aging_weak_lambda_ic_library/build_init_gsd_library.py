#!/usr/bin/env python3
"""Build a multi-frame init GSD for aging_weak_lambda from IC library trajectories.

Maps output replica index ``r`` to IC-library stream ``r % n_source`` and frame
``r // n_source``, matching the OpenMM-style ``--frame replica`` convention.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import gsd.hoomd
import numpy as np


def source_mapping(output_replica: int, n_source_replicas: int) -> tuple[int, int]:
    """Return (source_replica, source_frame) for an output replica index."""
    return output_replica % n_source_replicas, output_replica // n_source_replicas


def min_frames_required(n_output_replicas: int, n_source_replicas: int) -> int:
    """Minimum frames each source trajectory must contain."""
    if n_output_replicas <= 0:
        raise ValueError("n_output_replicas must be positive")
    if n_source_replicas <= 0:
        raise ValueError("n_source_replicas must be positive")
    return (n_output_replicas - 1) // n_source_replicas + 1


def add_cavity_particle_at_origin(snapshot) -> None:
    """Append a cavity photon ``L`` at the origin (q=0, switch-protocol IC)."""
    if "L" in snapshot.particles.types:
        raise ValueError("snapshot already contains a cavity particle")

    box = np.array(snapshot.configuration.box[:3], dtype=float)
    newpos = np.zeros(3, dtype=float)
    image_flags = np.floor((newpos + box / 2.0) / box).astype(int)
    wrapped_position = newpos - image_flags * box

    old_n = int(snapshot.particles.N)
    new_n = old_n + 1
    snapshot.particles.types = list(snapshot.particles.types) + ["L"]

    def _grow_array(name: str, default) -> None:
        current = getattr(snapshot.particles, name, None)
        if current is None:
            return
        arr = np.asarray(current)
        if arr.ndim == 1:
            grown = np.empty(new_n, dtype=arr.dtype)
            grown[:old_n] = arr
            grown[old_n] = default
        elif arr.ndim == 2:
            grown = np.empty((new_n, arr.shape[1]), dtype=arr.dtype)
            grown[:old_n] = arr
            grown[old_n] = default
        else:
            raise ValueError(f"unsupported array rank for particles.{name}")
        setattr(snapshot.particles, name, grown)

    _grow_array("typeid", 2)
    _grow_array("position", wrapped_position)
    _grow_array("velocity", np.zeros(3, dtype=float))
    _grow_array("charge", 0.0)
    _grow_array("mass", 1.0)
    _grow_array("diameter", 1.0)
    _grow_array("image", image_flags)
    _grow_array("orientation", np.array([1.0, 0.0, 0.0, 0.0], dtype=float))
    _grow_array("angmom", np.zeros(4, dtype=float))
    _grow_array("body", -1)
    _grow_array("moment_inertia", np.zeros(3, dtype=float))

    snapshot.particles.N = new_n


def prod_trajectory_path(base: Path, source_replica: int) -> Path:
    """Path to a library replica production GSD trajectory."""
    return (
        base
        / f"replicas/replica_{source_replica}"
        / "no_cavity"
        / f"prod-{source_replica}.gsd"
    )


def validate_source_trajectories(
    base: Path,
    n_source_replicas: int,
    n_output_replicas: int,
) -> dict[str, object]:
    """Ensure every source stream has enough frames for the requested mapping."""
    required_frames = min_frames_required(n_output_replicas, n_source_replicas)
    summary: dict[str, object] = {
        "required_frames_per_source": required_frames,
        "sources": {},
    }

    for source_replica in range(n_source_replicas):
        path = prod_trajectory_path(base, source_replica)
        if not path.is_file():
            raise FileNotFoundError(f"missing source trajectory: {path}")
        with gsd.hoomd.open(path, "r") as traj:
            n_frames = len(traj)
            if n_frames < required_frames:
                raise RuntimeError(
                    f"{path.name} for replica_{source_replica} has {n_frames} frames; "
                    f"need at least {required_frames} for {n_output_replicas} outputs"
                )
            elapsed_ps = float(traj[-1].log["Time/elapsed_ps"][0])
            particles_n = int(traj[0].particles.N)
        summary["sources"][str(source_replica)] = {
            "path": str(path),
            "frames": n_frames,
            "elapsed_ps": elapsed_ps,
            "particles_n": particles_n,
        }

    return summary


def build_init_gsd(
    base: Path,
    output_path: Path,
    n_output_replicas: int = 500,
    n_source_replicas: int = 8,
    add_cavity: bool = True,
) -> dict[str, object]:
    """Write a multi-frame init GSD and return a reproducibility manifest."""
    source_summary = validate_source_trajectories(
        base=base,
        n_source_replicas=n_source_replicas,
        n_output_replicas=n_output_replicas,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, int | float]] = []

    with gsd.hoomd.open(output_path, "w") as out_traj:
        for output_replica in range(n_output_replicas):
            source_replica, source_frame = source_mapping(
                output_replica, n_source_replicas
            )
            src_path = prod_trajectory_path(base, source_replica)
            with gsd.hoomd.open(src_path, "r") as src_traj:
                frame = copy.deepcopy(src_traj[source_frame])
            if add_cavity:
                add_cavity_particle_at_origin(frame)
            out_traj.append(frame)
            entries.append(
                {
                    "output_replica": output_replica,
                    "source_replica": source_replica,
                    "source_frame": source_frame,
                    "elapsed_ps": float(frame.log["Time/elapsed_ps"][0]),
                }
            )

    manifest: dict[str, object] = {
        "output_path": str(output_path),
        "n_output_replicas": n_output_replicas,
        "n_source_replicas": n_source_replicas,
        "add_cavity": add_cavity,
        "mapping": "source_replica = output_replica % n_source; "
        "source_frame = output_replica // n_source",
        "sources": source_summary["sources"],
        "required_frames_per_source": source_summary["required_frames_per_source"],
        "entries": entries,
    }

    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="IC library root (default: aging_weak_lambda_ic_library/)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output init GSD (default: <base>/init-500-from-ic-library.gsd)",
    )
    parser.add_argument("--n-output-replicas", type=int, default=500)
    parser.add_argument("--n-source-replicas", type=int, default=8)
    parser.add_argument(
        "--no-cavity",
        action="store_true",
        help="Keep molecular-only frames (production adds L at runtime)",
    )
    args = parser.parse_args()

    output = args.output or (args.base / "init-500-from-ic-library.gsd")
    manifest = build_init_gsd(
        base=args.base,
        output_path=output,
        n_output_replicas=args.n_output_replicas,
        n_source_replicas=args.n_source_replicas,
        add_cavity=not args.no_cavity,
    )
    print(f"wrote {output} ({manifest['n_output_replicas']} frames)")
    print(f"manifest: {output.with_suffix('.manifest.json')}")
    for source_replica, info in manifest["sources"].items():
        print(
            f"  source replica_{source_replica}: frames={info['frames']} "
            f"elapsed={info['elapsed_ps']:.1f} ps"
        )


if __name__ == "__main__":
    main()
