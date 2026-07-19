#!/usr/bin/env python3
"""Refresh IC checkpoints from in-progress prod GSD trajectories."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import gsd.hoomd


def refresh_checkpoints(
    base: Path,
    n_replicas: int,
    target_runtime_ps: float,
) -> dict[str, object]:
    """Copy the last prod frame to ic_checkpoints and compute remaining runtime."""
    replicas_root = base / "replicas"
    ic_dir = base / "ic_checkpoints"
    ic_dir.mkdir(parents=True, exist_ok=True)

    summary: dict[str, object] = {
        "target_runtime_ps": target_runtime_ps,
        "replicas": {},
    }

    for replica in range(n_replicas):
        prod_path = replicas_root / f"replica_{replica}" / "no_cavity" / f"prod-{replica}.gsd"
        ic_path = ic_dir / f"replica_{replica}.gsd"
        entry: dict[str, object] = {
            "prod_path": str(prod_path),
            "ic_path": str(ic_path),
        }

        if prod_path.is_file():
            with gsd.hoomd.open(prod_path, "r") as prod_traj:
                if len(prod_traj) == 0:
                    raise RuntimeError(f"{prod_path} exists but has zero frames")
                last_frame = prod_traj[-1]
                elapsed_ps = float(last_frame.log["Time/elapsed_ps"][0])
                n_frames = len(prod_traj)
            with gsd.hoomd.open(ic_path, "w") as ic_traj:
                ic_traj.append(last_frame)
            entry.update(
                {
                    "source": "prod",
                    "frames": n_frames,
                    "elapsed_ps": elapsed_ps,
                    "remaining_ps": max(0.0, target_runtime_ps - elapsed_ps),
                }
            )
        elif ic_path.is_file():
            with gsd.hoomd.open(ic_path, "r") as ic_traj:
                last_frame = ic_traj[-1]
                n_frames = len(ic_traj)
                elapsed_ps = float(last_frame.log["Time/elapsed_ps"][0])
            entry.update(
                {
                    "source": "existing_ic",
                    "frames": n_frames,
                    "elapsed_ps": elapsed_ps,
                    "remaining_ps": max(0.0, target_runtime_ps - elapsed_ps),
                }
            )
        else:
            raise FileNotFoundError(
                f"no prod or ic checkpoint found for replica {replica}"
            )

        summary["replicas"][str(replica)] = entry

    manifest_path = base / "logs" / "resume_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    parser.add_argument("--n-replicas", type=int, default=8)
    parser.add_argument("--target-runtime-ps", type=float, default=100000.0)
    args = parser.parse_args()

    summary = refresh_checkpoints(
        base=args.base,
        n_replicas=args.n_replicas,
        target_runtime_ps=args.target_runtime_ps,
    )
    for replica, entry in summary["replicas"].items():
        print(
            f"replica_{replica}: source={entry['source']} "
            f"frames={entry.get('frames')} elapsed={entry['elapsed_ps']:.1f} ps "
            f"remaining={entry['remaining_ps']:.1f} ps"
        )


if __name__ == "__main__":
    main()
