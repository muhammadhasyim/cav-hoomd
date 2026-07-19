#!/usr/bin/env python3
"""Write a molecular-only GSD by removing cavity (L) particles from each frame."""
from __future__ import annotations

import argparse
from pathlib import Path

import gsd.hoomd
import numpy as np


def strip_frame(frame):
    typeid = np.array(frame.particles.typeid)
    if "L" not in frame.particles.types:
        return frame
    l_id = list(frame.particles.types).index("L")
    keep = typeid != l_id
    if keep.all():
        return frame
    snap = gsd.hoomd.Frame()
    snap.configuration = frame.configuration
    snap.particles.types = list(frame.particles.types)
    snap.particles.N = int(np.sum(keep))
    for attr in (
        "typeid", "position", "charge", "mass", "diameter", "image",
        "velocity", "orientation", "angmom", "body", "moment_inertia",
    ):
        if hasattr(frame.particles, attr):
            val = getattr(frame.particles, attr)
            setattr(snap.particles, attr, val[keep])
    if hasattr(frame, "bonds") and frame.bonds.N:
        snap.bonds = frame.bonds
    return snap


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()

    with gsd.hoomd.open(args.input, "r") as fin:
        with gsd.hoomd.open(args.output, "w") as fout:
            for i, frame in enumerate(fin):
                out = strip_frame(frame)
                fout.append(out)
                if i == 0:
                    print(
                        f"frame0: N={out.particles.N} types={list(out.particles.types)}"
                    )
            print(f"Wrote {len(fin)} frames to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
