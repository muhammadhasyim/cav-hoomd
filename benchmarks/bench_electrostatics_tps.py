#!/usr/bin/env python3
"""Benchmark physics-only throughput for reaction field vs PPPM via 05_advanced_run."""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "examples" / "05_advanced_run.py"
GSD = REPO / "square_lambda0.025_diffeq" / "prod-0.gsd"


def bench(electrostatics: str, runtime_ps: float = 30.0, eps_rf: float = 0.0,
          coulomb_rcut: float = 15.0) -> tuple[float, float]:
    out = Path(f"benchmarks/results/tps_{electrostatics}")
    out.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(RUNNER),
        "--molecular-bath", "bussi", "--cavity-bath", "langevin",
        "--lambda-coupling", "0.03", "--coupling-type", "step",
        "--temperature", "100.0", "--frequency", "1560.0",
        "--runtime", str(runtime_ps), "--switch-time", "200.0",
        "--input-gsd", str(GSD), "--frame", "-1",
        "--device", "GPU", "--gpu-id", "0",
        "--fixed-timestep", "--timestep", "1.0",
        "--disable-gsd", "--disable-hdf5-output", "--disable-temp-tracker",
        "--console-output-period-ps", "100.0",
        "--replicas", "0", "--seed", "1",
        "--electrostatics", electrostatics,
        "--eps-rf", str(eps_rf),
        "--coulomb-rcut", str(coulomb_rcut),
    ]
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env["PYTHONPATH"] = f"/opt/hoomd-md:{env.get('PYTHONPATH', '')}"
    log = out / "run.log"
    subprocess.run(cmd, cwd=out, env=env, check=True,
                   stdout=log.open("w"), stderr=subprocess.STDOUT)
    text = log.read_text(encoding="utf-8")
    rows = re.findall(
        r"^\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", text, re.MULTILINE)
    if not rows:
        wt = re.search(r"Wall time: ([\d.]+) seconds", text)
        if wt:
            steps = runtime_ps * 1000
            tps = steps / float(wt.group(1))
            return tps, tps * 0.0864
        raise RuntimeError(f"No timing data in {log}")
    tps = float(rows[-1][1])
    ns_day = float(rows[-1][3])
    return tps, ns_day


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runtime-ps", type=float, default=30.0)
    args = ap.parse_args()
    results = {}
    for label, method, eps, rc in [
        ("pppm", "pppm", 0.0, 15.0),
        ("reaction_field", "reaction_field", 0.0, 15.0),
    ]:
        tps, ns = bench(method, args.runtime_ps, eps, rc)
        results[label] = {"steps_per_s": tps, "ns_per_day": ns}
        print(f"{label}: {tps:.1f} steps/s, {ns:.1f} ns/day")
    rf_tps = results["reaction_field"]["steps_per_s"]
    print(f"\nTarget 2500 steps/s: {'PASS' if rf_tps >= 2500 else 'FAIL'} ({rf_tps:.0f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
