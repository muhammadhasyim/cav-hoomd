#!/usr/bin/env python3
"""Split a packed manifest and emit part1/part2 SLURM scripts."""

from __future__ import annotations

import argparse
from pathlib import Path


def split_manifest(manifest: Path, part1: Path, part2: Path) -> tuple[int, int]:
    """Write first and second halves of a manifest TSV."""
    lines = manifest.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"empty manifest: {manifest}")
    midpoint = len(lines) // 2
    part1.write_text("\n".join(lines[:midpoint]) + "\n", encoding="utf-8")
    part2.write_text("\n".join(lines[midpoint:]) + "\n", encoding="utf-8")
    return len(lines[:midpoint]), len(lines[midpoint:])


def write_part_sbatch(
    *,
    template: Path,
    output: Path,
    manifest: Path,
    job_name: str,
    array_end: int,
) -> None:
    """Clone a packed sbatch file with a new manifest path and array range."""
    text = template.read_text(encoding="utf-8")
    lines = text.splitlines()
    out_lines: list[str] = []
    for line in lines:
        if line.startswith("#SBATCH --job-name="):
            out_lines.append(f"#SBATCH --job-name={job_name}")
        elif line.startswith("#SBATCH --array="):
            concurrent = line.split("%", 1)[1].rstrip("]")
            out_lines.append(f"#SBATCH --array=0-{array_end}%{concurrent}")
        elif "--manifest " in line:
            out_lines.append(f"    --manifest {manifest} \\")
        else:
            out_lines.append(line)
    output.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    output.chmod(0o755)


def _job_name_from_template(template: Path) -> str:
    """Extract ``#SBATCH --job-name=...`` from a packed sbatch template."""
    for line in template.read_text(encoding="utf-8").splitlines():
        if line.startswith("#SBATCH --job-name="):
            return line.split("=", 1)[1]
    raise ValueError(f"no job name in {template}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--part1-name", default="job_part1.sbatch")
    parser.add_argument("--part2-name", default="job_part2.sbatch")
    parser.add_argument(
        "--part1-job-name",
        default=None,
        help="Override part1 SLURM job name (default: <template>-part1)",
    )
    parser.add_argument(
        "--part2-job-name",
        default=None,
        help="Override part2 SLURM job name (default: <template>-part2)",
    )
    args = parser.parse_args()

    plan_dir = args.plan_dir
    manifest = plan_dir / "manifest.tsv"
    template = plan_dir / "job.sbatch"
    if not manifest.is_file() or not template.is_file():
        raise SystemExit(f"missing manifest or job.sbatch in {plan_dir}")

    base_job_name = _job_name_from_template(template)
    part1_job_name = args.part1_job_name or f"{base_job_name}-part1"
    part2_job_name = args.part2_job_name or f"{base_job_name}-part2"

    n1, n2 = split_manifest(
        manifest,
        plan_dir / "manifest_part1.tsv",
        plan_dir / "manifest_part2.tsv",
    )
    write_part_sbatch(
        template=template,
        output=plan_dir / args.part1_name,
        manifest=plan_dir / "manifest_part1.tsv",
        job_name=part1_job_name,
        array_end=n1 - 1,
    )
    write_part_sbatch(
        template=template,
        output=plan_dir / args.part2_name,
        manifest=plan_dir / "manifest_part2.tsv",
        job_name=part2_job_name,
        array_end=n2 - 1,
    )
    print(f"split manifest: part1={n1} part2={n2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
