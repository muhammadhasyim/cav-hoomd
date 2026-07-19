"""Tests for packed manifest splitting."""

from __future__ import annotations

from pathlib import Path

from examples.slurm.split_packed_manifest import split_manifest, write_part_sbatch


def test_split_manifest_divides_evenly(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(f"0p03\t0.03\t{i}\t" for i in range(4)) + "\n",
        encoding="utf-8",
    )
    part1 = tmp_path / "manifest_part1.tsv"
    part2 = tmp_path / "manifest_part2.tsv"

    n1, n2 = split_manifest(manifest, part1, part2)

    assert n1 == 2
    assert n2 == 2
    assert len(part1.read_text(encoding="utf-8").splitlines()) == 2
    assert len(part2.read_text(encoding="utf-8").splitlines()) == 2


def test_write_part_sbatch_rewrites_manifest_and_array(tmp_path: Path) -> None:
    template = tmp_path / "job.sbatch"
    template.write_text(
        "\n".join(
            [
                "#!/bin/bash",
                "#SBATCH --job-name=base-job",
                "#SBATCH --array=0-9%24",
                "python run.py \\",
                "    --manifest /path/manifest.tsv \\",
                "    --task-index ${SLURM_ARRAY_TASK_ID}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "manifest_part1.tsv"
    manifest.write_text("0p03\t0.03\t4\t\n", encoding="utf-8")
    output = tmp_path / "job_part1.sbatch"

    write_part_sbatch(
        template=template,
        output=output,
        manifest=manifest,
        job_name="base-job-part1",
        array_end=0,
    )

    text = output.read_text(encoding="utf-8")
    assert "#SBATCH --job-name=base-job-part1" in text
    assert "#SBATCH --array=0-0%24" in text
    assert f"--manifest {manifest} \\" in text
