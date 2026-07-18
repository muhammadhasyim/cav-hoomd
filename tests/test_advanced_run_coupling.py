"""Regression tests for the aging campaign coupling switch protocol."""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from hoomd.cavitymd.variants.basic import StepVariant


advanced_run = importlib.import_module("examples.05_advanced_run")


class FakeTimeTracker:
    """Minimal mutable material-time source for StepVariant evaluation."""

    elapsed_time = 0.0


def test_switch_time_defaults_advanced_runner_to_step_coupling() -> None:
    assert (
        advanced_run.resolve_coupling_variant_type(
            switch_time_ps=200.0,
            requested_type=None,
        )
        == "step"
    )


def test_no_switch_time_defaults_advanced_runner_to_constant_coupling() -> None:
    assert (
        advanced_run.resolve_coupling_variant_type(
            switch_time_ps=None,
            requested_type=None,
        )
        == "constant"
    )


def test_constant_coupling_is_rejected_when_switch_time_is_set() -> None:
    with pytest.raises(ValueError, match="requires step"):
        advanced_run.resolve_coupling_variant_type(
            switch_time_ps=200.0,
            requested_type="constant",
        )


def test_step_variant_is_zero_before_switch_and_target_afterward() -> None:
    tracker = FakeTimeTracker()
    coupling = StepVariant(
        target_value=0.03,
        switch_time_ps=200.0,
        time_tracker=tracker,
    )

    tracker.elapsed_time = 199.999
    assert coupling(0) == 0.0
    tracker.elapsed_time = 200.0
    assert coupling(0) == pytest.approx(0.03)


def test_run_single_experiment_passes_step_variant_to_simulation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    class FakeSimulation:
        def __init__(self, **kwargs: object) -> None:
            captured.update(kwargs)

        def run(self) -> int:
            return 0

    monkeypatch.setattr(advanced_run, "CavityMDSimulation", FakeSimulation)
    monkeypatch.chdir(tmp_path)

    success = advanced_run.run_single_experiment(
        molecular_thermo="bussi",
        cavity_thermo="langevin",
        finite_q=False,
        coupling=0.03,
        temperature=100.0,
        frequency=1560.0,
        replica=4,
        frame=-1,
        runtime_ps=2500.0,
        molecular_tau=5.0,
        cavity_tau=5.0,
        enable_fkt=True,
        fkt_kmag=1.0085,
        fkt_wavevectors=50,
        fkt_ref_interval=200.0,
        fkt_max_refs=13,
        switch_time_ps=200.0,
        input_gsd=str(tmp_path / "init-0.gsd"),
    )

    assert success
    assert captured["coupling_variant_type"] == "step"
    assert captured["switch_time_ps"] == 200.0
    assert captured["lambda_coupling"] == 0.03
    assert captured.get("couplstr") is None
