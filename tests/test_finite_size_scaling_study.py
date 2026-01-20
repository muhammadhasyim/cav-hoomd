"""Tests for finite_size_scaling_study CLI configuration."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_module():
    module_path = Path(__file__).resolve().parents[1] / "examples" / "finite_size_scaling_study.py"
    spec = importlib.util.spec_from_file_location("finite_size_scaling_study", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_build_parser_defaults():
    module = _load_module()
    parser = module.build_parser()
    args = parser.parse_args([])

    assert args.collective_coupling == 0.02
    assert args.n_values == "100,250,500,1000,2500,5000,10000"
    assert args.device == "GPU"
    assert args.molecular_thermostat == "bussi"
    assert args.cavity_thermostat == "langevin"
    assert args.molecular_thermostat_tau == 5.0
    assert args.cavity_thermostat_tau == 5.0
    assert args.gsd_output_period_ps == 1.0
    assert args.energy_output_period_ps == 0.1
    assert args.enable_hdf5_output is True
    assert args.enable_dipole_autocorr is True

