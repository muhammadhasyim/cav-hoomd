"""Regression scaffold for optimizing ``CavityForceComputeGPU``."""

from pathlib import Path

import hoomd
import numpy as np
import pytest

from hoomd.cavitymd import _cavitymd
from hoomd.cavitymd.forces import CavityForce


PROJECT_ROOT = Path(__file__).resolve().parents[2]
GPU_SOURCE = PROJECT_ROOT / "src" / "CavityForceComputeGPU.cc"
GPU_KERNEL_SOURCE = PROJECT_ROOT / "src" / "CavityForceComputeGPU.cu"
GPU_HEADER = PROJECT_ROOT / "src" / "CavityForceComputeGPU.h"
PARAMS_HEADER = PROJECT_ROOT / "src" / "CavityForceParams.h"


@pytest.fixture(scope="module")
def gpu_device():
    """Return an active HOOMD GPU device or skip when none is usable."""
    try:
        device = hoomd.device.GPU(notice_level=0)
    except RuntimeError as error:
        pytest.skip(f"HOOMD GPU unavailable: {error}")

    if not device._cpp_exec_conf.isCUDAEnabled():
        pytest.skip("HOOMD did not activate a GPU execution configuration")
    return device


def _function_body(source: str, signature: str) -> str:
    """Extract a C++ function body using balanced braces."""
    signature_start = source.index(signature)
    body_start = source.index("{", signature_start)
    depth = 0
    for index in range(body_start, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[body_start:index + 1]
    raise AssertionError(f"Unterminated function body for {signature}")


def _make_simulation(device):
    """Create a deterministic O/N molecular system with one L photon."""
    snapshot = hoomd.Snapshot(device.communicator)
    if snapshot.communicator.rank == 0:
        snapshot.configuration.box = [12.0, 12.0, 12.0, 0.0, 0.0, 0.0]
        snapshot.particles.N = 5
        snapshot.particles.types = ["O", "N", "L"]
        snapshot.particles.typeid[:] = [0, 1, 0, 1, 2]
        snapshot.particles.position[:] = [
            [3.0, 2.0, 0.0],
            [-2.0, -1.0, 0.5],
            [1.0, -3.0, -0.5],
            [-3.0, 3.0, 1.0],
            [0.4, -0.6, 0.8],
        ]
        snapshot.particles.charge[:] = [-0.70, 0.35, 0.20, 0.15, 0.0]
        snapshot.particles.mass[:] = [16.0, 14.0, 16.0, 14.0, 1.0]

    simulation = hoomd.Simulation(device=device, seed=42)
    simulation.create_state_from_snapshot(snapshot)
    cavity_force = CavityForce(
        kvector=[0.0, 0.0, 1.0],
        lambda_coupling=0.03,
        omegac=0.011,
        phmass=1.7,
    )
    simulation.operations.integrator = hoomd.md.Integrator(
        dt=1.0e-8,
        methods=[hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())],
        forces=[cavity_force],
    )
    simulation.operations.tuners.append(
        hoomd.tune.ParticleSorter(
            trigger=hoomd.trigger.Periodic(1),
            grid=8,
        )
    )
    simulation.run(0)
    return simulation, cavity_force


def _read_force_result(simulation, cavity_force):
    """Read force vectors by stable particle tag and all cavity energies."""
    with simulation.state.cpu_local_snapshot as snapshot:
        tags = np.array(snapshot.particles.tag, copy=True)
    with cavity_force.cpu_local_force_arrays as arrays:
        local_forces = np.array(arrays.force, copy=True)

    forces_by_tag = local_forces[np.argsort(tags)]
    energies = np.array(
        [
            cavity_force.harmonic_energy,
            cavity_force.coupling_energy,
            cavity_force.dipole_self_energy,
        ]
    )
    return tags, forces_by_tag, energies


def test_gpu_hot_path_avoids_host_work():
    """Require the optimized GPU hot path to remain entirely device-side."""
    gpu_source = GPU_SOURCE.read_text()
    kernel_source = GPU_KERNEL_SOURCE.read_text()
    compute_body = _function_body(
        gpu_source,
        "void CavityForceComputeGPU::computeForces",
    )
    sort_callback_body = _function_body(
        gpu_source,
        "void CavityForceComputeGPU::updatePhotonIndex",
    )
    launcher_body = _function_body(
        kernel_source,
        "hipError_t gpu_compute_cavity_forces",
    )

    forbidden_patterns = {
        "explicit device synchronization": (
            compute_body,
            "hipDeviceSynchronize",
        ),
        "host access to GPU arrays": (
            compute_body,
            "access_location::host",
        ),
        "host access in the particle-sort callback": (
            sort_callback_body,
            "access_location::host",
        ),
        "explicit synchronization in the particle-sort callback": (
            sort_callback_body,
            "hipDeviceSynchronize",
        ),
        "host-to-device copy in the hot path": (
            compute_body + sort_callback_body + launcher_body,
            "hipMemcpyHostToDevice",
        ),
        "device-to-host copy in the hot path": (
            sort_callback_body + launcher_body,
            "hipMemcpyDeviceToHost",
        ),
        "full photon search on every evaluation": (
            launcher_body,
            "gpu_find_photon_kernel",
        ),
        "non-deterministic atomic dipole accumulation": (
            kernel_source,
            "atomicAdd",
        ),
    }
    violations = [
        description
        for description, (body, pattern) in forbidden_patterns.items()
        if pattern in body
    ]

    assert not violations, (
        "GPU CavityForce hot path still performs forbidden work: "
        + ", ".join(violations)
    )
    assert PARAMS_HEADER.is_file()
    assert "struct cavity_force_params" not in kernel_source
    assert "Autotuner" not in GPU_HEADER.read_text()
    constructor_body = _function_body(
        gpu_source,
        "CavityForceComputeGPU::CavityForceComputeGPU",
    )
    assert "getNRanks() > 1" in constructor_body
    assert "single-rank" in constructor_body
    assert "getGlobalParticleNumberChangeSignal()" in gpu_source
    assert "handleParticleNumberChange" in gpu_source


def test_gpu_matches_cpu_before_and_after_particle_reorder(gpu_device):
    """Compare forces and three energies before and after a real reorder."""
    cpu_simulation, cpu_force = _make_simulation(
        hoomd.device.CPU(notice_level=0)
    )
    gpu_simulation, gpu_force = _make_simulation(gpu_device)

    assert cpu_force.implementation == "cpp"
    assert gpu_force.implementation == "cuda"

    cpu_tags_before, cpu_forces_before, cpu_energies_before = (
        _read_force_result(cpu_simulation, cpu_force)
    )
    gpu_tags_before, gpu_forces_before, gpu_energies_before = (
        _read_force_result(gpu_simulation, gpu_force)
    )

    np.testing.assert_allclose(
        gpu_forces_before,
        cpu_forces_before,
        rtol=1.0e-6,
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        gpu_energies_before,
        cpu_energies_before,
        rtol=1.0e-10,
        atol=1.0e-12,
    )

    cpu_simulation.run(1)
    gpu_simulation.run(1)
    cpu_tags_after, cpu_forces_after, cpu_energies_after = _read_force_result(
        cpu_simulation,
        cpu_force,
    )
    gpu_tags_after, gpu_forces_after, gpu_energies_after = _read_force_result(
        gpu_simulation,
        gpu_force,
    )

    assert not np.array_equal(cpu_tags_after, cpu_tags_before)
    assert not np.array_equal(gpu_tags_after, gpu_tags_before)
    np.testing.assert_allclose(
        gpu_forces_after,
        cpu_forces_after,
        rtol=1.0e-6,
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        gpu_energies_after,
        cpu_energies_after,
        rtol=1.0e-10,
        atol=1.0e-12,
    )


def test_force_matches_dse_derivative_for_nonunit_phmass(gpu_device):
    """Check the analytic molecular force for a non-unit photon mass."""
    positions = np.array(
        [
            [3.0, 2.0, 0.0],
            [-2.0, -1.0, 0.5],
            [1.0, -3.0, -0.5],
            [-3.0, 3.0, 1.0],
            [0.4, -0.6, 0.8],
        ]
    )
    charges = np.array([-0.70, 0.35, 0.20, 0.15, 0.0])
    lambda_coupling = 0.03
    omegac = 0.011
    phmass = 1.7
    dipole = np.sum(charges[:4, None] * positions[:4], axis=0)
    coupling = lambda_coupling * omegac
    displacement = positions[-1, :2] + (
        lambda_coupling / (phmass * omegac)
    ) * dipole[:2]
    expected_forces = np.zeros((5, 3))
    expected_forces[:4, :2] = (
        -coupling * charges[:4, None] * displacement[None, :]
    )
    spring_constant = phmass * omegac**2
    expected_forces[-1] = (
        -spring_constant * positions[-1]
        - coupling * np.array([dipole[0], dipole[1], 0.0])
    )

    for device in [hoomd.device.CPU(notice_level=0), gpu_device]:
        simulation, cavity_force = _make_simulation(device)
        _, actual_forces, _ = _read_force_result(simulation, cavity_force)
        np.testing.assert_allclose(
            actual_forces,
            expected_forces,
            rtol=2.0e-10,
            atol=2.0e-12,
        )


def test_gpu_large_non_aligned_reduction_matches_analytic_reference(
    gpu_device,
):
    """Validate deterministic reduction when the first pass exceeds 256 blocks."""
    molecule_count = 65_537
    particle_count = molecule_count + 1
    molecule_index = np.arange(molecule_count)
    positions = np.zeros((particle_count, 3), dtype=np.float64)
    positions[:molecule_count, 0] = (
        molecule_index % 257 - 128
    ) * 0.01
    positions[:molecule_count, 1] = (
        (molecule_index // 257) % 257 - 128
    ) * 0.008
    positions[-1] = [0.4, -0.6, 0.8]
    charges = np.zeros(particle_count, dtype=np.float64)
    charges[:molecule_count] = (
        molecule_index % 7 - 3
    ) * 1.0e-4

    snapshot = hoomd.Snapshot(gpu_device.communicator)
    if snapshot.communicator.rank == 0:
        snapshot.configuration.box = [
            20.0,
            20.0,
            20.0,
            0.0,
            0.0,
            0.0,
        ]
        snapshot.particles.N = particle_count
        snapshot.particles.types = ["O", "L"]
        snapshot.particles.typeid[:] = 0
        snapshot.particles.typeid[-1] = 1
        snapshot.particles.position[:] = positions
        snapshot.particles.charge[:] = charges
        snapshot.particles.mass[:] = 1.0

    lambda_coupling = 0.03
    omegac = 0.011
    phmass = 1.7
    simulation = hoomd.Simulation(device=gpu_device, seed=43)
    simulation.create_state_from_snapshot(snapshot)
    cavity_force = CavityForce(
        kvector=[0.0, 0.0, 1.0],
        lambda_coupling=lambda_coupling,
        omegac=omegac,
        phmass=phmass,
    )
    simulation.operations.integrator = hoomd.md.Integrator(
        dt=1.0e-8,
        methods=[hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())],
        forces=[cavity_force],
    )
    simulation.run(0)

    _, forces_by_tag, gpu_energies = _read_force_result(
        simulation,
        cavity_force,
    )
    dipole = np.sum(
        charges[:molecule_count, None] * positions[:molecule_count],
        axis=0,
        dtype=np.float64,
    )
    spring_constant = phmass * omegac**2
    coupling = lambda_coupling * omegac
    photon_coordinate = positions[-1]
    displacement = photon_coordinate[:2] + (
        lambda_coupling / (phmass * omegac)
    ) * dipole[:2]

    selected_tags = np.array([0, 1, 255, 256, 65_536])
    expected_forces = np.zeros((selected_tags.size + 1, 3))
    expected_forces[:-1, :2] = (
        -coupling
        * charges[selected_tags, None]
        * displacement[None, :]
    )
    expected_forces[-1] = (
        -spring_constant * photon_coordinate
        - coupling * np.array([dipole[0], dipole[1], 0.0])
    )
    actual_forces = np.vstack(
        [forces_by_tag[selected_tags], forces_by_tag[-1]]
    )
    expected_energies = np.array(
        [
            0.5 * spring_constant * np.dot(
                photon_coordinate,
                photon_coordinate,
            ),
            coupling * np.dot(dipole[:2], photon_coordinate[:2]),
            0.5
            * coupling**2
            / spring_constant
            * np.dot(dipole[:2], dipole[:2]),
        ]
    )

    np.testing.assert_allclose(
        actual_forces,
        expected_forces,
        rtol=2.0e-8,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        gpu_energies,
        expected_energies,
        rtol=2.0e-10,
        atol=2.0e-12,
    )


@pytest.mark.parametrize("photon_count", [0, 2])
def test_gpu_requires_exactly_one_photon(gpu_device, photon_count):
    """GPU construction rejects systems without exactly one L photon."""
    snapshot = hoomd.Snapshot(gpu_device.communicator)
    if snapshot.communicator.rank == 0:
        snapshot.configuration.box = [8.0, 8.0, 8.0, 0.0, 0.0, 0.0]
        snapshot.particles.N = 3
        snapshot.particles.types = ["O", "L"]
        snapshot.particles.typeid[:] = 0
        snapshot.particles.typeid[:photon_count] = 1
        snapshot.particles.position[:] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
        ]
        snapshot.particles.charge[:] = [-0.5, 0.25, 0.25]
        snapshot.particles.mass[:] = 1.0

    simulation = hoomd.Simulation(device=gpu_device, seed=44)
    simulation.create_state_from_snapshot(snapshot)
    with pytest.raises(
        RuntimeError,
        match=rf"requires exactly one.*found {photon_count}",
    ):
        _cavitymd.CavityForceComputeGPU(
            simulation.state._cpp_sys_def,
            0.011,
            hoomd.variant.Constant(0.03),
            1.7,
        )


@pytest.mark.parametrize("photon_count", [0, 2])
def test_public_gpu_force_propagates_photon_validation(
    gpu_device,
    photon_count,
):
    """Public GPU attachment must not hide photon validation failures."""
    snapshot = hoomd.Snapshot(gpu_device.communicator)
    if snapshot.communicator.rank == 0:
        snapshot.configuration.box = [8.0, 8.0, 8.0, 0.0, 0.0, 0.0]
        snapshot.particles.N = 3
        snapshot.particles.types = ["O", "L"]
        snapshot.particles.typeid[:] = 0
        snapshot.particles.typeid[:photon_count] = 1
        snapshot.particles.position[:] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
        ]
        snapshot.particles.charge[:] = [-0.5, 0.25, 0.25]
        snapshot.particles.mass[:] = 1.0

    simulation = hoomd.Simulation(device=gpu_device, seed=45)
    simulation.create_state_from_snapshot(snapshot)
    cavity_force = CavityForce(
        kvector=[0.0, 0.0, 1.0],
        lambda_coupling=0.03,
        omegac=0.011,
        phmass=1.7,
    )
    simulation.operations.integrator = hoomd.md.Integrator(
        dt=1.0e-8,
        methods=[hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())],
        forces=[cavity_force],
    )

    with pytest.raises(
        RuntimeError,
        match=rf"requires exactly one.*found {photon_count}",
    ):
        simulation.run(0)
    assert cavity_force.implementation == "cuda"
    assert cavity_force._force_impl is None


def _make_cpu_system_definition():
    """Build a minimal CPU system definition for direct constructor tests."""
    device = hoomd.device.CPU(notice_level=0)
    snapshot = hoomd.Snapshot(device.communicator)
    if snapshot.communicator.rank == 0:
        snapshot.configuration.box = [6.0, 6.0, 6.0, 0.0, 0.0, 0.0]
        snapshot.particles.N = 1
        snapshot.particles.types = ["L"]
        snapshot.particles.typeid[0] = 0
        snapshot.particles.mass[0] = 1.0
    simulation = hoomd.Simulation(device=device, seed=46)
    simulation.create_state_from_snapshot(snapshot)
    return simulation.state._cpp_sys_def


@pytest.mark.parametrize(
    ("omegac", "phmass", "message"),
    [
        (0.0, 1.0, "omegac"),
        (-0.01, 1.0, "omegac"),
        (np.nan, 1.0, "omegac"),
        (np.inf, 1.0, "omegac"),
        (0.011, 0.0, "phmass"),
        (0.011, -1.0, "phmass"),
        (0.011, np.nan, "phmass"),
        (0.011, np.inf, "phmass"),
    ],
)
def test_shared_constructor_rejects_invalid_physical_parameters(
    omegac,
    phmass,
    message,
):
    """CPU and GPU base construction rejects nonphysical parameters."""
    with pytest.raises(RuntimeError, match=message):
        _cavitymd.CavityForceCompute(
            _make_cpu_system_definition(),
            omegac,
            hoomd.variant.Constant(0.03),
            phmass,
        )


@pytest.mark.parametrize("lambda_coupling", [np.nan, np.inf, -np.inf])
def test_shared_constructor_rejects_nonfinite_lambda(lambda_coupling):
    """The sampled lambda variant value must be finite."""
    with pytest.raises(RuntimeError, match="lambda"):
        _cavitymd.CavityForceCompute(
            _make_cpu_system_definition(),
            0.011,
            hoomd.variant.Constant(lambda_coupling),
            1.7,
        )


def test_set_params_validates_and_recomputes_spring_constant():
    """Parameter updates validate first and keep K consistent."""
    force = _cavitymd.CavityForceCompute(
        _make_cpu_system_definition(),
        0.011,
        hoomd.variant.Constant(0.03),
        1.7,
    )
    force.setParams(0.02, hoomd.variant.Constant(0.04), 2.5)
    params = force.getParams()
    assert params["omegac"] == pytest.approx(0.02)
    assert params["lambda_coupling"] == pytest.approx(0.04)
    assert params["phmass"] == pytest.approx(2.5)
    assert params["K"] == pytest.approx(2.5 * 0.02**2)

    for omegac, coupling, phmass, message in [
        (0.0, 0.04, 2.5, "omegac"),
        (0.02, np.nan, 2.5, "lambda"),
        (0.02, 0.04, 0.0, "phmass"),
    ]:
        with pytest.raises(RuntimeError, match=message):
            force.setParams(
                omegac,
                hoomd.variant.Constant(coupling),
                phmass,
            )


def test_public_set_params_updates_numeric_and_variant_lambda():
    """Public updates keep Python and active native lambda state aligned."""
    _, cavity_force = _make_simulation(hoomd.device.CPU(notice_level=0))

    cavity_force.set_params(lambda_coupling=0.07)
    numeric_variant = cavity_force.lambda_coupling_variant
    assert isinstance(numeric_variant, hoomd.variant.Constant)
    assert numeric_variant.value == pytest.approx(0.07)
    assert cavity_force._param_dict["lambda_coupling"] is numeric_variant
    assert cavity_force._force_impl.getParams()[
        "lambda_coupling"
    ] == pytest.approx(0.07)

    ramp = hoomd.variant.Ramp(
        A=0.08,
        B=0.10,
        t_start=0,
        t_ramp=10,
    )
    cavity_force.set_params(lambda_coupling=ramp)
    assert cavity_force.lambda_coupling_variant is ramp
    assert cavity_force._param_dict["lambda_coupling"] is ramp
    assert cavity_force.uses_variant_coupling
    assert cavity_force._force_impl.getParams()[
        "lambda_coupling"
    ] == pytest.approx(0.08)


def _public_parameter_state(cavity_force):
    """Capture Python and native parameter state for rollback assertions."""
    return {
        "omegac": cavity_force.omegac,
        "phmass": cavity_force.phmass,
        "kvector_object": cavity_force._param_dict._dict["kvector"],
        "lambda_variant": cavity_force.lambda_coupling_variant,
        "dict_omegac": cavity_force._param_dict["omegac"],
        "dict_phmass": cavity_force._param_dict["phmass"],
        "dict_lambda": cavity_force._param_dict["lambda_coupling"],
        "native": dict(cavity_force._force_impl.getParams()),
    }


def _assert_public_parameter_state(cavity_force, expected):
    """Assert that a rejected update changed no public or native state."""
    assert cavity_force.omegac == expected["omegac"]
    assert cavity_force.phmass == expected["phmass"]
    assert (
        cavity_force._param_dict._dict["kvector"]
        is expected["kvector_object"]
    )
    assert cavity_force.lambda_coupling_variant is expected["lambda_variant"]
    assert cavity_force._param_dict["omegac"] == expected["dict_omegac"]
    assert cavity_force._param_dict["phmass"] == expected["dict_phmass"]
    assert (
        cavity_force._param_dict["lambda_coupling"]
        is expected["dict_lambda"]
    )
    assert dict(cavity_force._force_impl.getParams()) == expected["native"]


def test_public_set_params_is_transactional_on_rejection():
    """Rejected public updates preserve Python and native parameter state."""
    _, cavity_force = _make_simulation(hoomd.device.CPU(notice_level=0))
    expected = _public_parameter_state(cavity_force)

    rejected_updates = [
        ({"omegac": 0.0}, "omegac"),
        ({"phmass": -1.0}, "phmass"),
        ({"lambda_coupling": np.nan}, "lambda"),
        ({"omegac": 0.02, "unknown_parameter": 1.0}, "Unknown"),
    ]
    for update, message in rejected_updates:
        with pytest.raises((RuntimeError, ValueError), match=message):
            cavity_force.set_params(**update)
        _assert_public_parameter_state(cavity_force, expected)
