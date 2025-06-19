# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Test the BussiReservoir thermostat."""

import pytest
import hoomd
import numpy as np


def test_bussi_reservoir_basic():
    """Test basic functionality of BussiReservoir thermostat."""
    # Create a simple simulation
    device = hoomd.device.CPU()
    simulation = hoomd.Simulation(device=device, seed=42)
    
    # Create a simple system
    L = 10
    N = 100
    snapshot = hoomd.Snapshot(device.communicator)
    if snapshot.communicator.rank == 0:
        box = [L, L, L, 0, 0, 0]
        snapshot.configuration.box = box
        snapshot.particles.N = N
        snapshot.particles.position[:] = np.random.uniform(-L/2, L/2, (N, 3))
        snapshot.particles.velocity[:] = np.random.normal(0, 1, (N, 3))
        snapshot.particles.typeid[:] = 0
        snapshot.particles.types = ['A']
    
    simulation.create_state_from_snapshot(snapshot)
    
    # Create integrator with BussiReservoir thermostat
    import hoomd.bussi_reservoir as bussi_res
    
    bussi = bussi_res.thermostats.BussiReservoir(kT=1.5, tau=1.0)
    
    nve = hoomd.md.methods.ConstantVolume(
        filter=hoomd.filter.All(),
        thermostat=bussi
    )
    
    integrator = hoomd.md.Integrator(dt=0.001, methods=[nve])
    simulation.operations.integrator = integrator
    
    # Add LJ potential
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
    lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
    lj.r_cut[('A', 'A')] = 2.5
    simulation.operations.integrator.forces.append(lj)
    
    # Test that we can access thermostat properties
    assert bussi.kT == 1.5
    assert bussi.tau == 1.0
    
    # Run a short simulation
    simulation.run(0)  # Initialize
    
    # Check initial reservoir energies are zero
    assert bussi.total_reservoir_energy == 0.0
    assert bussi.reservoir_energy_translational == 0.0
    assert bussi.reservoir_energy_rotational == 0.0
    
    # Run simulation
    simulation.run(100)
    
    # Check that reservoir energies are being tracked
    # (We can't predict exact values, but they should be non-zero after rescaling)
    print(f"Total reservoir energy: {bussi.total_reservoir_energy}")
    print(f"Translational: {bussi.reservoir_energy_translational}")
    print(f"Rotational: {bussi.reservoir_energy_rotational}")
    
    # Test reset functionality
    bussi.reset_reservoir_energy()
    assert bussi.total_reservoir_energy == 0.0
    assert bussi.reservoir_energy_translational == 0.0
    assert bussi.reservoir_energy_rotational == 0.0


def test_bussi_reservoir_logging():
    """Test that BussiReservoir properties can be logged."""
    # Create a simple simulation (similar to above)
    device = hoomd.device.CPU()
    simulation = hoomd.Simulation(device=device, seed=42)
    
    L = 5
    N = 50
    snapshot = hoomd.Snapshot(device.communicator)
    if snapshot.communicator.rank == 0:
        box = [L, L, L, 0, 0, 0]
        snapshot.configuration.box = box
        snapshot.particles.N = N
        snapshot.particles.position[:] = np.random.uniform(-L/2, L/2, (N, 3))
        snapshot.particles.velocity[:] = np.random.normal(0, 1, (N, 3))
        snapshot.particles.typeid[:] = 0
        snapshot.particles.types = ['A']
    
    simulation.create_state_from_snapshot(snapshot)
    
    # Create integrator with BussiReservoir thermostat
    import hoomd.bussi_reservoir as bussi_res
    
    bussi = bussi_res.thermostats.BussiReservoir(kT=2.0, tau=0.5)
    
    nve = hoomd.md.methods.ConstantVolume(
        filter=hoomd.filter.All(),
        thermostat=bussi
    )
    
    integrator = hoomd.md.Integrator(dt=0.002, methods=[nve])
    simulation.operations.integrator = integrator
    
    # Add LJ potential
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
    lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
    lj.r_cut[('A', 'A')] = 2.5
    simulation.operations.integrator.forces.append(lj)
    
    # Create logger
    logger = hoomd.logging.Logger()
    logger.add(bussi, quantities=[
        'total_reservoir_energy',
        'reservoir_energy_translational', 
        'reservoir_energy_rotational',
        'instantaneous_reservoir_total'
    ])
    
    # Test that logging works
    simulation.run(0)
    log_data = logger.log()
    
    # Check that all logged quantities are present
    assert ('BussiReservoir', 'total_reservoir_energy') in log_data
    assert ('BussiReservoir', 'reservoir_energy_translational') in log_data
    assert ('BussiReservoir', 'reservoir_energy_rotational') in log_data
    assert ('BussiReservoir', 'instantaneous_reservoir_total') in log_data


if __name__ == '__main__':
    test_bussi_reservoir_basic()
    test_bussi_reservoir_logging()
    print("All tests passed!") 