// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef HOOMD_BUSSI_RESERVOIR_THERMOSTAT_H
#define HOOMD_BUSSI_RESERVOIR_THERMOSTAT_H

#include "Thermostat.h"
#include <hoomd/RandomNumbers.h>
#include <hoomd/RNGIdentifiers.h>

namespace hoomd::md
{

/** Extended Bussi thermostat that tracks reservoir energy.

    This class extends the standard Bussi stochastic velocity rescaling thermostat
    to track the energy that is dumped into or taken from the thermal reservoir.
    
    The reservoir energy is calculated as the difference in kinetic energy before
    and after velocity rescaling at each timestep.
*/
class BussiReservoirThermostat : public Thermostat
{
public:
    /** Construct the thermostat.

        @param T Temperature set point over time.
        @param group Group of particles this thermostat is applied to.
        @param thermo Use to compute the thermodynamic properties of the group.
        @param sysdef Used to access the simulation seed and MPI communicator.
        @param tau Thermostat time constant.
    */
    BussiReservoirThermostat(std::shared_ptr<Variant> T,
                           std::shared_ptr<ParticleGroup> group,
                           std::shared_ptr<ComputeThermo> thermo,
                           std::shared_ptr<SystemDefinition> sysdef,
                           Scalar tau)
        : Thermostat(T, group, thermo, sysdef), m_tau(tau)
    {
        resetReservoirEnergy();
    }

    std::array<Scalar, 2> getRescalingFactorsOne(uint64_t timestep, Scalar deltaT) override
    {
        if (deltaT == 0.0)
        {
            return {1.0, 1.0};
        }

        m_thermo->compute(timestep);

        const auto translational_dof = m_thermo->getTranslationalDOF();
        const auto rotational_dof = m_thermo->getRotationalDOF();
        const auto translational_kinetic_energy = m_thermo->getTranslationalKineticEnergy();
        const auto rotational_kinetic_energy = m_thermo->getRotationalKineticEnergy();
        
        if ((translational_dof != 0 && translational_kinetic_energy == 0)
            || (rotational_dof != 0 && rotational_kinetic_energy == 0))
        {
            throw std::runtime_error("Bussi thermostat requires non-zero initial momenta.");
        }

        unsigned int instance_id = 0;
        if (m_group->getNumMembersGlobal() > 0)
            instance_id = m_group->getMemberTag(0);
        RandomGenerator rng(Seed(RNGIdentifier::BussiThermostat, timestep, m_sysdef->getSeed()),
                           instance_id);

        const auto set_T = m_T->operator()(timestep);

        // Compute rescaling factors
        const auto translational_factor = compute_rescale_factor(translational_kinetic_energy,
                                                               translational_dof,
                                                               deltaT,
                                                               set_T,
                                                               rng);
        const auto rotational_factor = compute_rescale_factor(rotational_kinetic_energy, 
                                                            rotational_dof, 
                                                            deltaT, 
                                                            set_T, 
                                                            rng);

        // Calculate energy dumped to/taken from reservoir
        // Reservoir energy = KE_old - KE_new = KE_old * (1 - alpha^2)
        // Positive means energy dumped to reservoir, negative means taken from reservoir
        Scalar delta_trans = translational_kinetic_energy * (1.0 - translational_factor * translational_factor);
        Scalar delta_rot = rotational_kinetic_energy * (1.0 - rotational_factor * rotational_factor);
        
        // Accumulate reservoir energies
        m_reservoir_energy_translational += delta_trans;
        m_reservoir_energy_rotational += delta_rot;
        
        // Track instantaneous values for logging
        m_instantaneous_reservoir_translational = delta_trans;
        m_instantaneous_reservoir_rotational = delta_rot;

        return {translational_factor, rotational_factor};
    }

    /// Get the thermostat time constant.
    Scalar getTau() const
    {
        return m_tau;
    }

    /// Set the thermostat time constant.
    void setTau(Scalar tau)
    {
        m_tau = tau;
    }

    /// Get cumulative reservoir energy from translational degrees of freedom.
    Scalar getReservoirEnergyTranslational() const
    {
        return m_reservoir_energy_translational;
    }

    /// Get cumulative reservoir energy from rotational degrees of freedom.
    Scalar getReservoirEnergyRotational() const
    {
        return m_reservoir_energy_rotational;
    }

    /// Get total cumulative reservoir energy.
    Scalar getTotalReservoirEnergy() const
    {
        return m_reservoir_energy_translational + m_reservoir_energy_rotational;
    }

    /// Get instantaneous reservoir energy change from translational DOF (last timestep).
    Scalar getInstantaneousReservoirTranslational() const
    {
        return m_instantaneous_reservoir_translational;
    }

    /// Get instantaneous reservoir energy change from rotational DOF (last timestep).
    Scalar getInstantaneousReservoirRotational() const
    {
        return m_instantaneous_reservoir_rotational;
    }

    /// Get total instantaneous reservoir energy change (last timestep).
    Scalar getInstantaneousReservoirTotal() const
    {
        return m_instantaneous_reservoir_translational + m_instantaneous_reservoir_rotational;
    }

    /// Reset all reservoir energy counters to zero.
    void resetReservoirEnergy()
    {
        m_reservoir_energy_translational = 0.0;
        m_reservoir_energy_rotational = 0.0;
        m_instantaneous_reservoir_translational = 0.0;
        m_instantaneous_reservoir_rotational = 0.0;
    }

protected:
    /// Thermostat time constant
    Scalar m_tau;
    
    /// Cumulative reservoir energies
    Scalar m_reservoir_energy_translational = 0.0;
    Scalar m_reservoir_energy_rotational = 0.0;
    
    /// Instantaneous reservoir energy changes (for logging)
    Scalar m_instantaneous_reservoir_translational = 0.0;
    Scalar m_instantaneous_reservoir_rotational = 0.0;

    /** Compute the rescaling factor (now with correct sign determination)

        @param K kinetic energy
        @param degrees_of_freedom Number of degrees of freedom with this energy.
        @param deltaT Time step size.
        @param set_T Temperature set point.
        @param rng Random number generator.
    **/
    Scalar compute_rescale_factor(Scalar K,
                                double degrees_of_freedom,
                                Scalar deltaT,
                                Scalar set_T,
                                RandomGenerator& rng)
    {
        if (degrees_of_freedom == 0)
            return Scalar(1.0);

        double time_decay_factor = 0.0;
        if (m_tau != 0.0)
        {
            time_decay_factor = exp(-deltaT / m_tau);
        }

        NormalDistribution<double> normal(1.0);
        Scalar r_normal_one = normal(rng);

        GammaDistribution<double> gamma((degrees_of_freedom - 1.0) / Scalar(2.0), Scalar(1.0));
        double r_gamma = 0.0;
        if (degrees_of_freedom > 1.0)
        {
            r_gamma = 2.0 * gamma(rng);
        }

        double v = set_T / 2.0 / K;
        double term1 = v * (1.0 - time_decay_factor) * (r_gamma + r_normal_one * r_normal_one);
        double term2 = 2.0 * r_normal_one * sqrt(v * (1.0 - time_decay_factor) * time_decay_factor);

        // Calculate α² (always positive)
        double alpha_squared = time_decay_factor + term1 + term2;
        double alpha_magnitude = sqrt(alpha_squared);
        
        // Determine sign according to equation (A8) from Bussi et al. 2009
        // sign[α(t)] = sign[R(t) + √(cNf*K*(t)/((1-c)K̄*))]
        double c = time_decay_factor;  // c = exp(-Δt/τ)
        double K_bar = set_T * degrees_of_freedom / 2.0;  // average kinetic energy
        double sign_term = r_normal_one + sqrt(c * degrees_of_freedom * K / ((1.0 - c) * K_bar));
        
        // Apply the correct sign
        if (sign_term >= 0.0)
        {
            return Scalar(alpha_magnitude);
        }
        else
        {
            return Scalar(-alpha_magnitude);
        }
    }
};

} // namespace hoomd::md

#endif // HOOMD_BUSSI_RESERVOIR_THERMOSTAT_H 