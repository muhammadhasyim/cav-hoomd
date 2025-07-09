// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*! \file CavityForceComputeGPU.cc
    \brief GPU implementation of cavity force computation
*/

#include "CavityForceComputeGPU.h"
#include "hoomd/Autotuner.h"

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

#include <stdexcept>
#include <iostream>

namespace hoomd
{
namespace cavitymd
{

// External GPU kernel declaration
extern "C" hipError_t gpu_compute_cavity_forces_driver(
    Scalar4* d_force,
    const Scalar4* d_pos,
    const Scalar4* d_vel,
    const Scalar* d_charge,
    const int3* d_image,
    const BoxDim& box,
    unsigned int N,
    unsigned int L_typeid,
    Scalar current_couplstr, // Now pass evaluated coupling strength
    Scalar omegac,
    Scalar phmass,
    Scalar damping_ratio,
    Scalar* d_temp_energy,
    int* d_photon_idx,
    Scalar3* d_dipole_result,
    unsigned int block_size);

/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param omegac Cavity frequency in atomic units
    \param couplstr Coupling strength variant in atomic units
    \param phmass Photon mass (default 1.0)
    \param damping_ratio Damping ratio (default 0.0)
*/
CavityForceComputeGPU::CavityForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                                             Scalar omegac,
                                             std::shared_ptr<Variant> couplstr,
                                             Scalar phmass,
                                             Scalar damping_ratio)
    : CavityForceCompute(sysdef, omegac, couplstr, phmass, damping_ratio)
{
    m_exec_conf->msg->notice(5) << "Constructing CavityForceComputeGPU" << std::endl;
    
    // Initialize autotuners
    m_tuner_photon_dipole.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(m_exec_conf)},
                                                  m_exec_conf,
                                                  "cavity_photon_dipole"));
    m_tuner_forces.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(m_exec_conf)},
                                          m_exec_conf,
                                          "cavity_forces"));
    
    // Add autotuners to base class list
    m_autotuners.insert(m_autotuners.end(), {m_tuner_photon_dipole, m_tuner_forces});
    
    // Allocate GPU memory
    GPUArray<Scalar> temp_energy(3, m_exec_conf);
    m_temp_energy.swap(temp_energy);
    
    GPUArray<int> photon_idx_gpu(1, m_exec_conf);
    m_photon_idx_gpu.swap(photon_idx_gpu);
    
    GPUArray<Scalar3> dipole_gpu(1, m_exec_conf);
    m_dipole_gpu.swap(dipole_gpu);
    
    // Display initial coupling strength
    Scalar initial_coupling = couplstr ? (*couplstr)(0) : 0.0;
    
    std::cout << "CavityForceComputeGPU initialized: "
              << "omegac=" << omegac << " a.u., "
              << "couplstr=" << initial_coupling << " a.u. (at t=0), "
              << "K=" << m_params.K << " a.u., "
              << "damping_ratio=" << damping_ratio << std::endl;
}

CavityForceComputeGPU::~CavityForceComputeGPU()
{
    m_exec_conf->msg->notice(5) << "Destroying CavityForceComputeGPU" << std::endl;
}

/*! \param timestep Current time step
    This function computes the cavity-molecule interaction forces on the GPU
*/
void CavityForceComputeGPU::computeForces(uint64_t timestep)
{
    // Access the particle data arrays
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar> d_charge(m_pdata->getCharges(), access_location::device, access_mode::read);
    ArrayHandle<int3> d_image(m_pdata->getImages(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar> d_virial(m_virial, access_location::device, access_mode::overwrite);
    
    // Zero virial
    hipMemset(d_virial.data, 0, sizeof(Scalar) * m_virial.getNumElements());
    
    const BoxDim& box = m_pdata->getGlobalBox();
    unsigned int N = m_pdata->getN();
    
    if (N == 0) return;
    
    // Get current coupling strength at this timestep
    Scalar current_couplstr = m_params.couplstr ? (*m_params.couplstr)(timestep) : Scalar(0.0);
    
    // Access GPU memory
    ArrayHandle<Scalar> d_temp_energy(m_temp_energy, access_location::device, access_mode::overwrite);
    ArrayHandle<int> d_photon_idx(m_photon_idx_gpu, access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar3> d_dipole_result(m_dipole_gpu, access_location::device, access_mode::overwrite);
    
    // Get L_typeid
    unsigned int L_typeid = m_pdata->getTypeByName("L");
    
    // Start autotuning and determine optimal block size
    m_tuner_forces->begin();
    unsigned int block_size = m_tuner_forces->getParam()[0];
    
    try {
        // Call GPU driver function with evaluated coupling strength
        hipError_t error = gpu_compute_cavity_forces_driver(
            d_force.data,
            d_pos.data,
            d_vel.data,
            d_charge.data,
            d_image.data,
            box,
            N,
            L_typeid,
            current_couplstr, // Pass evaluated coupling strength
            m_params.omegac,
            m_params.phmass,
            m_params.damping_ratio,
            d_temp_energy.data,
            d_photon_idx.data,
            d_dipole_result.data,
            block_size);
        
        if (error != hipSuccess) {
            throw std::runtime_error("GPU cavity force computation failed");
        }
        
        // Update autotuner
        m_tuner_forces->end();
        
        // Copy energy results back to host
        ArrayHandle<Scalar> h_temp_energy(m_temp_energy, access_location::host, access_mode::read);
        m_harmonic_energy = h_temp_energy.data[0];
        m_coupling_energy = h_temp_energy.data[1];
        m_dipole_self_energy = h_temp_energy.data[2];
        
    } catch (const std::exception& e) {
        m_exec_conf->msg->error() << "GPU cavity force computation error: " << e.what() << std::endl;
        throw;
    }
}

namespace detail
{
void export_CavityForceComputeGPU(pybind11::module& m)
{
    pybind11::class_<CavityForceComputeGPU, CavityForceCompute, std::shared_ptr<CavityForceComputeGPU>>(
        m, "CavityForceComputeGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar, std::shared_ptr<Variant>, Scalar, Scalar>(),
             pybind11::arg("sysdef"), pybind11::arg("omegac"), pybind11::arg("couplstr"), 
             pybind11::arg("phmass") = 1.0, pybind11::arg("damping_ratio") = 0.0);
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd 