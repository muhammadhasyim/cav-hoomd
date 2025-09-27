// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "PerturbationForceComputeGPU.h"

/*! \file PerturbationForceComputeGPU.cc
    \brief Contains code for the PerturbationForceComputeGPU class
*/

namespace hoomd
{
namespace cavitymd
{

// Declare the GPU kernel driver
hipError_t gpu_compute_perturbation_forces(Scalar4* d_force,
                                          const Scalar4* d_pos,
                                          const int3* d_image,
                                          const BoxDim& box,
                                          unsigned int N,
                                          const perturbation_force_params& params,
                                          unsigned int L_typeid,
                                          Scalar* d_energy_result,
                                          unsigned int block_size);

/*! \param sysdef SystemDefinition containing the ParticleData
    \param kvector Wave vector for perturbation
    \param amplitude Perturbation amplitude h₀
    \param sign +1 for (+) clone, -1 for (-) clone
*/
PerturbationForceComputeGPU::PerturbationForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                                                         Scalar3 kvector,
                                                         Scalar amplitude,
                                                         Scalar sign)
    : PerturbationForceCompute(sysdef, kvector, amplitude, sign)
{
    // Initialize GPU arrays
    m_temp_energy = GPUArray<Scalar>(1, m_exec_conf);
    
    // Initialize autotuner
    m_tuner.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(m_exec_conf)},
                                   m_exec_conf, "perturbation_force"));
    
    m_exec_conf->msg->notice(5) << "Constructing PerturbationForceComputeGPU" << std::endl;
}

PerturbationForceComputeGPU::~PerturbationForceComputeGPU()
{
    m_exec_conf->msg->notice(5) << "Destroying PerturbationForceComputeGPU" << std::endl;
}

/*! \param timestep Current time step
    GPU implementation of perturbation force computation
*/
void PerturbationForceComputeGPU::computeForces(uint64_t timestep)
{
    if (!m_params.enabled)
    {
        // Zero forces when disabled
        ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::overwrite);
        hipMemset(d_force.data, 0, m_force.getNumElements() * sizeof(Scalar4));
        m_perturbation_energy = Scalar(0.0);
        return;
    }
    
    // Access particle data on GPU
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<int3> d_image(m_pdata->getImages(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar> d_energy(m_temp_energy, access_location::device, access_mode::overwrite);
    
    unsigned int N = m_pdata->getN();
    BoxDim box = m_pdata->getGlobalBox();
    
    // Get the typeid for cavity particle 'L' to exclude it
    unsigned int L_typeid;
    try {
        L_typeid = m_pdata->getTypeByName("L");
    } catch (...) {
        // No cavity particle type, use invalid type
        L_typeid = UINT_MAX;
    }
    
    // Auto-tune block size
    m_tuner->begin();
    unsigned int block_size = m_tuner->getParam()[0];
    
    // Launch GPU kernel
    hipError_t error = gpu_compute_perturbation_forces(d_force.data,
                                                       d_pos.data,
                                                       d_image.data,
                                                       box,
                                                       N,
                                                       m_params,
                                                       L_typeid,
                                                       d_energy.data,
                                                       block_size);
    
    if (error != hipSuccess)
    {
        m_exec_conf->msg->error() << "Error launching perturbation force kernel: " 
                                  << hipGetErrorString(error) << std::endl;
        throw std::runtime_error("GPU kernel launch failed");
    }
    
    m_tuner->end();
    
    // Copy energy result back to host
    ArrayHandle<Scalar> h_energy(m_temp_energy, access_location::host, access_mode::read);
    m_perturbation_energy = h_energy.data[0];
    
    if (m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
}

namespace detail
{
void export_PerturbationForceComputeGPU(pybind11::module& m)
{
    pybind11::class_<PerturbationForceComputeGPU, PerturbationForceCompute, std::shared_ptr<PerturbationForceComputeGPU>>(
        m, "PerturbationForceComputeGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar3, Scalar, Scalar>());
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd
