// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "DipoleResponseForceComputeGPU.h"
#include "DipoleResponseForceComputeGPU.cuh"

namespace py = pybind11;

/*! \file DipoleResponseForceComputeGPU.cc
    \brief Contains code for the DipoleResponseForceComputeGPU class
*/

namespace hoomd
{
namespace cavitymd
{

/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param field_vector Electric field direction vector
    \param field_strength Electric field strength E₀
    \param sign +1 for (+) clone, -1 for (-) clone
    \param exclude_cavity Whether to exclude cavity particles
*/
DipoleResponseForceComputeGPU::DipoleResponseForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                                                             Scalar3 field_vector,
                                                             Scalar field_strength,
                                                             Scalar sign,
                                                             bool exclude_cavity)
    : DipoleResponseForceCompute(sysdef, field_vector, field_strength, sign, exclude_cavity)
{
    m_exec_conf->msg->notice(5) << "Constructing DipoleResponseForceComputeGPU" << std::endl;
    
    // Initialize autotuner with new API  
    std::vector<std::vector<unsigned int>> dimension_ranges;
    dimension_ranges.push_back({32, 64, 128, 256, 512, 1024});  // Block size options
    
    m_tuner.reset(new Autotuner<1>(dimension_ranges, 
                                   m_exec_conf, 
                                   "dipole_response_force", 
                                   5));
    
    // Allocate GPU arrays (will be resized as needed)
    GPUArray<Scalar> energy_sum(1, m_exec_conf);
    m_energy_sum.swap(energy_sum);
    
    std::cout << "DipoleResponseForceComputeGPU initialized successfully" << std::endl;
}

void DipoleResponseForceComputeGPU::setAutotunerParams(bool enable, unsigned int period)
{
    // Note: The new Autotuner API doesn't have setEnabled/setPeriod methods
    // Autotuning is enabled by default and controlled via constructor parameters
}

/*! \param timestep Current time step
*/
void DipoleResponseForceComputeGPU::computeForces(uint64_t timestep)
{
    // If electric field is disabled, use CPU implementation to zero forces
    if (!m_params.enabled)
    {
        DipoleResponseForceCompute::computeForces(timestep);
        return;
    }
    
    // Start profiling (GPU implementation)
    // Note: Profiling code removed for simplicity
    
    // Access particle data on GPU
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar> d_charge(m_pdata->getCharges(), access_location::device, access_mode::read);
    ArrayHandle<int3> d_image(m_pdata->getImages(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::overwrite);
    
    unsigned int N = m_pdata->getN();
    const BoxDim& box = m_pdata->getBox();
    
    // Find cavity particle type ID
    unsigned int L_typeid = m_pdata->getTypeByName("L");
    
    // Resize and allocate memory for energy reduction on GPU
    unsigned int num_blocks = (N + 1023) / 1024;  // Rough estimate, will be updated by tuner
    m_energy_sum.resize(num_blocks);
    ArrayHandle<Scalar> d_energy_sum(m_energy_sum, access_location::device, access_mode::overwrite);
    
    // Run the autotuner to find optimal block size
    m_tuner->begin();
    unsigned int block_size = m_tuner->getParam()[0];  // Get first dimension
    num_blocks = (N + block_size - 1) / block_size;
    
    // Ensure we have enough space for energy reduction
    if (m_energy_sum.getNumElements() < num_blocks)
    {
        m_energy_sum.resize(num_blocks);
    }
    
    // Launch GPU kernel
    gpu_compute_dipole_response_forces(d_force.data,
                                      d_pos.data,
                                      d_charge.data,
                                      d_image.data,
                                      box,
                                      N,
                                      m_params,
                                      L_typeid,
                                      d_energy_sum.data,
                                      block_size,
                                      num_blocks);
    
    // Check for errors
    if (m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    
    m_tuner->end();
    
    // Reduce energy sum on GPU
    ArrayHandle<Scalar> h_energy_sum(m_energy_sum, access_location::host, access_mode::read);
    m_electric_field_energy = Scalar(0.0);
    for (unsigned int i = 0; i < num_blocks; i++)
    {
        m_electric_field_energy += h_energy_sum.data[i];
    }
    
    // End profiling (GPU implementation)
    // Note: Profiling code removed for simplicity
}

namespace detail
{
void export_DipoleResponseForceComputeGPU(pybind11::module& m)
{
    pybind11::class_<DipoleResponseForceComputeGPU, DipoleResponseForceCompute, std::shared_ptr<DipoleResponseForceComputeGPU>>(
        m, "DipoleResponseForceComputeGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar3, Scalar, Scalar, bool>());
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd
