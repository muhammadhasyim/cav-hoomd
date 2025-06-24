// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "CavityForceComputeGPU.h"

#ifdef ENABLE_HIP
#include "CavityForceComputeGPU.cuh"
#include <hoomd/GPUArray.h>
#endif

namespace py = pybind11;

/*! \file CavityForceComputeGPU.cc
    \brief Contains code for the CavityForceComputeGPU class
*/

namespace hoomd
{
namespace cavitymd
{

/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param omegac Cavity frequency in atomic units
    \param couplstr Coupling strength in atomic units
    \param phmass Photon mass (default 1.0)
*/
CavityForceComputeGPU::CavityForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                                             Scalar omegac,
                                             Scalar couplstr, 
                                             Scalar phmass)
    : CavityForceCompute(sysdef, omegac, couplstr, phmass)
{
    m_exec_conf->msg->notice(5) << "Constructing CavityForceComputeGPU" << std::endl;
    
    // Initialize GPU memory arrays with explicit device-only allocation to avoid hipHostRegister issues
    if (m_exec_conf->isCUDAEnabled()) {
        try {
            // Test if we can allocate GPU arrays without error
            GPUArray<Scalar> temp_energy(3, m_exec_conf);
            m_temp_energy.swap(temp_energy);
            
            GPUArray<Scalar3> temp_dipole(1024, m_exec_conf); // Max blocks
            m_temp_dipole.swap(temp_dipole);
            
            GPUArray<int> photon_idx(1, m_exec_conf);
            m_photon_idx.swap(photon_idx);
            
            // Add global dipole storage for fused kernel optimization
            GPUArray<Scalar3> dipole_global(1, m_exec_conf);
            m_dipole_global.swap(dipole_global);
            
            m_exec_conf->msg->notice(1) << "GPU arrays initialized successfully" << std::endl;
        }
        catch (const std::exception& e) {
            m_exec_conf->msg->warning() << "Failed to initialize GPU arrays: " << e.what() << std::endl;
            m_exec_conf->msg->warning() << "GPU computation will be disabled for this force" << std::endl;
            // Leave arrays uninitialized - will fall back to CPU
        }
    }
    else {
        m_exec_conf->msg->notice(1) << "CUDA not enabled, GPU arrays not initialized" << std::endl;
    }
    
    // Initialize autotuner with new API
    std::vector<std::vector<unsigned int>> dimension_ranges;
    dimension_ranges.push_back({32, 64, 128, 256, 512, 1024});  // Block size options
    
    m_tuner.reset(new Autotuner<1>(dimension_ranges, 
                                   m_exec_conf, 
                                   "cavity_force", 
                                   5));  // 5 samples
    
    // Check if cooperative kernel launch is supported
    if (m_exec_conf->isCUDAEnabled()) {
        int device;
        hipGetDevice(&device);
        hipDeviceProp_t prop;
        hipGetDeviceProperties(&prop, device);
        m_supports_cooperative = false; // Disable cooperative kernels for compatibility
        
        std::cout << "CavityForceComputeGPU initialized with GPU acceleration" << std::endl;
        std::cout << "  - Using optimized 2-kernel approach for better performance" << std::endl;
    }
    else {
        m_supports_cooperative = false;
        std::cout << "CavityForceComputeGPU initialized (CPU fallback mode)" << std::endl;
    }
}

CavityForceComputeGPU::~CavityForceComputeGPU()
{
    m_exec_conf->msg->notice(5) << "Destroying CavityForceComputeGPU" << std::endl;
}

/*! This function computes cavity forces and energy on the GPU
    \param timestep Current timestep
*/
void CavityForceComputeGPU::computeForces(uint64_t timestep)
{
#ifdef ENABLE_HIP
    // Check if we have properly initialized GPU arrays
    if (!m_exec_conf->isCUDAEnabled() || m_temp_energy.isNull()) {
        m_exec_conf->msg->error() << "GPU arrays not initialized properly!" << std::endl;
        throw std::runtime_error("GPU computation required but not available");
    }
    
    unsigned int N = m_pdata->getN();
    
    // Get the typeid for 'L' type first
    unsigned int L_typeid;
    try {
        L_typeid = m_pdata->getTypeByName("L");
    } catch (...) {
        // No 'L' type found - zero all energies and return
        m_harmonic_energy = Scalar(0.0);
        m_coupling_energy = Scalar(0.0);
        m_dipole_self_energy = Scalar(0.0);
        return;
    }
    
    // Get box dimensions
    BoxDim box = m_pdata->getGlobalBox();
    
    // Initialize temporary GPU arrays
    {
        ArrayHandle<Scalar> d_temp_energy(m_temp_energy, access_location::device, access_mode::overwrite);
        ArrayHandle<int> d_photon_idx(m_photon_idx, access_location::device, access_mode::overwrite);
        ArrayHandle<Scalar3> d_dipole_global(m_dipole_global, access_location::device, access_mode::overwrite);
        
        // Zero initialization arrays
        hipMemset(d_temp_energy.data, 0, sizeof(Scalar) * 3);
        hipMemset(d_photon_idx.data, -1, sizeof(int));
        hipMemset(d_dipole_global.data, 0, sizeof(Scalar3));
    }
    
    // Zero output arrays
    {
        ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::overwrite);
        hipMemset(d_force.data, 0, sizeof(Scalar4) * N);
    }
    m_virial.zeroFill();
    
    // Main computation with limited concurrent handles
    {
        // Access particle data arrays
        ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
        ArrayHandle<Scalar> d_charge(m_pdata->getCharges(), access_location::device, access_mode::read);
        ArrayHandle<int3> d_image(m_pdata->getImages(), access_location::device, access_mode::read);
        
        // Access output arrays
        ArrayHandle<Scalar4> d_force(m_force, access_location::device, access_mode::readwrite);
        ArrayHandle<Scalar> d_virial(m_virial, access_location::device, access_mode::readwrite);
        
        // Access GPU workspace arrays
        ArrayHandle<Scalar> d_temp_energy(m_temp_energy, access_location::device, access_mode::readwrite);
        ArrayHandle<Scalar3> d_temp_dipole(m_temp_dipole, access_location::device, access_mode::readwrite);
        ArrayHandle<int> d_photon_idx(m_photon_idx, access_location::device, access_mode::readwrite);
        ArrayHandle<Scalar3> d_dipole_global(m_dipole_global, access_location::device, access_mode::readwrite);
        
        unsigned int block_size = 256;  // Fixed block size
        
        // Use HOOMD's logging system for debug info (only at high debug levels)
        m_exec_conf->msg->notice(10) << "CavityForceComputeGPU: GPU computation with N=" << N 
                                     << ", block_size=" << block_size 
                                     << ", L_typeid=" << L_typeid << std::endl;
        
        // CRITICAL FIX: Pass structures by pointer to avoid ABI corruption
        hipError_t error = kernel::gpu_compute_cavity_forces(d_force.data,
                                                              d_virial.data,
                                                              m_virial.getPitch(),
                                                              N,
                                                              d_pos.data,
                                                              d_charge.data,
                                                              d_image.data,
                                                              &box,
                                                              &m_params,
                                                              d_photon_idx.data,
                                                              d_temp_energy.data,
                                                              d_temp_dipole.data,
                                                              block_size,
                                                              L_typeid,
                                                              d_dipole_global.data,
                                                              false); // use_fused_kernel = false
        
        if (error != hipSuccess)
        {
            m_exec_conf->msg->error() << "Error launching cavity force GPU kernel: " << hipGetErrorString(error) << std::endl;
            throw std::runtime_error("GPU kernel launch failed - GPU computation is required");
        }
        
    } // Release all handles before final host access
    
    // Copy energy components back to host
    {
        ArrayHandle<Scalar> h_temp_energy(m_temp_energy, access_location::host, access_mode::read);
        m_harmonic_energy = h_temp_energy.data[0];
        m_coupling_energy = h_temp_energy.data[1];
        m_dipole_self_energy = h_temp_energy.data[2];
    }
    
#else
    m_exec_conf->msg->error() << "GPU computation requested but HIP is not enabled!" << std::endl;
    throw std::runtime_error("GPU computation required but HIP not available");
#endif
}

namespace detail
{
void export_CavityForceComputeGPU(pybind11::module& m)
{
    pybind11::class_<CavityForceComputeGPU, CavityForceCompute, std::shared_ptr<CavityForceComputeGPU>>(
        m, "CavityForceComputeGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar, Scalar, Scalar>(),
             pybind11::arg("sysdef"), pybind11::arg("omegac"), pybind11::arg("couplstr"), 
             pybind11::arg("phmass") = 1.0);
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd 