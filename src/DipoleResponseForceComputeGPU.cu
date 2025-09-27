/*! \file DipoleResponseForceComputeGPU.cu
    \brief Implements GPU kernel for electric field force computation for dipole moment FDR measurements
*/

#include "hoomd/ParticleData.cuh"
#include "hoomd/HOOMDMath.h"

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

#include <assert.h>

namespace hoomd
{
namespace cavitymd
{

//! Parameters for dipole response force computation  
struct dipole_response_force_params
{
    Scalar3 field_vector;           //!< Electric field direction (normalized)
    Scalar field_strength;          //!< E₀ electric field strength
    Scalar sign;                   //!< +1 for (+) clone, -1 for (-) clone
    bool enabled;                  //!< Whether electric field is active
    bool exclude_cavity;           //!< Whether to exclude cavity particles
    Scalar3 force_per_unit_charge; //!< Precomputed force per unit charge
};

namespace kernel
{

//! Optimized dipole response force computation kernel
__global__ void gpu_compute_dipole_response_forces(Scalar4* d_force,
                                                   const Scalar4* d_pos,
                                                   const Scalar* d_charge,
                                                   const int3* d_image,
                                                   const BoxDim box,
                                                   unsigned int N,
                                                   const dipole_response_force_params params,
                                                   unsigned int L_typeid,
                                                   Scalar* d_energy_result,
                                                   unsigned int block_size,
                                                   unsigned int num_blocks)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Shared memory for energy reduction
    extern __shared__ Scalar sdata[];
    
    // Initialize shared memory
    unsigned int tid = threadIdx.x;
    sdata[tid] = Scalar(0.0);
    
    // Process particles assigned to this thread
    if (idx < N)
    {
        // Get particle type
        int type = __scalar_as_int(d_pos[idx].w);
        
        // Check if we should exclude cavity particles
        bool skip_particle = params.exclude_cavity && ((unsigned int)type == L_typeid);
        
        if (!skip_particle)
        {
            // Get particle charge
            Scalar q_i = d_charge[idx];
            
            // Only compute force for charged particles
            if (fabsf(q_i) > Scalar(1e-12))
            {
                // Apply force: F⃗ᵢ = qᵢ * (sign * E₀ * Ê)
                d_force[idx].x = q_i * params.force_per_unit_charge.x;
                d_force[idx].y = q_i * params.force_per_unit_charge.y;
                d_force[idx].z = q_i * params.force_per_unit_charge.z;
                d_force[idx].w = Scalar(0.0);  // No single-particle potential energy
                
                // Calculate unwrapped position for energy computation
                Scalar3 unwrapped_pos;
                Scalar3 L = box.getL();
                unwrapped_pos.x = d_pos[idx].x + d_image[idx].x * L.x;
                unwrapped_pos.y = d_pos[idx].y + d_image[idx].y * L.y;
                unwrapped_pos.z = d_pos[idx].z + d_image[idx].z * L.z;
                
                // Compute electric field energy contribution: U = -sign * E₀ * qᵢ (Ê · r⃗ᵢ)
                Scalar field_dot_r = params.field_vector.x * unwrapped_pos.x +
                                     params.field_vector.y * unwrapped_pos.y +
                                     params.field_vector.z * unwrapped_pos.z;
                Scalar energy_contribution = -params.sign * params.field_strength * q_i * field_dot_r;
                
                // Store energy contribution in shared memory for reduction
                sdata[tid] = energy_contribution;
            }
            else
            {
                // No force on uncharged particles
                d_force[idx].x = Scalar(0.0);
                d_force[idx].y = Scalar(0.0);
                d_force[idx].z = Scalar(0.0);
                d_force[idx].w = Scalar(0.0);
            }
        }
        else
        {
            // Skip cavity particle - set force to zero
            d_force[idx].x = Scalar(0.0);
            d_force[idx].y = Scalar(0.0);
            d_force[idx].z = Scalar(0.0);
            d_force[idx].w = Scalar(0.0);
        }
    }
    
    __syncthreads();
    
    // Perform block-wise energy reduction
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    // Write block result to global memory
    if (tid == 0)
    {
        d_energy_result[blockIdx.x] = sdata[0];
    }
}

} // end namespace kernel

//! Kernel launcher for dipole response force computation
hipError_t gpu_compute_dipole_response_forces(Scalar4* d_force,
                                             const Scalar4* d_pos,
                                             const Scalar* d_charge,
                                             const int3* d_image,
                                             const BoxDim box,
                                             unsigned int N,
                                             const dipole_response_force_params& params,
                                             unsigned int L_typeid,
                                             Scalar* d_energy_result,
                                             unsigned int block_size,
                                             unsigned int num_blocks)
{
    // Calculate shared memory size
    unsigned int shared_mem_size = block_size * sizeof(Scalar);
    
    // Launch kernel
    hipLaunchKernelGGL(kernel::gpu_compute_dipole_response_forces,
                       dim3(num_blocks),
                       dim3(block_size),
                       shared_mem_size,
                       0,
                       d_force,
                       d_pos,
                       d_charge,
                       d_image,
                       box,
                       N,
                       params,
                       L_typeid,
                       d_energy_result,
                       block_size,
                       num_blocks);
    
    return hipGetLastError();
}

} // end namespace cavitymd
} // end namespace hoomd
