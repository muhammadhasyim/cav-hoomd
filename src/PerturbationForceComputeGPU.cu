/*! \file PerturbationForceComputeGPU.cu
    \brief Implements GPU kernel for perturbation force computation for FDR measurements
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

//! Parameters for perturbation force computation  
struct perturbation_force_params
{
    Scalar3 kvector;        //!< Wave vector for perturbation
    Scalar amplitude;       //!< h₀ perturbation amplitude
    Scalar sign;           //!< +1 for (+) clone, -1 for (-) clone
    bool enabled;          //!< Whether perturbation is active
    Scalar k_magnitude;     //!< |k| for reference
    Scalar k_dot_k;         //!< |k|² for optimizations
};

namespace kernel
{

//! Optimized perturbation force computation kernel
__global__ void gpu_compute_perturbation_forces(Scalar4* d_force,
                                               const Scalar4* d_pos,
                                               const int3* d_image,
                                               const BoxDim box,
                                               unsigned int N,
                                               const perturbation_force_params params,
                                               unsigned int L_typeid,
                                               Scalar* d_energy_result)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < N && params.enabled)
    {
        Scalar4 pos = d_pos[idx];
        int type = __scalar_as_int(pos.w);
        
        // Skip cavity particle (type 'L')
        if (type != (int)L_typeid)
        {
            // Unwrap position
            int3 img = d_image[idx];
            Scalar3 L = box.getL();
            Scalar3 unwrapped_pos;
            unwrapped_pos.x = pos.x + Scalar(img.x) * L.x;
            unwrapped_pos.y = pos.y + Scalar(img.y) * L.y;
            unwrapped_pos.z = pos.z + Scalar(img.z) * L.z;
            
            // Compute k·r
            Scalar k_dot_r = params.kvector.x * unwrapped_pos.x +
                            params.kvector.y * unwrapped_pos.y +
                            params.kvector.z * unwrapped_pos.z;
            
            // Compute trigonometric functions
            Scalar sin_kr = sin(k_dot_r);
            
            // Perturbation force: F = sign * h₀ * k * sin(k·r)
            Scalar force_magnitude = params.sign * params.amplitude * sin_kr;
            
            // Apply force in k direction
            Scalar4 force;
            force.x = force_magnitude * params.kvector.x;
            force.y = force_magnitude * params.kvector.y;
            force.z = force_magnitude * params.kvector.z;
            force.w = 0.0;  // No particle potential energy
            
            d_force[idx] = force;
        }
        else
        {
            // Zero force for cavity particle
            d_force[idx] = make_scalar4(0, 0, 0, 0);
        }
    }
    else if (idx < N)
    {
        // Zero force when disabled or out of bounds
        d_force[idx] = make_scalar4(0, 0, 0, 0);
    }
    
    // Set energy to zero for now (simplified implementation)
    if (idx == 0)
    {
        *d_energy_result = 0.0;
    }
}


} // end namespace kernel

//! Kernel driver for perturbation force computation
hipError_t gpu_compute_perturbation_forces(Scalar4* d_force,
                                          const Scalar4* d_pos,
                                          const int3* d_image,
                                          const BoxDim& box,
                                          unsigned int N,
                                          const perturbation_force_params& params,
                                          unsigned int L_typeid,
                                          Scalar* d_energy_result,
                                          unsigned int block_size)
{
    if (!params.enabled)
    {
        // Zero all forces when disabled
        hipMemset(d_force, 0, N * sizeof(Scalar4));
        hipMemset(d_energy_result, 0, sizeof(Scalar));
        return hipSuccess;
    }
    
    // Calculate grid dimensions
    unsigned int num_blocks = (N + block_size - 1) / block_size;
    dim3 grid(num_blocks, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // Use provided temporary energy storage (allocated in the GPU compute class)
    // d_energy_result points to already allocated memory
    
    // For now, just compute forces and set energy to zero
    // Launch perturbation force kernel (simplified version)
    unsigned int shared_bytes = block_size * sizeof(Scalar);
    kernel::gpu_compute_perturbation_forces<<<grid, threads, shared_bytes>>>(
        d_force, d_pos, d_image, box, N, params, L_typeid, d_energy_result);
    
    // Check for kernel errors
    hipError_t error = hipGetLastError();
    if (error != hipSuccess)
        return error;
    
    return hipGetLastError();
}

} // end namespace cavitymd
} // end namespace hoomd
