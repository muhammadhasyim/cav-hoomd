/*! \file CavityForceComputeGPU.cu
    \brief Implements GPU kernel for cavity force computation
    
    DEBUG FLAGS:
    - Define CAVITY_DEBUG_VERBOSE during compilation to enable verbose GPU kernel output
    - C++ debug messages use HOOMD's notice level 9 (enable with msg.setNoticeLevel(9))
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

//! Parameters for cavity force computation  
struct cavity_force_params
{
    Scalar omegac;     //!< Cavity frequency in atomic units
    Scalar couplstr;   //!< Coupling strength in atomic units  
    Scalar K;          //!< Spring constant (phmass * omegac^2)
    Scalar phmass;     //!< Photon mass
};

namespace kernel
{

//! Optimized single-pass kernel that finds photon and computes dipole
__global__ void gpu_compute_photon_and_dipole(const Scalar4* d_pos,
                                               const Scalar* d_charge,
                                               const int3* d_image,
                                               const BoxDim box,
                                               unsigned int N,
                                               unsigned int L_typeid,
                                               int* d_photon_idx,
                                               Scalar3* d_dipole_result)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    
    // Use external shared memory for dynamic allocation
    extern __shared__ Scalar3 sdata[];
    __shared__ int shared_photon_idx;
    
    // Initialize shared memory
    if (tid == 0) {
        shared_photon_idx = -1;
    }
    
    // Bounds check for shared memory access
    if (tid < blockDim.x) {
        sdata[tid] = make_scalar3(0, 0, 0);
    }
    __syncthreads();
    
    // Each thread processes one particle
    Scalar3 local_dipole = make_scalar3(0, 0, 0);
    int local_photon_candidate = -1;
    
    if (idx < N) {
        Scalar4 pos = d_pos[idx];
        int type = __scalar_as_int(pos.w);
        
        // Check if this is the photon particle
        if (type == (int)L_typeid) {
            local_photon_candidate = (int)idx;
        } else {
            // Compute dipole contribution from molecular particles
            int3 img = d_image[idx];
            Scalar charge = d_charge[idx];
            
            // Unwrap position
            Scalar3 L = box.getL();
            Scalar3 unwrapped_pos;
            unwrapped_pos.x = pos.x + Scalar(img.x) * L.x;
            unwrapped_pos.y = pos.y + Scalar(img.y) * L.y;
            unwrapped_pos.z = pos.z + Scalar(img.z) * L.z;
            
            // Accumulate charge * position (only x,y components needed for cavity)
            local_dipole.x = charge * unwrapped_pos.x;
            local_dipole.y = charge * unwrapped_pos.y;
            local_dipole.z = 0;
        }
    }
    
    // Store photon candidate if found
    if (local_photon_candidate != -1) {
        atomicExch(&shared_photon_idx, local_photon_candidate);
    }
    
    // Load dipole into shared memory for reduction
    if (tid < blockDim.x) {
        sdata[tid] = local_dipole;
    }
    __syncthreads();
    
    // Tree reduction for dipole calculation
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s && tid + s < blockDim.x) {
            sdata[tid].x += sdata[tid + s].x;
            sdata[tid].y += sdata[tid + s].y;
        }
        __syncthreads();
    }
    
    // Store block-level results
    if (tid == 0) {
        // Use atomics to accumulate dipole across blocks
        atomicAdd(&d_dipole_result[0].x, sdata[0].x);
        atomicAdd(&d_dipole_result[0].y, sdata[0].y);
        
        // Store photon index if found
        if (shared_photon_idx != -1) {
            atomicExch(d_photon_idx, shared_photon_idx);
        }
    }
}

//! Optimized force computation kernel
__global__ void gpu_compute_forces_optimized(Scalar4* d_force,
                                              const Scalar4* d_pos,
                                              const Scalar* d_charge,
                                              const int3* d_image,
                                              const BoxDim box,
                                              unsigned int N,
                                              int photon_idx,
                                              const Scalar3 dipole_xy,
                                              const cavity_force_params params,
                                              Scalar* d_temp_energy,
                                              unsigned int L_typeid)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx >= N || photon_idx < 0)
        return;
    
    // Get photon position (unwrapped) - computed once per block
    __shared__ Scalar3 shared_q_photon;
    __shared__ Scalar3 shared_q_photon_xy;
    
    if (threadIdx.x == 0) {
        Scalar4 photon_pos_data = d_pos[photon_idx];
        int3 photon_img = d_image[photon_idx];
        Scalar3 L = box.getL();
        
        shared_q_photon.x = photon_pos_data.x + Scalar(photon_img.x) * L.x;
        shared_q_photon.y = photon_pos_data.y + Scalar(photon_img.y) * L.y;
        shared_q_photon.z = photon_pos_data.z + Scalar(photon_img.z) * L.z;
        
        shared_q_photon_xy = make_scalar3(shared_q_photon.x, shared_q_photon.y, 0);
    }
    __syncthreads();
    
    // Compute energy components (only in first thread)
    if (idx == 0) {
        Scalar harmonic_energy = Scalar(0.5) * params.K * 
                               (shared_q_photon.x * shared_q_photon.x + 
                                shared_q_photon.y * shared_q_photon.y + 
                                shared_q_photon.z * shared_q_photon.z);
        Scalar coupling_energy = params.couplstr * 
                               (dipole_xy.x * shared_q_photon_xy.x + dipole_xy.y * shared_q_photon_xy.y);
        Scalar dipole_self_energy = Scalar(0.5) * (params.couplstr * params.couplstr / params.K) *
                                  (dipole_xy.x * dipole_xy.x + dipole_xy.y * dipole_xy.y);
        
        d_temp_energy[0] = harmonic_energy;
        d_temp_energy[1] = coupling_energy;
        d_temp_energy[2] = dipole_self_energy;
        
        Scalar total_energy = harmonic_energy + coupling_energy + dipole_self_energy;
        // DO NOT assign energy to particle potential energy - prevents double-counting
        d_force[photon_idx].w = 0.0;
    }
    
    // Compute forces for this particle
    if (idx < N) {
        int type = __scalar_as_int(d_pos[idx].w);
        
        if (type != (int)L_typeid) { // Molecular particle
            Scalar charge = d_charge[idx];
            
            // Dq = q_photon_xy + (g/K) * dipole_xy
            Scalar3 Dq;
            Dq.x = shared_q_photon_xy.x + (params.couplstr / params.K) * dipole_xy.x;
            Dq.y = shared_q_photon_xy.y + (params.couplstr / params.K) * dipole_xy.y;
            Dq.z = 0;
            
            // Force = -g * charge * Dq
            Scalar3 force;
            force.x = -params.couplstr * charge * Dq.x;
            force.y = -params.couplstr * charge * Dq.y;
            force.z = Scalar(0.0);
            
            d_force[idx].x = force.x;
            d_force[idx].y = force.y;
            d_force[idx].z = force.z;
            
        } else if ((int)idx == photon_idx) { // Cavity particle
            // Force = -K * q_photon - g * dipole_xy
            Scalar3 photon_force;
            photon_force.x = -params.K * shared_q_photon.x - params.couplstr * dipole_xy.x;
            photon_force.y = -params.K * shared_q_photon.y - params.couplstr * dipole_xy.y;
            photon_force.z = -params.K * shared_q_photon.z;
            
            d_force[idx].x = photon_force.x;
            d_force[idx].y = photon_force.y;
            d_force[idx].z = photon_force.z;
        }
    }
}

//! GPU kernel for computing cavity forces
__global__ void gpu_compute_cavity_force_kernel_legacy(Scalar4* d_force,
                                                     const Scalar4* d_pos,
                                                     const Scalar* d_charge,
                                                     const int3* d_image,
                                                     const BoxDim box,
                                                     const unsigned int N,
                                                     const Scalar3 q_photon_xy,
                                                     const Scalar3 dipole_xy,
                                                     const cavity_force_params params,
                                                     unsigned int L_typeid)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Check bounds
    if (idx >= N) return;
    
    // Share photon position
    __shared__ Scalar3 shared_q_photon_xy;
    __shared__ Scalar3 shared_q_photon;
    if (threadIdx.x == 0) {
        shared_q_photon_xy = q_photon_xy;
        
        // For photon force calculation, we need full q_photon
        // Find photon particle and get its position
        for (unsigned int i = 0; i < N; i++) {
            if (__scalar_as_int(d_pos[i].w) == (int)L_typeid) {
                shared_q_photon.x = d_pos[i].x;
                shared_q_photon.y = d_pos[i].y;
                shared_q_photon.z = d_pos[i].z;
                break;
            }
        }
    }
    __syncthreads();
    
    if (idx < N) {
        int type = __scalar_as_int(d_pos[idx].w);
        
        if (type != (int)L_typeid) { // Molecular particle
            Scalar charge = d_charge[idx];
            
            // Dq = q_photon_xy + (g/K) * dipole_xy
            Scalar3 Dq;
            Dq.x = shared_q_photon_xy.x + (params.couplstr / params.K) * dipole_xy.x;
            Dq.y = shared_q_photon_xy.y + (params.couplstr / params.K) * dipole_xy.y;
            Dq.z = 0;
            
            // Force = -g * charge * Dq
            Scalar3 force;
            force.x = -params.couplstr * charge * Dq.x;
            force.y = -params.couplstr * charge * Dq.y;
            force.z = Scalar(0.0);
            
            d_force[idx].x = force.x;
            d_force[idx].y = force.y;
            d_force[idx].z = force.z;
            
        } else { // Photon particle (type == L_typeid)
            // Force = -K * q_photon - g * dipole_xy
            Scalar3 photon_force;
            photon_force.x = -params.K * shared_q_photon.x - params.couplstr * dipole_xy.x;
            photon_force.y = -params.K * shared_q_photon.y - params.couplstr * dipole_xy.y;
            photon_force.z = -params.K * shared_q_photon.z;
            
            d_force[idx].x = photon_force.x;
            d_force[idx].y = photon_force.y;
            d_force[idx].z = photon_force.z;
        }
    }
}

//! GPU kernel to find the photon particle
__global__ void gpu_find_photon_kernel(const Scalar4* d_pos,
                                        const unsigned int N,
                                        int* d_photon_idx,
                                        unsigned int L_typeid)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx == 0) {
        // Initialize to -1 (not found)
        *d_photon_idx = -1;
        
        // Search for photon particle
        for (unsigned int i = 0; i < N; i++) {
            int particle_type = __scalar_as_int(d_pos[i].w);
            
            if (particle_type == (int)L_typeid) {
                *d_photon_idx = i;
#ifdef CAVITY_DEBUG_VERBOSE
                printf("GPU: Found photon particle at index %d\n", i);
#endif
                break;
            }
        }
        
        if (*d_photon_idx == -1) {
#ifdef CAVITY_DEBUG_VERBOSE
            printf("GPU: No photon particle found! Searched %d particles for L_typeid %d\n", N, (int)L_typeid);
#endif
        }
    }
}

//! GPU kernel to reduce dipole moment using shared memory
__global__ void gpu_compute_dipole_kernel(const Scalar4* d_pos,
                                          const Scalar* d_charge,
                                          const int3* d_image,
                                          const BoxDim box,
                                          const unsigned int N,
                                          const int photon_idx,
                                          Scalar3* d_temp_dipole,
                                          unsigned int L_typeid)
{
    extern __shared__ Scalar3 sdata[];
    
    unsigned int tid = threadIdx.x;
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Initialize shared memory
    sdata[tid] = make_scalar3(0.0, 0.0, 0.0);
    
    // Compute partial dipole for this thread
    if (idx < N && (int)idx != photon_idx && __scalar_as_int(d_pos[idx].w) != (int)L_typeid) {
        Scalar charge = d_charge[idx];
        Scalar4 pos_data = d_pos[idx];
        int3 img = d_image[idx];
        
        // Get box dimensions
        Scalar3 L = box.getL();
        
        // Compute unwrapped position
        Scalar3 unwrapped_pos;
        unwrapped_pos.x = pos_data.x + Scalar(img.x) * L.x;
        unwrapped_pos.y = pos_data.y + Scalar(img.y) * L.y;
        unwrapped_pos.z = pos_data.z + Scalar(img.z) * L.z;
        
        // Accumulate dipole contribution
        sdata[tid].x = charge * unwrapped_pos.x;
        sdata[tid].y = charge * unwrapped_pos.y;
        sdata[tid].z = charge * unwrapped_pos.z;
    }
    
    __syncthreads();
    
    // Reduction in shared memory
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid].x += sdata[tid + s].x;
            sdata[tid].y += sdata[tid + s].y;
            sdata[tid].z += sdata[tid + s].z;
        }
        __syncthreads();
    }
    
    // Write result for this block to global memory
    if (tid == 0) {
        d_temp_dipole[blockIdx.x] = sdata[0];
    }
}

//! New main GPU kernel for cavity force computation
__global__ void gpu_compute_cavity_force_kernel(Scalar4* d_force,
                                                 const Scalar4* d_pos,
                                                 const Scalar* d_charge,
                                                 const int3* d_image,
                                                 const BoxDim box,
                                                 const unsigned int N,
                                                 const int photon_idx,
                                                 const Scalar3 dipole_xy,
                                                 const cavity_force_params params,
                                                 Scalar* d_temp_energy,
                                                 unsigned int L_typeid)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // BOUNDS CHECK: Prevent illegal memory access
    if (idx >= N || N == 0 || idx >= 100000 || photon_idx < 0 || photon_idx >= (int)N) return;
    
    // Get photon position (unwrapped) - shared across all threads
    __shared__ Scalar3 shared_q_photon;
    __shared__ Scalar3 shared_q_photon_xy;
    
    if (threadIdx.x == 0) {
        Scalar4 photon_pos_data = d_pos[photon_idx];
        int3 photon_img = d_image[photon_idx];
        Scalar3 L = box.getL();
        
        shared_q_photon.x = photon_pos_data.x + Scalar(photon_img.x) * L.x;
        shared_q_photon.y = photon_pos_data.y + Scalar(photon_img.y) * L.y;
        shared_q_photon.z = photon_pos_data.z + Scalar(photon_img.z) * L.z;
        
        shared_q_photon_xy = make_scalar3(shared_q_photon.x, shared_q_photon.y, 0);
    }
    __syncthreads();
    
    if ((int)idx != photon_idx)
    {
        // Molecular particle
        Scalar charge = d_charge[idx];
        
        // Correct formula: Dq = q_photon_xy + (g/K) * dipole_xy
        Scalar3 Dq;
        Dq.x = shared_q_photon_xy.x + (params.couplstr / params.K) * dipole_xy.x;
        Dq.y = shared_q_photon_xy.y + (params.couplstr / params.K) * dipole_xy.y;
        Dq.z = 0;
        
        Scalar3 force;
        force.x = -params.couplstr * charge * Dq.x;
        force.y = -params.couplstr * charge * Dq.y;
        force.z = Scalar(0.0);
        
        d_force[idx].x = force.x;
        d_force[idx].y = force.y;
        d_force[idx].z = force.z;
    }
    else if ((int)idx == photon_idx)
    {
        // Photon particle
#ifdef CAVITY_DEBUG_VERBOSE
        printf("GPU: Computing photon force for particle %d\n", idx);
#endif
        
        Scalar3 photon_force;
        photon_force.x = -params.K * shared_q_photon.x - params.couplstr * dipole_xy.x;
        photon_force.y = -params.K * shared_q_photon.y - params.couplstr * dipole_xy.y;
        photon_force.z = -params.K * shared_q_photon.z;
        
        d_force[idx].x = photon_force.x;
        d_force[idx].y = photon_force.y;
        d_force[idx].z = photon_force.z;
        
        // Energy computation - write energies every time since they should be the same
        {
            Scalar harmonic_energy = Scalar(0.5) * params.K * (shared_q_photon.x * shared_q_photon.x + shared_q_photon.y * shared_q_photon.y + shared_q_photon.z * shared_q_photon.z);
            Scalar coupling_energy = params.couplstr * (dipole_xy.x * shared_q_photon_xy.x + dipole_xy.y * shared_q_photon_xy.y);
            Scalar dipole_self_energy = Scalar(0.5) * params.couplstr * params.couplstr / params.K * (dipole_xy.x * dipole_xy.x + dipole_xy.y * dipole_xy.y);
            
            // Always write energies since they should be the same regardless of particle indexing
#ifdef CAVITY_DEBUG_VERBOSE
            printf("GPU: ENERGY WRITE H=%.6f, C=%.6f, D=%.6f (photon_idx=%d)\n", harmonic_energy, coupling_energy, dipole_self_energy, idx);
#endif
            d_temp_energy[0] = harmonic_energy;
            d_temp_energy[1] = coupling_energy;
            d_temp_energy[2] = dipole_self_energy;
            
            // DO NOT store energy in particle potential energy - prevents double-counting
            Scalar total_energy = harmonic_energy + coupling_energy + dipole_self_energy;
            d_force[idx].w = 0.0;
        }
    }
}

//! Legacy dipole reduction kernel
__global__ void gpu_reduce_dipole_kernel(const Scalar3* d_temp_dipole,
                                          Scalar3* d_dipole_result,
                                          unsigned int num_blocks)
{
    extern __shared__ Scalar3 sdata[];
    
    unsigned int tid = threadIdx.x;
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Initialize shared memory
    sdata[tid] = make_scalar3(0.0, 0.0, 0.0);
    
    // Load data from global memory
    if (idx < num_blocks) {
        sdata[tid] = d_temp_dipole[idx];
    }
    
    __syncthreads();
    
    // Reduction in shared memory
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s && tid + s < num_blocks) {
            sdata[tid].x += sdata[tid + s].x;
            sdata[tid].y += sdata[tid + s].y;
            sdata[tid].z += sdata[tid + s].z;
        }
        __syncthreads();
    }
    
    // Write result for this block to global memory
    if (tid == 0) {
        d_dipole_result[blockIdx.x] = sdata[0];
    }
}

//! GPU function to compute cavity forces
hipError_t gpu_compute_cavity_force(Scalar4* d_force,
                                    const Scalar4* d_pos,
                                    const Scalar* d_charge,
                                    const int3* d_image,
                                    const BoxDim* box,
                                    const unsigned int N,
                                    const cavity_force_params& params,
                                    Scalar* d_temp_energy,
                                    int* d_photon_idx,
                                    Scalar3* d_temp_dipole,
                                    Scalar3* d_dipole_global,
                                    unsigned int L_typeid,
                                    unsigned int& block_size)
{
    // Argument validation
    if (d_force == nullptr || d_pos == nullptr ||
        d_charge == nullptr || d_image == nullptr || box == nullptr ||
        d_temp_energy == nullptr || d_photon_idx == nullptr ||
        d_temp_dipole == nullptr || d_dipole_global == nullptr)
    {
        return hipErrorInvalidValue;
    }
    
    if (N == 0) {
        return hipSuccess;
    }
    
    // Find the photon particle
    {
        unsigned int find_block_size = 1; // Single thread kernel
        unsigned int find_grid_size = 1;
        
        hipLaunchKernelGGL(gpu_find_photon_kernel, 
                          dim3(find_grid_size), dim3(find_block_size), 0, 0,
                          d_pos, N, d_photon_idx, L_typeid);
        
        hipError_t error = hipGetLastError();
        if (error != hipSuccess) return error;
    }
    
    // Get photon index from device
    int photon_idx_host;
    hipError_t error = hipMemcpy(&photon_idx_host, d_photon_idx, sizeof(int), hipMemcpyDeviceToHost);
    if (error != hipSuccess) return error;
    
    // If no photon found, skip computation
    if (photon_idx_host == -1) {
        // Zero forces
        hipMemset(d_force, 0, N * sizeof(Scalar4));
        
        // Zero energies
        Scalar zero_energies[3] = {0.0, 0.0, 0.0};
        hipMemcpy(d_temp_energy, zero_energies, 3 * sizeof(Scalar), hipMemcpyHostToDevice);
        
        return hipSuccess;
    }
    
    // Compute dipole moment
    {
        unsigned int dipole_block_size = block_size;
        unsigned int dipole_grid_size = (N + dipole_block_size - 1) / dipole_block_size;
        unsigned int dipole_shmem_size = dipole_block_size * sizeof(Scalar3);
        
        hipLaunchKernelGGL(gpu_compute_dipole_kernel, 
                          dim3(dipole_grid_size), dim3(dipole_block_size), dipole_shmem_size, 0,
                          d_pos, d_charge, d_image, *box, N, photon_idx_host, d_temp_dipole, L_typeid);
        
        error = hipGetLastError();
        if (error != hipSuccess) return error;
        
        // Reduce partial dipoles
        if (dipole_grid_size > 1) {
            unsigned int reduce_block_size = min(dipole_grid_size, 256u);
            unsigned int reduce_grid_size = 1;
            unsigned int reduce_shmem_size = reduce_block_size * sizeof(Scalar3);
            
            hipLaunchKernelGGL(gpu_reduce_dipole_kernel, 
                              dim3(reduce_grid_size), dim3(reduce_block_size), reduce_shmem_size, 0,
                              d_temp_dipole, d_dipole_global, dipole_grid_size);
            
            error = hipGetLastError();
            if (error != hipSuccess) return error;
        } else {
            // Single block result
            hipMemcpy(d_dipole_global, d_temp_dipole, sizeof(Scalar3), hipMemcpyDeviceToDevice);
        }
    }
    
    // Get dipole result from device
    Scalar3 dipole_result;
    error = hipMemcpy(&dipole_result, d_dipole_global, sizeof(Scalar3), hipMemcpyDeviceToHost);
    if (error != hipSuccess) return error;
    
    // Compute forces
    {
        // IMPORTANT: Don't reset the atomic flag here - let energies persist across multiple kernel launches per timestep
        // The flag will be reset at the beginning of each new timestep in the host code
        
        unsigned int force_grid_size = (N + block_size - 1) / block_size;
        
        hipLaunchKernelGGL(gpu_compute_cavity_force_kernel, 
                          dim3(force_grid_size), dim3(block_size), 0, 0,
                          d_force, d_pos, d_charge, d_image, *box, N, photon_idx_host,
                          dipole_result, params, d_temp_energy, L_typeid);
        
        error = hipGetLastError();
        if (error != hipSuccess) return error;
    }
    
    return hipSuccess;
}

} // end namespace kernel
} // end namespace cavitymd
} // end namespace hoomd 