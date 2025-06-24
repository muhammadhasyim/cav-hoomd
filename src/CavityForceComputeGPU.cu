// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "hoomd/ParticleData.cuh"
#include "hoomd/HOOMDMath.h"

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

#include <assert.h>

/*! \file CavityForceComputeGPU.cu
    \brief Defines GPU kernel code for computing cavity forces on the GPU
*/

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
        int type = (int)pos.w;
        
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
                                              Scalar* d_virial,
                                              const size_t virial_pitch,
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
        d_force[photon_idx].w = total_energy;
    }
    
    // Compute forces for this particle
    if (idx < N) {
        int type = (int)d_pos[idx].w;
        
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
            
            // Compute virial contribution
            Scalar4 pos_data = d_pos[idx];
            int3 img = d_image[idx];
            Scalar3 L = box.getL();
            
            Scalar3 unwrapped_pos;
            unwrapped_pos.x = pos_data.x + Scalar(img.x) * L.x;
            unwrapped_pos.y = pos_data.y + Scalar(img.y) * L.y;
            unwrapped_pos.z = pos_data.z + Scalar(img.z) * L.z;
            
            // Virial = r_i · F_i
            Scalar virial[6];
            virial[0] = unwrapped_pos.x * force.x; // xx
            virial[1] = unwrapped_pos.x * force.y; // xy  
            virial[2] = unwrapped_pos.x * force.z; // xz
            virial[3] = unwrapped_pos.y * force.y; // yy
            virial[4] = unwrapped_pos.y * force.z; // yz
            virial[5] = unwrapped_pos.z * force.z; // zz
            
            for (unsigned int j = 0; j < 6; j++) {
                d_virial[j * virial_pitch + idx] += virial[j];
            }
            
        } else if ((int)idx == photon_idx) { // Cavity particle
            // Force = -K * q_photon - g * dipole_xy
            Scalar3 photon_force;
            photon_force.x = -params.K * shared_q_photon.x - params.couplstr * dipole_xy.x;
            photon_force.y = -params.K * shared_q_photon.y - params.couplstr * dipole_xy.y;
            photon_force.z = -params.K * shared_q_photon.z;
            
            d_force[idx].x = photon_force.x;
            d_force[idx].y = photon_force.y;
            d_force[idx].z = photon_force.z;
            
            // Compute virial contribution for cavity particle
            Scalar virial[6];
            virial[0] = shared_q_photon.x * photon_force.x; // xx
            virial[1] = shared_q_photon.x * photon_force.y; // xy  
            virial[2] = shared_q_photon.x * photon_force.z; // xz
            virial[3] = shared_q_photon.y * photon_force.y; // yy
            virial[4] = shared_q_photon.y * photon_force.z; // yz
            virial[5] = shared_q_photon.z * photon_force.z; // zz
            
            for (unsigned int j = 0; j < 6; j++) {
                d_virial[j * virial_pitch + idx] += virial[j];
            }
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
    
    // BOUNDS CHECK: Prevent illegal memory access
    if (idx >= N || N == 0 || idx >= 100000) return;
    
    // CRITICAL FIX: Add null pointer check for device arrays
    if (d_pos == nullptr || d_photon_idx == nullptr) return;
    
    // CRITICAL FIX: In HOOMD, particle type is stored as int in .w component
    // Don't use __float_as_int() or __scalar_as_int() - just cast directly
    int particle_type = (int)d_pos[idx].w;
    
    if (particle_type == (int)L_typeid)
    {
        atomicExch(d_photon_idx, (int)idx);
    }
}

//! GPU kernel to compute the dipole moment of molecular particles
__global__ void gpu_compute_dipole_kernel(const Scalar4* d_pos,
                                           const Scalar* d_charge,
                                           const int3* d_image,
                                           const BoxDim box,
                                           const unsigned int N,
                                           const int photon_idx,
                                           Scalar3* d_temp_dipole)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // BOUNDS CHECK: Prevent illegal memory access
    if (idx >= N || N == 0 || idx >= 100000 || photon_idx < 0 || photon_idx >= (int)N) return;
    
    // Use external shared memory for dynamic allocation
    extern __shared__ Scalar3 sdata[];
    
    Scalar3 dipole = make_scalar3(0, 0, 0);
    
    if ((int)idx != photon_idx)
    {
        Scalar4 pos_data = d_pos[idx];
        Scalar charge = d_charge[idx];
        int3 img = d_image[idx];
        
        Scalar3 L = box.getL();
        Scalar3 unwrapped_pos;
        unwrapped_pos.x = pos_data.x + Scalar(img.x) * L.x;
        unwrapped_pos.y = pos_data.y + Scalar(img.y) * L.y;
        unwrapped_pos.z = pos_data.z + Scalar(img.z) * L.z;
        
        dipole.x = charge * unwrapped_pos.x;
        dipole.y = charge * unwrapped_pos.y;
        dipole.z = charge * unwrapped_pos.z;
    }
    
    // Store in shared memory for reduction
    if (threadIdx.x < blockDim.x) {
        sdata[threadIdx.x] = dipole;
    }
    __syncthreads();
    
    // Perform reduction in shared memory
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (threadIdx.x < s && threadIdx.x + s < blockDim.x)
        {
            sdata[threadIdx.x].x += sdata[threadIdx.x + s].x;
            sdata[threadIdx.x].y += sdata[threadIdx.x + s].y;
            sdata[threadIdx.x].z += sdata[threadIdx.x + s].z;
        }
        __syncthreads();
    }
    
    // Write block result to global memory
    if (threadIdx.x == 0)
    {
        d_temp_dipole[blockIdx.x] = sdata[0];
    }
}

//! GPU kernel to compute cavity forces
__global__ void gpu_compute_cavity_force_kernel(Scalar4* d_force,
                                                 Scalar* d_virial,
                                                 const size_t virial_pitch,
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
    
    if ((int)idx != photon_idx)
    {
        // Molecular particle
        Scalar charge = d_charge[idx];
        
        Scalar3 Dq;
        Dq.x = dipole_xy.x / params.K;
        Dq.y = dipole_xy.y / params.K;
        Dq.z = 0;
        
        Scalar3 force;
        force.x = -params.couplstr * charge * Dq.x;
        force.y = -params.couplstr * charge * Dq.y;
        force.z = Scalar(0.0);
        
        d_force[idx].x = force.x;
        d_force[idx].y = force.y;
        d_force[idx].z = force.z;
        
        Scalar4 pos_data = d_pos[idx];
        int3 img = d_image[idx];
        
        Scalar3 L = box.getL();
        Scalar3 unwrapped_pos;
        unwrapped_pos.x = pos_data.x + Scalar(img.x) * L.x;
        unwrapped_pos.y = pos_data.y + Scalar(img.y) * L.y;
        unwrapped_pos.z = pos_data.z + Scalar(img.z) * L.z;
        
        Scalar virial[6];
        virial[0] = unwrapped_pos.x * force.x;
        virial[1] = unwrapped_pos.x * force.y;
        virial[2] = unwrapped_pos.x * force.z;
        virial[3] = unwrapped_pos.y * force.y;
        virial[4] = unwrapped_pos.y * force.z;
        virial[5] = unwrapped_pos.z * force.z;
        
        for (unsigned int j = 0; j < 6; j++)
        {
            d_virial[j * virial_pitch + idx] += virial[j];
        }
    }
    else if ((int)idx == photon_idx)
    {
        // Photon particle
        Scalar4 photon_pos = d_pos[photon_idx];
        Scalar3 q_photon = make_scalar3(photon_pos.x, photon_pos.y, photon_pos.z);
        
        Scalar3 photon_force;
        photon_force.x = -params.K * q_photon.x - params.couplstr * dipole_xy.x;
        photon_force.y = -params.K * q_photon.y - params.couplstr * dipole_xy.y;
        photon_force.z = -params.K * q_photon.z;
        
        d_force[idx].x = photon_force.x;
        d_force[idx].y = photon_force.y;
        d_force[idx].z = photon_force.z;
        
        Scalar virial[6];
        virial[0] = q_photon.x * photon_force.x;
        virial[1] = q_photon.x * photon_force.y;
        virial[2] = q_photon.x * photon_force.z;
        virial[3] = q_photon.y * photon_force.y;
        virial[4] = q_photon.y * photon_force.z;
        virial[5] = q_photon.z * photon_force.z;
        
        for (unsigned int j = 0; j < 6; j++)
        {
            d_virial[j * virial_pitch + idx] += virial[j];
        }
        
        // Energy computation
        if (threadIdx.x == 0 && blockIdx.x == 0)
        {
            Scalar harmonic_energy = Scalar(0.5) * params.K * (q_photon.x * q_photon.x + q_photon.y * q_photon.y + q_photon.z * q_photon.z);
            Scalar coupling_energy = params.couplstr * (dipole_xy.x * q_photon.x + dipole_xy.y * q_photon.y);
            Scalar dipole_self_energy = Scalar(0.5) * params.couplstr * params.couplstr / params.K * (dipole_xy.x * dipole_xy.x + dipole_xy.y * dipole_xy.y);
            
            d_temp_energy[0] = harmonic_energy;
            d_temp_energy[1] = coupling_energy;
            d_temp_energy[2] = dipole_self_energy;
        }
    }
}

//! Legacy dipole reduction kernel
__global__ void gpu_reduce_dipole_kernel(const Scalar3* d_temp_dipole,
                                          Scalar3* d_dipole_result,
                                          unsigned int num_blocks)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Use external shared memory for dynamic allocation
    extern __shared__ Scalar3 sdata[];
    
    Scalar3 dipole = make_scalar3(0, 0, 0);
    if (idx < num_blocks)
    {
        dipole = d_temp_dipole[idx];
    }
    
    if (threadIdx.x < blockDim.x) {
        sdata[threadIdx.x] = dipole;
    }
    __syncthreads();
    
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (threadIdx.x < s && threadIdx.x + s < blockDim.x)
        {
            sdata[threadIdx.x].x += sdata[threadIdx.x + s].x;
            sdata[threadIdx.x].y += sdata[threadIdx.x + s].y;
            sdata[threadIdx.x].z += sdata[threadIdx.x + s].z;
        }
        __syncthreads();
    }
    
    if (threadIdx.x == 0)
    {
        d_dipole_result[0] = sdata[0];
    }
}

//! Optimized driver function
hipError_t gpu_compute_cavity_forces(Scalar4* d_force,
                                      Scalar* d_virial,
                                      const size_t virial_pitch,
                                      const unsigned int N,
                                      const Scalar4* d_pos,
                                      const Scalar* d_charge,
                                      const int3* d_image,
                                      const BoxDim* box,
                                      const cavity_force_params* params,
                                      int* d_photon_idx,
                                      Scalar* d_temp_energy,
                                      Scalar3* d_temp_dipole,
                                      const unsigned int block_size,
                                      unsigned int L_typeid,
                                      Scalar3* d_dipole_global,
                                      bool use_fused_kernel)
{
    // Validate all device pointers
    if (d_force == nullptr || d_virial == nullptr || d_pos == nullptr || 
        d_charge == nullptr || d_image == nullptr || d_photon_idx == nullptr ||
        d_temp_energy == nullptr || d_temp_dipole == nullptr || box == nullptr || params == nullptr) {
        return hipErrorInvalidValue;
    }
    
    if (N == 0) {
        return hipSuccess;
    }
    
    unsigned int num_blocks = (N + block_size - 1) / block_size;
    
    if (use_fused_kernel && d_dipole_global != nullptr) {
        // Use optimized 2-kernel approach (find photon+dipole, then forces)
        hipMemset(d_photon_idx, -1, sizeof(int));
        hipMemset(d_dipole_global, 0, sizeof(Scalar3));
        hipMemset(d_force, 0, N * sizeof(Scalar4));
        
        // Kernel 1: Find photon and compute dipole in one pass
        size_t shared_mem_size = block_size * sizeof(Scalar3);
        
        hipLaunchKernelGGL(gpu_compute_photon_and_dipole, dim3(num_blocks), dim3(block_size), shared_mem_size, 0,
                           d_pos, d_charge, d_image, *box, N, L_typeid, d_photon_idx, d_dipole_global);
        
        hipError_t error = hipGetLastError();
        if (error != hipSuccess) {
            return error;
        }
        
        // Copy photon index to check if found
        int photon_idx_host;
        hipMemcpy(&photon_idx_host, d_photon_idx, sizeof(int), hipMemcpyDeviceToHost);
        
        if (photon_idx_host < 0) {
            hipMemset(d_temp_energy, 0, 3 * sizeof(Scalar));
            return hipSuccess;
        }
        
        // Copy dipole result for kernel parameter (minimal transfer)
        Scalar3 dipole_host;
        hipMemcpy(&dipole_host, d_dipole_global, sizeof(Scalar3), hipMemcpyDeviceToHost);
        Scalar3 dipole_xy = make_scalar3(dipole_host.x, dipole_host.y, 0);
        
        // Kernel 2: Compute forces using optimized kernel
        hipLaunchKernelGGL(gpu_compute_forces_optimized, dim3(num_blocks), dim3(block_size), 0, 0,
                           d_force, d_virial, virial_pitch, d_pos, d_charge, d_image, *box, N,
                           photon_idx_host, dipole_xy, *params, d_temp_energy, L_typeid);
        
        error = hipGetLastError();
        if (error != hipSuccess) {
            return error;
        }
        
        return hipSuccess;
    }
    
    // Fallback to legacy multi-kernel approach
    hipMemset(d_photon_idx, -1, sizeof(int));
    
    hipLaunchKernelGGL(gpu_find_photon_kernel, dim3(num_blocks), dim3(block_size), 0, 0,
                       d_pos, N, d_photon_idx, L_typeid);
    
    hipError_t error = hipGetLastError();
    if (error != hipSuccess) {
        return error;
    }
    
    int photon_idx_host;
    error = hipMemcpy(&photon_idx_host, d_photon_idx, sizeof(int), hipMemcpyDeviceToHost);
    if (error != hipSuccess) {
        return error;
    }
    
    if (photon_idx_host < 0)
    {
        hipMemset(d_force, 0, N * sizeof(Scalar4));
        hipMemset(d_temp_energy, 0, 3 * sizeof(Scalar));
        return hipSuccess;
    }
    
    size_t shared_mem_size = block_size * sizeof(Scalar3);
    
    hipLaunchKernelGGL(gpu_compute_dipole_kernel, dim3(num_blocks), dim3(block_size), shared_mem_size, 0,
                       d_pos, d_charge, d_image, *box, N, photon_idx_host, d_temp_dipole);
    
    error = hipGetLastError();
    if (error != hipSuccess) {
        return error;
    }
    
    unsigned int reduce_blocks = (num_blocks + block_size - 1) / block_size;
    if (reduce_blocks > 1)
    {
        hipLaunchKernelGGL(gpu_reduce_dipole_kernel, dim3(reduce_blocks), dim3(block_size), shared_mem_size, 0,
                           d_temp_dipole, d_temp_dipole, num_blocks);
        
        error = hipGetLastError();
        if (error != hipSuccess) {
            return error;
        }
    }
    
    Scalar3 dipole_host;
    error = hipMemcpy(&dipole_host, d_temp_dipole, sizeof(Scalar3), hipMemcpyDeviceToHost);
    if (error != hipSuccess) {
        return error;
    }
    
    Scalar3 dipole_xy = make_scalar3(dipole_host.x, dipole_host.y, 0);
    
    hipLaunchKernelGGL(gpu_compute_cavity_force_kernel, dim3(num_blocks), dim3(block_size), 0, 0,
                       d_force, d_virial, virial_pitch, d_pos, d_charge, d_image, *box, N, photon_idx_host,
                       dipole_xy, *params, d_temp_energy, L_typeid);
    
    error = hipGetLastError();
    if (error != hipSuccess) {
        return error;
    }
    
    return hipSuccess;
}

} // end namespace kernel
} // end namespace cavitymd
} // end namespace hoomd 