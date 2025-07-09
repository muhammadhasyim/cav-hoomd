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
    Scalar omegac;        //!< Cavity frequency in atomic units
    Scalar couplstr;      //!< Coupling strength in atomic units (evaluated at current timestep)
    Scalar K;             //!< Spring constant (phmass * omegac^2)
    Scalar phmass;        //!< Photon mass
    Scalar damping_ratio; //!< Damping ratio (dimensionless)
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
    // Define block size for grid-stride loop
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Shared memory for reduction 
    __shared__ Scalar3 shared_dipole[256];  // Adjust size as needed
    __shared__ int shared_photon_idx[256];
    
    // Initialize shared memory
    if (threadIdx.x < 256) {
        shared_dipole[threadIdx.x] = make_scalar3(0, 0, 0);
        shared_photon_idx[threadIdx.x] = -1;
    }
    __syncthreads();
    
    // Initialize local variables
    Scalar3 local_dipole = make_scalar3(0, 0, 0);
    int local_photon_idx = -1;
    
    // Grid-stride loop to handle all particles
    for (unsigned int i = idx; i < N; i += blockDim.x * gridDim.x) {
        int type = __scalar_as_int(d_pos[i].w);
        
        if (type == (int)L_typeid) {
            // Found photon particle
            local_photon_idx = (int)i;
        } else {
            // Molecular particle - compute unwrapped position and add to dipole
            Scalar charge = d_charge[i];
            Scalar3 L = box.getL();
            int3 img = d_image[i];
            
            Scalar3 unwrapped_pos;
            unwrapped_pos.x = d_pos[i].x + Scalar(img.x) * L.x;
            unwrapped_pos.y = d_pos[i].y + Scalar(img.y) * L.y;
            unwrapped_pos.z = d_pos[i].z + Scalar(img.z) * L.z;
            
            local_dipole.x += charge * unwrapped_pos.x;
            local_dipole.y += charge * unwrapped_pos.y;
            local_dipole.z += charge * unwrapped_pos.z;
        }
    }
    
    // Store results in shared memory
    shared_dipole[threadIdx.x] = local_dipole;
    shared_photon_idx[threadIdx.x] = local_photon_idx;
    __syncthreads();
    
    // Reduction for dipole moment
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) {
            shared_dipole[threadIdx.x].x += shared_dipole[threadIdx.x + stride].x;
            shared_dipole[threadIdx.x].y += shared_dipole[threadIdx.x + stride].y;
            shared_dipole[threadIdx.x].z += shared_dipole[threadIdx.x + stride].z;
            
            // For photon index, choose the valid one (not -1)
            if (shared_photon_idx[threadIdx.x] == -1 && shared_photon_idx[threadIdx.x + stride] != -1) {
                shared_photon_idx[threadIdx.x] = shared_photon_idx[threadIdx.x + stride];
            }
        }
        __syncthreads();
    }
    
    // Write results to global memory
    if (threadIdx.x == 0) {
        // Use atomic operations for global reduction
        atomicAdd(&d_dipole_result->x, shared_dipole[0].x);
        atomicAdd(&d_dipole_result->y, shared_dipole[0].y);
        atomicAdd(&d_dipole_result->z, shared_dipole[0].z);
        
        if (shared_photon_idx[0] != -1) {
            *d_photon_idx = shared_photon_idx[0];
        }
    }
}

//! GPU kernel for computing cavity forces - minimal version
__global__ void gpu_compute_cavity_forces_minimal(Scalar4* d_force,
                                                   const Scalar4* d_pos,
                                                   const Scalar4* d_vel,
                                                   const Scalar* d_charge,
                                                   const int3* d_image,
                                                   const BoxDim box,
                                                   unsigned int N,
                                                   unsigned int L_typeid,
                                                   int photon_idx,
                                                   Scalar3 dipole_xy,
                                                   Scalar* d_temp_energy,
                                                   cavity_force_params params)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Check bounds
    if (idx >= N) return;
    
    // Share photon position
    __shared__ Scalar3 shared_q_photon_xy;
    __shared__ Scalar3 shared_q_photon;
    if (threadIdx.x == 0) {
        shared_q_photon_xy = dipole_xy;
        
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
            
        } else if ((int)idx == photon_idx) { // Cavity particle
            // Force = -K * q_photon - g * dipole_xy - gamma * velocity
            // where gamma = 2 * damping_ratio * sqrt(K)
            Scalar gamma = 2.0 * params.damping_ratio * sqrt(params.K);
            Scalar4 photon_velocity = d_vel[photon_idx];
            Scalar3 photon_force;
            photon_force.x = -params.K * shared_q_photon.x - params.couplstr * dipole_xy.x - gamma * photon_velocity.x;
            photon_force.y = -params.K * shared_q_photon.y - params.couplstr * dipole_xy.y - gamma * photon_velocity.y;
            photon_force.z = -params.K * shared_q_photon.z - gamma * photon_velocity.z;
            
            d_force[idx].x = photon_force.x;
            d_force[idx].y = photon_force.y;
            d_force[idx].z = photon_force.z;
        }
    }
}

//! GPU kernel for computing cavity forces
__global__ void gpu_compute_cavity_forces(Scalar4* d_force,
                                           const Scalar4* d_pos,
                                           const Scalar4* d_vel,
                                           const Scalar* d_charge,
                                           const int3* d_image,
                                           const BoxDim box,
                                           unsigned int N,
                                           unsigned int L_typeid,
                                           int photon_idx,
                                           Scalar3 q_photon_xy,
                                           Scalar3 dipole_xy,
                                           Scalar* d_temp_energy,
                                           cavity_force_params params)
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
        
        // Get photon velocity and add dissipation term
        // gamma = 2 * damping_ratio * sqrt(K)
        Scalar gamma = 2.0 * params.damping_ratio * sqrt(params.K);
        Scalar4 photon_velocity = d_vel[photon_idx];
        Scalar3 photon_force;
        photon_force.x = -params.K * shared_q_photon.x - params.couplstr * dipole_xy.x - gamma * photon_velocity.x;
        photon_force.y = -params.K * shared_q_photon.y - params.couplstr * dipole_xy.y - gamma * photon_velocity.y;
        photon_force.z = -params.K * shared_q_photon.z - gamma * photon_velocity.z;
        
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
            d_force[idx].w = 0.0;
        }
    }
}

//! Legacy dipole reduction kernel
__global__ void gpu_compute_dipole_reduction(const Scalar4* d_pos,
                                              const Scalar* d_charge,
                                              const int3* d_image,
                                              const BoxDim box,
                                              unsigned int N,
                                              unsigned int L_typeid,
                                              Scalar3* d_dipole_result)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Shared memory for reduction
    __shared__ Scalar3 shared_dipole[256];
    
    // Initialize shared memory
    Scalar3 local_dipole = make_scalar3(0, 0, 0);
    
    // Grid-stride loop
    for (unsigned int i = idx; i < N; i += blockDim.x * gridDim.x) {
        int type = __scalar_as_int(d_pos[i].w);
        if (type != (int)L_typeid) { // Not cavity particle
            Scalar charge = d_charge[i];
            Scalar3 L = box.getL();
            int3 img = d_image[i];
            
            // Compute unwrapped position
            Scalar3 unwrapped_pos;
            unwrapped_pos.x = d_pos[i].x + Scalar(img.x) * L.x;
            unwrapped_pos.y = d_pos[i].y + Scalar(img.y) * L.y;
            unwrapped_pos.z = d_pos[i].z + Scalar(img.z) * L.z;
            
            // Add to dipole moment
            local_dipole.x += charge * unwrapped_pos.x;
            local_dipole.y += charge * unwrapped_pos.y;
            local_dipole.z += charge * unwrapped_pos.z;
        }
    }
    
    // Store in shared memory for reduction
    shared_dipole[threadIdx.x] = local_dipole;
    __syncthreads();
    
    // Block reduction
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) {
            shared_dipole[threadIdx.x].x += shared_dipole[threadIdx.x + stride].x;
            shared_dipole[threadIdx.x].y += shared_dipole[threadIdx.x + stride].y;
            shared_dipole[threadIdx.x].z += shared_dipole[threadIdx.x + stride].z;
        }
        __syncthreads();
    }
    
    // Write block result to global memory
    if (threadIdx.x == 0) {
        atomicAdd(&d_dipole_result->x, shared_dipole[0].x);
        atomicAdd(&d_dipole_result->y, shared_dipole[0].y);
        atomicAdd(&d_dipole_result->z, shared_dipole[0].z);
    }
}

} // end namespace kernel
} // end namespace cavitymd
} // end namespace hoomd 