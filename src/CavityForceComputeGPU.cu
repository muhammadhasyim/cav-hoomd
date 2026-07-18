/*! \file CavityForceComputeGPU.cu
    \brief Implements device-side cavity force and energy computation.
*/

#include "CavityForceParams.h"
#include "hoomd/ParticleData.cuh"

namespace hoomd
{
namespace cavitymd
{
namespace kernel
{

constexpr unsigned int cavity_block_size = 256;
constexpr unsigned int device_not_local = 0xffffffffu;

//! Refresh the local photon index without leaving the device.
__global__ void gpu_update_photon_index_kernel(unsigned int* d_photon_idx,
                                               const unsigned int* d_rtag,
                                               const unsigned int photon_tag,
                                               const bool has_photon)
    {
    if (threadIdx.x == 0 && blockIdx.x == 0)
        {
        d_photon_idx[0] = has_photon ? d_rtag[photon_tag] : device_not_local;
        }
    }

//! Compute one deterministic dipole partial sum per particle block.
__global__ void gpu_compute_dipole_partials(const Scalar4* d_pos,
                                            const Scalar* d_charge,
                                            const int3* d_image,
                                            const BoxDim box,
                                            const unsigned int N,
                                            const unsigned int* d_photon_idx,
                                            Scalar3* d_partial)
    {
    extern __shared__ Scalar3 shared_dipole[];

    const unsigned int thread_idx = threadIdx.x;
    const unsigned int particle_idx = blockIdx.x * blockDim.x + thread_idx;
    const unsigned int photon_idx = d_photon_idx[0];
    Scalar3 contribution = make_scalar3(Scalar(0.0), Scalar(0.0), Scalar(0.0));

    if (particle_idx < N && particle_idx != photon_idx)
        {
        const Scalar4 position = d_pos[particle_idx];
        const int3 image = d_image[particle_idx];
        const Scalar3 box_length = box.getL();
        const Scalar charge = d_charge[particle_idx];
        contribution.x
            = charge * (position.x + Scalar(image.x) * box_length.x);
        contribution.y
            = charge * (position.y + Scalar(image.y) * box_length.y);
        contribution.z
            = charge * (position.z + Scalar(image.z) * box_length.z);
        }

    shared_dipole[thread_idx] = contribution;
    __syncthreads();

    for (unsigned int offset = blockDim.x / 2; offset > 0; offset >>= 1)
        {
        if (thread_idx < offset)
            {
            shared_dipole[thread_idx].x += shared_dipole[thread_idx + offset].x;
            shared_dipole[thread_idx].y += shared_dipole[thread_idx + offset].y;
            shared_dipole[thread_idx].z += shared_dipole[thread_idx + offset].z;
            }
        __syncthreads();
        }

    if (thread_idx == 0)
        {
        d_partial[blockIdx.x] = shared_dipole[0];
        }
    }

//! Reduce arbitrary many partials into one partial per block.
__global__ void gpu_reduce_dipole_partials(const Scalar3* d_input,
                                           Scalar3* d_output,
                                           const unsigned int input_count)
    {
    extern __shared__ Scalar3 shared_dipole[];

    const unsigned int thread_idx = threadIdx.x;
    const unsigned int input_idx = blockIdx.x * blockDim.x + thread_idx;
    Scalar3 value = make_scalar3(Scalar(0.0), Scalar(0.0), Scalar(0.0));
    if (input_idx < input_count)
        {
        value = d_input[input_idx];
        }
    shared_dipole[thread_idx] = value;
    __syncthreads();

    for (unsigned int offset = blockDim.x / 2; offset > 0; offset >>= 1)
        {
        if (thread_idx < offset)
            {
            shared_dipole[thread_idx].x += shared_dipole[thread_idx + offset].x;
            shared_dipole[thread_idx].y += shared_dipole[thread_idx + offset].y;
            shared_dipole[thread_idx].z += shared_dipole[thread_idx + offset].z;
            }
        __syncthreads();
        }

    if (thread_idx == 0)
        {
        d_output[blockIdx.x] = shared_dipole[0];
        }
    }

//! Compute every force component and all cavity energies from the final dipole.
__global__ void gpu_compute_force_and_energy(Scalar4* d_force,
                                             const Scalar4* d_pos,
                                             const Scalar* d_charge,
                                             const int3* d_image,
                                             const BoxDim box,
                                             const cavity_force_params params,
                                             const Scalar3* d_dipole,
                                             Scalar* d_energy,
                                             const unsigned int N,
                                             const unsigned int* d_photon_idx)
    {
    const unsigned int particle_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (particle_idx >= N)
        {
        return;
        }

    const unsigned int photon_idx = d_photon_idx[0];
    if (photon_idx >= N)
        {
        d_force[particle_idx]
            = make_scalar4(Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0));
        if (particle_idx == 0)
            {
            d_energy[0] = Scalar(0.0);
            d_energy[1] = Scalar(0.0);
            d_energy[2] = Scalar(0.0);
            }
        return;
        }

    const Scalar3 dipole = d_dipole[0];
    const Scalar4 photon_position = d_pos[photon_idx];
    const int3 photon_image = d_image[photon_idx];
    const Scalar3 box_length = box.getL();
    const Scalar3 photon_coordinate
        = make_scalar3(photon_position.x + Scalar(photon_image.x) * box_length.x,
                       photon_position.y + Scalar(photon_image.y) * box_length.y,
                       photon_position.z + Scalar(photon_image.z) * box_length.z);
    const Scalar coupling = params.lambda_coupling * params.omegac;

    if (particle_idx == 0)
        {
        d_energy[0]
            = Scalar(0.5) * params.K
              * (photon_coordinate.x * photon_coordinate.x
                 + photon_coordinate.y * photon_coordinate.y
                 + photon_coordinate.z * photon_coordinate.z);
        d_energy[1]
            = coupling
              * (dipole.x * photon_coordinate.x
                 + dipole.y * photon_coordinate.y);
        d_energy[2]
            = Scalar(0.5) * coupling * coupling / params.K
              * (dipole.x * dipole.x + dipole.y * dipole.y);
        }

    Scalar4 force = make_scalar4(Scalar(0.0),
                                 Scalar(0.0),
                                 Scalar(0.0),
                                 Scalar(0.0));
    if (particle_idx == photon_idx)
        {
        force.x = -params.K * photon_coordinate.x - coupling * dipole.x;
        force.y = -params.K * photon_coordinate.y - coupling * dipole.y;
        force.z = -params.K * photon_coordinate.z;
        }
    else
        {
        const Scalar charge = d_charge[particle_idx];
        const Scalar dipole_scale
            = params.lambda_coupling / (params.phmass * params.omegac);
        force.x = -coupling * charge
                  * (photon_coordinate.x + dipole_scale * dipole.x);
        force.y = -coupling * charge
                  * (photon_coordinate.y + dipole_scale * dipole.y);
        }

    d_force[particle_idx] = force;
    }

hipError_t gpu_update_photon_index(unsigned int* d_photon_idx,
                                   const unsigned int* d_rtag,
                                   unsigned int photon_tag,
                                   bool has_photon)
    {
    hipLaunchKernelGGL(gpu_update_photon_index_kernel,
                      dim3(1),
                      dim3(1),
                      0,
                      0,
                      d_photon_idx,
                      d_rtag,
                      photon_tag,
                      has_photon);
    return hipGetLastError();
    }

hipError_t gpu_compute_cavity_forces(Scalar4* d_force,
                                     const Scalar4* d_pos,
                                     const Scalar* d_charge,
                                     const int3* d_image,
                                     const BoxDim box,
                                     const cavity_force_params params,
                                     Scalar* d_energy,
                                     const unsigned int* d_photon_idx,
                                     Scalar3* d_partial_a,
                                     Scalar3* d_partial_b,
                                     unsigned int N,
                                     unsigned int partial_capacity)
    {
    if (N == 0)
        {
        return hipMemset(d_energy, 0, 3 * sizeof(Scalar));
        }

    unsigned int partial_count
        = (N + cavity_block_size - 1) / cavity_block_size;
    if (partial_count > partial_capacity)
        {
        return hipErrorInvalidValue;
        }

    constexpr size_t shared_bytes = cavity_block_size * sizeof(Scalar3);
    hipLaunchKernelGGL(gpu_compute_dipole_partials,
                      dim3(partial_count),
                      dim3(cavity_block_size),
                      shared_bytes,
                      0,
                      d_pos,
                      d_charge,
                      d_image,
                      box,
                      N,
                      d_photon_idx,
                      d_partial_a);
    hipError_t error = hipGetLastError();
    if (error != hipSuccess)
        {
        return error;
        }

    const Scalar3* reduction_input = d_partial_a;
    Scalar3* reduction_output = d_partial_b;
    while (partial_count > 1)
        {
        const unsigned int output_count
            = (partial_count + cavity_block_size - 1) / cavity_block_size;
        hipLaunchKernelGGL(gpu_reduce_dipole_partials,
                          dim3(output_count),
                          dim3(cavity_block_size),
                          shared_bytes,
                          0,
                          reduction_input,
                          reduction_output,
                          partial_count);
        error = hipGetLastError();
        if (error != hipSuccess)
            {
            return error;
            }

        partial_count = output_count;
        reduction_input = reduction_output;
        reduction_output
            = reduction_output == d_partial_a ? d_partial_b : d_partial_a;
        }

    const unsigned int force_grid_size
        = (N + cavity_block_size - 1) / cavity_block_size;
    hipLaunchKernelGGL(gpu_compute_force_and_energy,
                      dim3(force_grid_size),
                      dim3(cavity_block_size),
                      0,
                      0,
                      d_force,
                      d_pos,
                      d_charge,
                      d_image,
                      box,
                      params,
                      reduction_input,
                      d_energy,
                      N,
                      d_photon_idx);
    return hipGetLastError();
    }

} // end namespace kernel
} // end namespace cavitymd
} // end namespace hoomd
