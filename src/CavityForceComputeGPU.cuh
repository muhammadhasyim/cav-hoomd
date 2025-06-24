// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CAVITY_FORCE_COMPUTE_GPU_CUH__
#define __CAVITY_FORCE_COMPUTE_GPU_CUH__

#include "CavityForceCompute.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/ParticleData.cuh"

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

/*! \file CavityForceComputeGPU.cuh
    \brief Declares GPU kernel code for computing cavity forces on the GPU
*/

namespace hoomd
{
namespace cavitymd
{
namespace kernel
{

//! GPU kernel driver for computing cavity forces
#ifdef ENABLE_HIP
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
                                      Scalar3* d_dipole_global = nullptr,
                                      bool use_fused_kernel = false);
#endif

} // end namespace kernel
} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_COMPUTE_GPU_CUH__ 