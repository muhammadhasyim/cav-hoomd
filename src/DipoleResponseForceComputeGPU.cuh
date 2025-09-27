// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __DIPOLE_RESPONSE_FORCE_COMPUTE_GPU_CUH__
#define __DIPOLE_RESPONSE_FORCE_COMPUTE_GPU_CUH__

#include "DipoleResponseForceCompute.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/ParticleData.cuh"

#ifdef ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

/*! \file DipoleResponseForceComputeGPU.cuh
    \brief Declares GPU kernel code for computing dipole response forces
*/

namespace hoomd
{
namespace cavitymd
{

//! Kernel launcher for dipole response force computation on GPU
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
                                             unsigned int num_blocks);

} // end namespace cavitymd
} // end namespace hoomd

#endif
