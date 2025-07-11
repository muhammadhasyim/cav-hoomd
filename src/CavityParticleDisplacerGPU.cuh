// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"
#include "hoomd/VectorMath.h"

namespace hoomd
{
namespace cavitymd
{
namespace kernel
{

// Kernel to sum partial dipole moments
__global__ void gpu_sum_dipole_kernel(
    Scalar3* dipole_partial_sum,
    unsigned int num_blocks);

// Kernel to displace the cavity particle
__global__ void gpu_displace_cavity_particle_kernel(
    unsigned int N,
    unsigned int L_typeid,
    Scalar g,
    Scalar K,
    Scalar4* pos,
    const Scalar* charge,
    int3* image,
    const BoxDim box);

// Wrapper function to call the displacement kernel
cudaError_t gpu_cavity_displace(
    unsigned int N,
    unsigned int L_typeid,
    Scalar g,
    Scalar K,
    Scalar4* pos,
    const Scalar* charge,
    int3* image,
    const BoxDim& box);

} // end namespace kernel
} // end namespace cavitymd
} // end namespace hoomd 