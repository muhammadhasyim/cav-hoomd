// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CAVITY_FORCE_COMPUTE_GPU_H__
#define __CAVITY_FORCE_COMPUTE_GPU_H__

#include "CavityForceCompute.h"
#include "hoomd/Autotuner.h"

/*! \file CavityForceComputeGPU.h
    \brief Declares the CavityForceComputeGPU class
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

namespace hoomd
{
namespace cavitymd
{

//! Implements cavity force computation on the GPU
/*! Uses CUDA kernels to compute cavity-molecule interaction forces on the GPU.
    The coupling strength can be time-varying via the Variant interface.
    
    \ingroup computes
*/
class PYBIND11_EXPORT CavityForceComputeGPU : public CavityForceCompute
{
public:
    //! Constructs the compute
    CavityForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                          Scalar omegac,
                          std::shared_ptr<Variant> couplstr,
                          Scalar phmass = Scalar(1.0),
                          Scalar damping_ratio = Scalar(0.0));

    //! Destructor
    virtual ~CavityForceComputeGPU();

protected:
    //! Actually compute the forces on the GPU
    virtual void computeForces(uint64_t timestep);

    std::shared_ptr<Autotuner<1>> m_tuner_photon_dipole; //!< Autotuner for photon/dipole kernel
    std::shared_ptr<Autotuner<1>> m_tuner_forces;        //!< Autotuner for force kernel
    
    GPUArray<Scalar> m_temp_energy;    //!< Temporary storage for energy values on GPU
    GPUArray<int> m_photon_idx_gpu;    //!< GPU storage for photon index
    GPUArray<Scalar3> m_dipole_gpu;    //!< GPU storage for dipole moment
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_COMPUTE_GPU_H__ 