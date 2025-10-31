// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CAVITY_FORCE_COMPUTE_GPU_H__
#define __CAVITY_FORCE_COMPUTE_GPU_H__

#include "CavityForceCompute.h"
#include "hoomd/Autotuner.h"

/*! \file CavityForceComputeGPU.h
    \brief Declares a class for computing cavity forces on the GPU
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

namespace hoomd
{
namespace cavitymd
{

//! GPU accelerated cavity force computation
/*! Implements the cavity-molecule interaction force on the GPU
    \ingroup computes
*/
class PYBIND11_EXPORT CavityForceComputeGPU : public CavityForceCompute
{
public:
    //! Constructs the compute
    CavityForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                          Scalar omegac,
                          std::shared_ptr<Variant> lambda_coupling,
                          Scalar phmass = Scalar(1.0));

    //! Destructor
    virtual ~CavityForceComputeGPU();

protected:
    //! Actually compute the forces (GPU implementation)
    virtual void computeForces(uint64_t timestep);

private:
    GPUArray<Scalar> m_temp_energy;     //!< Temporary storage for energy components
    GPUArray<Scalar3> m_temp_dipole;    //!< Temporary storage for dipole calculation
    GPUArray<int> m_photon_idx;         //!< Storage for photon particle index
    GPUArray<Scalar3> m_dipole_global;  //!< Global dipole result storage
    
    std::unique_ptr<Autotuner<1>> m_tuner; //!< Autotuner for kernel optimization
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_COMPUTE_GPU_H__ 