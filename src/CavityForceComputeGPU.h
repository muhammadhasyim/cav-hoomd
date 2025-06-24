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

//! Computes cavity-molecule interaction forces on the GPU
/*! GPU implementation of cavity force computation using CUDA kernels.
    
    \ingroup computes
*/
class PYBIND11_EXPORT CavityForceComputeGPU : public CavityForceCompute
{
public:
    //! Constructs the compute
    CavityForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                          Scalar omegac,
                          Scalar couplstr, 
                          Scalar phmass = Scalar(1.0));

    //! Destructor
    virtual ~CavityForceComputeGPU();

protected:
    //! Actually compute the forces on the GPU
    virtual void computeForces(uint64_t timestep);
    
private:
    std::shared_ptr<Autotuner<1>> m_tuner;  //!< Autotuner for block size
    
    // Device memory for temporary storage
    GPUArray<Scalar> m_temp_energy;     //!< Temporary energy storage on GPU
    GPUArray<Scalar3> m_temp_dipole;    //!< Temporary dipole storage on GPU  
    GPUArray<int> m_photon_idx;         //!< Photon particle index on GPU
    GPUArray<Scalar3> m_dipole_global;  //!< Global dipole storage for fused kernel
    
    bool m_supports_cooperative;        //!< Whether device supports cooperative kernel launches
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_COMPUTE_GPU_H__ 