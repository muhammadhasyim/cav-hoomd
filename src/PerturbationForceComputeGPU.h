// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PERTURBATION_FORCE_COMPUTE_GPU_H__
#define __PERTURBATION_FORCE_COMPUTE_GPU_H__

#include "PerturbationForceCompute.h"
#include "hoomd/Autotuner.h"

/*! \file PerturbationForceComputeGPU.h
    \brief Declares a class for computing perturbation forces on the GPU
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

namespace hoomd
{
namespace cavitymd
{

//! GPU accelerated perturbation force computation
/*! Implements the external perturbation force on the GPU for FDR measurements
    \ingroup computes
*/
class PYBIND11_EXPORT PerturbationForceComputeGPU : public PerturbationForceCompute
{
public:
    //! Constructs the compute
    PerturbationForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                               Scalar3 kvector,
                               Scalar amplitude,
                               Scalar sign = Scalar(1.0));

    //! Destructor
    virtual ~PerturbationForceComputeGPU();

protected:
    //! Actually compute the forces (GPU implementation)
    virtual void computeForces(uint64_t timestep);

private:
    GPUArray<Scalar> m_temp_energy;     //!< Temporary storage for energy reduction
    
    std::unique_ptr<Autotuner<1>> m_tuner; //!< Autotuner for kernel optimization
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __PERTURBATION_FORCE_COMPUTE_GPU_H__
