// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __DIPOLE_RESPONSE_FORCE_COMPUTE_GPU_H__
#define __DIPOLE_RESPONSE_FORCE_COMPUTE_GPU_H__

#include "DipoleResponseForceCompute.h"
#include "hoomd/Autotuner.h"

/*! \file DipoleResponseForceComputeGPU.h
    \brief Declares a class for computing electric field forces on the GPU for dipole moment FDR measurements
*/

namespace hoomd
{
namespace cavitymd
{

//! GPU implementation of DipoleResponseForceCompute
/*! GPU accelerated implementation of electric field forces for dipole moment FDR measurements.
    
    This class inherits from DipoleResponseForceCompute and overrides computeForces()
    to use GPU kernels for efficient force computation.
    
    \ingroup computes
*/
class PYBIND11_EXPORT DipoleResponseForceComputeGPU : public DipoleResponseForceCompute
{
public:
    //! Constructs the compute
    DipoleResponseForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                                  Scalar3 field_vector,
                                  Scalar field_strength,
                                  Scalar sign = Scalar(1.0),
                                  bool exclude_cavity = true);

    //! Destructor
    virtual ~DipoleResponseForceComputeGPU() { }

    //! Set autotuner parameters
    virtual void setAutotunerParams(bool enable, unsigned int period);

protected:
    //! Actually compute the forces on the GPU
    virtual void computeForces(uint64_t timestep);

private:
    std::unique_ptr<Autotuner<1>> m_tuner; //!< Autotuner for block size optimization
    GPUArray<Scalar> m_energy_sum;     //!< Partial energy sums for reduction
};

//! Exports the DipoleResponseForceComputeGPU class to python
void export_DipoleResponseForceComputeGPU(pybind11::module& m);

} // end namespace cavitymd
} // end namespace hoomd

#endif
