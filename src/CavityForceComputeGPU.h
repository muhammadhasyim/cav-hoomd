// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CAVITY_FORCE_COMPUTE_GPU_H__
#define __CAVITY_FORCE_COMPUTE_GPU_H__

#include "CavityForceCompute.h"

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

    Exactly one particle of type L must exist. Particle-count changes trigger
    revalidation, but changing particle types without changing the particle
    count is unsupported after attachment.

    Multi-rank execution is unsupported because the dipole reduction is local
    to one GPU.

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

    //! Get harmonic energy, copying the device values only when needed
    Scalar getHarmonicEnergy() override;

    //! Get coupling energy, copying the device values only when needed
    Scalar getCouplingEnergy() override;

    //! Get dipole self-energy, copying the device values only when needed
    Scalar getDipoleSelfEnergy() override;

protected:
    //! Actually compute the forces (GPU implementation)
    void computeForces(uint64_t timestep) override;

private:
    //! Locate the photon once and cache its stable tag and current local index
    void initializePhoton();

    //! Refresh the local photon index after ParticleData reorders particles
    void updatePhotonIndex();

    //! Revalidate the unique photon after global particle-count changes
    void handleParticleNumberChange();

    //! Transfer and cache the three device energy values
    void updateEnergyCache();

    GPUArray<Scalar> m_device_energy; //!< Three energy components on the device
    GPUArray<unsigned int> m_device_photon_idx; //!< Current photon index on the device
    GPUArray<Scalar3> m_partial_dipole_a; //!< First ping-pong reduction buffer
    GPUArray<Scalar3> m_partial_dipole_b; //!< Second ping-pong reduction buffer

    unsigned int m_photon_tag; //!< Stable tag of the L photon
    bool m_has_photon; //!< True when an L photon was found at construction
    bool m_energy_cache_valid; //!< True after the current energies reach the host
    size_t m_reduction_capacity; //!< Allocated partial dipole elements
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_COMPUTE_GPU_H__ 