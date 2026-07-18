// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CAVITY_FORCE_COMPUTE_H__
#define __CAVITY_FORCE_COMPUTE_H__

#include "CavityForceParams.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include <memory>

/*! \file CavityForceCompute.h
    \brief Declares a class for computing cavity-molecule interaction forces
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>
#include "hoomd/Variant.h"

namespace hoomd
{
namespace cavitymd
{

//! Computes cavity-molecule interaction forces
/*! Implements the force from the cavity Hamiltonian:
    H = (1/2) * K * q² + g * q · d + (g²/2K) * d²
    
    where q is the cavity mode position, d is the molecular dipole moment,
    g is the coupling strength, and K is the cavity spring constant.
    
    The cavity particle must have type name 'L'.
    Only x,y components of the cavity mode and dipole are used.
    
    \ingroup computes
*/
class PYBIND11_EXPORT CavityForceCompute : public ForceCompute
{
public:
    //! Constructs the compute
    CavityForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                       Scalar omegac,
                       std::shared_ptr<Variant> lambda_coupling,
                       Scalar phmass = Scalar(1.0));

    //! Destructor
    virtual ~CavityForceCompute();

    //! Set parameters
    void setParams(Scalar omegac, std::shared_ptr<Variant> lambda_coupling, Scalar phmass = Scalar(1.0));
    
    //! Get parameters as dictionary
    pybind11::dict getParams();
    
    //! Get harmonic energy component
    virtual Scalar getHarmonicEnergy();
    
    //! Get coupling energy component  
    virtual Scalar getCouplingEnergy();
    
    //! Get dipole self-energy component
    virtual Scalar getDipoleSelfEnergy();

protected:
    //! Validate physical parameters before changing force state
    static void validateParams(Scalar omegac,
                               Scalar lambda_coupling,
                               Scalar phmass);

    //! Actually compute the forces
    virtual void computeForces(uint64_t timestep);
    
    //! Find the photon particle index (-1 if not found)
    int findPhotonParticle(const Scalar4* pos_data, unsigned int N);
    
    //! Compute unwrapped positions
    void computeUnwrappedPositions(std::vector<vec3<Scalar>>& unwrapped_pos,
                                   const Scalar4* pos,
                                   const int3* image,
                                   const BoxDim& box,
                                   unsigned int N);
    
    //! Compute molecular dipole moment
    vec3<Scalar> computeDipoleMoment(const std::vector<vec3<Scalar>>& unwrapped_pos,
                                     const Scalar* charge,
                                     unsigned int N,
                                     int photon_idx);
    

    cavity_force_params m_params;  //!< Force parameters

    //!< Dimensionless coupling parameter (lambda) variant
    std::shared_ptr<Variant> m_lambda_coupling;

    
    //! Energy components (now protected for GPU access)
    Scalar m_harmonic_energy;      //!< (1/2) * K * q²
    Scalar m_coupling_energy;      //!< g * q · d
    Scalar m_dipole_self_energy;   //!< (g²/2K) * d²
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_COMPUTE_H__ 