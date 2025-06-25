// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CAVITY_FORCE_COMPUTE_H__
#define __CAVITY_FORCE_COMPUTE_H__

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

namespace hoomd
{
namespace cavitymd
{

//! Parameters for cavity force computation
struct cavity_force_params
{
    Scalar omegac;     //!< Cavity frequency in atomic units
    Scalar couplstr;   //!< Coupling strength in atomic units  
    Scalar K;          //!< Spring constant (phmass * omegac^2)
    Scalar phmass;     //!< Photon mass
    
#ifndef __HIPCC__
    cavity_force_params() : omegac(0.), couplstr(0.), K(0.), phmass(1.) {}
    
    cavity_force_params(Scalar _omegac, Scalar _couplstr, Scalar _phmass) 
        : omegac(_omegac), couplstr(_couplstr), phmass(_phmass)
    {
        K = phmass * omegac * omegac;
    }
    
    pybind11::dict asDict()
    {
        pybind11::dict v;
        v["omegac"] = omegac;
        v["couplstr"] = couplstr;
        v["K"] = K;
        v["phmass"] = phmass;
        return v;
    }
#endif
} __attribute__((aligned(16)));

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
                       Scalar couplstr, 
                       Scalar phmass = Scalar(1.0));

    //! Destructor
    virtual ~CavityForceCompute();

    //! Set parameters
    void setParams(Scalar omegac, Scalar couplstr, Scalar phmass = Scalar(1.0));
    
    //! Get parameters as dictionary
    pybind11::dict getParams();
    
    //! Get harmonic energy component
    Scalar getHarmonicEnergy();
    
    //! Get coupling energy component  
    Scalar getCouplingEnergy();
    
    //! Get dipole self-energy component
    Scalar getDipoleSelfEnergy();

protected:
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
    
    //! Energy components (now protected for GPU access)
    Scalar m_harmonic_energy;      //!< (1/2) * K * q²
    Scalar m_coupling_energy;      //!< g * q · d
    Scalar m_dipole_self_energy;   //!< (g²/2K) * d²
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_COMPUTE_H__ 