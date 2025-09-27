// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PERTURBATION_FORCE_COMPUTE_H__
#define __PERTURBATION_FORCE_COMPUTE_H__

#include "hoomd/ForceCompute.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include <memory>

/*! \file PerturbationForceCompute.h
    \brief Declares a class for computing external perturbation forces for FDR measurements
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

//! Parameters for perturbation force computation
struct perturbation_force_params
{
    Scalar3 kvector;        //!< Wave vector for perturbation (user-specified)
    Scalar amplitude;       //!< h₀ perturbation amplitude
    Scalar sign;           //!< +1 for (+) clone, -1 for (-) clone
    bool enabled;          //!< Whether perturbation is active
    
    // Computed values for efficiency
    Scalar k_magnitude;     //!< |k| for reference
    Scalar k_dot_k;         //!< |k|² for optimizations
    
#ifndef __HIPCC__
    perturbation_force_params() : amplitude(0.), sign(1.), enabled(false), k_magnitude(0.), k_dot_k(0.) 
    {
        kvector = make_scalar3(0, 0, 0);
    }
    
    perturbation_force_params(Scalar3 _kvector, Scalar _amplitude, Scalar _sign = Scalar(1.0)) 
        : kvector(_kvector), amplitude(_amplitude), sign(_sign), enabled(true)
    {
        k_magnitude = sqrt(kvector.x * kvector.x + kvector.y * kvector.y + kvector.z * kvector.z);
        k_dot_k = k_magnitude * k_magnitude;
    }
    
    pybind11::dict asDict()
    {
        pybind11::dict v;
        v["kvector_x"] = kvector.x;
        v["kvector_y"] = kvector.y;
        v["kvector_z"] = kvector.z;
        v["amplitude"] = amplitude;
        v["sign"] = sign;
        v["enabled"] = enabled;
        v["k_magnitude"] = k_magnitude;
        return v;
    }
#endif
} __attribute__((aligned(16)));

//! Computes external perturbation forces for FDR measurements
/*! Implements the external perturbation force:
    F_pert = sign * h₀ * k * sin(k·r)
    
    from the external potential:
    U_ext = -sign * h₀ * cos(k·r)
    
    where k is the user-specified wave vector, h₀ is the perturbation amplitude,
    and sign is +1 for (+) clone or -1 for (-) clone.
    
    The perturbation is applied only to molecular particles (not cavity particle).
    
    \ingroup computes
*/
class PYBIND11_EXPORT PerturbationForceCompute : public ForceCompute
{
public:
    //! Constructs the compute
    PerturbationForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                            Scalar3 kvector,
                            Scalar amplitude,
                            Scalar sign = Scalar(1.0));

    //! Destructor
    virtual ~PerturbationForceCompute();

    //! Set parameters
    void setParams(Scalar3 kvector, Scalar amplitude, Scalar sign = Scalar(1.0));
    
    //! Get parameters as dictionary
    pybind11::dict getParams();
    
    //! Enable/disable perturbation
    void setEnabled(bool enabled);
    
    //! Get total perturbation energy
    Scalar getPerturbationEnergy();

protected:
    //! Actually compute the forces
    virtual void computeForces(uint64_t timestep);
    
    //! Compute unwrapped positions
    void computeUnwrappedPositions(std::vector<vec3<Scalar>>& unwrapped_pos,
                                   const Scalar4* pos,
                                   const int3* image,
                                   const BoxDim& box,
                                   unsigned int N);

    perturbation_force_params m_params;  //!< Force parameters
    
    //! Energy component
    Scalar m_perturbation_energy;        //!< Total perturbation energy
};

} // end namespace cavitymd
} // end namespace hoomd

#endif // __PERTURBATION_FORCE_COMPUTE_H__
