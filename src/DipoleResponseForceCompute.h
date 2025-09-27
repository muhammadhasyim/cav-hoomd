// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __DIPOLE_RESPONSE_FORCE_COMPUTE_H__
#define __DIPOLE_RESPONSE_FORCE_COMPUTE_H__

#include "hoomd/ForceCompute.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include <memory>

/*! \file DipoleResponseForceCompute.h
    \brief Declares a class for computing electric field forces for dipole moment FDR measurements
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

//! Parameters for dipole response force computation
struct dipole_response_force_params
{
    Scalar3 field_vector;      //!< Electric field direction (normalized)
    Scalar field_strength;     //!< E₀ electric field strength
    Scalar sign;              //!< +1 for (+) clone, -1 for (-) clone
    bool enabled;             //!< Whether electric field is active
    bool exclude_cavity;      //!< Whether to exclude cavity particles
    
    // Computed values for efficiency
    Scalar3 force_per_unit_charge;  //!< sign * E₀ * field_vector
    
#ifndef __HIPCC__
    dipole_response_force_params() : field_strength(0.), sign(1.), enabled(false), exclude_cavity(true)
    {
        field_vector = make_scalar3(0, 0, 1);  // Default z-direction
        force_per_unit_charge = make_scalar3(0, 0, 0);
    }
    
    dipole_response_force_params(Scalar3 _field_vector, Scalar _field_strength, Scalar _sign = Scalar(1.0), bool _exclude_cavity = true) 
        : field_vector(_field_vector), field_strength(_field_strength), sign(_sign), enabled(true), exclude_cavity(_exclude_cavity)
    {
        // Normalize field vector
        Scalar magnitude = sqrt(field_vector.x * field_vector.x + field_vector.y * field_vector.y + field_vector.z * field_vector.z);
        if (magnitude > Scalar(1e-12))
        {
            field_vector.x /= magnitude;
            field_vector.y /= magnitude;
            field_vector.z /= magnitude;
        }
        else
        {
            field_vector = make_scalar3(0, 0, 1);  // Fallback to z-direction
        }
        
        // Compute force per unit charge: F⃗/q = sign * E₀ * Ê
        Scalar force_magnitude = sign * field_strength;
        force_per_unit_charge.x = force_magnitude * field_vector.x;
        force_per_unit_charge.y = force_magnitude * field_vector.y;
        force_per_unit_charge.z = force_magnitude * field_vector.z;
    }
    
    pybind11::dict asDict()
    {
        pybind11::dict v;
        v["field_vector_x"] = field_vector.x;
        v["field_vector_y"] = field_vector.y;
        v["field_vector_z"] = field_vector.z;
        v["field_strength"] = field_strength;
        v["sign"] = sign;
        v["enabled"] = enabled;
        v["exclude_cavity"] = exclude_cavity;
        v["force_per_unit_charge_x"] = force_per_unit_charge.x;
        v["force_per_unit_charge_y"] = force_per_unit_charge.y;
        v["force_per_unit_charge_z"] = force_per_unit_charge.z;
        return v;
    }
#endif
} __attribute__((aligned(16)));

//! Computes electric field forces for dipole moment FDR measurements
/*! Implements the electric field force:
    F⃗ᵢ = sign * E₀ * qᵢ * Ê
    
    from the electric field potential:
    U = -sign * E₀ * Σᵢ qᵢ (Ê · r⃗ᵢ)
    
    where Ê is the normalized electric field direction, E₀ is the field strength,
    qᵢ is the charge on particle i, and sign is +1 for (+) clone or -1 for (-) clone.
    
    This force couples to the total dipole moment μ⃗ = Σᵢ qᵢ r⃗ᵢ and enables measurement
    of dipole moment susceptibility for FDR analysis.
    
    The electric field is applied to all charged particles, with an option to exclude
    cavity particles.
    
    \ingroup computes
*/
class PYBIND11_EXPORT DipoleResponseForceCompute : public ForceCompute
{
public:
    //! Constructs the compute
    DipoleResponseForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                              Scalar3 field_vector,
                              Scalar field_strength,
                              Scalar sign = Scalar(1.0),
                              bool exclude_cavity = true);

    //! Destructor
    virtual ~DipoleResponseForceCompute();

    //! Set parameters
    void setParams(Scalar3 field_vector, Scalar field_strength, Scalar sign = Scalar(1.0), bool exclude_cavity = true);
    
    //! Get parameters as dictionary
    pybind11::dict getParams();
    
    //! Enable/disable electric field
    void setEnabled(bool enabled);
    
    //! Check if electric field is enabled
    bool getEnabled() const { return m_params.enabled; }
    
    //! Get total electric field energy
    Scalar getElectricFieldEnergy() const { return m_electric_field_energy; }

protected:
    //! Actually compute the forces
    virtual void computeForces(uint64_t timestep);

protected:
    dipole_response_force_params m_params;  //!< Parameters for force computation
    Scalar m_electric_field_energy;         //!< Total electric field energy
};

//! Exports the DipoleResponseForceCompute class to python
void export_DipoleResponseForceCompute(pybind11::module& m);

} // end namespace cavitymd
} // end namespace hoomd

#endif
