// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "DipoleResponseForceCompute.h"
#include <stdexcept>
#include <iostream>
#include <cmath>

namespace py = pybind11;

/*! \file DipoleResponseForceCompute.cc
    \brief Contains code for the DipoleResponseForceCompute class
*/

namespace hoomd
{
namespace cavitymd
{

/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param field_vector Electric field direction vector
    \param field_strength Electric field strength E₀
    \param sign +1 for (+) clone, -1 for (-) clone
    \param exclude_cavity Whether to exclude cavity particles
*/
DipoleResponseForceCompute::DipoleResponseForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                                                       Scalar3 field_vector,
                                                       Scalar field_strength,
                                                       Scalar sign,
                                                       bool exclude_cavity)
    : ForceCompute(sysdef)
{
    m_exec_conf->msg->notice(5) << "Constructing DipoleResponseForceCompute" << std::endl;
    
    // Initialize parameters
    m_params = dipole_response_force_params(field_vector, field_strength, sign, exclude_cavity);
    
    // Initialize energy
    m_electric_field_energy = Scalar(0.0);
    
    std::cout << "DipoleResponseForceCompute initialized: "
              << "field_vector=[" << m_params.field_vector.x << ", " << m_params.field_vector.y << ", " << m_params.field_vector.z << "], "
              << "field_strength=" << m_params.field_strength << ", "
              << "sign=" << m_params.sign << ", "
              << "exclude_cavity=" << (m_params.exclude_cavity ? "true" : "false") << std::endl;
}

DipoleResponseForceCompute::~DipoleResponseForceCompute()
{
    m_exec_conf->msg->notice(5) << "Destroying DipoleResponseForceCompute" << std::endl;
}

void DipoleResponseForceCompute::setParams(Scalar3 field_vector, Scalar field_strength, Scalar sign, bool exclude_cavity)
{
    m_params = dipole_response_force_params(field_vector, field_strength, sign, exclude_cavity);
}

pybind11::dict DipoleResponseForceCompute::getParams()
{
    return m_params.asDict();
}

void DipoleResponseForceCompute::setEnabled(bool enabled)
{
    m_params.enabled = enabled;
}

/*! \param timestep Current time step
*/
void DipoleResponseForceCompute::computeForces(uint64_t timestep)
{
    // If electric field is disabled, set all forces to zero and return
    if (!m_params.enabled)
    {
        // Access force data
        ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
        unsigned int N = m_pdata->getN();
        
        // Zero all forces
        for (unsigned int i = 0; i < N; i++)
        {
            h_force.data[i].x = Scalar(0.0);
            h_force.data[i].y = Scalar(0.0);
            h_force.data[i].z = Scalar(0.0);
            h_force.data[i].w = Scalar(0.0);
        }
        
        m_electric_field_energy = Scalar(0.0);
        return;
    }
    
    // Start profiling (if available)
    // Note: Profiling code removed for simplicity
    
    // Access particle data
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_charge(m_pdata->getCharges(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
    
    // Get unwrapped positions for energy calculation
    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::read);
    const BoxDim& box = m_pdata->getBox();
    
    unsigned int N = m_pdata->getN();
    
    // Find cavity particle type ID (assumed to be type with name 'L')
    unsigned int L_typeid = m_pdata->getTypeByName("L");
    
    // Initialize energy accumulator
    m_electric_field_energy = Scalar(0.0);
    
    // Compute forces on all particles
    for (unsigned int i = 0; i < N; i++)
    {
        // Get particle type
        int type = __scalar_as_int(h_pos.data[i].w);
        
        // Skip cavity particle if requested
        if (m_params.exclude_cavity && (unsigned int)type == L_typeid)
        {
            h_force.data[i].x = Scalar(0.0);
            h_force.data[i].y = Scalar(0.0);
            h_force.data[i].z = Scalar(0.0);
            h_force.data[i].w = Scalar(0.0);
            continue;
        }
        
        // Get particle charge
        Scalar q_i = h_charge.data[i];
        
        // Only apply force to charged particles
        if (fabs(q_i) > Scalar(1e-12))
        {
            // Apply force: F⃗ᵢ = qᵢ * (sign * E₀ * Ê)
            h_force.data[i].x = q_i * m_params.force_per_unit_charge.x;
            h_force.data[i].y = q_i * m_params.force_per_unit_charge.y;
            h_force.data[i].z = q_i * m_params.force_per_unit_charge.z;
            h_force.data[i].w = Scalar(0.0);  // No single-particle potential energy
            
            // Calculate unwrapped position for energy computation
            Scalar3 unwrapped_pos = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
            Scalar3 L = box.getL();
            unwrapped_pos.x += h_image.data[i].x * L.x;
            unwrapped_pos.y += h_image.data[i].y * L.y;
            unwrapped_pos.z += h_image.data[i].z * L.z;
            
            // Accumulate electric field energy: U = -sign * E₀ * Σᵢ qᵢ (Ê · r⃗ᵢ)
            Scalar field_dot_r = m_params.field_vector.x * unwrapped_pos.x +
                                 m_params.field_vector.y * unwrapped_pos.y +
                                 m_params.field_vector.z * unwrapped_pos.z;
            m_electric_field_energy += -m_params.sign * m_params.field_strength * q_i * field_dot_r;
        }
        else
        {
            // No force on uncharged particles
            h_force.data[i].x = Scalar(0.0);
            h_force.data[i].y = Scalar(0.0);
            h_force.data[i].z = Scalar(0.0);
            h_force.data[i].w = Scalar(0.0);
        }
    }
    
    // End profiling (if available)
    // Note: Profiling code removed for simplicity
}

namespace detail
{
void export_DipoleResponseForceCompute(pybind11::module& m)
{
    pybind11::class_<DipoleResponseForceCompute, ForceCompute, std::shared_ptr<DipoleResponseForceCompute>>(
        m, "DipoleResponseForceCompute")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar3, Scalar, Scalar, bool>())
        .def("setParams", &DipoleResponseForceCompute::setParams)
        .def("getParams", &DipoleResponseForceCompute::getParams)
        .def("setEnabled", &DipoleResponseForceCompute::setEnabled)
        .def("getEnabled", &DipoleResponseForceCompute::getEnabled)
        .def("getElectricFieldEnergy", &DipoleResponseForceCompute::getElectricFieldEnergy);
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd
