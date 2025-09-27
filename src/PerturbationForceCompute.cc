// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "PerturbationForceCompute.h"
#include <stdexcept>
#include <iostream>
#include <cmath>

namespace py = pybind11;

/*! \file PerturbationForceCompute.cc
    \brief Contains code for the PerturbationForceCompute class
*/

namespace hoomd
{
namespace cavitymd
{

/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param kvector Wave vector for perturbation
    \param amplitude Perturbation amplitude h₀
    \param sign +1 for (+) clone, -1 for (-) clone
*/
PerturbationForceCompute::PerturbationForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                                                   Scalar3 kvector,
                                                   Scalar amplitude,
                                                   Scalar sign)
    : ForceCompute(sysdef)
{
    m_exec_conf->msg->notice(5) << "Constructing PerturbationForceCompute" << std::endl;
    
    // Initialize parameters
    m_params = perturbation_force_params(kvector, amplitude, sign);
    
    // Initialize energy
    m_perturbation_energy = Scalar(0.0);
    
    std::cout << "PerturbationForceCompute initialized: "
              << "kvector=[" << m_params.kvector.x << ", " << m_params.kvector.y << ", " << m_params.kvector.z << "], "
              << "amplitude=" << m_params.amplitude << ", "
              << "sign=" << m_params.sign << ", "
              << "|k|=" << m_params.k_magnitude << std::endl;
}

PerturbationForceCompute::~PerturbationForceCompute()
{
    m_exec_conf->msg->notice(5) << "Destroying PerturbationForceCompute" << std::endl;
}

void PerturbationForceCompute::setParams(Scalar3 kvector, Scalar amplitude, Scalar sign)
{
    m_params = perturbation_force_params(kvector, amplitude, sign);
}

pybind11::dict PerturbationForceCompute::getParams()
{
    return m_params.asDict();
}

void PerturbationForceCompute::setEnabled(bool enabled)
{
    m_params.enabled = enabled;
    if (enabled)
        std::cout << "PerturbationForce enabled" << std::endl;
    else
        std::cout << "PerturbationForce disabled" << std::endl;
}

Scalar PerturbationForceCompute::getPerturbationEnergy()
{
    return m_perturbation_energy;
}

void PerturbationForceCompute::computeUnwrappedPositions(std::vector<vec3<Scalar>>& unwrapped_pos,
                                                        const Scalar4* pos,
                                                        const int3* image,
                                                        const BoxDim& box,
                                                        unsigned int N)
{
    unwrapped_pos.resize(N);
    Scalar3 L = box.getL();
    
    for (unsigned int i = 0; i < N; i++)
    {
        unwrapped_pos[i].x = pos[i].x + Scalar(image[i].x) * L.x;
        unwrapped_pos[i].y = pos[i].y + Scalar(image[i].y) * L.y;
        unwrapped_pos[i].z = pos[i].z + Scalar(image[i].z) * L.z;
    }
}

/*! \param timestep Current time step
    Compute perturbation forces: F = sign * h₀ * k * sin(k·r)
*/
void PerturbationForceCompute::computeForces(uint64_t timestep)
{
    if (!m_params.enabled)
        return;
        
    // Access particle data
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
    
    // Zero forces first
    memset((void*)h_force.data, 0, sizeof(Scalar4) * m_force.getNumElements());
    
    unsigned int N = m_pdata->getN();
    BoxDim box = m_pdata->getGlobalBox();
    
    // Get unwrapped positions
    std::vector<vec3<Scalar>> unwrapped_pos;
    computeUnwrappedPositions(unwrapped_pos, h_pos.data, h_image.data, box, N);
    
    // Get the typeid for cavity particle 'L' to exclude it
    unsigned int L_typeid;
    try {
        L_typeid = m_pdata->getTypeByName("L");
    } catch (...) {
        // No cavity particle type, that's fine
        L_typeid = UINT_MAX;  // Invalid type
    }
    
    // Initialize energy accumulator
    m_perturbation_energy = Scalar(0.0);
    
    // Compute forces on all particles (except cavity particle)
    for (unsigned int i = 0; i < N; i++)
    {
        int type = __scalar_as_int(h_pos.data[i].w);
        
        // Skip cavity particle
        if (type == (int)L_typeid)
            continue;
            
        // Compute k·r
        Scalar k_dot_r = m_params.kvector.x * unwrapped_pos[i].x +
                         m_params.kvector.y * unwrapped_pos[i].y +
                         m_params.kvector.z * unwrapped_pos[i].z;
        
        // Compute force: F = sign * h₀ * k * sin(k·r)
        Scalar sin_kr = sin(k_dot_r);
        Scalar cos_kr = cos(k_dot_r);
        Scalar force_magnitude = m_params.sign * m_params.amplitude * sin_kr;
        
        // Apply force in k direction
        h_force.data[i].x = force_magnitude * m_params.kvector.x;
        h_force.data[i].y = force_magnitude * m_params.kvector.y;
        h_force.data[i].z = force_magnitude * m_params.kvector.z;
        h_force.data[i].w = 0.0;  // No particle potential energy contribution
        
        // Accumulate total perturbation energy: U = -sign * h₀ * cos(k·r)
        m_perturbation_energy += -m_params.sign * m_params.amplitude * cos_kr;
    }
}

namespace detail
{
void export_PerturbationForceCompute(pybind11::module& m)
{
    pybind11::class_<PerturbationForceCompute, ForceCompute, std::shared_ptr<PerturbationForceCompute>>(
        m, "PerturbationForceCompute")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar3, Scalar, Scalar>())
        .def("setParams", &PerturbationForceCompute::setParams)
        .def("getParams", &PerturbationForceCompute::getParams)
        .def("setEnabled", &PerturbationForceCompute::setEnabled)
        .def("getPerturbationEnergy", &PerturbationForceCompute::getPerturbationEnergy);
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd
