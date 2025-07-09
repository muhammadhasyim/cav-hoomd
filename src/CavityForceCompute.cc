// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "CavityForceCompute.h"
#include <stdexcept>
#include <iostream>

namespace py = pybind11;

/*! \file CavityForceCompute.cc
    \brief Contains code for the CavityForceCompute class
*/

namespace hoomd
{
namespace cavitymd
{

/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param omegac Cavity frequency in atomic units
    \param couplstr Coupling strength variant in atomic units
    \param phmass Photon mass (default 1.0)
    \param damping_ratio Damping ratio (default 0.0)
*/
CavityForceCompute::CavityForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                                       Scalar omegac,
                                       std::shared_ptr<Variant> couplstr, 
                                       Scalar phmass,
                                       Scalar damping_ratio)
    : ForceCompute(sysdef), m_params(omegac, couplstr, phmass, damping_ratio)
{
    m_exec_conf->msg->notice(5) << "Constructing CavityForceCompute" << std::endl;
    
    // Initialize energy components
    m_harmonic_energy = Scalar(0.0);
    m_coupling_energy = Scalar(0.0);
    m_dipole_self_energy = Scalar(0.0);
    
    // Compute effective gamma from damping ratio
    Scalar gamma = 2.0 * damping_ratio * sqrt(m_params.K);
    
    // Display initial coupling strength
    Scalar initial_coupling = couplstr ? (*couplstr)(0) : 0.0;
    
    std::cout << "CavityForceCompute initialized: "
              << "omegac=" << omegac << " a.u., "
              << "couplstr=" << initial_coupling << " a.u. (at t=0), "  
              << "K=" << m_params.K << " a.u., "
              << "damping_ratio=" << damping_ratio << ", "
              << "gamma=" << gamma << " a.u." << std::endl;
}

CavityForceCompute::~CavityForceCompute()
{
    m_exec_conf->msg->notice(5) << "Destroying CavityForceCompute" << std::endl;
}

void CavityForceCompute::setParams(Scalar omegac, std::shared_ptr<Variant> couplstr, Scalar phmass, Scalar damping_ratio)
{
    m_params.omegac = omegac;
    m_params.couplstr = couplstr;
    m_params.phmass = phmass;
    m_params.damping_ratio = damping_ratio;
    m_params.K = phmass * omegac * omegac;
}

py::dict CavityForceCompute::getParams()
{
    return m_params.asDict();
}

Scalar CavityForceCompute::getHarmonicEnergy()
{
    return m_harmonic_energy;
}

Scalar CavityForceCompute::getCouplingEnergy()
{
    return m_coupling_energy;
}

Scalar CavityForceCompute::getDipoleSelfEnergy()
{
    return m_dipole_self_energy;
}

int CavityForceCompute::findPhotonParticle(const Scalar4* pos_data, unsigned int N)
{
    unsigned int L_typeid = m_pdata->getTypeByName("L");
    
    for (unsigned int i = 0; i < N; i++)
    {
        int type = __scalar_as_int(pos_data[i].w);
        if (type == (int)L_typeid)
        {
            return (int)i;
        }
    }
    return -1; // Not found
}

void CavityForceCompute::computeUnwrappedPositions(std::vector<vec3<Scalar>>& unwrapped_pos,
                                                   const Scalar4* pos,
                                                   const int3* image,
                                                   const BoxDim& box,
                                                   unsigned int N)
{
    Scalar3 box_L = box.getL();
    vec3<Scalar> L(box_L.x, box_L.y, box_L.z);
    unwrapped_pos.resize(N);
    
    for (unsigned int i = 0; i < N; i++)
    {
        unwrapped_pos[i].x = pos[i].x + Scalar(image[i].x) * L.x;
        unwrapped_pos[i].y = pos[i].y + Scalar(image[i].y) * L.y;
        unwrapped_pos[i].z = pos[i].z + Scalar(image[i].z) * L.z;
    }
}

vec3<Scalar> CavityForceCompute::computeDipoleMoment(const std::vector<vec3<Scalar>>& unwrapped_pos,
                                                     const Scalar* charge,
                                                     unsigned int N,
                                                     int photon_idx)
{
    vec3<Scalar> dipole = vec3<Scalar>(0, 0, 0);
    
    for (unsigned int i = 0; i < N; i++)
    {
        if ((int)i != photon_idx) // Skip photon particle
        {
            dipole.x += charge[i] * unwrapped_pos[i].x;
            dipole.y += charge[i] * unwrapped_pos[i].y;
            dipole.z += charge[i] * unwrapped_pos[i].z;
        }
    }
    
    return dipole;
}

/*! \param timestep Current time step
    This function computes the cavity-molecule interaction forces
*/
void CavityForceCompute::computeForces(uint64_t timestep)
{
    // access the particle data arrays
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_charge(m_pdata->getCharges(), access_location::host, access_mode::read);
    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar> h_virial(m_virial, access_location::host, access_mode::overwrite);

    // Zero virial
    memset((void*)h_virial.data, 0, sizeof(Scalar) * m_virial.getNumElements());

    const BoxDim& box = m_pdata->getGlobalBox();
    unsigned int N = m_pdata->getN();
    
    // Find photon particle
    int photon_idx = findPhotonParticle(h_pos.data, N);
    if (photon_idx == -1)
    {
        m_exec_conf->msg->error() << "Cavity particle with type 'L' not found!" << std::endl;
        throw std::runtime_error("Cavity particle not found");
    }
    
    // Get current coupling strength at this timestep
    Scalar current_couplstr = m_params.couplstr ? (*m_params.couplstr)(timestep) : Scalar(0.0);
    
    // Compute unwrapped positions
    std::vector<vec3<Scalar>> unwrapped_pos;
    computeUnwrappedPositions(unwrapped_pos, h_pos.data, h_image.data, box, N);
    
    // Compute dipole moment (excluding photon)
    vec3<Scalar> dipole = computeDipoleMoment(unwrapped_pos, h_charge.data, N, photon_idx);
    
    // Get photon position (only x,y components used)
    vec3<Scalar> q_photon = unwrapped_pos[photon_idx];
    vec3<Scalar> q_photon_xy = vec3<Scalar>(q_photon.x, q_photon.y, 0);
    vec3<Scalar> dipole_xy = vec3<Scalar>(dipole.x, dipole.y, 0);
    
    // Compute energy contributions using current coupling strength
    m_harmonic_energy = Scalar(0.5) * m_params.K * dot(q_photon, q_photon);
    m_coupling_energy = current_couplstr * dot(dipole_xy, q_photon_xy);
    m_dipole_self_energy = Scalar(0.5) * (current_couplstr * current_couplstr / m_params.K) * dot(dipole_xy, dipole_xy);
    
    // DO NOT assign energy to particle potential energy - prevents double-counting
    // Energy is accessed directly through force object methods
    h_force.data[photon_idx].w = 0.0;
    
    // Compute forces on molecular particles using current coupling strength
    vec3<Scalar> Dq = q_photon_xy + (current_couplstr / m_params.K) * dipole_xy;
    
    // Get the typeid for 'L' type
    unsigned int L_typeid = m_pdata->getTypeByName("L");
    
    for (unsigned int i = 0; i < N; i++)
    {
        int type = __scalar_as_int(h_pos.data[i].w);
        if (type != (int)L_typeid) // Not cavity particle
        {
            Scalar charge = h_charge.data[i];
            vec3<Scalar> force = -current_couplstr * charge * Dq;
            
            h_force.data[i].x = force.x;
            h_force.data[i].y = force.y;
            h_force.data[i].z = Scalar(0.0); // Zero z-component
        }
    }
    
    // Force on photon particle: F = -K * q - g * dipole_xy - gamma * velocity
    // where gamma = 2 * damping_ratio * sqrt(K)
    Scalar gamma = 2.0 * m_params.damping_ratio * sqrt(m_params.K);
    vec3<Scalar> photon_velocity = vec3<Scalar>(h_vel.data[photon_idx].x, h_vel.data[photon_idx].y, h_vel.data[photon_idx].z);
    vec3<Scalar> photon_force = -m_params.K * q_photon - current_couplstr * dipole_xy - gamma * photon_velocity;
    
    h_force.data[photon_idx].x = photon_force.x;
    h_force.data[photon_idx].y = photon_force.y;
    h_force.data[photon_idx].z = photon_force.z;
    h_force.data[photon_idx].w = 0.0; // No potential energy contribution to particle
}

namespace detail
{
void export_CavityForceCompute(pybind11::module& m)
{
    pybind11::class_<CavityForceCompute, ForceCompute, std::shared_ptr<CavityForceCompute>>(
        m, "CavityForceCompute")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar, std::shared_ptr<Variant>, Scalar, Scalar>(),
             pybind11::arg("sysdef"), pybind11::arg("omegac"), pybind11::arg("couplstr"), 
             pybind11::arg("phmass") = 1.0, pybind11::arg("damping_ratio") = 0.0)
        .def("setParams", &CavityForceCompute::setParams)
        .def("getParams", &CavityForceCompute::getParams)
        .def("getHarmonicEnergy", &CavityForceCompute::getHarmonicEnergy)
        .def("getCouplingEnergy", &CavityForceCompute::getCouplingEnergy)
        .def("getDipoleSelfEnergy", &CavityForceCompute::getDipoleSelfEnergy);
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd 