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
    \param couplstr Coupling strength in atomic units
    \param phmass Photon mass (default 1.0)
*/
CavityForceCompute::CavityForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                                       Scalar omegac,
                                       Scalar couplstr, 
                                       Scalar phmass)
    : ForceCompute(sysdef), m_params(omegac, couplstr, phmass)
{
    m_exec_conf->msg->notice(5) << "Constructing CavityForceCompute" << std::endl;
    
    // Initialize energy components
    m_harmonic_energy = Scalar(0.0);
    m_coupling_energy = Scalar(0.0);
    m_dipole_self_energy = Scalar(0.0);
    
    std::cout << "CavityForceCompute initialized: "
              << "omegac=" << omegac << " a.u., "
              << "couplstr=" << couplstr << " a.u., "  
              << "K=" << m_params.K << " a.u." << std::endl;
}

CavityForceCompute::~CavityForceCompute()
{
    m_exec_conf->msg->notice(5) << "Destroying CavityForceCompute" << std::endl;
}

void CavityForceCompute::setParams(Scalar omegac, Scalar couplstr, Scalar phmass)
{
    m_params = cavity_force_params(omegac, couplstr, phmass);
}

pybind11::dict CavityForceCompute::getParams()
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
    // Find photon particle with type name 'L'
    std::shared_ptr<ParticleData> pdata = m_pdata;
    
    // Get the type mapping to find the typeid for 'L'
    unsigned int L_typeid = pdata->getTypeByName("L");
    
    for (unsigned int i = 0; i < N; i++)
    {
        if (__scalar_as_int(pos_data[i].w) == (int)L_typeid)
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
    Scalar3 L_scalar3 = box.getL();
    vec3<Scalar> L = vec3<Scalar>(L_scalar3.x, L_scalar3.y, L_scalar3.z);
    unwrapped_pos.resize(N);
    
    for (unsigned int i = 0; i < N; i++)
    {
        vec3<Scalar> p = vec3<Scalar>(pos[i]);
        int3 img = image[i];
        
        // Unwrap position
        unwrapped_pos[i].x = p.x + Scalar(img.x) * L.x;
        unwrapped_pos[i].y = p.y + Scalar(img.y) * L.y;
        unwrapped_pos[i].z = p.z + Scalar(img.z) * L.z;
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
            dipole += charge[i] * unwrapped_pos[i];
        }
    }
    
    return dipole;
}

/*! This function computes cavity forces and energy
    \param timestep Current timestep
*/
void CavityForceCompute::computeForces(uint64_t timestep)
{
    // Access particle data
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar> h_charge(m_pdata->getCharges(), access_location::host, access_mode::read);
    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::read);
    
    unsigned int N = m_pdata->getN();
    
    // Zero force arrays
    memset(h_force.data, 0, sizeof(Scalar4) * N);
    
    // Find photon particle
    int photon_idx = findPhotonParticle(h_pos.data, N);
    if (photon_idx == -1)
    {
        // No photon particle found - zero all energies
        m_harmonic_energy = Scalar(0.0);
        m_coupling_energy = Scalar(0.0);
        m_dipole_self_energy = Scalar(0.0);
        return;
    }
    
    // Get box dimensions
    BoxDim box = m_pdata->getGlobalBox();
    
    // Compute unwrapped positions
    std::vector<vec3<Scalar>> unwrapped_pos;
    computeUnwrappedPositions(unwrapped_pos, h_pos.data, h_image.data, box, N);
    
    // Compute molecular dipole moment  
    vec3<Scalar> dipole = computeDipoleMoment(unwrapped_pos, h_charge.data, N, photon_idx);
    
    // Get photon position (only x,y components used)
    vec3<Scalar> q_photon = unwrapped_pos[photon_idx];
    vec3<Scalar> q_photon_xy = vec3<Scalar>(q_photon.x, q_photon.y, 0);
    vec3<Scalar> dipole_xy = vec3<Scalar>(dipole.x, dipole.y, 0);
    
    // Compute energy contributions
    m_harmonic_energy = Scalar(0.5) * m_params.K * dot(q_photon, q_photon);
    m_coupling_energy = m_params.couplstr * dot(dipole_xy, q_photon_xy);
    m_dipole_self_energy = Scalar(0.5) * (m_params.couplstr * m_params.couplstr / m_params.K) * dot(dipole_xy, dipole_xy);
    
    Scalar total_cavity_energy = m_harmonic_energy + m_coupling_energy + m_dipole_self_energy;
    
    // DO NOT assign energy to particle potential energy - prevents double-counting
    // Energy is accessed directly through force object methods
    h_force.data[photon_idx].w = 0.0;
    
    // Compute forces on molecular particles
    vec3<Scalar> Dq = q_photon_xy + (m_params.couplstr / m_params.K) * dipole_xy;
    
    // Get the typeid for 'L' type
    unsigned int L_typeid = m_pdata->getTypeByName("L");
    
    for (unsigned int i = 0; i < N; i++)
    {
        int type = __scalar_as_int(h_pos.data[i].w);
        if (type != (int)L_typeid) // Not cavity particle
        {
            Scalar charge = h_charge.data[i];
            vec3<Scalar> force = -m_params.couplstr * charge * Dq;
            
            h_force.data[i].x = force.x;
            h_force.data[i].y = force.y;
            h_force.data[i].z = Scalar(0.0); // Zero z-component
        }
    }
    
    // Force on photon particle
    vec3<Scalar> photon_force = -m_params.K * q_photon - m_params.couplstr * dipole_xy;
    
    h_force.data[photon_idx].x = photon_force.x;
    h_force.data[photon_idx].y = photon_force.y;
    h_force.data[photon_idx].z = photon_force.z;
}

namespace detail
{
void export_CavityForceCompute(pybind11::module& m)
{
    pybind11::class_<CavityForceCompute, ForceCompute, std::shared_ptr<CavityForceCompute>>(
        m, "CavityForceCompute")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar, Scalar, Scalar>(),
             pybind11::arg("sysdef"), pybind11::arg("omegac"), pybind11::arg("couplstr"), 
             pybind11::arg("phmass") = 1.0)
        .def("setParams", &CavityForceCompute::setParams)
        .def("getParams", &CavityForceCompute::getParams)
        .def("getHarmonicEnergy", &CavityForceCompute::getHarmonicEnergy)
        .def("getCouplingEnergy", &CavityForceCompute::getCouplingEnergy)
        .def("getDipoleSelfEnergy", &CavityForceCompute::getDipoleSelfEnergy);
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd 