// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*! \file CavityForceComputeGPU.cc
    \brief Implements CavityForceComputeGPU class
    
    DEBUG LOGGING:
    - GPU debug messages use HOOMD notice level 9 (set msg.setNoticeLevel(9) to see them)
    - GPU kernel verbose messages require CAVITY_DEBUG_VERBOSE compile flag
*/

#include "CavityForceComputeGPU.h"

#include <algorithm>

#ifdef ENABLE_HIP
#include "CavityForceComputeGPU.cuh"
#include <hoomd/GPUArray.h>
#endif

namespace py = pybind11;

namespace hoomd
{
namespace cavitymd
{

/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param omegac Cavity frequency in atomic units
    \param lambda_coupling Dimensionless coupling parameter (lambda)
    \param phmass Photon mass (default 1.0)
*/
CavityForceComputeGPU::CavityForceComputeGPU(std::shared_ptr<SystemDefinition> sysdef,
                                             Scalar omegac,
                                             std::shared_ptr<Variant> lambda_coupling,
                                             Scalar phmass)
    : CavityForceCompute(sysdef, omegac, lambda_coupling, phmass),
      m_photon_tag(NOT_LOCAL),
      m_has_photon(false),
      m_energy_cache_valid(false),
      m_reduction_capacity(1)
{
    m_exec_conf->msg->notice(5) << "Constructing CavityForceComputeGPU" << std::endl;

    if (!m_exec_conf->isCUDAEnabled())
        {
        throw std::runtime_error("CavityForceComputeGPU requires a GPU execution configuration");
        }
    if (m_exec_conf->getNRanks() > 1)
        {
        throw std::runtime_error(
            "CavityForceComputeGPU supports only single-rank execution");
        }

    GPUArray<Scalar> device_energy(3, m_exec_conf);
    m_device_energy.swap(device_energy);
    GPUArray<unsigned int> device_photon_idx(1, m_exec_conf);
    m_device_photon_idx.swap(device_photon_idx);
    GPUArray<Scalar3> partial_dipole_a(1, m_exec_conf);
    m_partial_dipole_a.swap(partial_dipole_a);
    GPUArray<Scalar3> partial_dipole_b(1, m_exec_conf);
    m_partial_dipole_b.swap(partial_dipole_b);

    initializePhoton();
    updatePhotonIndex();
    m_pdata->getParticleSortSignal()
        .template connect<CavityForceComputeGPU,
                          &CavityForceComputeGPU::updatePhotonIndex>(this);
    m_pdata->getGlobalParticleNumberChangeSignal()
        .template connect<CavityForceComputeGPU,
                          &CavityForceComputeGPU::handleParticleNumberChange>(this);
}

CavityForceComputeGPU::~CavityForceComputeGPU()
{
    m_pdata->getParticleSortSignal()
        .template disconnect<CavityForceComputeGPU,
                             &CavityForceComputeGPU::updatePhotonIndex>(this);
    m_pdata->getGlobalParticleNumberChangeSignal()
        .template disconnect<CavityForceComputeGPU,
                             &CavityForceComputeGPU::handleParticleNumberChange>(this);
    m_exec_conf->msg->notice(5) << "Destroying CavityForceComputeGPU" << std::endl;
}

void CavityForceComputeGPU::initializePhoton()
{
    m_photon_tag = NOT_LOCAL;
    m_has_photon = false;
    unsigned int photon_count = 0;
    unsigned int L_typeid;
    try
        {
        L_typeid = m_pdata->getTypeByName("L");
        }
    catch (const std::runtime_error&)
        {
        throw std::runtime_error(
            "CavityForceComputeGPU requires exactly one particle of type L; found 0");
        }

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(),
                               access_location::host,
                               access_mode::read);
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(),
                                    access_location::host,
                                    access_mode::read);
    const unsigned int N = m_pdata->getN();
    for (unsigned int idx = 0; idx < N; ++idx)
        {
        if (__scalar_as_int(h_pos.data[idx].w) == static_cast<int>(L_typeid))
            {
            m_photon_tag = h_tag.data[idx];
            ++photon_count;
            }
        }
    if (photon_count != 1)
        {
        throw std::runtime_error(
            "CavityForceComputeGPU requires exactly one particle of type L; found "
            + std::to_string(photon_count));
        }
    m_has_photon = true;
}

void CavityForceComputeGPU::updatePhotonIndex()
{
    ArrayHandle<unsigned int> d_rtag(m_pdata->getRTags(),
                                     access_location::device,
                                     access_mode::read);
    ArrayHandle<unsigned int> d_photon_idx(m_device_photon_idx,
                                           access_location::device,
                                           access_mode::overwrite);
    const hipError_t error = kernel::gpu_update_photon_index(d_photon_idx.data,
                                                             d_rtag.data,
                                                             m_photon_tag,
                                                             m_has_photon);
    if (error != hipSuccess)
        {
        throw std::runtime_error(std::string("Error refreshing GPU photon index: ")
                                 + hipGetErrorString(error));
        }
}

void CavityForceComputeGPU::handleParticleNumberChange()
{
    initializePhoton();
    updatePhotonIndex();
}

void CavityForceComputeGPU::updateEnergyCache()
{
    if (m_energy_cache_valid)
        {
        return;
        }

    ArrayHandle<Scalar> h_energy(m_device_energy,
                                 access_location::host,
                                 access_mode::read);
    m_harmonic_energy = h_energy.data[0];
    m_coupling_energy = h_energy.data[1];
    m_dipole_self_energy = h_energy.data[2];
    m_energy_cache_valid = true;
}

Scalar CavityForceComputeGPU::getHarmonicEnergy()
{
    updateEnergyCache();
    return m_harmonic_energy;
}

Scalar CavityForceComputeGPU::getCouplingEnergy()
{
    updateEnergyCache();
    return m_coupling_energy;
}

Scalar CavityForceComputeGPU::getDipoleSelfEnergy()
{
    updateEnergyCache();
    return m_dipole_self_energy;
}

/*! This function computes cavity forces and energy on the GPU
    \param timestep Current timestep
*/
void CavityForceComputeGPU::computeForces(uint64_t timestep)
{
#ifdef ENABLE_HIP
    const Scalar lambda_coupling = (*m_lambda_coupling)(timestep);
    validateParams(m_params.omegac, lambda_coupling, m_params.phmass);
    m_params.lambda_coupling = lambda_coupling;
    m_energy_cache_valid = false;

    const unsigned int N = m_pdata->getN();
    const BoxDim box = m_pdata->getGlobalBox();
    constexpr unsigned int block_size = 256;
    const size_t required_capacity
        = std::max<size_t>(1, (N + block_size - 1) / block_size);
    if (required_capacity > m_reduction_capacity)
        {
        m_partial_dipole_a.resize(required_capacity);
        m_partial_dipole_b.resize(required_capacity);
        m_reduction_capacity = required_capacity;
        }

    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(),
                               access_location::device,
                               access_mode::read);
    ArrayHandle<Scalar> d_charge(m_pdata->getCharges(),
                                 access_location::device,
                                 access_mode::read);
    ArrayHandle<int3> d_image(m_pdata->getImages(),
                              access_location::device,
                              access_mode::read);
    ArrayHandle<Scalar4> d_force(m_force,
                                 access_location::device,
                                 access_mode::overwrite);
    ArrayHandle<Scalar> d_energy(m_device_energy,
                                 access_location::device,
                                 access_mode::overwrite);
    ArrayHandle<unsigned int> d_photon_idx(m_device_photon_idx,
                                            access_location::device,
                                            access_mode::read);
    ArrayHandle<Scalar3> d_partial_a(m_partial_dipole_a,
                                     access_location::device,
                                     access_mode::overwrite);
    ArrayHandle<Scalar3> d_partial_b(m_partial_dipole_b,
                                     access_location::device,
                                     access_mode::overwrite);

    const hipError_t error
        = kernel::gpu_compute_cavity_forces(d_force.data,
                                            d_pos.data,
                                            d_charge.data,
                                            d_image.data,
                                            box,
                                            m_params,
                                            d_energy.data,
                                            d_photon_idx.data,
                                            d_partial_a.data,
                                            d_partial_b.data,
                                            N,
                                            static_cast<unsigned int>(
                                                m_reduction_capacity));
    if (error != hipSuccess)
        {
        throw std::runtime_error(std::string("Error launching cavity force GPU kernel: ")
                                 + hipGetErrorString(error));
        }
#else
    m_exec_conf->msg->error() << "GPU computation requested but HIP is not enabled!" << std::endl;
    throw std::runtime_error("GPU computation required but HIP not available");
#endif
}

namespace detail
{
void export_CavityForceComputeGPU(pybind11::module& m)
{
    pybind11::class_<CavityForceComputeGPU, CavityForceCompute, std::shared_ptr<CavityForceComputeGPU>>(
        m, "CavityForceComputeGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar, std::shared_ptr<Variant>, Scalar>(),
             pybind11::arg("sysdef"), pybind11::arg("omegac"), pybind11::arg("lambda_coupling"), 
             pybind11::arg("phmass") = 1.0);
}
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd 