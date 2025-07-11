// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "CavityParticleDisplacer.h"
#include "hoomd/VectorMath.h"
#include <vector>

namespace hoomd
{
namespace cavitymd
{

CavityParticleDisplacer::CavityParticleDisplacer(std::shared_ptr<SystemDefinition> sysdef,
                                             std::shared_ptr<Trigger> trigger,
                                             std::shared_ptr<Variant> couplstr,
                                             Scalar omegac,
                                             Scalar phmass)
    : Updater(sysdef, trigger),
      m_couplstr(couplstr),
      m_omegac(omegac),
      m_phmass(phmass),
      m_has_run(false)
    {
    m_exec_conf->msg->notice(5) << "Constructing CavityParticleDisplacer" << std::endl;
    m_K = m_phmass * m_omegac * m_omegac;
    }

CavityParticleDisplacer::~CavityParticleDisplacer()
    {
    m_exec_conf->msg->notice(5) << "Destroying CavityParticleDisplacer" << std::endl;
    }

int CavityParticleDisplacer::findPhotonParticle(const Scalar4* pos_data, unsigned int N)
{
    // Find photon particle with type name 'L'
    std::shared_ptr<ParticleData> pdata = m_pdata;
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

void CavityParticleDisplacer::computeUnwrappedPositions(std::vector<vec3<Scalar>>& unwrapped_pos,
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
        unwrapped_pos[i].x = p.x + Scalar(img.x) * L.x;
        unwrapped_pos[i].y = p.y + Scalar(img.y) * L.y;
        unwrapped_pos[i].z = p.z + Scalar(img.z) * L.z;
        }
}

vec3<Scalar> CavityParticleDisplacer::computeDipoleMoment(const std::vector<vec3<Scalar>>& unwrapped_pos,
                                                     const Scalar* charge,
                                                     unsigned int N,
                                                     int photon_idx)
{
    vec3<Scalar> dipole = vec3<Scalar>(0, 0, 0);

    for (unsigned int i = 0; i < N; i++)
        {
        if ((int)i != photon_idx)
            {
            dipole += charge[i] * unwrapped_pos[i];
            }
        }

    return dipole;
}

void CavityParticleDisplacer::update(uint64_t timestep)
    {
    m_exec_conf->msg->notice(1) << "CavityParticleDisplacer::update running at timestep " << timestep << std::endl;

    // Get current coupling strength
    Scalar g_current = (*m_couplstr)(timestep);
    
    // Check if we've already processed the switch or if coupling is still zero
    if (m_has_run || g_current == 0.0) {
        // Either already displaced, or coupling is still zero - nothing to do
        return;
    }
    
    // If we reach here, g_current > 0.0 and we haven't run yet
    // This means the coupling has just switched from 0 to non-zero
    // Time to perform the displacement!
    
    m_exec_conf->msg->notice(1) << "Coupling has switched ON! Performing finite-q displacement..." << std::endl;
    m_exec_conf->msg->notice(1) << "  Current coupling g: " << g_current << std::endl;

    // Get particle data
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_charge(m_pdata->getCharges(), access_location::host, access_mode::read);
    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::readwrite);

    unsigned int N = m_pdata->getN();

    // Find photon particle
    int photon_idx = findPhotonParticle(h_pos.data, N);
    if (photon_idx == -1)
        {
        m_exec_conf->msg->warning() << "cavitymd.update: No cavity particle found, skipping displacement." << std::endl;
        m_has_run = true;
        return;
        }

    // Get box dimensions
    BoxDim box = m_pdata->getGlobalBox();
    Scalar3 L = box.getL();

    // Compute unwrapped positions
    std::vector<vec3<Scalar>> unwrapped_pos;
    computeUnwrappedPositions(unwrapped_pos, h_pos.data, h_image.data, box, N);

    // Compute molecular dipole moment
    vec3<Scalar> dipole = computeDipoleMoment(unwrapped_pos, h_charge.data, N, photon_idx);
    vec3<Scalar> dipole_xy(dipole.x, dipole.y, 0);

    // Calculate new equilibrium position using current coupling
    vec3<Scalar> q_eq = - (g_current / m_K) * dipole_xy;

    // Get current cavity position for comparison
    vec3<Scalar> q_photon_old = unwrapped_pos[photon_idx];

    m_exec_conf->msg->notice(1) << "CavityParticleDisplacer displacement details:" << std::endl;
    m_exec_conf->msg->notice(1) << "  Timestep: " << timestep << std::endl;
    m_exec_conf->msg->notice(1) << "  Coupling g: " << g_current << std::endl;
    m_exec_conf->msg->notice(1) << "  Dipole moment (xy): " << dipole_xy.x << ", " << dipole_xy.y << std::endl;
    m_exec_conf->msg->notice(1) << "  Old cavity position: " << q_photon_old.x << ", " << q_photon_old.y << ", " << q_photon_old.z << std::endl;
    m_exec_conf->msg->notice(1) << "  New equilibrium position (xy only): " << q_eq.x << ", " << q_eq.y << ", " << q_eq.z << std::endl;

    // Create new position: use equilibrium x,y coordinates but preserve original z coordinate
    vec3<Scalar> new_pos_unwrapped;
    new_pos_unwrapped.x = q_eq.x;  // New equilibrium x
    new_pos_unwrapped.y = q_eq.y;  // New equilibrium y
    new_pos_unwrapped.z = q_photon_old.z;  // Preserve original z
    
    m_exec_conf->msg->notice(1) << "  Final new position (z preserved): " << new_pos_unwrapped.x << ", " << new_pos_unwrapped.y << ", " << new_pos_unwrapped.z << std::endl;

    // Wrap the new position back into the box and find the new image flags
    int3 new_image;
    new_image.x = rint(new_pos_unwrapped.x / L.x);
    new_image.y = rint(new_pos_unwrapped.y / L.y);
    new_image.z = rint(new_pos_unwrapped.z / L.z);

    vec3<Scalar> new_pos_wrapped;
    new_pos_wrapped.x = new_pos_unwrapped.x - Scalar(new_image.x) * L.x;
    new_pos_wrapped.y = new_pos_unwrapped.y - Scalar(new_image.y) * L.y;
    new_pos_wrapped.z = new_pos_unwrapped.z - Scalar(new_image.z) * L.z;
    
    // Update the particle position
    h_pos.data[photon_idx].x = new_pos_wrapped.x;
    h_pos.data[photon_idx].y = new_pos_wrapped.y;
    h_pos.data[photon_idx].z = new_pos_wrapped.z;
    h_image.data[photon_idx] = new_image;

    // For finite-q displacement, we only modify position, never velocities
    // The cavity particle velocity should remain unchanged
    
    m_exec_conf->msg->notice(1) << "Cavity particle displaced (position only, velocity unchanged)" << std::endl;
    m_exec_conf->msg->notice(1) << "  Current velocity: " << h_vel.data[photon_idx].x << ", " << h_vel.data[photon_idx].y << ", " << h_vel.data[photon_idx].z << std::endl;
    
    // Calculate displacement magnitude in x,y plane only (z is preserved)
    Scalar displacement_xy = sqrt((new_pos_unwrapped.x - q_photon_old.x)*(new_pos_unwrapped.x - q_photon_old.x) + 
                                  (new_pos_unwrapped.y - q_photon_old.y)*(new_pos_unwrapped.y - q_photon_old.y));
    m_exec_conf->msg->notice(1) << "  Displacement magnitude (xy plane): " << displacement_xy << std::endl;
    m_exec_conf->msg->notice(1) << "  Z coordinate preserved (no displacement)" << std::endl;
    
    // CRITICAL: Close array handles BEFORE calling synchronization functions
    // This ensures all data is properly committed to the particle data structures
    h_pos.~ArrayHandle();
    h_vel.~ArrayHandle();
    h_charge.~ArrayHandle();
    h_image.~ArrayHandle();
    
    // CRITICAL: Notify HOOMD-blue that particle data has been modified
    // This invalidates all cached force calculations and neighbor lists
    m_pdata->notifyParticleSort();
    
    // CRITICAL: Force immediate synchronization of all data structures
    // This ensures the new position is immediately visible to all force calculations
    if (m_exec_conf->isCUDAEnabled()) {
        // GPU case: ensure device synchronization
        #ifdef ENABLE_HIP
        hipDeviceSynchronize();
        #endif
    }
    
    m_exec_conf->msg->notice(1) << "System state synchronized after position update" << std::endl;
    m_exec_conf->msg->notice(1) << "All cached data invalidated to ensure energy consistency" << std::endl;
    
    // Mark as completed so we don't run again
    m_has_run = true;
    }

namespace detail
{
void export_CavityParticleDisplacer(pybind11::module& m)
    {
    pybind11::class_<CavityParticleDisplacer, Updater, std::shared_ptr<CavityParticleDisplacer>>(
        m, "CavityParticleDisplacer")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<Trigger>, std::shared_ptr<Variant>, Scalar, Scalar>(),
             pybind11::arg("sysdef"),
             pybind11::arg("trigger"),
             pybind11::arg("couplstr"),
             pybind11::arg("omegac"),
             pybind11::arg("phmass") = 1.0);
    }
} // end namespace detail

} // end namespace cavitymd
} // end namespace hoomd 