// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include "hoomd/Updater.h"
#include "hoomd/Variant.h"
#include <memory>

/*! \file CavityParticleDisplacer.h
    \brief Declares an updater that displaces the cavity particle.
*/

namespace hoomd
{
namespace cavitymd
{

class PYBIND11_EXPORT CavityParticleDisplacer : public Updater
    {
    public:
        //! Constructs the updater
        CavityParticleDisplacer(std::shared_ptr<SystemDefinition> sysdef,
                                std::shared_ptr<Trigger> trigger,
                                std::shared_ptr<Variant> couplstr,
                                Scalar omegac,
                                Scalar phmass);

        //! Destructor
        virtual ~CavityParticleDisplacer();

        //! Performs the displacement
        virtual void update(uint64_t timestep);

    protected:
        std::shared_ptr<Variant> m_couplstr;  // The coupling strength variant
        Scalar m_omegac;                      // Cavity frequency
        Scalar m_phmass;                      // Photon mass
        Scalar m_K;                           // Harmonic force constant K = m * omega^2
        bool m_has_run;                       // Flag to ensure the displacement happens only once

    private:
        // Helper functions to find the cavity particle and perform calculations
        int findPhotonParticle(const Scalar4* pos_data, unsigned int N);
        void computeUnwrappedPositions(std::vector<vec3<Scalar>>& unwrapped_pos,
                                       const Scalar4* pos,
                                       const int3* image,
                                       const BoxDim& box,
                                       unsigned int N);
        vec3<Scalar> computeDipoleMoment(const std::vector<vec3<Scalar>>& unwrapped_pos,
                                         const Scalar* charge,
                                         unsigned int N,
                                         int photon_idx);
    };

} // end namespace cavitymd
} // end namespace hoomd 