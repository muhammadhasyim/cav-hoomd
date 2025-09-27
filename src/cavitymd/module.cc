// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "CavityForceCompute.h"
#include "CavityParticleDisplacer.h"
#include "DipoleResponseForceCompute.h"
// #include "PerturbationForceCompute.h"

#ifdef ENABLE_HIP
#include "CavityForceComputeGPU.h"
#include "DipoleResponseForceComputeGPU.h"
// #include "PerturbationForceComputeGPU.h"
#endif

namespace py = pybind11;

/*! \file module.cc
    \brief Export classes to python
*/

using namespace hoomd;
using namespace hoomd::cavitymd;

// Forward declarations
namespace hoomd { namespace cavitymd { namespace detail { void export_CavityForceCompute(pybind11::module& m); } } }
namespace hoomd { namespace cavitymd { namespace detail { void export_DipoleResponseForceCompute(pybind11::module& m); } } }
// namespace hoomd { namespace cavitymd { namespace detail { void export_PerturbationForceCompute(pybind11::module& m); } } }
#ifdef ENABLE_HIP
namespace hoomd { namespace cavitymd { namespace detail { void export_CavityForceComputeGPU(pybind11::module& m); } } }
namespace hoomd { namespace cavitymd { namespace detail { void export_DipoleResponseForceComputeGPU(pybind11::module& m); } } }
// namespace hoomd { namespace cavitymd { namespace detail { void export_PerturbationForceComputeGPU(pybind11::module& m); } } }
#endif
namespace hoomd { namespace cavitymd { namespace detail { void export_CavityParticleDisplacer(pybind11::module& m); } } }

// specify the python module. Note that the name must exactly match the filename of the output
// library
PYBIND11_MODULE(_cavitymd, m)
{
    hoomd::cavitymd::detail::export_CavityForceCompute(m);
    hoomd::cavitymd::detail::export_DipoleResponseForceCompute(m);
    // Temporarily disable PerturbationForceCompute to fix symbol issues
    // hoomd::cavitymd::detail::export_PerturbationForceCompute(m);

#ifdef ENABLE_HIP
    hoomd::cavitymd::detail::export_CavityForceComputeGPU(m);
    hoomd::cavitymd::detail::export_DipoleResponseForceComputeGPU(m);
    // hoomd::cavitymd::detail::export_PerturbationForceComputeGPU(m);
#endif

    hoomd::cavitymd::detail::export_CavityParticleDisplacer(m);
} 