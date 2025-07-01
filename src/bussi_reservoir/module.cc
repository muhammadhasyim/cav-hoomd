// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "../BussiReservoirThermostat.h"

namespace py = pybind11;

/*! \file module.cc
    \brief Export classes to python
*/

using namespace hoomd;
using namespace hoomd::md;

// Forward declaration
namespace hoomd { namespace md { void export_BussiReservoirThermostat(pybind11::module& m); } }

// specify the python module. Note that the name must exactly match the filename of the output
// library
PYBIND11_MODULE(_bussi_reservoir, m)
{
    export_BussiReservoirThermostat(m);
} 