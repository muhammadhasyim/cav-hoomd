// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "BussiReservoirThermostat.h"
#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace md
    {

void export_BussiReservoirThermostat(pybind11::module& m);

// Set the name of the python module to bussi_reservoir
PYBIND11_MODULE(_bussi_reservoir, m)
    {
        export_BussiReservoirThermostat(m);

#ifdef ENABLE_HIP
        // TODO: Add GPU support if needed
#endif
    }

    } // end namespace md
    } // end namespace hoomd
