// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "BussiReservoirThermostat.h"

namespace hoomd::md
{

void export_BussiReservoirThermostat(pybind11::module& m)
{
    // Export with proper inheritance from Thermostat base class
    pybind11::class_<BussiReservoirThermostat, Thermostat, std::shared_ptr<BussiReservoirThermostat>>(
        m, "BussiReservoirThermostat")
        .def(pybind11::init<std::shared_ptr<Variant>,
                            std::shared_ptr<ParticleGroup>,
                            std::shared_ptr<ComputeThermo>,
                            std::shared_ptr<SystemDefinition>,
                            Scalar>())
        .def_property("tau", &BussiReservoirThermostat::getTau, &BussiReservoirThermostat::setTau)
        // Base class methods that need to be exposed
        .def_property("kT", &BussiReservoirThermostat::getT, &BussiReservoirThermostat::setT)
        // Derived class specific methods
        .def("getReservoirEnergyTranslational", &BussiReservoirThermostat::getReservoirEnergyTranslational)
        .def("getReservoirEnergyRotational", &BussiReservoirThermostat::getReservoirEnergyRotational)
        .def("getTotalReservoirEnergy", &BussiReservoirThermostat::getTotalReservoirEnergy)
        .def("getInstantaneousReservoirTranslational", &BussiReservoirThermostat::getInstantaneousReservoirTranslational)
        .def("getInstantaneousReservoirRotational", &BussiReservoirThermostat::getInstantaneousReservoirRotational)
        .def("getInstantaneousReservoirTotal", &BussiReservoirThermostat::getInstantaneousReservoirTotal)
        .def("resetReservoirEnergy", &BussiReservoirThermostat::resetReservoirEnergy);
}

} // namespace hoomd::md 