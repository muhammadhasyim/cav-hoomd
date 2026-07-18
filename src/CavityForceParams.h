// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __CAVITY_FORCE_PARAMS_H__
#define __CAVITY_FORCE_PARAMS_H__

#include "hoomd/HOOMDMath.h"

#ifndef __HIPCC__
#include <pybind11/pybind11.h>
#endif

namespace hoomd
{
namespace cavitymd
{

//! Parameters shared by host and device cavity force code.
struct cavity_force_params
    {
    Scalar omegac;           //!< Cavity frequency in atomic units
    Scalar lambda_coupling;  //!< Dimensionless coupling parameter
    Scalar K;                //!< Spring constant (phmass * omegac^2)
    Scalar phmass;           //!< Photon mass

#ifndef __HIPCC__
    cavity_force_params() : omegac(0.), lambda_coupling(0.), K(0.), phmass(1.) { }

    cavity_force_params(Scalar _omegac, Scalar _lambda_coupling, Scalar _phmass)
        : omegac(_omegac), lambda_coupling(_lambda_coupling), phmass(_phmass)
        {
        K = phmass * omegac * omegac;
        }

    pybind11::dict asDict()
        {
        pybind11::dict values;
        values["omegac"] = omegac;
        values["lambda_coupling"] = lambda_coupling;
        values["K"] = K;
        values["phmass"] = phmass;
        return values;
        }
#endif
    } __attribute__((aligned(16)));

} // end namespace cavitymd
} // end namespace hoomd

#endif // __CAVITY_FORCE_PARAMS_H__
