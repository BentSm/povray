//******************************************************************************
///
/// @file core/shape/fractal/distestimator.h
///
/// This module contains prototypes and definitions for `distestimator.cpp'.
///
/// @copyright
/// @parblock
///
/// Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
/// Copyright 1991-2016 Persistence of Vision Raytracer Pty. Ltd.
///
/// POV-Ray is free software: you can redistribute it and/or modify
/// it under the terms of the GNU Affero General Public License as
/// published by the Free Software Foundation, either version 3 of the
/// License, or (at your option) any later version.
///
/// POV-Ray is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU Affero General Public License for more details.
///
/// You should have received a copy of the GNU Affero General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.
///
/// ----------------------------------------------------------------------------
///
/// POV-Ray is based on the popular DKB raytracer version 2.12.
/// DKBTrace was originally written by David K. Buck.
/// DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
///
/// @endparblock
///
//******************************************************************************

#ifndef POVRAY_CORE_FRACTAL_DISTESTIMATOR_H
#define POVRAY_CORE_FRACTAL_DISTESTIMATOR_H

#include "core/coretypes.h"

#include "base/pov_err.h"

#include "core/math/vector.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/types.h"

namespace pov
{

/* This used to be Fractal_Tolerance from hcmplx.cpp.  I am using the value from
   hcmplx.cpp (1e-8), which was different from the value of Fractal_Tolerance in
   fractal.cpp (1e-7).  I have assumed that the two are truly meant to be separate
   constants. */
const DBL kDistanceEstimatorTolerance = 1e-8;

typedef DBL EstimatorFunc(const FractalRules *pRules, DBL norm, int iters, const Vector4d& direction,
                          const Fractal *pFractal, FractalIterData *pIterData);

struct DistanceEstimator
{
    EstimatorFunc *pEstim;
    EstimatorType eType;
};

namespace estimators
{

extern const DistanceEstimator kNone;
extern const DistanceEstimator kNewton;
extern const DistanceEstimator kNewtonOrig;

static inline const DistanceEstimator& BadEstimator()
{
    throw POV_EXCEPTION_STRING("Unsupported distance estimator for fractal type.");
}

}

}

#endif
