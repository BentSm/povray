//******************************************************************************
///
/// @file core/shape/fractal/distestimator.cpp
///
/// This module implements distance estimators for use in fractals.
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

// configcore.h must always be the first POV file included in core *.cpp files (pulls in platform config)
#include "core/configcore.h"

#include "core/shape/fractal/distestimator.h"

#include "core/shape/fractal/quaternion.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

namespace estimators
{

DBL EstimatorNone(const FractalRules *, DBL, int, const Vector4d&, const Fractal *pFractal, FractalIterData *)
{
    return pFractal->Precision;
}

DBL EstimatorNewton(const FractalRules *pRules, DBL norm, int iters, const Vector4d& direction,
                    const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL step, fValue;

    step = pRules->CalcDirDeriv(direction, iters, pFractal, pIterData);

    step *= (step < 0 ? -1 : 1);

    if (step < kDerivativeTolerance)
    {
        step = kDerivativeTolerance;
    }

    return (norm - pFractal->Exit_Value) / step;
}

DBL EstimatorNewtonOrig(const FractalRules *pRules, DBL norm, int iters, const Vector4d& direction,
                        const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL step, fValue;

    step = pRules->CalcDirDeriv(direction, iters, pFractal, pIterData);

    if (step < -kDerivativeTolerance)
    {
        step *= -1;

        fValue = norm - pFractal->Exit_Value;

        if (fValue > pFractal->Trust_Vector[iters] * step)
        {
            return pFractal->Precision_Vector[iters];
        }
        else if (fValue > pFractal->Precision_Vector[iters] * step)
        {
            return fValue / step;
        }
    }

    return pFractal->Precision;
}

const DistanceEstimator kNone = { EstimatorNone, kNoEstimator };
const DistanceEstimator kNewton = { AddLimiting<EstimatorNewton>, kNewtonEstimator };
const DistanceEstimator kNewtonOrig = { EstimatorNewtonOrig, kOrigNewtonEstimator };

}

}
