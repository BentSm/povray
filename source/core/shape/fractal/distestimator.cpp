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

DBL EstimatorNone(const FractalRules *, DBL, int, const Vector3d&, const Fractal *pFractal, FractalIterData *)
{
    return pFractal->Precision;
}

DBL EstimatorNewton(const FractalRules *pRules, DBL norm, int iters, const Vector3d& direction,
                    const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL step, fValue, trustAmt;
    Vector3d normal;

    trustAmt = pow(pFractal->Jump_Decay, iters) * pFractal->Jump_Max * pFractal->Precision;
    if (trustAmt < pFractal->Precision * pFractal->Jump_Min)
    {
        return pFractal->Precision;
    }

    step = pRules->CalcDirDeriv(direction, iters, pFractal, pIterData);

    step *= (step < 0 ? -2 : 2);

    fValue = norm - pFractal->Exit_Value;

    if (fValue > trustAmt * step)
    {
        return trustAmt;
    }
    else if (fValue > pFractal->Precision * step)
    {
        return fValue / step;
    }

    return pFractal->Precision;
}

DBL EstimatorNewtonOrig(const FractalRules *pRules, DBL norm, int iters, const Vector3d& direction,
                        const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL step, fValue;
    Vector3d normal;

    step = pRules->CalcDirDeriv(direction, iters, pFractal, pIterData);

    if (step < -kDistanceEstimatorTolerance)
    {
        step *= -2;

        fValue = norm - pFractal->Exit_Value;

        if (fValue > pFractal->Precision * pFractal->Jump_Max * step)
        {
            return pFractal->Precision;
        }
        else if (fValue > pFractal->Precision * step)
        {
            return fValue / step;
        }
    }

    return pFractal->Precision;
}

DBL EstimatorSpecialOrig_QuatSqr(const FractalRules *rules, DBL norm, int iters, const Vector3d& direction,
                                 const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL tmp, nProd, pow;
    int j;

    QuaternionSqrFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionSqrFractalRules::AuxIterData *>(pIterData->auxIter.data());

    tmp = dot(pFractal->SliceNorm, direction);

    nProd = 1.0 + tmp * tmp;

    pow = 1.0 / 2.0;

    for (j = 0; j < iters; ++j)
    {
        nProd *= pAuxIterStack[j].sNorm;
        pow /= 2.0;
    }

    return pow / sqrt(nProd) * log(norm);
}

DBL EstimatorSpecialOrig_QuatCube(const FractalRules *rules, DBL norm, int iters, const Vector3d& direction,
                                  const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL tmp, nProd, pow;
    int j;

    QuaternionCubeFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionCubeFractalRules::AuxIterData *>(pIterData->auxIter.data());

    tmp = dot(pFractal->SliceNorm, direction);

    nProd = 1.0 + tmp * tmp;

    pow = 1.0 / 3.0;

    for (j = 0; j < iters; ++j)
    {
        nProd *= pAuxIterStack[j].sNorm;
        pow /= 3.0;
    }

    return pow / sqrt(nProd) * log(norm);
}

}

}
