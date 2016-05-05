//******************************************************************************
///
/// @file core/shape/fractal/distestimator.h
///
/// This module implements distance estimators for use in fractals.
///
/// @copyright
/// @parblock
///
/// Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
/// Copyright 1991-2015 Persistence of Vision Raytracer Pty. Ltd.
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

class DistEstimatorNone
{
public:
    template <class Rules>
    static inline DBL Call(const Rules *, DBL, int, const Vector3d&, const Fractal *pFractal, void *)
    {
        return pFractal->Precision;
    }

    static const EstimatorType eType = kNoEstimator;
};

class DistEstimatorNewton
{
public:
    template <class Rules>
    static inline DBL Call(const Rules *pRules, DBL norm, int iters, const Vector3d& direction,
                           const Fractal *pFractal, void *pIterData)
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

        if (step > kDistanceEstimatorTolerance * 2)
        {
            fValue = norm - pFractal->Exit_Value;

            if (fValue > trustAmt * step)
            {
                return trustAmt;
            }
            else if (fValue > pFractal->Precision * step)
            {
                return fValue / step;
            }
        }

        return pFractal->Precision;
    }

    static const EstimatorType eType = kNewtonEstimator;
};

class DistEstimatorNewtonOrig
{
public:
    template <class Rules>
    static inline DBL Call(const Rules *pRules, DBL norm, int iters, const Vector3d& direction,
                           const Fractal *pFractal, void *pIterData)
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

    static const EstimatorType eType = kOrigNewtonEstimator;
};

class DistEstimatorSpecialOrig_QuatSqr
{
public:
    template <class Rules>
    static inline DBL Call(const Rules *rules, DBL norm, int iters, const Vector3d& direction,
                           const Fractal *pFractal, void *pIterData)
    {
        DBL tmp, nProd, pow;
        int j;

        typename Rules::IterationData *pIterStack =
            reinterpret_cast<typename Rules::IterationData *>(pIterData);

        tmp = dot(pFractal->SliceNorm, direction);

        nProd = 1.0 + tmp * tmp;

        pow = 1.0 / 2.0;

        for (j = 0; j < iters; ++j)
        {
            nProd *= pIterStack[j].sNorm;
            pow /= 2.0;
        }

        return pow / sqrt(nProd) * log(norm);
    }

    static const EstimatorType eType = kOrigSpecialEstimators;
};

class DistEstimatorSpecialOrig_QuatCube
{
public:
    template <class Rules>
    static inline DBL Call(const Rules *rules, DBL norm, int iters, const Vector3d& direction,
                           const Fractal *pFractal, void *pIterData)
    {
        DBL tmp, nProd, pow;
        int j;

        typename Rules::IterationData *pIterStack =
            reinterpret_cast<typename Rules::IterationData *>(pIterData);

        tmp = dot(pFractal->SliceNorm, direction);

        nProd = 1.0 + tmp * tmp;

        pow = 1.0 / 3.0;

        for (j = 0; j < iters; ++j)
        {
            nProd *= pIterStack[j].sNorm;
            pow /= 3.0;
        }

        return pow / sqrt(nProd) * log(norm);
    }

    static const EstimatorType eType = kOrigSpecialEstimators;
};

}

#endif
