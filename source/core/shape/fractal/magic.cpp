//******************************************************************************
///
/// @file core/shape/fractal/magic.cpp
///
/// This module contains the generic implementation of FractalRules subclasses.
///
/// (This was originally named for the templating 'magic' that it used, but
/// could equally well be described as the 'magic' that makes everything work!)
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

#include "core/shape/fractal/magic.h"

#include "base/pov_err.h"

#include "core/math/complexfn.h"
#include "core/math/vector.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/distestimator.h"
#include "core/shape/fractal/space.h"
#include "core/shape/fractal/util.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

const DistanceEstimator& MagicRulesBase::
GetEstimatorFromType(EstimatorType estimatorType, EstimatorType defaultEstimator, EstimatorType legacyEstimator,
                     const DistanceEstimator& (*ExtraEstimators)(EstimatorType eType))
{
    EstimatorType tgtType = estimatorType;
    if (tgtType == kDefaultEstimator)
        tgtType = defaultEstimator;
    if (tgtType == kLegacyEstimator)
        tgtType = legacyEstimator;
    switch (tgtType)
    {
    case kInvalidEstimator:
        return estimators::BadEstimator();

    case kNoEstimator:
        return estimators::kNone;

    case kNewtonEstimator:
        return estimators::kNewton;

    case kOrigNewtonEstimator:
        return estimators::kNewtonOrig;

    default:
        if (ExtraEstimators != NULL)
            return (*ExtraEstimators)(tgtType);
        else
            return estimators::BadEstimator();
    }
}

int MagicFractalRules::
Iterate(const Vector4d& iPoint, const Fractal *pFractal, const Vector4d& direction, DBL *pDist,
        FractalIterData *pIterData) const
{
    Vector4d *pIterStack = static_cast<Vector4d *>(pIterData->mainIter.data());
    Vector4d v = iPoint;
    DBL norm, exitValue;
    int i;

    pIterStack[0] = v;

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        norm = v.lengthSqr();

        if (norm > exitValue)
        {
            if (pDist != NULL)
            {
                *pDist = (*(mEstimator.pEstim))(this, norm, i, direction, pFractal, pIterData);
            }

            return i;
        }

        IterateCalc(v, norm, i, pFractal, pIterData);

        pIterStack[i+1] = v;
    }

    norm = v.lengthSqr();

    if (norm > exitValue)
    {
        if (pDist != NULL)
        {
            *pDist = (*(mEstimator.pEstim))(this, norm, pFractal->Num_Iterations, direction, pFractal, pIterData);
        }

        return pFractal->Num_Iterations;
    }

    if (pDist != NULL)
    {
        *pDist = pFractal->Precision;
    }

    return pFractal->Num_Iterations + 1;
}

DBL MagicFractalRules::
CalcDirDeriv(const Vector4d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Vector4d *pIterStack = static_cast<Vector4d *>(pIterData->mainIter.data());
    Vector4d d = pIterStack[nMax];
    DBL mult = 1.0;
    int i;

    for (i = nMax - 1; i >= 0; i--)
    {
        GradientCalc(d, pIterStack[i], i, mult, pFractal, pIterData);
    }

    return mult * dot(d, dir);
}

void MagicFractalRules::
CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    Vector4d *pTIterStack = static_cast<Vector4d *>(pTIterData->mainIter.data()),
        *pPIterStack = (pPIterData != NULL ? static_cast<Vector4d *>(pPIterData->mainIter.data()) : NULL);
    Vector4d d;
    DBL mult, dist;
    int newNMax = nMax, i;
    bool foundDisc = false;

    if (pPIterData != NULL && pFractal->Discontinuity_Test > 0)
    {
        for (i = 0; i < nMax; i++)
        {
            if (DiscontinuityCheck(d, dist, pTIterStack[i], pPIterStack[i],
                                   i, pFractal, pTIterData, pPIterData))
            {
                foundDisc = true;
                newNMax = i;
                break;
            }
        }
    }

    if (!foundDisc)
    {
        d = pTIterStack[nMax];
    }

    for (i = newNMax - 1; i >= 0; i--)
    {
        GradientCalc(d, pTIterStack[i], i, mult, pFractal, pTIterData);
    }

    rResult[X] = dot(d, mSpace4D.transformedX());
    rResult[Y] = dot(d, mSpace4D.transformedY());
    rResult[Z] = dot(d, mSpace4D.transformedZ());
}

bool MagicFractalRules::
DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}

}

