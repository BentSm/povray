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

int MagicQuaternionFractalRules::
Iterate(const Vector4d& iPoint, const Fractal *pFractal, const Vector4d& direction, DBL *pDist,
        FractalIterData *pIterData) const
{
    int i;
    Vector4d v = iPoint;
    DBL norm, exitValue;

    Vector4d *pIterStack = static_cast<Vector4d *>(pIterData->mainIter.data());

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

DBL MagicQuaternionFractalRules::
CalcDirDeriv(const Vector4d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Vector4d d = dir;
    DBL mult = 1.0;
    int i;

    Vector4d *pIterStack = static_cast<Vector4d *>(pIterData->mainIter.data());

    for (i = 0; i < nMax; i++)
    {
        DirDerivCalc(d, pIterStack[i], i, mult, false, pFractal, pIterData);
    }

    return mult * dot(d, pIterStack[nMax]);
}

void MagicQuaternionFractalRules::
CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    Vector4d nX = mSpace4D.transformedX(), nY = mSpace4D.transformedY(), nZ = mSpace4D.transformedZ();
    int i;
    DBL mult = 1.0;

    Vector4d *pTIterStack = static_cast<Vector4d *>(pTIterData->mainIter.data()),
        *pPIterStack = (pPIterData != NULL ? static_cast<Vector4d *>(pPIterData->mainIter.data()) : NULL);

    for (i = 0; i < nMax; i++)
    {
        if (pPIterData != NULL && pFractal->Discontinuity_Test > 0)
        {
            Vector4d d;
            DBL dist;
            if (DiscontinuityCheck(d, dist, pTIterStack[i], pPIterStack[i],
                                   i, pFractal, pTIterData, pPIterData))
            {
                rResult[X] = mult * dot(nX, d);
                rResult[Y] = mult * dot(nY, d);
                rResult[Z] = mult * dot(nZ, d);
                return;
            }
        }

        DirDerivCalc(nX, pTIterStack[i], i, mult, false, pFractal, pTIterData);
        DirDerivCalc(nY, pTIterStack[i], i, mult, true, pFractal, pTIterData);
        DirDerivCalc(nZ, pTIterStack[i], i, mult, true, pFractal, pTIterData);
    }

    rResult[X] = mult * dot(nX, pTIterStack[nMax]);
    rResult[Y] = mult * dot(nY, pTIterStack[nMax]);
    rResult[Z] = mult * dot(nZ, pTIterStack[nMax]);
}

bool MagicQuaternionFractalRules::
DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}


int MagicHypercomplexFractalRules::
Iterate(const Vector4d& iPoint, const Fractal *pFractal, const Vector4d& direction, DBL *pDist,
        FractalIterData *pIterData) const
{
    int i;
    Vector4d d = iPoint;
    DBL norm, exitValue;

    Vector4d *pIterStack = static_cast<Vector4d *>(pIterData->mainIter.data());

    pIterStack[0] = d;

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        norm = 0.5 * d.lengthSqr();

        if (norm > exitValue)
        {
            if (pDist != NULL)
            {
                *pDist = (*(mEstimator.pEstim))(this, norm, i, direction, pFractal, pIterData);
            }

            return i;
        }

        IterateCalc(d, norm, i, pFractal, pIterData);

        pIterStack[i+1] = d;

    }

    norm = 0.5 * d.lengthSqr();

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

DBL MagicHypercomplexFractalRules::
CalcDirDeriv(const Vector4d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Vector4d d = dir;
    int i;
    DBL mult = 0.5;

    Vector4d *pIterStack = static_cast<Vector4d *>(pIterData->mainIter.data());

    for (i = 0; i < nMax; i++)
    {
        DerivCalc(d, pIterStack[i], i, mult, pFractal, pIterData);
    }

    return mult * dot(d, pIterStack[nMax]);
}

void MagicHypercomplexFractalRules::
CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    Vector4d d(1.0, 0.0, 1.0, 0.0);
    int i;
    DBL mult = 0.5;

    /*
     * The fact that the functions used for iterating are well-behaved in the hypercomplexes
     * allows for great simplification of computations here.  In particular, the existence of a
     * (sort of) derivative that behaves more-or-less uniformly (e.g., d(i*f)/dz=i*(df/dz)) is of
     * much use.  This is not necessarily the case for a general hypercomplex function, though,
     * and is a problem for many of even the most basic quaternionic functions (e.g., z^2).
     */

    Vector4d *pTIterStack = static_cast<Vector4d *>(pTIterData->mainIter.data()),
        *pPIterStack = (pPIterData != NULL ? static_cast<Vector4d *>(pPIterData->mainIter.data()) : NULL);

    for (i = 0; i < nMax; i++)
    {
        if (pPIterData != NULL && pFractal->Discontinuity_Test > 0)
        {
            Vector4d dc;
            DBL dist;
            if (DiscontinuityCheck(dc, dist, pTIterStack[i], pPIterStack[i],
                                   i, pFractal, pTIterData, pPIterData))
            {
                d[Y] *= -1.0;
                d[W] *= -1.0;

                complex_fn::Mult(AsComplex(d, 0), AsComplex(d, 0), AsComplex(dc, 0));
                complex_fn::Mult(AsComplex(d, 1), AsComplex(d, 1), AsComplex(dc, 1));

                rResult[X] = mult * dot(d, mSpace4D.transformedX());
                rResult[Y] = mult * dot(d, mSpace4D.transformedY());
                rResult[Z] = mult * dot(d, mSpace4D.transformedZ());

                return;
            }
        }

        DerivCalc(d, pTIterStack[i], i, mult, pFractal, pTIterData);
    }

    d[Y] *= -1.0;
    d[W] *= -1.0;

    complex_fn::Mult(AsComplex(d, 0), AsComplex(d, 0), AsComplex(pTIterStack[nMax], 0));
    complex_fn::Mult(AsComplex(d, 1), AsComplex(d, 1), AsComplex(pTIterStack[nMax], 1));

    rResult[X] = mult * dot(d, mSpace4D.transformedX());
    rResult[Y] = mult * dot(d, mSpace4D.transformedY());
    rResult[Z] = mult * dot(d, mSpace4D.transformedZ());
}

bool MagicHypercomplexFractalRules::
DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}

}
