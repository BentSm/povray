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
        return BadEstimator();

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
            return BadEstimator();
    }
}

// This was adapted from the Sphere::Intersect code
bool MagicRulesBase::
Bound(const BasicRay& ray, const Fractal *pFractal, DBL *pDepthMin, DBL *pDepthMax) const
{
    DBL OCSquared, t_Closest_Approach, Half_Chord, t_Half_Chord_Squared, OCWComp, DirSquared, DirWComp;
    Vector3d Center_To_Origin;

    Center_To_Origin = ray.Origin - pFractal->Center;

    OCWComp = pFractal->SliceDistNorm - dot(Center_To_Origin, pFractal->SliceNorm);

    DirWComp = -dot(ray.Direction, pFractal->SliceNorm);

    OCSquared = Center_To_Origin.lengthSqr() + Sqr(OCWComp);

    DirSquared = ray.Direction.lengthSqr() + Sqr(DirWComp);

    t_Closest_Approach = -(dot(Center_To_Origin, ray.Direction) + OCWComp * DirWComp) / DirSquared;

    if ((OCSquared >= pFractal->Radius_Squared) && (t_Closest_Approach < EPSILON))
        return false;

    t_Half_Chord_Squared = (pFractal->Radius_Squared - OCSquared) / DirSquared + Sqr(t_Closest_Approach);

    if (t_Half_Chord_Squared > EPSILON)
    {
        Half_Chord = sqrt(t_Half_Chord_Squared);

        *pDepthMin = t_Closest_Approach - Half_Chord;
        *pDepthMax = t_Closest_Approach + Half_Chord;

        return true;
    }

    return false;
}


int MagicQuaternionFractalRules::
Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction, DBL *pDist,
        FractalIterData *pIterData) const
{
    int i;
    VECTOR_4D v = {iPoint[X], iPoint[Y], iPoint[Z], pFractal->SliceDistNorm - dot(pFractal->SliceNorm, iPoint)};
    DBL norm, exitValue;

    VECTOR_4D *pIterStack = static_cast<VECTOR_4D *>(pIterData->mainIter.data());

    Assign_Vector_4D(pIterStack[0], v);

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        V4D_Dot(norm, v, v);

        if (norm > exitValue)
        {
            if (pDist != NULL)
            {
                *pDist = (*(mEstimator.pEstim))(this, norm, i, direction, pFractal, pIterData);
            }

            return i;
        }

        IterateCalc(v, norm, i, pFractal, pIterData);

        Assign_Vector_4D(pIterStack[i+1], v);

    }

    V4D_Dot(norm, v, v);

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

void MagicQuaternionFractalRules::
CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    VECTOR_4D nX = {1.0, 0.0, 0.0, -pFractal->SliceNorm[X]},
              nY = {0.0, 1.0, 0.0, -pFractal->SliceNorm[Y]},
              nZ = {0.0, 0.0, 1.0, -pFractal->SliceNorm[Z]};
    int i;

    VECTOR_4D *pTIterStack = static_cast<VECTOR_4D *>(pTIterData->mainIter.data()),
        *pPIterStack = (pPIterData != NULL ? static_cast<VECTOR_4D *>(pPIterData->mainIter.data()) : NULL);

    for (i = 0; i < nMax; i++)
    {
        if (pPIterData != NULL && pFractal->Discontinuity_Test > 0)
        {
            VECTOR_4D d;
            DBL dist;
            if (DiscontinuityCheck(d, dist, pTIterStack[i], pPIterStack[i],
                                   i, pFractal, pTIterData, pPIterData))
            {
                V4D_Dot(rResult[X], nX, d);
                V4D_Dot(rResult[Y], nY, d);
                V4D_Dot(rResult[Z], nZ, d);
                return;
            }
        }

        DirDerivCalc(nX, pTIterStack[i], i, false, pFractal, pTIterData);
        DirDerivCalc(nY, pTIterStack[i], i, true, pFractal, pTIterData);
        DirDerivCalc(nZ, pTIterStack[i], i, true, pFractal, pTIterData);
    }

    V4D_Dot(rResult[X], nX, pTIterStack[nMax]);
    V4D_Dot(rResult[Y], nY, pTIterStack[nMax]);
    V4D_Dot(rResult[Z], nZ, pTIterStack[nMax]);
}

DBL MagicQuaternionFractalRules::
CalcDirDeriv(const Vector3d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const
{
    VECTOR_4D d = {dir[X], dir[Y], dir[Z], -dot(pFractal->SliceNorm, dir)};
    int i;
    DBL res;

    VECTOR_4D *pIterStack = static_cast<VECTOR_4D *>(pIterData->mainIter.data());

    for (i = 0; i < nMax; i++)
    {
        DirDerivCalc(d, pIterStack[i], i, false, pFractal, pIterData);
    }

    V4D_Dot(res, d, pIterStack[nMax]);
    return res;
}

bool MagicQuaternionFractalRules::
DiscontinuityCheck(VECTOR_4D& rD, DBL& rDist, const VECTOR_4D& t, const VECTOR_4D& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}


int MagicHypercomplexFractalRules::
Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction, DBL *pDist,
        FractalIterData *pIterData) const
{
    int i;
    Duplex d;
    DBL norm, exitValue;
    VECTOR_4D v = {iPoint[X], iPoint[Y], iPoint[Z], pFractal->SliceDistNorm - dot(pFractal->SliceNorm, iPoint)};

    Duplex *pIterStack = static_cast<Duplex *>(pIterData->mainIter.data());

    ComputeDuplexFromHypercomplex(d, v);

    AssignDuplex(pIterStack[0], d);

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        norm = 0.5 * (d[0].x * d[0].x + d[0].y * d[0].y + d[1].x * d[1].x + d[1].y * d[1].y);

        if (norm > exitValue)
        {
            if (pDist != NULL)
            {
                *pDist = (*(mEstimator.pEstim))(this, norm, i, direction, pFractal, pIterData);
            }

            return i;
        }

        IterateCalc(d, norm, i, pFractal, pIterData);

        AssignDuplex(pIterStack[i+1], d);

    }

    norm = 0.5 * (d[0].x * d[0].x + d[0].y * d[0].y + d[1].x * d[1].x + d[1].y * d[1].y);

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

void MagicHypercomplexFractalRules::
CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    int i;
    Duplex d = {{1.0, 0.0}, {1.0, 0.0}};
    VECTOR_4D v;

    /*
     * The fact that the functions used for iterating are well-behaved in the hypercomplexes
     * allows for great simplification of computations here.  In particular, the existence of a
     * (sort of) derivative that behaves more-or-less uniformly (e.g., d(i*f)/dz=i*(df/dz)) is of
     * much use.  This is not necessarily the case for a general hypercomplex function, though,
     * and is a problem for many of even the most basic quaternionic functions (e.g., z^2).
     */

    Duplex *pTIterStack = static_cast<Duplex *>(pTIterData->mainIter.data()),
        *pPIterStack = (pPIterData != NULL ? static_cast<Duplex *>(pPIterData->mainIter.data()) : NULL);

    for (i = 0; i < nMax; i++)
    {
        if (pPIterData != NULL && pFractal->Discontinuity_Test > 0)
        {
            Duplex dc;
            DBL dist;
            if (DiscontinuityCheck(dc, dist, pTIterStack[i], pPIterStack[i],
                                   i, pFractal, pTIterData, pPIterData))
            {
                d[0].y *= -1.0;
                d[1].y *= -1.0;

                complex_fn::Mult(d[0], d[0], dc[0]);
                complex_fn::Mult(d[1], d[1], dc[1]);

                ComputeHypercomplexFromDuplex(v, d);

                rResult[X] = (v[X] - v[W] * pFractal->SliceNorm[X]);
                rResult[Y] = (v[Y] - v[W] * pFractal->SliceNorm[Y]);
                rResult[Z] = (v[Z] - v[W] * pFractal->SliceNorm[Z]);

                return;
            }
        }

        DerivCalc(d, pTIterStack[i], i, pFractal, pTIterData);
    }

    d[0].y *= -1.0;
    d[1].y *= -1.0;

    complex_fn::Mult(d[0], d[0], pTIterStack[nMax][0]);
    complex_fn::Mult(d[1], d[1], pTIterStack[nMax][1]);

    ComputeHypercomplexFromDuplex(v, d);

    rResult[X] = (v[X] - v[W] * pFractal->SliceNorm[X]);
    rResult[Y] = (v[Y] - v[W] * pFractal->SliceNorm[Y]);
    rResult[Z] = (v[Z] - v[W] * pFractal->SliceNorm[Z]);
}

bool MagicHypercomplexFractalRules::
DiscontinuityCheck(Duplex& rD, DBL& rDist, const Duplex& t, const Duplex& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}

}

