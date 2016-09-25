//******************************************************************************
///
/// @file core/shape/fractal/magicimpl.h
///
/// This module contains the generic implementation of FractalRules subclasses.
///
/// (This was originally named for the templating 'magic' that it uses, but
/// could equally well be described as the 'magic' that makes everything work!)
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

#ifndef POVRAY_CORE_FRACTAL_MAGICIMPL_H
#define POVRAY_CORE_FRACTAL_MAGICIMPL_H

#include "core/coretypes.h"
#include "core/shape/fractal/magic.h"

#include "base/pov_err.h"

#include "core/math/complexfn.h"
#include "core/math/vector.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/util.h"

namespace pov
{

template <template <class> class RulesClass, class Estimator, class BaseRules>
void MagicFractalRulesBase<RulesClass, Estimator, BaseRules>::
CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    static_cast<const RulesClass<Estimator> *>(this)->NonVirtualCalcNormalWDisc(rResult, nMax, pFractal, pTIterData, pPIterData);
}

// This was adapted from the Sphere::Intersect code
template <template <class> class RulesClass, class Estimator, class BaseRules>
bool MagicFractalRulesBase<RulesClass, Estimator, BaseRules>::
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
        return(false);

    t_Half_Chord_Squared = (pFractal->Radius_Squared - OCSquared) / DirSquared + Sqr(t_Closest_Approach);

    if (t_Half_Chord_Squared > EPSILON)
    {
        Half_Chord = sqrt(t_Half_Chord_Squared);

        *pDepthMin = t_Closest_Approach - Half_Chord;
        *pDepthMax = t_Closest_Approach + Half_Chord;

        return(true);
    }

    return(false);
}

template <template <class> class RulesClass, class Estimator, class BaseRules>
int MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>::
Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction, DBL *pDist,
        void *pIterData) const
{
    int i;
    VECTOR_4D v = {iPoint[X], iPoint[Y], iPoint[Z], pFractal->SliceDistNorm - dot(pFractal->SliceNorm, iPoint)};
    DBL norm, exitValue;

    typename RulesClass<Estimator>::IterationData *pIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pIterData);

    Assign_Vector_4D(pIterStack[0].point, v);

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        V4D_Dot(norm, v, v);

        if (norm > exitValue)
        {
            if (pDist != NULL)
            {
                *pDist = Estimator::Call(static_cast<const RulesClass<Estimator> *>(this), norm, i, direction,
                                         pFractal, pIterData);
            }

            return i;
        }

        static_cast<const RulesClass<Estimator> *>(this)->
            IterateCalc(v, norm, i, pFractal, pIterData);

        Assign_Vector_4D(pIterStack[i+1].point, v);

    }

    V4D_Dot(norm, v, v);

    if (norm > exitValue)
    {
        if (pDist != NULL)
        {
            *pDist = Estimator::Call(static_cast<const RulesClass<Estimator> *>(this), norm,
                                     pFractal->Num_Iterations, direction, pFractal, pIterData);
        }

        return pFractal->Num_Iterations;
    }

    if (pDist != NULL)
    {
        *pDist = pFractal->Precision;
    }

    return pFractal->Num_Iterations + 1;

}

template <template <class> class RulesClass, class Estimator, class BaseRules>
template <bool disc>
void MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>::
TemplatedCalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    VECTOR_4D nX = {1.0, 0.0, 0.0, -pFractal->SliceNorm[X]},
              nY = {0.0, 1.0, 0.0, -pFractal->SliceNorm[Y]},
              nZ = {0.0, 0.0, 1.0, -pFractal->SliceNorm[Z]};
    int i;

    typename RulesClass<Estimator>::IterationData *pTIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pTIterData),
        *pPIterStack = reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pPIterData);

    for (i = 0; i < nMax; i++)
    {
        if (disc && pFractal->Discontinuity_Test > 0)
        {
            VECTOR_4D d;
            DBL dist;
            if (static_cast<const RulesClass<Estimator> *>(this)->
                DiscontinuityCheck(d, dist, pTIterStack[i].point, pPIterStack[i].point,
                                   i, pFractal, pTIterData, pPIterData))
            {
                V4D_Dot(rResult[X], nX, d);
                V4D_Dot(rResult[Y], nY, d);
                V4D_Dot(rResult[Z], nZ, d);
                return;
            }
        }

        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(nX, pTIterStack[i].point, i, false, pFractal, pTIterData);
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(nY, pTIterStack[i].point, i, true, pFractal, pTIterData);
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(nZ, pTIterStack[i].point, i, true, pFractal, pTIterData);
    }

    V4D_Dot(rResult[X], nX, pTIterStack[nMax].point);
    V4D_Dot(rResult[Y], nY, pTIterStack[nMax].point);
    V4D_Dot(rResult[Z], nZ, pTIterStack[nMax].point);
}

template <template <class> class RulesClass, class Estimator, class BaseRules>
DBL MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>::
CalcDirDeriv(const Vector3d& dir, int nMax, const Fractal *pFractal, void *pIterData) const
{
    VECTOR_4D d = {dir[X], dir[Y], dir[Z], -dot(pFractal->SliceNorm, dir)};
    int i;
    DBL res;

    typename RulesClass<Estimator>::IterationData *pIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pIterData);

    for (i = 0; i < nMax; i++)
    {
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(d, pIterStack[i].point, i, false, pFractal, pIterData);
    }

    V4D_Dot(res, d, pIterStack[nMax].point);
    return res;
}

template <template <class> class RulesClass, class Estimator, class BaseRules>
bool MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>::
DiscontinuityCheck(VECTOR_4D& rD, DBL& rDist, const VECTOR_4D& t, const VECTOR_4D& p,
                   int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}


template <template <class> class RulesClass, class Estimator, class BaseRules>
int MagicHypercomplexFractalRules<RulesClass, Estimator, BaseRules>::
Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction, DBL *pDist,
        void *pIterData) const
{
    int i;
    Duplex d;
    DBL norm, exitValue;
    VECTOR_4D v = {iPoint[X], iPoint[Y], iPoint[Z], pFractal->SliceDistNorm - dot(pFractal->SliceNorm, iPoint)};

    typename RulesClass<Estimator>::IterationData *pIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pIterData);

    ComputeDuplexFromHypercomplex(d, v);

    AssignDuplex(pIterStack[0].point, d);

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        norm = 0.5 * (d[0].x * d[0].x + d[0].y * d[0].y + d[1].x * d[1].x + d[1].y * d[1].y);

        if (norm > exitValue)
        {
            if (pDist != NULL)
            {
                *pDist = Estimator::Call(static_cast<const RulesClass<Estimator> *>(this), norm, i, direction,
                                         pFractal, pIterData);
            }

            return i;
        }

        static_cast<const RulesClass<Estimator> *>(this)->IterateCalc(d, norm, i, pFractal, pIterData);

        AssignDuplex(pIterStack[i+1].point, d);

    }

    norm = 0.5 * (d[0].x * d[0].x + d[0].y * d[0].y + d[1].x * d[1].x + d[1].y * d[1].y);

    if (norm > exitValue)
    {
        if (pDist != NULL)
        {
            *pDist = Estimator::Call(static_cast<const RulesClass<Estimator> *>(this), norm,
                                     pFractal->Num_Iterations, direction, pFractal, pIterData);
        }

        return pFractal->Num_Iterations;
    }

    if (pDist != NULL)
    {
        *pDist = pFractal->Precision;
    }

    return pFractal->Num_Iterations + 1;

}

template <template <class> class RulesClass, class Estimator, class BaseRules>
template <bool disc>
void MagicHypercomplexFractalRules<RulesClass, Estimator, BaseRules>::
TemplatedCalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
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

    typename RulesClass<Estimator>::IterationData *pTIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pTIterData),
        *pPIterStack = reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pPIterData);


    for (i = 0; i < nMax; i++)
    {
        if (disc && pFractal->Discontinuity_Test > 0)
        {
            Duplex dc;
            DBL dist;
            if (static_cast<const RulesClass<Estimator> *>(this)->
                DiscontinuityCheck(dc, dist, pTIterStack[i].point, pPIterStack[i].point,
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
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDerivCalc(d, pTIterStack[i].point, i, pFractal, pTIterData);
    }

    d[0].y *= -1.0;
    d[1].y *= -1.0;

    complex_fn::Mult(d[0], d[0], pTIterStack[nMax].point[0]);
    complex_fn::Mult(d[1], d[1], pTIterStack[nMax].point[1]);

    ComputeHypercomplexFromDuplex(v, d);

    rResult[X] = (v[X] - v[W] * pFractal->SliceNorm[X]);
    rResult[Y] = (v[Y] - v[W] * pFractal->SliceNorm[Y]);
    rResult[Z] = (v[Z] - v[W] * pFractal->SliceNorm[Z]);
}

template <template <class> class RulesClass, class Estimator, class BaseRules>
bool MagicHypercomplexFractalRules<RulesClass, Estimator, BaseRules>::
DiscontinuityCheck(Duplex& rD, DBL& rDist, const Duplex& t, const Duplex& p,
                   int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}

}

#endif
