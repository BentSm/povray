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
    DBL x, y, z, w;
    DBL norm, exitValue;

    typename RulesClass<Estimator>::IterationData *pIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pIterData);

    x = pIterStack[0].point[X] = iPoint[X];
    y = pIterStack[0].point[Y] = iPoint[Y];
    z = pIterStack[0].point[Z] = iPoint[Z];
    w = pIterStack[0].point[W] = pFractal->SliceDistNorm -
                                 pFractal->SliceNorm[X] * x -
                                 pFractal->SliceNorm[Y] * y -
                                 pFractal->SliceNorm[Z] * z;

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        norm = x * x + y * y + z * z + w * w;

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
            IterateCalc(x, y, z, w, norm, i, pFractal, pIterData);

        pIterStack[i+1].point[X] = x;
        pIterStack[i+1].point[Y] = y;
        pIterStack[i+1].point[Z] = z;
        pIterStack[i+1].point[W] = w;

    }

    norm = x * x + y * y + z * z + w * w;

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
    DBL n11 = 1.0, n12 = 0.0, n13 = 0.0, n14 = -pFractal->SliceNorm[X],
	n21 = 0.0, n22 = 1.0, n23 = 0.0, n24 = -pFractal->SliceNorm[Y],
	n31 = 0.0, n32 = 0.0, n33 = 1.0, n34 = -pFractal->SliceNorm[Z];
    int i;

    typename RulesClass<Estimator>::IterationData *pTIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pTIterData),
        *pPIterStack = reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pPIterData);


    for (i = 0; i < nMax; i++)
    {
        if (disc && pFractal->Discontinuity_Test > 0)
        {
            DBL dx, dy, dz, dw, dist;
            if (static_cast<const RulesClass<Estimator> *>(this)->
                DiscontinuityCheck(dx, dy, dz, dw, dist,
                                   pTIterStack[i].point[X], pTIterStack[i].point[Y],
                                   pTIterStack[i].point[Z], pTIterStack[i].point[W],
                                   pPIterStack[i].point[X], pPIterStack[i].point[Y],
                                   pPIterStack[i].point[Z], pPIterStack[i].point[W],
                                   i, pFractal, pTIterData, pPIterData))
            {
                rResult[X] = n11 * dx + n12 * dy + n13 * dz + n14 * dw;
                rResult[Y] = n21 * dx + n22 * dy + n23 * dz + n24 * dw;
                rResult[Z] = n31 * dx + n32 * dy + n33 * dz + n34 * dw;
                return;
            }
        }

        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(n11, n12, n13, n14,
                              pTIterStack[i].point[X], pTIterStack[i].point[Y],
                              pTIterStack[i].point[Z], pTIterStack[i].point[W],
                              i, false, pFractal, pTIterData);
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(n21, n22, n23, n24,
                              pTIterStack[i].point[X], pTIterStack[i].point[Y],
                              pTIterStack[i].point[Z], pTIterStack[i].point[W],
                              i, true, pFractal, pTIterData);
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(n31, n32, n33, n34,
                              pTIterStack[i].point[X], pTIterStack[i].point[Y],
                              pTIterStack[i].point[Z], pTIterStack[i].point[W],
                              i, true, pFractal, pTIterData);
    }

    rResult[X] = n11 * pTIterStack[nMax].point[X] + n12 * pTIterStack[nMax].point[Y] +
        n13 * pTIterStack[nMax].point[Z] + n14 * pTIterStack[nMax].point[W];
    rResult[Y] = n21 * pTIterStack[nMax].point[X] + n22 * pTIterStack[nMax].point[Y] +
        n23 * pTIterStack[nMax].point[Z] + n24 * pTIterStack[nMax].point[W];
    rResult[Z] = n31 * pTIterStack[nMax].point[X] + n32 * pTIterStack[nMax].point[Y] +
        n33 * pTIterStack[nMax].point[Z] + n34 * pTIterStack[nMax].point[W];
}

template <template <class> class RulesClass, class Estimator, class BaseRules>
DBL MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>::
CalcDirDeriv(const Vector3d& dir, int nMax, const Fractal *pFractal, void *pIterData) const
{
    DBL dx = dir[X], dy = dir[Y], dz = dir[Z], dw = -dot(pFractal->SliceNorm, dir);
    int i;

    typename RulesClass<Estimator>::IterationData *pIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pIterData);

    for (i = 0; i < nMax; i++)
    {
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDirDerivCalc(dx, dy, dz, dw,
                              pIterStack[i].point[X], pIterStack[i].point[Y],
                              pIterStack[i].point[Z], pIterStack[i].point[W],
                              i, false, pFractal, pIterData);
    }

    return dx * pIterStack[nMax].point[X] + dy * pIterStack[nMax].point[Y] +
        dz * pIterStack[nMax].point[Z] + dw * pIterStack[nMax].point[W];
}

template <template <class> class RulesClass, class Estimator, class BaseRules>
bool MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>::
DiscontinuityCheck(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL& rDist,
                   DBL tx, DBL ty, DBL tz, DBL tw, DBL px, DBL py, DBL pz, DBL pw,
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
    Complex c0, c1;
    DBL norm, exitValue;

    typename RulesClass<Estimator>::IterationData *pIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pIterData);

    ComputeDuplexFromHypercomplex(c0, c1, iPoint[X], iPoint[Y], iPoint[Z],
                                  pFractal->SliceDistNorm -
                                  pFractal->SliceNorm[X] * iPoint[X] -
                                  pFractal->SliceNorm[Y] * iPoint[Y] -
                                  pFractal->SliceNorm[Z] * iPoint[Z]);

    pIterStack[0].point[0] = c0;
    pIterStack[0].point[1] = c1;

    exitValue = pFractal->Exit_Value;

    for (i = 0; i < pFractal->Num_Iterations; i++)
    {
        norm = 0.5 * (complex_fn::Norm(c0) + complex_fn::Norm(c1));

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
            IterateCalc(c0, c1, norm, i, pFractal, pIterData);

        pIterStack[i+1].point[0] = c0;
        pIterStack[i+1].point[1] = c1;

    }

    norm = 0.5 * (complex_fn::Norm(c0) + complex_fn::Norm(c1));

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
    Complex c0, c1;
    DBL x, y, z, w;

    /*
     * The fact that the functions used for iterating are well-behaved in the hypercomplexes
     * allows for great simplification of computations here.  In particular, the existence of a
     * (sort of) derivative that behaves more-or-less uniformly (e.g., d(i*f)/dz=i*(df/dz)) is of
     * much use.  This is not necessarily the case for a general hypercomplex function, though,
     * and is a problem for many of even the most basic quaternionic functions (e.g., z^2).
     */

    c0 = c1 = 1.0;

    typename RulesClass<Estimator>::IterationData *pTIterStack =
        reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pTIterData),
        *pPIterStack = reinterpret_cast<typename RulesClass<Estimator>::IterationData *>(pPIterData);


    for (i = 0; i < nMax; i++)
    {
        if (disc && pFractal->Discontinuity_Test > 0)
        {
            Complex dc0, dc1;
            DBL dist;
            if (static_cast<const RulesClass<Estimator> *>(this)->
                DiscontinuityCheck(dc0, dc1, dist,
                                   pTIterStack[i].point[0], pTIterStack[i].point[1],
                                   pPIterStack[i].point[0], pPIterStack[i].point[1],
                                   i, pFractal, pTIterData, pPIterData))
            {
                c0.imag(-c0.imag());
                c1.imag(-c1.imag());

                complex_fn::Mult(c0, c0, dc0);
                complex_fn::Mult(c1, c1, dc1);

                ComputeHypercomplexFromDuplex(x, y, z, w, c0, c1);

                rResult[X] = (x - w * pFractal->SliceNorm[X]);
                rResult[Y] = (y - w * pFractal->SliceNorm[Y]);
                rResult[Z] = (z - w * pFractal->SliceNorm[Z]);

                return;
            }
        }
        static_cast<const RulesClass<Estimator> *>(this)->
            ApplyDerivCalc(c0, c1, pTIterStack[i].point[0], pTIterStack[i].point[1], i, pFractal,
                           pTIterData);
    }

    c0.imag(-c0.imag());
    c1.imag(-c1.imag());

    complex_fn::Mult(c0, c0, pTIterStack[nMax].point[0]);
    complex_fn::Mult(c1, c1, pTIterStack[nMax].point[1]);

    ComputeHypercomplexFromDuplex(x, y, z, w, c0, c1);

    rResult[X] = (x - w * pFractal->SliceNorm[X]);
    rResult[Y] = (y - w * pFractal->SliceNorm[Y]);
    rResult[Z] = (z - w * pFractal->SliceNorm[Z]);
}

template <template <class> class RulesClass, class Estimator, class BaseRules>
bool MagicHypercomplexFractalRules<RulesClass, Estimator, BaseRules>::
DiscontinuityCheck(Complex& rD0, Complex& rD1, DBL& rDist,
                   const Complex& t0, const Complex& t1, const Complex& p0, const Complex& p1,
                   int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    throw POV_EXCEPTION_STRING("Discontinuity detection not supported for this fractal type.");
}

}

#endif
