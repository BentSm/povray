//******************************************************************************
///
/// @file core/shape/fractal/space.cpp
///
/// This module contains the implementation for handling 4-D fractal spaces.
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

#include "core/shape/fractal/space.h"

#include "base/pov_err.h"

#include "core/math/vector.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/util.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

FractalSpace::FractalSpace(FractalTransformMethod transformMethod, FractalAlgebra algebra, const Vector4d& slice, DBL sliceDist)
{
    mTransformMethod = transformMethod;
    mAlgebra = algebra;
    mSliceActual = slice;
    mSliceDistActual = sliceDist;

    switch (mTransformMethod)
    {
    case kTransformProjection:
        mSliceNorm3d = Vector3d(slice / slice.w());
        mSliceDistNorm = sliceDist / slice.w();
        mTrX = Vector4d(1.0, 0.0, 0.0, -mSliceNorm3d.x());
        mTrY = Vector4d(0.0, 1.0, 0.0, -mSliceNorm3d.y());
        mTrZ = Vector4d(0.0, 0.0, 1.0, -mSliceNorm3d.z());
        mTr0 = Vector4d(0.0, 0.0, 0.0, mSliceDistNorm);
        break;
    case kTransformIsometric:
        mTrX = Vector4d(-slice.y(), slice.x(), -slice.w(), slice.z());
        mTrY = Vector4d(-slice.z(), slice.w(), slice.x(), -slice.y());
        mTrZ = Vector4d(-slice.w(), -slice.z(), slice.y(), slice.x());
        mTr0 = slice * sliceDist;
        break;
    default:
        throw POV_EXCEPTION_STRING("Unknown fractal mapping type.");
    }

    switch (mAlgebra)
    {
    case kQuaternion:
        break;
    case kHypercomplex:
        mTrX = DuplexFromHypercomplex(mTrX);
        mTrY = DuplexFromHypercomplex(mTrY);
        mTrZ = DuplexFromHypercomplex(mTrZ);
        mTr0 = DuplexFromHypercomplex(mTr0);
        break;
    default:
        throw POV_EXCEPTION_STRING("Unknown algebra for fractal mapping.");
    }
}

const Vector4d FractalSpace::
TransformTo4D(const Vector3d& point) const
{
    return mTrX * point.x() + mTrY * point.y() + mTrZ * point.z() + mTr0;
}

const Vector4d FractalSpace::
TransformDirTo4D(const Vector3d& dir) const
{
    return mTrX * dir.x() + mTrY * dir.y() + mTrZ * dir.z();
}

// This was adapted from the Sphere::Intersect code
bool FractalSpace::
Bound(const BasicRay& ray, const Fractal *pFractal, DBL *pDepthMin, DBL *pDepthMax) const
{
    DBL closestApproach, halfChord, halfChordSquared, ctrToOrigSquared, ctrToOrigWComp, dirSquared, dirWComp;
    Vector3d centerToOrigin;

    centerToOrigin = ray.Origin - pFractal->Center;

    if (mTransformMethod == kTransformIsometric)
    {
        ctrToOrigSquared = centerToOrigin.lengthSqr() + Sqr(mSliceDistActual);
        dirSquared = ray.Direction.lengthSqr();

        closestApproach = -dot(centerToOrigin, ray.Direction) / dirSquared;
    }
    else
    {
        ctrToOrigWComp = mSliceDistNorm - dot(centerToOrigin, mSliceNorm3d);
        dirWComp = -dot(ray.Direction, mSliceNorm3d);

        ctrToOrigSquared = centerToOrigin.lengthSqr() + Sqr(ctrToOrigWComp);
        dirSquared = ray.Direction.lengthSqr() + Sqr(dirWComp);

        closestApproach = -(dot(centerToOrigin, ray.Direction) + ctrToOrigWComp * dirWComp) / dirSquared;
    }

    if ((ctrToOrigSquared >= pFractal->Radius_Squared) && (closestApproach < EPSILON))
        return false;

    halfChordSquared = (pFractal->Radius_Squared - ctrToOrigSquared) / dirSquared + Sqr(closestApproach);

    if (halfChordSquared > EPSILON)
    {
        halfChord = sqrt(halfChordSquared);

        *pDepthMin = closestApproach - halfChord;
        *pDepthMax = closestApproach + halfChord;

        return true;
    }

    return false;
}

bool FractalSpace::
Compute_BBox(BoundingBox& rBBox, const Fractal *pFractal) const
{
    DBL dx, dy, dz, x0, y0, z0, r, sliceRadSquared;

    sliceRadSquared = pFractal->Radius_Squared - Sqr(mSliceDistActual);

    if (sliceRadSquared < 0)
        return false;

    switch (mTransformMethod)
    {
    case kTransformProjection:

        dx = sqrt(sliceRadSquared * (1 - Sqr(mSliceActual.x())));
        x0 = mSliceActual.x() * mSliceDistActual - dx;

        dy = sqrt(sliceRadSquared * (1 - Sqr(mSliceActual.y())));
        y0 = mSliceActual.y() * mSliceDistActual - dy;

        dz = sqrt(sliceRadSquared * (1 - Sqr(mSliceActual.z())));
        z0 = mSliceActual.z() * mSliceDistActual - dz;

        Make_BBox(rBBox, x0, y0, z0, 2.0 * dx, 2.0 * dy, 2.0 * dz);

        break;

    case kTransformIsometric:

        r = sqrt(sliceRadSquared);
        Make_BBox(rBBox, -r, -r, -r, 2.0 * r, 2.0 * r, 2.0 * r);

        break;

    default:
        throw POV_EXCEPTION_STRING("Unknown fractal mapping type.");
    }

    return true;
}

}

