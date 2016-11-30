//******************************************************************************
///
/// @file core/shape/fractal/space.h
///
/// This module contains prototypes for the handling of 4-D fractal spaces.
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

#ifndef POVRAY_CORE_FRACTAL_SPACE_H
#define POVRAY_CORE_FRACTAL_SPACE_H

#include "core/coretypes.h"
#include "core/math/vector.h"
#include "core/shape/fractal/types.h"

namespace pov
{

class FractalSpace
{
public:
    FractalSpace(FractalTransformMethod transformMethod, FractalAlgebra algebra, const Vector4d& slice, DBL sliceDist);
    const Vector4d& transformedX() const { return mTrX; }
    const Vector4d& transformedY() const { return mTrY; }
    const Vector4d& transformedZ() const { return mTrZ; }
    const Vector4d& transformed0() const { return mTr0; }

    const Vector4d TransformTo4D(const Vector3d& point) const;
    const Vector4d TransformDirTo4D(const Vector3d& dir) const;
    bool Bound(const BasicRay& ray, const Fractal *pFractal, DBL *pDepthMin, DBL *pDepthMax) const;
    bool Compute_BBox(BoundingBox& rBBox, const Fractal *pFractal) const;

protected:
    FractalTransformMethod mTransformMethod;
    FractalAlgebra mAlgebra;
    Vector4d mTrX, mTrY, mTrZ, mTr0, mSliceActual;
    Vector3d mSliceNorm3d;
    DBL mSliceDistActual, mSliceDistNorm;

};

}

#endif
