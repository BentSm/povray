//******************************************************************************
///
/// @file core/shape/fractal/util.h
///
/// This module contains miscellaneous fractal-related code.  If you can think
/// of a good home for any (or all) of this, feel free to move things around.
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

#ifndef POVRAY_CORE_FRACTAL_UTIL_H
#define POVRAY_CORE_FRACTAL_UTIL_H

#include "core/coretypes.h"

namespace pov
{

static inline void ComputeDuplexFromHypercomplex(Complex& rc0, Complex& rc1, DBL x, DBL y, DBL z, DBL w)
{
    rc0.x = x - w;
    rc0.y = y + z;
    rc1.x = x + w;
    rc1.y = y - z;
}

static inline void ComputeHypercomplexFromDuplex(DBL& rx, DBL& ry, DBL& rz, DBL& rw,
                                                 const Complex& c0, const Complex& c1)
{
    rx = .5 * (c0.x + c1.x);
    ry = .5 * (c0.y + c1.y);
    rz = .5 * (c0.y - c1.y);
    rw = .5 * (c1.x - c0.x);
}

}

#endif
