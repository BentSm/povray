//******************************************************************************
///
/// @file core/shape/fractal/cxx03compat.h
///
/// This module contains C++03-compatible replacements for initializer-lists.
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

#ifndef POVRAY_CORE_FRACTAL_CXX03COMPAT_H
#define POVRAY_CORE_FRACTAL_CXX03COMPAT_H

#include "core/coretypes.h"

#include "core/shape/fractal/types.h"

namespace pov
{

/* Helper functions to permit inline creation of certain const structs
   (and const arrays).  These are written in such a way as to permit simple
   replacement for C++11 or higher. */
/* C++11 initializer-lists are nice.  We're targeting C++03, though. */
#ifdef FRACTAL_USE_CXX11

#define CREATE_GENERAL(...) {__VA_ARGS__}

#define CreateComplex (Complex)CREATE_GENERAL
#define CreateFuncType (FractalFuncType)CREATE_GENERAL
#define CreateRulesInfo (FractalRulesInfo)CREATE_GENERAL

typedef Complex Duplex[2];
typedef DBL Quaternion[4];

#define INIT_DUPLEX(var, v0, v1) var{v0, v1}
#define INIT_QUATERNION(var, qx, qy, qz, qw) var{qx, qy, qz, qw}

#else

#include <cstdarg>

static inline const Complex CreateComplex(DBL x, DBL y)
{
    Complex c = { x, y };
    return c;
}

static inline const FractalFuncType CreateFuncType(FractalAlgebra algebra, FractalFunc_FuncType type, FractalFunc_VariantType variant)
{
    FractalFuncType f = { algebra, type, variant };
    return f;
}

static inline const FractalRulesInfo CreateRulesInfo(const FractalFuncType& funcType, EstimatorType estimatorType,
                                                     DiscontinuitySupportLevel discontinuitySupport, int iterationDataSize)
{
    FractalRulesInfo f = { funcType, estimatorType, discontinuitySupport, iterationDataSize };
    return f;
}

template <class T, unsigned N>
class IArray
{
public:
    IArray(const T t0, ...)
    {
        va_list vals;
        int i = 0;
        va_start(vals, t0);
        value[0] = t0;
        for (i = 1; i < N; i++)
            value[i] = va_arg(vals, const T);
        va_end(vals);
    }
    T& operator[](unsigned k) { return value[k]; }
    const T& operator[](unsigned k) const { return value[k]; }

private:
    T value[N];
};

typedef IArray<Complex, 2> Duplex;
typedef IArray<DBL, 4> Quaternion;

#define INIT_DUPLEX(var, v0, v1) var(v0, v1)
#define INIT_QUATERNION(var, qx, qy, qz, qw) var(qx, qy, qz, qw)

#endif

}

#endif
