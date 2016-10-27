//******************************************************************************
///
/// @file core/shape/fractal/util.h
///
/// This module contains miscellaneous fractal-related code.
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

#include <cstdarg>
#include <set>

#include "base/pov_err.h"

#include "core/shape/fractal/types.h"

namespace pov
{

static inline void ComputeDuplexFromHypercomplex(Duplex& rd, const VECTOR_4D& h)
{
    rd[0].x = h[X] - h[W];
    rd[0].y = h[Y] + h[Z];
    rd[1].x = h[X] + h[W];
    rd[1].y = h[Y] - h[Z];
}

static inline void ComputeHypercomplexFromDuplex(VECTOR_4D& rh, const Duplex &d)
{
    rh[X] = .5 * (d[0].x + d[1].x);
    rh[Y] = .5 * (d[0].y + d[1].y);
    rh[Z] = .5 * (d[0].y - d[1].y);
    rh[W] = .5 * (d[1].x - d[0].x);
}

static inline void AssignDuplex(Duplex& rD, const Duplex& s)
{
    rD[0] = s[0];
    rD[1] = s[1];
}

static inline const DistanceEstimator& BadEstimator()
{
    throw POV_EXCEPTION_STRING("Unsupported distance estimator for fractal type.");
}

template <typename T>
static inline int GetADataSize()
{
    return (typeid(NilData) == typeid(T) ? 0 : sizeof(T));
}

template <class Rules>
static inline const FractalDataSizes GetDataSizes()
{
    FractalDataSizes sizes = { GetADataSize<typename Rules::FixedData>(),
                               GetADataSize<typename Rules::MainIterData>(),
                               GetADataSize<typename Rules::AuxIterData>() };
    return sizes;
}

template <int n, typename T>
static inline const std::set<T> CreateSet(T t0, ...)
{
    std::set<T> s;
    int i;
    va_list tList;
    T t;
    va_start(tList, t0);
    s.insert(t0);
    for (i = 1; i < n; i++)
    {
        t = va_arg(tList, T);
        s.insert(t);
    }
    va_end(tList);
    return s;
}

/* Helper functions to permit inline creation of certain const structs
   (and const arrays).  These are written in such a way as to permit simple
   replacement for C++11 or higher. */
#ifndef FRACTAL_USE_CXX11

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

static inline const FractalRulesInfo CreateRulesInfo(const FractalFuncType& funcType, DiscontinuitySupportLevel discontinuitySupport,
                                                     const FractalDataSizes& sizes, EstimatorType estimatorType)
{
    FractalRulesInfo f = { funcType, discontinuitySupport, sizes, estimatorType };
    return f;
}

#define INIT_DUPLEX(var, v0, v1) var(v0, v1)
#define INIT_VECTOR_4D(var, qx, qy, qz, qw) var(qx, qy, qz, qw)

#else

/* The easy part. */

#define CREATE_GENERAL(...) {__VA_ARGS__}

#define CreateComplex (Complex)CREATE_GENERAL
#define CreateFuncType (FractalFuncType)CREATE_GENERAL
#define CreateRulesInfo (FractalRulesInfo)CREATE_GENERAL

#define INIT_DUPLEX(var, v0, v1) var{v0, v1}
#define INIT_VECTOR_4D(var, qx, qy, qz, qw) var{qx, qy, qz, qw}

#endif

}

#endif
