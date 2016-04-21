//******************************************************************************
///
/// @file core/shape/fractal/func.h
///
/// This module contains all defines, typedefs, and prototypes for `func.cpp`.
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

#ifndef POVRAY_CORE_FRACTAL_FUNC_H
#define POVRAY_CORE_FRACTAL_FUNC_H

#include "core/coretypes.h"

#include <map>
#include <set>

#include "core/math/complexfn.h"
#include "core/shape/fractal/cxx03compat.h"
#include "core/shape/fractal/types.h"

namespace pov
{

struct FractalFunc
{
    FractalFuncType type;
    complex_fn::FuncType *pFunc, *pDeriv;
    complex_fn::DiscontinuityTestFn *pDisc;
};

const std::map<FractalFuncType, FractalFunc>& AllFractalFuncs();

const FractalFunc& FractalFuncForType(const FractalFuncType& funcType);

inline DiscontinuitySupportLevel DiscontinuitySupport_Func(const FractalFunc& f)
{
    if (f.pDisc == NULL)
        return kDiscontinuityNotImplemented;
    else if (f.pDisc == &(complex_fn::False_DTest))
        return kDiscontinuityUnneeded;
    else
        return kDiscontinuitySupported;
}

inline const std::set<FractalFuncType> Func_FuncTypeSet(FractalAlgebra algebra)
{
    std::set<FractalFuncType> s;
    std::map<FractalFuncType, FractalFunc>::const_iterator f = AllFractalFuncs().begin();
    while (f != AllFractalFuncs().end())
    {
        s.insert(CreateFuncType(algebra, f->first.type, f->first.variant));
        f++;
    }
    return s;
}

}

#endif

