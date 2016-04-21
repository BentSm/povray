//******************************************************************************
///
/// @file core/shape/fractal/dispatch.h
///
/// Handles dispatching rule-creation requests to the appropriate class.
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

#ifndef POVRAY_CORE_FRACTAL_DISPATCH_H
#define POVRAY_CORE_FRACTAL_DISPATCH_H

#include "core/coretypes.h"

#include <cstdarg>
#include <map>
#include <set>

#include "core/shape/fractal/types.h"

namespace pov
{

class RulesDispatch
{
public:
    RulesDispatch(const FractalFuncType& fType, int priority = 0);
    RulesDispatch(const std::set<FractalFuncType> fTypes, int priority = 0);

    static FractalRulesPtr CreateNew(const FractalConstructorData& data);

protected:
    virtual FractalRulesPtr MakeRules(const FractalConstructorData& data) const = 0;

    const int mPriority;

    static std::map<FractalFuncType, const RulesDispatch *>& DispatchMap()
    {
        static std::map<FractalFuncType, const RulesDispatch *> theMap;
        return theMap;
    }
};

static inline const std::set<FractalFuncType> CreateFuncTypeSet(unsigned n, ...)
{
    std::set<FractalFuncType> s;
    int i;
    va_list typeList;
    FractalFuncType f;
    va_start(typeList, n);
    for (i = 0; i < n; i++)
    {
        f = va_arg(typeList, const FractalFuncType);
        s.insert(f);
    }
    va_end(typeList);
    return s;
}

void InitDispatch();

}

#endif
