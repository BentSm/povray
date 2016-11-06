//******************************************************************************
///
/// @file core/shape/fractal/dispatch.cpp
///
/// This module is the main dispatch class.  Simple, but important.
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
#include "core/shape/fractal/dispatch.h"

#ifdef DISPATCH_DEBUG
#include <iostream>
#endif

#include "core/shape/fractal/hypercomplex.h"
#include "core/shape/fractal/quaternion.h"

#include "base/pov_err.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

static inline void DumpDispatchMap(std::map<FractalFuncType, const RulesDispatch *> dMap)
{
#ifdef DISPATCH_DEBUG
    std::map<FractalFuncType, const RulesDispatch *>::const_iterator dispatchIter;
    std::cerr << "Start of dispatch map" << std::endl;
    for (dispatchIter = dMap.begin(); dispatchIter != dMap.end(); dispatchIter++)
        std::cerr << (int)dispatchIter->first.algebra << " " << (int)dispatchIter->first.type
                  << " " << (int)dispatchIter->first.variant << std::endl;
    std::cerr << "End of dispatch map" << std::endl;
#endif
}

RulesDispatch::RulesDispatch(CreatorFunc *func, const FractalFuncType& funcType, int priority) : mCreatorFunc(func), mPriority(priority)
{
    std::map<FractalFuncType, const RulesDispatch *>::const_iterator dispatchIter = DispatchMap().find(funcType);
    if (dispatchIter == DispatchMap().end() || dispatchIter->second->mPriority <= priority)
        DispatchMap()[funcType] = this;
    DumpDispatchMap(DispatchMap());
}

RulesDispatch::RulesDispatch(CreatorFunc *func, const std::set<FractalFuncType> funcTypes, int priority) : mCreatorFunc(func), mPriority(priority)
{
    std::set<FractalFuncType>::const_iterator c = funcTypes.begin();
    std::map<FractalFuncType, const RulesDispatch *>::const_iterator dispatchIter;
    while (c != funcTypes.end())
    {
        dispatchIter = DispatchMap().find(*c);
        if (dispatchIter == DispatchMap().end() || dispatchIter->second->mPriority <= priority)
            DispatchMap()[*c] = this;
        c++;
    }
    DumpDispatchMap(DispatchMap());
}

FractalRulesPtr RulesDispatch::CreateNew(const FractalConstructorData& data)
{
    std::map<FractalFuncType, const RulesDispatch *>::const_iterator dispatchIter;
    dispatchIter = DispatchMap().find(data.funcType);
    if (dispatchIter == DispatchMap().end())
    {
        DumpDispatchMap(DispatchMap());
        throw POV_EXCEPTION_STRING("Algebra/function type/variant combination unknown in fractal.");
    }
    else
        return (*(dispatchIter->second->mCreatorFunc))(data);
}

void InitDispatch() {
    HypercomplexDispatchInit();
    QuaternionDispatchInit();
}

}
