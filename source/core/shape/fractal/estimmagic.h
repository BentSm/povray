//******************************************************************************
///
/// @file core/shape/fractal/estimmagic.h
///
/// 'Magic' for dispatching to the estimator-specialized templates.  Also
/// contains the default estimator setup.
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

#ifndef POVRAY_CORE_FRACTAL_ESTIMMAGIC_H
#define POVRAY_CORE_FRACTAL_ESTIMMAGIC_H

#include "core/coretypes.h"

#include "core/math/vector.h"
#include "core/shape/fractal/types.h"
#include "core/shape/fractal/dispatch.h"

namespace pov
{

namespace magic
{

/* Additional estimator templates */

class ExtraEstim_None
{
public:
    template <template <class> class Rules>
    class Invoke
    {
    public:
        static FractalRulesPtr ExtraEstimators(const FractalConstructorData& data, EstimatorType actualType)
        {
            throw POV_EXCEPTION_STRING("Unsupported distance estimator for fractal type.");
        }
    };
};

template <EstimatorType eType, class EstimatorClass, class NextExtraEstim = ExtraEstim_None>
class ExtraEstim
{
public:
    template <template <class> class Rules>
    class Invoke
    {
    public:
        static FractalRulesPtr ExtraEstimators(const FractalConstructorData& data, EstimatorType actualType)
        {
            if (actualType == eType)
                return static_cast<FractalRulesPtr>(new Rules<EstimatorClass>(data));
            else
                return NextExtraEstim::template Invoke<Rules>::ExtraEstimators(data, actualType);
        }
    };
};

/* General estimator templates */

template <EstimatorType defaultEstim = kNewtonEstimator, EstimatorType legacyEstim = kOrigNewtonEstimator,
          class ExtraEstim = ExtraEstim_None>
class Estim_Standard
{
public:
    template <template <class> class Rules>
    class Invoke
    {
    public:
        static FractalRulesPtr CreateWithEstimator(const FractalConstructorData& data)
        {
            EstimatorType actualType = data.estimatorType;
            if (actualType == kDefaultEstimator)
                actualType = defaultEstim;
            else if (actualType == kLegacyEstimator)
                actualType = legacyEstim;
            switch (actualType)
            {
            case kNoEstimator:
                return static_cast<FractalRulesPtr>(new Rules<DistEstimatorNone>(data));
                break;
            case kNewtonEstimator:
                return static_cast<FractalRulesPtr>(new Rules<DistEstimatorNewton>(data));
                break;
            case kOrigNewtonEstimator:
                return static_cast<FractalRulesPtr>(new Rules<DistEstimatorNewtonOrig>(data));
                break;
            default:
                return ExtraEstim::template Invoke<Rules>::ExtraEstimators(data, actualType);
            }
        }
    };
};

template <class EstimatorClass>
class Estim_OrigSpecial
{
public:
    template <template <class> class Rules>
    class Invoke
    {
    public:
        static FractalRulesPtr CreateWithEstimator(const FractalConstructorData& data)
        {
            return Estim_Standard<kOrigSpecialEstimators, kOrigSpecialEstimators,
                ExtraEstim<kOrigSpecialEstimators, EstimatorClass> >::template Invoke<Rules>::
                CreateWithEstimator(data);
        }
    };
};

}

/* RulesDispatch subclass */

template <template <class> class Rules, class Estim = magic::Estim_Standard<> >
class MagicRulesDispatch : public RulesDispatch
{
public:
    MagicRulesDispatch(const FractalFuncType& fType, int priority = 0) : RulesDispatch(fType, priority) {}
    MagicRulesDispatch(const std::set<FractalFuncType> fTypes, int priority = 0) : RulesDispatch(fTypes, priority) {}

protected:
    virtual FractalRulesPtr MakeRules(const FractalConstructorData& data) const
    {
        return Estim::template Invoke<Rules>::CreateWithEstimator(data);
    }
};

}

#endif
