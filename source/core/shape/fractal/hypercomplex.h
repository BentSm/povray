//******************************************************************************
///
/// @file core/shape/fractal/hypercomplex.h
///
/// This module contains all defines, typedefs, and prototypes for `hypercomplex.cpp`.
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

#ifndef POVRAY_CORE_FRACTAL_HYPERCOMPLEX_H
#define POVRAY_CORE_FRACTAL_HYPERCOMPLEX_H

#include "core/coretypes.h"
#include "core/math/vector.h"
#include "core/shape/fractal/cxx03compat.h"
#include "core/shape/fractal/func.h"
#include "core/shape/fractal/magic.h"
#include "core/shape/fractal/types.h"

namespace pov
{

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

/*****************************************************************************
* Global typedefs
******************************************************************************/

/*****************************************************************************
* Global variables
******************************************************************************/

/*****************************************************************************
* Global functions
******************************************************************************/

template <template <class> class RulesClass, class Estimator, class BaseRules = FractalRules>
class HypercomplexFractalRulesBase : public MagicHypercomplexFractalRules<RulesClass, Estimator, BaseRules>
{
public:
    typedef struct {
        Complex point[2];
    } IterationData;

    HypercomplexFractalRulesBase(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport) :
        MagicHypercomplexFractalRules<RulesClass, Estimator, BaseRules>(data, discontinuitySupport) {}
};

template <class Estimator>
class HypercomplexSqrFractalRules : public HypercomplexFractalRulesBase<HypercomplexSqrFractalRules, Estimator>
{
public:
    HypercomplexSqrFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase<pov::HypercomplexSqrFractalRules, Estimator>(data, kDiscontinuityUnneeded) {}
    inline void IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
                               int iter, const Fractal *pFractal, void *pIterData) const;

protected:
    using MagicHypercomplexFractalRules<pov::HypercomplexSqrFractalRules, Estimator>::mDuplexJuliaParm;
};

template <class Estimator>
class HypercomplexCubeFractalRules : public HypercomplexFractalRulesBase<HypercomplexCubeFractalRules, Estimator>
{
public:
    HypercomplexCubeFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase<pov::HypercomplexCubeFractalRules, Estimator>(data, kDiscontinuityUnneeded) {}
    inline void IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
                               int iter, const Fractal *pFractal, void *pIterData) const;

protected:
    using MagicHypercomplexFractalRules<pov::HypercomplexCubeFractalRules, Estimator>::mDuplexJuliaParm;
};

template <class Estimator>
class HypercomplexRecipFractalRules : public HypercomplexFractalRulesBase<HypercomplexRecipFractalRules, Estimator>
{
public:
    HypercomplexRecipFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase<pov::HypercomplexRecipFractalRules, Estimator>(data, kDiscontinuityUnneeded) {}
    inline void IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
                               int iter, const Fractal *pFractal, void *pIterData) const;

protected:
    using MagicHypercomplexFractalRules<pov::HypercomplexRecipFractalRules, Estimator>::mDuplexJuliaParm;
};

template <class Estimator>
class HypercomplexFuncFractalRules : public HypercomplexFractalRulesBase<HypercomplexFuncFractalRules, Estimator>
{
public:
    HypercomplexFuncFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase<pov::HypercomplexFuncFractalRules,
                                     Estimator>(data, DiscontinuitySupport_Func(FractalFuncForType(data.funcType))),
        mFunc(FractalFuncForType(data.funcType)), mExponent(data.exponent) {}
    inline void IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
                               int iter, const Fractal *pFractal, void *pIterData) const;

    bool DiscontinuityCheck(Complex& rD0, Complex& rD1, DBL& rDist,
                            const Complex& t0, const Complex& t1, const Complex& p0, const Complex& p1,
                            int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const;

protected:
    const FractalFunc mFunc;
    const Complex mExponent;
    using MagicHypercomplexFractalRules<pov::HypercomplexFuncFractalRules, Estimator>::mDuplexJuliaParm;
};

void HypercomplexDispatchInit();

}

#endif // POVRAY_CORE_FRACTAL_HYPERCOMPLEX_H
