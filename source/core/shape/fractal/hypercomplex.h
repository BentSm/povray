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

#ifndef POVRAY_CORE_FRACTAL_HYPERCOMPLEX_H
#define POVRAY_CORE_FRACTAL_HYPERCOMPLEX_H

#include "core/coretypes.h"
#include "core/math/vector.h"
#include "core/shape/fractal/func.h"
#include "core/shape/fractal/magic.h"
#include "core/shape/fractal/types.h"
#include "core/shape/fractal/util.h"

namespace pov
{


class HypercomplexFractalRulesBase : public MagicFractalRules
{
public:
    HypercomplexFractalRulesBase(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                                 const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        mDuplexJuliaParm(DuplexFromHypercomplex(data.juliaParm)),
        MagicFractalRules(data, discontinuitySupport, sizes, estimator) { }

protected:
    const Vector4d mDuplexJuliaParm;

};

class HypercomplexSqrFractalRules : public HypercomplexFractalRulesBase
{
public:
    HypercomplexSqrFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase(data, kDiscontinuityUnneeded,
                                     GetDataSizes<HypercomplexSqrFractalRules>(),
                                     GetEstimatorFromType(data.estimatorType)) { }

    virtual void IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;
};

class HypercomplexCubeFractalRules : public HypercomplexFractalRulesBase
{
public:
    HypercomplexCubeFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase(data, kDiscontinuityUnneeded,
                                     GetDataSizes<HypercomplexCubeFractalRules>(),
                                     GetEstimatorFromType(data.estimatorType)) { }

    virtual void IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;
};

class HypercomplexRecipFractalRules : public HypercomplexFractalRulesBase
{
public:
    HypercomplexRecipFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase(data, kDiscontinuityUnneeded,
                                     GetDataSizes<HypercomplexRecipFractalRules>(),
                                     GetEstimatorFromType(data.estimatorType)) { }

    virtual void IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;
};

class HypercomplexFuncFractalRules : public HypercomplexFractalRulesBase
{
public:
    HypercomplexFuncFractalRules(const FractalConstructorData& data) :
        HypercomplexFractalRulesBase(data, DiscontinuitySupport_Func(FractalFuncForType(data.funcType)),
                                     GetDataSizes<HypercomplexFuncFractalRules>(),
                                     GetEstimatorFromType(data.estimatorType)),
        mFunc(FractalFuncForType(data.funcType)), mExponent(data.exponent) { }

    virtual void IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;

    virtual bool DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;

protected:
    const FractalFunc mFunc;
    const IComplex mExponent;
};

void HypercomplexDispatchInit();

}

#endif // POVRAY_CORE_FRACTAL_HYPERCOMPLEX_H
