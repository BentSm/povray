//******************************************************************************
///
/// @file core/shape/fractal/quaternion.h
///
/// This module contains all defines, typedefs, and prototypes for `quaternion.cpp`.
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

#ifndef POVRAY_CORE_FRACTAL_QUATERNION_H
#define POVRAY_CORE_FRACTAL_QUATERNION_H

#include "core/coretypes.h"
#include "core/math/vector.h"
#include "core/shape/fractal/func.h"
#include "core/shape/fractal/magic.h"
#include "core/shape/fractal/types.h"
#include "core/shape/fractal/util.h"

namespace pov
{

struct BasicRay;
class Fractal;

template <template <class> class RulesClass, class Estimator, class BaseRules = FractalRules>
class QuaternionFractalRulesBase : public MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>
{
public:
    typedef struct {
        Quaternion point;
        DBL sNorm;
    } IterationData;

    QuaternionFractalRulesBase(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport) :
        MagicQuaternionFractalRules<RulesClass, Estimator, BaseRules>(data, discontinuitySupport) {}
};

template <class Estimator>
class QuaternionSqrFractalRules : public QuaternionFractalRulesBase<QuaternionSqrFractalRules, Estimator>
{
public:
    QuaternionSqrFractalRules(const FractalConstructorData& data) :
        QuaternionFractalRulesBase<pov::QuaternionSqrFractalRules, Estimator>(data, kDiscontinuityUnneeded) {}
    inline void IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const;

protected:
    using MagicQuaternionFractalRules<pov::QuaternionSqrFractalRules, Estimator>::mJuliaParm;
};

template <class Estimator>
class QuaternionCubeFractalRules : public QuaternionFractalRulesBase<QuaternionCubeFractalRules, Estimator>
{
public:
    typedef struct {
        Quaternion point;
        DBL sNorm, cVal;
    } IterationData;

    QuaternionCubeFractalRules(const FractalConstructorData& data) :
        QuaternionFractalRulesBase<pov::QuaternionCubeFractalRules, Estimator>(data, kDiscontinuityUnneeded) {}
    inline void IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const;

protected:
    using MagicQuaternionFractalRules<pov::QuaternionCubeFractalRules, Estimator>::mJuliaParm;
};

template <class Estimator>
class QuaternionRecipFractalRules : public QuaternionFractalRulesBase<QuaternionRecipFractalRules, Estimator>
{
public:
    QuaternionRecipFractalRules(const FractalConstructorData& data) :
        QuaternionFractalRulesBase<pov::QuaternionRecipFractalRules, Estimator>(data, kDiscontinuityUnneeded) {}
    inline void IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const;

protected:
    using MagicQuaternionFractalRules<pov::QuaternionRecipFractalRules, Estimator>::mJuliaParm;
};

template <class Estimator>
class QuaternionFuncFractalRules : public QuaternionFractalRulesBase<QuaternionFuncFractalRules, Estimator>
{
public:
    typedef struct {
        Quaternion point;
        DBL nNorm, normFVal;
    } IterationData;

    QuaternionFuncFractalRules(const FractalConstructorData& data) :
        QuaternionFractalRulesBase<pov::QuaternionFuncFractalRules,
                                   Estimator>(data, DiscontinuitySupport_Func(FractalFuncForType(data.funcType))),
        mFunc(FractalFuncForType(data.funcType)), mExponent(data.exponent) {}
    inline void IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const;

    bool DiscontinuityCheck(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL& rDist,
                            DBL tx, DBL ty, DBL tz, DBL tw, DBL px, DBL py, DBL pz, DBL pw,
                            int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const;

protected:
    const FractalFunc mFunc;
    const Complex mExponent;
    using MagicQuaternionFractalRules<pov::QuaternionFuncFractalRules, Estimator>::mJuliaParm;
};

template <class Estimator>
class QuaternionPwrFractalRules : public QuaternionFractalRulesBase<QuaternionPwrFractalRules, Estimator>
{
public:
    typedef struct {
        Quaternion point;
        DBL nNorm[2], normFVal[2];
        Complex expVal, lg[2];
    } IterationData;

    QuaternionPwrFractalRules(const FractalConstructorData& data) :
        QuaternionFractalRulesBase<pov::QuaternionPwrFractalRules, Estimator>(data, kDiscontinuitySupported),
        mExponent(data.exponent), mExponentConj(CreateComplex(data.exponent.x, -data.exponent.y)) {}

    inline void IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const;
    inline void ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const;

    bool DiscontinuityCheck(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL& rDist,
                            DBL tx, DBL ty, DBL tz, DBL tw, DBL px, DBL py, DBL pz, DBL pw,
                            int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const;

protected:
    const Complex mExponent, mExponentConj;
    using MagicQuaternionFractalRules<pov::QuaternionPwrFractalRules, Estimator>::mJuliaParm;
};

void QuaternionDispatchInit();

}

#endif // POVRAY_CORE_FRACTAL_QUATERNION_H
