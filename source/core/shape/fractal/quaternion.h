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

class QuaternionSqrFractalRules : public MagicQuaternionFractalRules
{
public:
    typedef struct {
        DBL sNorm;
    } AuxIterData;

    QuaternionSqrFractalRules(const FractalConstructorData& data) :
        MagicQuaternionFractalRules(data, kDiscontinuityUnneeded,
                                    GetDataSizes<QuaternionSqrFractalRules>(),
                                    GetEstimatorFromType(data.estimatorType, kOrigSpecialEstimators,
                                                         kOrigSpecialEstimators, QuaternionSqrFractalRules::ExtraEstimators)) { }

    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& v, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;

protected:
    static const DistanceEstimator& ExtraEstimators(EstimatorType tgtType);
};

class QuaternionCubeFractalRules : public MagicQuaternionFractalRules
{
public:
    typedef struct {
        DBL sNorm, cVal;
    } AuxIterData;

    QuaternionCubeFractalRules(const FractalConstructorData& data) :
        MagicQuaternionFractalRules(data, kDiscontinuityUnneeded,
                                    GetDataSizes<QuaternionCubeFractalRules>(),
                                    GetEstimatorFromType(data.estimatorType, kOrigSpecialEstimators,
                                                         kOrigSpecialEstimators, QuaternionCubeFractalRules::ExtraEstimators)) { }

    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& v, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;

protected:
    static const DistanceEstimator& ExtraEstimators(EstimatorType tgtType);
};

class QuaternionRecipFractalRules : public MagicQuaternionFractalRules
{
public:
    typedef struct {
        DBL sNorm;
    } AuxIterData;

    QuaternionRecipFractalRules(const FractalConstructorData& data) :
        MagicQuaternionFractalRules(data, kDiscontinuityUnneeded,
                                    GetDataSizes<QuaternionRecipFractalRules>(),
                                    GetEstimatorFromType(data.estimatorType)) { }

    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& v, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;
};

class QuaternionFuncFractalRules : public MagicQuaternionFractalRules
{
public:
    typedef struct {
        Vector4d point;
        DBL nNorm, normFVal;
    } AuxIterData;

    QuaternionFuncFractalRules(const FractalConstructorData& data) :
        MagicQuaternionFractalRules(data, DiscontinuitySupport_Func(FractalFuncForType(data.funcType)),
                                    GetDataSizes<QuaternionFuncFractalRules>(),
                                    GetEstimatorFromType(data.estimatorType)),
        mFunc(FractalFuncForType(data.funcType)), mExponent(data.exponent) { }

    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& v, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;

    virtual bool DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;

protected:
    const FractalFunc mFunc;
    const IComplex mExponent;
};

class QuaternionPwrFractalRules : public MagicQuaternionFractalRules
{
public:
    typedef struct {
        DBL nNorm[2], normFVal[2];
        Complex expVal, lg[2];
    } AuxIterData;

    QuaternionPwrFractalRules(const FractalConstructorData& data) :
        MagicQuaternionFractalRules(data, kDiscontinuitySupported,
                                    GetDataSizes<QuaternionPwrFractalRules>(),
                                    GetEstimatorFromType(data.estimatorType)),
        mExponent(data.exponent), mExponentConj(data.exponent[X], -data.exponent[Y]) { }

    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void GradientCalc(Vector4d& rD, const Vector4d& v, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const;

    virtual bool DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;

protected:
    const IComplex mExponent, mExponentConj;
};

void QuaternionDispatchInit();

namespace estimators
{

extern const DistanceEstimator kSpecialOrig_QuatSqr;
extern const DistanceEstimator kSpecialOrig_QuatCube;

}

}

#endif // POVRAY_CORE_FRACTAL_QUATERNION_H
