//******************************************************************************
///
/// @file core/shape/fractal/quaternion.cpp
///
/// This module implements Quaternion algebra julia fractals.
///
/// @author Pascal Massimino (original code)
/// @author Tim Wegner (revisions and updates for POV-Ray 3.x)
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
#include "core/shape/fractal/quaternion.h"

#include "core/math/complexfn.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/dispatch.h"
#include "core/shape/fractal/distestimator.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

#define CalcGenDerivs(n1,n2,n3,n4,normFVal,ny,nz,nw,tmp2)               \
    tmpx = (n1) * tmp2.x - tmp2.y *                                     \
        ((n2) * ny + (n3) * nz + (n4) * nw);				\
    tmpy = (n2) * (1.0 - ny * ny) * normFVal + ny *                     \
        ((n1) * tmp2.y + (n2) * ny * tmp2.x +				\
         (n3) * nz * (tmp2.x - normFVal) +                              \
         (n4) * nw * (tmp2.x - normFVal));                              \
    tmpz = (n3) * (1.0 - nz * nz) * normFVal + nz *                     \
        ((n1) * tmp2.y + (n3) * nz * tmp2.x +				\
         (n2) * ny * (tmp2.x - normFVal) +                              \
         (n4) * nw * (tmp2.x - normFVal));                              \
    tmpw = (n4) * (1.0 - nw * nw) * normFVal + nw *                     \
        ((n1) * tmp2.y + (n4) * nw * tmp2.x +				\
         (n2) * ny * (tmp2.x - normFVal) +                              \
         (n3) * nz * (tmp2.x - normFVal));                              \
                                                                        \
    (n1) = tmpx; (n2) = tmpy; (n3) = tmpz; (n4) = tmpw

#define AltGenDerivs(n1,n2,n3,n4,tmp2)                                  \
    (n1) *= tmp2.x; (n2) *= tmp2.x; (n3) *= tmp2.x; (n4) *= tmp2.x

#define ComponentComplexMult(n1,n2,c0)          \
    tmpx = (n1) * (c0).x - (n2) * (c0).y;       \
    (n2) = (n2) * (c0).x + (n1) * (c0).y;       \
    (n1) = tmpx

void QuaternionSqrFractalRules::
IterateCalc(VECTOR_4D &rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    DBL tmp;

    QuaternionSqrFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionSqrFractalRules::AuxIterData *>(pIterData->auxIter.data());

    pAuxIterStack[iter].sNorm = norm;

    tmp = 2.0 * rV[X];

    rV[X] = tmp * rV[X] + mJuliaParm[X] - norm;
    rV[Y] = tmp * rV[Y] + mJuliaParm[Y];
    rV[Z] = tmp * rV[Z] + mJuliaParm[Z];
    rV[W] = tmp * rV[W] + mJuliaParm[W];
}

void QuaternionSqrFractalRules::
DirDerivCalc(VECTOR_4D &rD, const VECTOR_4D& v, int iter, bool samePoint, const Fractal *pFractal, FractalIterData *pIterData) const
{
    DBL tmp;

    tmp = rD[X] * v[X] - rD[Y] * v[Y] - rD[Z] * v[Z] - rD[W] * v[W];

    rD[Y] = rD[X] * v[Y] + v[X] * rD[Y];
    rD[Z] = rD[X] * v[Z] + v[X] * rD[Z];
    rD[W] = rD[X] * v[W] + v[X] * rD[W];
    // Order matters here!
    rD[X] = tmp;
}

void QuaternionCubeFractalRules::
IterateCalc(VECTOR_4D& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    DBL cVal;

    QuaternionCubeFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionCubeFractalRules::AuxIterData *>(pIterData->auxIter.data());

    pAuxIterStack[iter].sNorm = norm;

    pAuxIterStack[iter].cVal = cVal = 4.0 * rV[X] * rV[X] - norm;

    rV[X] = (cVal - 2.0 * norm) * rV[X] + mJuliaParm[X];
    rV[Y] = cVal * rV[Y] + mJuliaParm[Y];
    rV[Z] = cVal * rV[Z] + mJuliaParm[Z];
    rV[W] = cVal * rV[W] + mJuliaParm[W];
}

void QuaternionCubeFractalRules::
DirDerivCalc(VECTOR_4D &rD, const VECTOR_4D& v, int iter, bool samePoint, const Fractal *pFractal, FractalIterData *pIterData) const
{
    static DBL tmp1, cVal, norm;
    DBL tmp2;

    if (!samePoint)
    {
        QuaternionCubeFractalRules::AuxIterData *pAuxIterStack =
            static_cast<QuaternionCubeFractalRules::AuxIterData *>(pIterData->auxIter.data());

        norm = pAuxIterStack[iter].sNorm;
        cVal = pAuxIterStack[iter].cVal;
        tmp1 = 2.0 * norm + cVal;
    }

    tmp2 = 2.0 * (3.0 * rD[X] * v[X] - rD[Y] * v[Y] - rD[Z] * v[Z] - rD[W] * v[W]);

    rD[X] = 3.0 * (-rD[X] * tmp1 + v[X] * tmp2);
    rD[Y] = rD[Y] * cVal + v[Y] * tmp2;
    rD[Z] = rD[Z] * cVal + v[Z] * tmp2;
    rD[W] = rD[W] * cVal + v[W] * tmp2;
}

void QuaternionRecipFractalRules::
IterateCalc(VECTOR_4D &rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    QuaternionRecipFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionRecipFractalRules::AuxIterData *>(pIterData->auxIter.data());

    pAuxIterStack[iter].sNorm = norm;

    rV[X] = rV[X] / norm + mJuliaParm[X];
    rV[Y] = -rV[Y] / norm + mJuliaParm[Y];
    rV[Z] = -rV[Z] / norm + mJuliaParm[Z];
    rV[W] = -rV[W] / norm + mJuliaParm[W];
}

void QuaternionRecipFractalRules::
DirDerivCalc(VECTOR_4D &rD, const VECTOR_4D& v, int iter, bool samePoint, const Fractal *pFractal, FractalIterData *pIterData) const
{
    static DBL norm;
    DBL tmp;

    if (!samePoint)
    {
        QuaternionRecipFractalRules::AuxIterData *pAuxIterStack =
            static_cast<QuaternionRecipFractalRules::AuxIterData *>(pIterData->auxIter.data());

        norm = pAuxIterStack[iter].sNorm;
    }

    tmp = 2.0 * (rD[X] * v[X] + rD[Y] * v[Y] + rD[Z] * v[Z] + rD[W] * v[W]) / norm;

    rD[X] = (rD[X] - tmp * v[X]) / norm;
    rD[Y] = (tmp * v[Y] - rD[Y]) / norm;
    rD[Z] = (tmp * v[Z] - rD[Z]) / norm;
    rD[W] = (tmp * v[W] - rD[W]) / norm;
}

void QuaternionFuncFractalRules::
IterateCalc(VECTOR_4D& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Complex tmp1, tmp2;
    DBL nNorm, normFVal;

    QuaternionFuncFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionFuncFractalRules::AuxIterData *>(pIterData->auxIter.data());

    pAuxIterStack[iter].nNorm = nNorm = sqrt(rV[Y] * rV[Y] + rV[Z] * rV[Z] + rV[W] * rV[W]);

    tmp1.x = rV[X];
    tmp1.y = nNorm;

    (*(mFunc.pFunc))(tmp2, tmp1, mExponent);

    if (nNorm == 0.0)
    {
        pAuxIterStack[iter].normFVal = normFVal = tmp2.y;
        if (normFVal != 0.0)
            return;
    }
    else
    {
        pAuxIterStack[iter].normFVal = normFVal = tmp2.y / nNorm;
    }

    rV[X] = tmp2.x + mJuliaParm[X];
    rV[Y] = rV[Y] * normFVal + mJuliaParm[Y];
    rV[Z] = rV[Z] * normFVal + mJuliaParm[Z];
    rV[W] = rV[W] * normFVal + mJuliaParm[W];
}

void QuaternionFuncFractalRules::
DirDerivCalc(VECTOR_4D& rD, const VECTOR_4D& v, int iter, bool samePoint, const Fractal *pFractal, FractalIterData *pIterData) const
{
    static DBL ny, nz, nw, nNorm, normFVal;
    static Complex tmp2;
    DBL tmpx, tmpy, tmpz, tmpw;

    if (!samePoint)
    {
        QuaternionFuncFractalRules::AuxIterData *pAuxIterStack =
            static_cast<QuaternionFuncFractalRules::AuxIterData *>(pIterData->auxIter.data());

        Complex tmp1;

        nNorm = pAuxIterStack[iter].nNorm;
        normFVal = pAuxIterStack[iter].normFVal;

        if (nNorm == 0.0 && normFVal != 0.0)
            return;

        tmp1.x = v[X];
        tmp1.y = nNorm;

        (*(mFunc.pDeriv))(tmp2, tmp1, mExponent);

        if (nNorm != 0.0)
        {
            ny = v[Y] / nNorm;
            nz = v[Z] / nNorm;
            nw = v[W] / nNorm;
        }
    }
    else if (nNorm == 0.0 && normFVal != 0.0)
        return;

    if (nNorm == 0.0)
    {
        AltGenDerivs(rD[X], rD[Y], rD[Z], rD[W], tmp2);
    }
    else
    {
        CalcGenDerivs(rD[X], rD[Y], rD[Z], rD[W], normFVal, ny, nz, nw, tmp2);
    }
}

bool QuaternionFuncFractalRules::
DiscontinuityCheck(VECTOR_4D& rD, DBL& rDist, const VECTOR_4D& t, const VECTOR_4D& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    Complex tmp, tPt, pPt;
    DBL dist, nc;

    QuaternionFuncFractalRules::AuxIterData *pTAuxIterStack =
        static_cast<QuaternionFuncFractalRules::AuxIterData *>(pTIterData->auxIter.data()),
        *pPAuxIterStack = static_cast<QuaternionFuncFractalRules::AuxIterData *>(pPIterData->auxIter.data());

    tPt.x = t[X];
    tPt.y = pTAuxIterStack[iter].nNorm;

    pPt.x = p[X];
    pPt.y = pPAuxIterStack[iter].nNorm;

    if ((*(mFunc.pDisc))(tmp, dist, tPt, pPt, mExponent))
    {
        rD[X] = tmp.x;
        if (pTAuxIterStack[iter].nNorm == 0.0)
        {
            rD[Y] = tmp.y;
            rD[Z] = rD[W] = 0.0;
        }
        else
        {
            nc = tmp.y / pTAuxIterStack[iter].nNorm;
            rD[Y] = t[Y] * nc;
            rD[Z] = t[Z] * nc;
            rD[W] = t[W] * nc;
        }
        rDist = dist;
        return true;
    }
    else
        return false;
}

void QuaternionPwrFractalRules::
IterateCalc(VECTOR_4D& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Complex tmp, lg0, lg1;
    DBL nNorm1, nNorm2, normFVal1, normFVal2;

    QuaternionPwrFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionPwrFractalRules::AuxIterData *>(pIterData->auxIter.data());

    pAuxIterStack[iter].nNorm[0] = nNorm1 = sqrt(rV[Y] * rV[Y] + rV[Z] * rV[Z] + rV[W] * rV[W]);

    tmp.x = rV[X];
    tmp.y = nNorm1;

    complex_fn::Ln(tmp, tmp);

    if (nNorm1 == 0.0)
    {
        pAuxIterStack[iter].normFVal[0] = normFVal1 = tmp.y;
    }
    else
    {
        pAuxIterStack[iter].normFVal[0] = normFVal1 = tmp.y / nNorm1;
    }

    lg0.x = tmp.x;
    lg0.y = rV[Y] * normFVal1;

    lg1.x = rV[Z] * normFVal1;
    lg1.y = rV[W] * normFVal1;

    complex_fn::Mult(lg0, lg0, mExponent);
    if (mInfo.funcType.variant == kVar_Left)
    {
        complex_fn::Mult(lg1, lg1, mExponent);
    }
    else
    {
        complex_fn::Mult(lg1, lg1, mExponentConj);
    }

    pAuxIterStack[iter].lg[0] = lg0;
    pAuxIterStack[iter].lg[1] = lg1;

    pAuxIterStack[iter].nNorm[1] = nNorm2 = sqrt(lg0.y * lg0.y + lg1.x * lg1.x + lg1.y * lg1.y);

    tmp.x = lg0.x;
    tmp.y = nNorm2;

    complex_fn::Exp(tmp, tmp);

    pAuxIterStack[iter].expVal = tmp;

    if (nNorm2 == 0.0)
    {
        pAuxIterStack[iter].normFVal[1] = normFVal2 = tmp.y;
    }
    else
    {
        pAuxIterStack[iter].normFVal[1] = normFVal2 = tmp.y / nNorm2;
    }

    rV[X] = tmp.x + mJuliaParm[X];
    rV[Y] = lg0.y * normFVal2 + mJuliaParm[Y];
    rV[Z] = lg1.x * normFVal2 + mJuliaParm[Z];
    rV[W] = lg1.y * normFVal2 + mJuliaParm[W];
}

void QuaternionPwrFractalRules::
DirDerivCalc(VECTOR_4D& rD, const VECTOR_4D& v, int iter, bool samePoint, const Fractal *pFractal, FractalIterData *pIterData) const
{
    static DBL ny1, ny2, nz1, nz2, nw1, nw2, nNorm1, nNorm2, normFVal1, normFVal2;
    static Complex tmp21, tmp22;
    DBL tmpx, tmpy, tmpz, tmpw;

    if (!samePoint)
    {
        Complex tmp1;

        QuaternionPwrFractalRules::AuxIterData *pAuxIterStack =
            static_cast<QuaternionPwrFractalRules::AuxIterData *>(pIterData->auxIter.data());

        nNorm1 = pAuxIterStack[iter].nNorm[0];
        normFVal1 = pAuxIterStack[iter].normFVal[0];

        tmp1.x = v[X];
        tmp1.y = nNorm1;

        complex_fn::Recip(tmp21, tmp1);

        nNorm2 = pAuxIterStack[iter].nNorm[1];
        normFVal2 = pAuxIterStack[iter].normFVal[1];

        tmp22 = pAuxIterStack[iter].expVal;

        if (nNorm1 != 0.0)
        {
            ny1 = v[Y] / nNorm1;
            nz1 = v[Z] / nNorm1;
            nw1 = v[W] / nNorm1;
        }

        if (nNorm2 != 0.0)
        {
            ny2 = pAuxIterStack[iter].lg[0].y / nNorm2;
            nz2 = pAuxIterStack[iter].lg[1].x / nNorm2;
            nw2 = pAuxIterStack[iter].lg[1].y / nNorm2;
        }
    }

    if (nNorm1 == 0.0)
    {
        AltGenDerivs(rD[X], rD[Y], rD[Z], rD[W], tmp21);
    }
    else
    {
        CalcGenDerivs(rD[X], rD[Y], rD[Z], rD[W], normFVal1, ny1, nz1, nw1, tmp21);
    }

    ComponentComplexMult(rD[X], rD[Y], mExponent);
    if (mInfo.funcType.variant == kVar_Left)
    {
        ComponentComplexMult(rD[Z], rD[W], mExponent);
    }
    else
    {
        ComponentComplexMult(rD[Z], rD[W], mExponentConj);
    }

    if (nNorm2 == 0.0)
    {
        AltGenDerivs(rD[X], rD[Y], rD[Z], rD[W], tmp22);
    }
    else
    {
        CalcGenDerivs(rD[X], rD[Y], rD[Z], rD[W], normFVal2, ny2, nz2, nw2, tmp22);
    }
}

bool QuaternionPwrFractalRules::
DiscontinuityCheck(VECTOR_4D& rD, DBL& rDist, const VECTOR_4D& t, const VECTOR_4D& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    Complex tmp, tPt, pPt;
    DBL dist, nc;

    QuaternionPwrFractalRules::AuxIterData *pTAuxIterStack =
        static_cast<QuaternionPwrFractalRules::AuxIterData *>(pTIterData->auxIter.data()),
        *pPAuxIterStack = static_cast<QuaternionPwrFractalRules::AuxIterData *>(pPIterData->auxIter.data());

    tPt.x = t[X];
    tPt.y = pTAuxIterStack[iter].nNorm[0];

    pPt.x = p[X];
    pPt.y = pPAuxIterStack[iter].nNorm[0];

    if (complex_fn::NegReal_DTest(tmp, dist, tPt, pPt, mExponent))
    {
        rD[X] = tmp.x;
        if (pTAuxIterStack[iter].nNorm[0] == 0.0)
        {
            rD[Y] = tmp.y;
            rD[Z] = rD[W] = 0.0;
        }
        else
        {
            nc = tmp.y / pTAuxIterStack[iter].nNorm[0];
            rD[Y] = t[Y] * nc;
            rD[Z] = t[Z] * nc;
            rD[W] = t[W] * nc;
        }
        rDist = dist;
        return true;
    }
    else
        return false;
}

namespace estimators
{


DBL EstimatorSpecialOrig_QuatSqr(const FractalRules *rules, DBL norm, int iters, const Vector3d& direction,
                                 const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL tmp, nProd, pow;
    int j;

    QuaternionSqrFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionSqrFractalRules::AuxIterData *>(pIterData->auxIter.data());

    tmp = dot(pFractal->SliceNorm, direction);

    nProd = 1.0 + tmp * tmp;

    pow = 1.0 / 2.0;

    for (j = 0; j < iters; ++j)
    {
        nProd *= pAuxIterStack[j].sNorm;
        pow /= 2.0;
    }

    return pow / sqrt(nProd) * log(norm);
}

DBL EstimatorSpecialOrig_QuatCube(const FractalRules *rules, DBL norm, int iters, const Vector3d& direction,
                                  const Fractal *pFractal, FractalIterData *pIterData)
{
    DBL tmp, nProd, pow;
    int j;

    QuaternionCubeFractalRules::AuxIterData *pAuxIterStack =
        static_cast<QuaternionCubeFractalRules::AuxIterData *>(pIterData->auxIter.data());

    tmp = dot(pFractal->SliceNorm, direction);

    nProd = 1.0 + tmp * tmp;

    pow = 1.0 / 3.0;

    for (j = 0; j < iters; ++j)
    {
        nProd *= pAuxIterStack[j].sNorm;
        pow /= 3.0;
    }

    return pow / sqrt(nProd) * log(norm);
}

const DistanceEstimator kSpecialOrig_QuatSqr = { EstimatorSpecialOrig_QuatSqr, kOrigSpecialEstimators };
const DistanceEstimator kSpecialOrig_QuatCube = { EstimatorSpecialOrig_QuatCube, kOrigSpecialEstimators };

}

const DistanceEstimator& QuaternionSqrFractalRules::
ExtraEstimators(EstimatorType tgtType)
{
    switch (tgtType)
    {
    case kOrigSpecialEstimators:
        return estimators::kSpecialOrig_QuatSqr;
    default:
        return estimators::BadEstimator();
    }
}

const DistanceEstimator& QuaternionCubeFractalRules::
ExtraEstimators(EstimatorType tgtType)
{
    switch (tgtType)
    {
    case kOrigSpecialEstimators:
        return estimators::kSpecialOrig_QuatCube;
    default:
        return estimators::BadEstimator();
    }
}

void QuaternionDispatchInit() {
    static const RulesDispatch QuatSqrDispatch(MakeCreatorFunc<QuaternionSqrFractalRules>,
                                               CreateFuncType(kQuaternion, kFunc_Sqr, kVar_Normal));

    static const RulesDispatch QuatCubeDispatch(MakeCreatorFunc<QuaternionCubeFractalRules>,
                                                CreateFuncType(kQuaternion, kFunc_Cube, kVar_Normal));

    static const RulesDispatch QuatRecipDispatch(MakeCreatorFunc<QuaternionRecipFractalRules>,
                                                 CreateFuncType(kQuaternion, kFunc_Reciprocal, kVar_Normal));

    static const RulesDispatch QuatFuncDispatch(MakeCreatorFunc<QuaternionFuncFractalRules>,
                                                Func_FuncTypeSet(kQuaternion), -1);

    static const RulesDispatch QuatPwrDispatch(MakeCreatorFunc<QuaternionPwrFractalRules>,
                                               CreateSet<2>(CreateFuncType(kQuaternion, kFunc_Pwr, kVar_Left),
                                                            CreateFuncType(kQuaternion, kFunc_Pwr, kVar_Right)));
}

}
