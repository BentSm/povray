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

// configcore.h must always be the first POV file included in core *.cpp files (pulls in platform config)
#include "core/configcore.h"
#include "core/shape/fractal/quaternion.h"

#include "core/math/complexfn.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/dispatch.h"
#include "core/shape/fractal/distestimator.h"
#include "core/shape/fractal/estimmagic.h"
#include "core/shape/fractal/magicimpl.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

#define CalcGenDerivs(n1,n2,n3,n4,normFVal,ny,nz,nw,tmp2)               \
    tmpx = (n1) * tmp2.real() - tmp2.imag() *                           \
        ((n2) * ny + (n3) * nz + (n4) * nw);				\
    tmpy = (n2) * (1.0 - ny * ny) * normFVal + ny *                     \
        ((n1) * tmp2.imag() + (n2) * ny * tmp2.real() +                 \
         (n3) * nz * (tmp2.real() - normFVal) +                         \
         (n4) * nw * (tmp2.real() - normFVal));                         \
    tmpz = (n3) * (1.0 - nz * nz) * normFVal + nz *                     \
        ((n1) * tmp2.imag() + (n3) * nz * tmp2.real() +                 \
         (n2) * ny * (tmp2.real() - normFVal) +                         \
         (n4) * nw * (tmp2.real() - normFVal));                         \
    tmpw = (n4) * (1.0 - nw * nw) * normFVal + nw *                     \
        ((n1) * tmp2.imag() + (n4) * nw * tmp2.real() +                 \
         (n2) * ny * (tmp2.real() - normFVal) +                         \
         (n3) * nz * (tmp2.real() - normFVal));                         \
                                                                        \
    (n1) = tmpx; (n2) = tmpy; (n3) = tmpz; (n4) = tmpw

#define AltGenDerivs(n1,n2,n3,n4,tmp2)                                  \
    (n1) *= tmp2.real(); (n2) *= tmp2.real();                           \
    (n3) *= tmp2.real(); (n4) *= tmp2.real()

#define ComponentComplexMult(n1,n2,c0)                  \
    tmpx = (n1) * (c0).real() - (n2) * (c0).imag();     \
    (n2) = (n2) * (c0).real() + (n1) * (c0).imag();     \
    (n1) = tmpx

template <class Estimator>
inline void QuaternionSqrFractalRules<Estimator>::
IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    DBL tmp;

    typename QuaternionSqrFractalRules::IterationData *pIterStack =
        reinterpret_cast<typename QuaternionSqrFractalRules::IterationData *>(pIterData);

    pIterStack[iter].sNorm = norm;

    tmp = 2.0 * rX;

    rY = tmp * rY + mJuliaParm[Y];
    rZ = tmp * rZ + mJuliaParm[Z];
    rW = tmp * rW + mJuliaParm[W];
    rX = tmp * rX + mJuliaParm[X] - norm;
}

template <class Estimator>
inline void QuaternionSqrFractalRules<Estimator>::
ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const
{
    DBL tmp;

    tmp = rDx * x - rDy * y - rDz * z - rDw * w;
    rDy = rDx * y + x * rDy;
    rDz = rDx * z + x * rDz;
    rDw = rDx * w + x * rDw;
    rDx = tmp;
}

template <class Estimator>
inline void QuaternionCubeFractalRules<Estimator>::
IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    DBL cVal;

    typename QuaternionCubeFractalRules::IterationData *pIterStack =
        reinterpret_cast<typename QuaternionCubeFractalRules::IterationData *>(pIterData);

    pIterStack[iter].sNorm = norm;

    pIterStack[iter].cVal = cVal = 4.0 * rX * rX - norm;

    rX = (cVal - 2.0 * norm) * rX + mJuliaParm[X];
    rY = cVal * rY + mJuliaParm[Y];
    rZ = cVal * rZ + mJuliaParm[Z];
    rW = cVal * rW + mJuliaParm[W];
}

template <class Estimator>
inline void QuaternionCubeFractalRules<Estimator>::
ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const
{
    static DBL tmp1, cVal, norm;
    DBL tmp2;

    if (!samePoint)
    {
        typename QuaternionCubeFractalRules::IterationData *pIterStack =
            reinterpret_cast<typename QuaternionCubeFractalRules::IterationData *>(pIterData);

        norm = pIterStack[iter].sNorm;
        cVal = pIterStack[iter].cVal;
        tmp1 = 2.0 * norm + cVal;
    }

    tmp2 = 2.0 * (3.0 * rDx * x - rDy * y - rDz * z - rDw * w);
    rDx = 3.0 * (-rDx * tmp1 + x * tmp2);
    rDy = rDy * cVal + y * tmp2;
    rDz = rDz * cVal + z * tmp2;
    rDw = rDw * cVal + w * tmp2;
}

template <class Estimator>
inline void QuaternionRecipFractalRules<Estimator>::
IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    DBL tmp;

    typename QuaternionRecipFractalRules::IterationData *pIterStack =
        reinterpret_cast<typename QuaternionRecipFractalRules::IterationData *>(pIterData);

    pIterStack[iter].sNorm = norm;

    rX = rX / norm + mJuliaParm[X];
    rY = -rY / norm + mJuliaParm[Y];
    rZ = -rZ / norm + mJuliaParm[Z];
    rW = -rW / norm + mJuliaParm[W];
}

template <class Estimator>
inline void QuaternionRecipFractalRules<Estimator>::
ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const
{
    static DBL norm;
    DBL tmp;

    if (!samePoint)
    {
        typename QuaternionRecipFractalRules::IterationData *pIterStack =
            reinterpret_cast<typename QuaternionRecipFractalRules::IterationData *>(pIterData);

        norm = pIterStack[iter].sNorm;
    }

    tmp = 2.0 * (rDx * x + rDy * y + rDz * z + rDw * w) / norm;
    rDx = (rDx - tmp * x) / norm;
    rDy = (tmp * y - rDy) / norm;
    rDz = (tmp * z - rDz) / norm;
    rDw = (tmp * w - rDw) / norm;
}

template <class Estimator>
inline void QuaternionFuncFractalRules<Estimator>::
IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    Complex tmp1, tmp2;
    DBL nNorm, normFVal;

    typename QuaternionFuncFractalRules::IterationData *pIterStack =
        reinterpret_cast<typename QuaternionFuncFractalRules::IterationData *>(pIterData);

    pIterStack[iter].nNorm = nNorm = sqrt(rY * rY + rZ * rZ + rW * rW);

    tmp1 = Complex(rX, nNorm);

    (*(mFunc.pFunc))(tmp2, tmp1, mExponent);

    if (nNorm == 0.0)
    {
        pIterStack[iter].normFVal = normFVal = tmp2.imag();
        if (normFVal != 0.0)
            return;
    }
    else
    {
        pIterStack[iter].normFVal = normFVal = tmp2.imag() / nNorm;
    }

    rX = tmp2.real() + mJuliaParm[X];
    rY = rY * normFVal + mJuliaParm[Y];
    rZ = rZ * normFVal + mJuliaParm[Z];
    rW = rW * normFVal + mJuliaParm[W];
}

template <class Estimator>
inline void QuaternionFuncFractalRules<Estimator>::
ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const
{
    static DBL ny, nz, nw, nNorm, normFVal;
    static Complex tmp2;
    DBL tmpx, tmpy, tmpz, tmpw;

    if (!samePoint)
    {
        typename QuaternionFuncFractalRules::IterationData *pIterStack =
            reinterpret_cast<typename QuaternionFuncFractalRules::IterationData *>(pIterData);

        Complex tmp1;

        nNorm = pIterStack[iter].nNorm;
        normFVal = pIterStack[iter].normFVal;

        if (nNorm == 0.0 && normFVal != 0.0)
            return;

        tmp1 = Complex(x, nNorm);

        (*(mFunc.pDeriv))(tmp2, tmp1, mExponent);

        if (nNorm != 0.0)
        {
            ny = y / nNorm;
            nz = z / nNorm;
            nw = w / nNorm;
        }
    }
    else if (nNorm == 0.0 && normFVal != 0.0)
        return;

    if (nNorm == 0.0)
    {
        AltGenDerivs(rDx, rDy, rDz, rDw, tmp2);
    }
    else
    {
        CalcGenDerivs(rDx, rDy, rDz, rDw, normFVal, ny, nz, nw, tmp2);
    }
}

template <class Estimator>
bool QuaternionFuncFractalRules<Estimator>::
DiscontinuityCheck(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL& rDist,
                   DBL tx, DBL ty, DBL tz, DBL tw, DBL px, DBL py, DBL pz, DBL pw,
                   int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    Complex tmp, tPt, pPt;
    DBL dist, nc;

    typename QuaternionFuncFractalRules::IterationData *pTIterStack =
        reinterpret_cast<typename QuaternionFuncFractalRules::IterationData *>(pTIterData),
        *pPIterStack = reinterpret_cast<typename QuaternionFuncFractalRules::IterationData *>(pPIterData);

    tPt = Complex(tx, pTIterStack[iter].nNorm);

    pPt = Complex(px, pPIterStack[iter].nNorm);

    if ((*(mFunc.pDisc))(tmp, dist, tPt, pPt, mExponent))
    {
        rDx = tmp.real();
        if (pTIterStack[iter].nNorm == 0.0)
        {
            rDy = tmp.imag();
            rDz = rDw = 0.0;
        }
        else
        {
            nc = tmp.imag() / pTIterStack[iter].nNorm;
            rDy = ty * nc;
            rDz = tz * nc;
            rDw = tw * nc;
        }
        rDist = dist;
        return true;
    }
    else
        return false;
}

template <class Estimator>
inline void QuaternionPwrFractalRules<Estimator>::
IterateCalc(DBL& rX, DBL& rY, DBL& rZ, DBL& rW, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    Complex tmp1, tmp2, lg0, lg1;
    DBL nNorm1, nNorm2, normFVal1, normFVal2;

    typename QuaternionPwrFractalRules::IterationData *pIterStack =
        reinterpret_cast<typename QuaternionPwrFractalRules::IterationData *>(pIterData);

    pIterStack[iter].nNorm[0] = nNorm1 = sqrt(rY * rY + rZ * rZ + rW * rW);

    tmp1 = Complex(rX, nNorm1);

    complex_fn::Ln(tmp1, tmp1);

    if (nNorm1 == 0.0)
    {
        pIterStack[iter].normFVal[0] = normFVal1 = tmp1.imag();
    }
    else
    {
        pIterStack[iter].normFVal[0] = normFVal1 = tmp1.imag() / nNorm1;
    }

    lg0 = Complex(tmp1.real(), rY * normFVal1);

    lg1 = Complex(rZ * normFVal1, rW * normFVal1);

    complex_fn::Mult(lg0, lg0, mExponent);
    if (InfoRulesBase<>::mInfo.funcType.variant == kVar_Left)
    {
        complex_fn::Mult(lg1, lg1, mExponent);
    }
    else
    {
        complex_fn::Mult(lg1, lg1, mExponentConj);
    }

    pIterStack[iter].lg[0] = lg0;
    pIterStack[iter].lg[1] = lg1;

    pIterStack[iter].nNorm[1] = nNorm2 = sqrt(complex_fn::Norm(lg0) + complex_fn::Norm(lg1));

    tmp1 = Complex(lg0.real(), nNorm2);

    complex_fn::Exp(tmp1, tmp1);

    pIterStack[iter].expVal = tmp1;

    if (nNorm2 == 0.0)
    {
        pIterStack[iter].normFVal[1] = normFVal2 = tmp1.imag();
    }
    else
    {
        pIterStack[iter].normFVal[1] = normFVal2 = tmp1.imag() / nNorm2;
    }

    rX = tmp1.real() + mJuliaParm[X];
    rY = lg0.imag() * normFVal2 + mJuliaParm[Y];
    rZ = lg1.real() * normFVal2 + mJuliaParm[Z];
    rW = lg1.imag() * normFVal2 + mJuliaParm[W];
}

template <class Estimator>
inline void QuaternionPwrFractalRules<Estimator>::
ApplyDirDerivCalc(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL x, DBL y, DBL z, DBL w,
                  int iter, bool samePoint, const Fractal *pFractal, void *pIterData) const
{
    static DBL ny1, ny2, nz1, nz2, nw1, nw2, nNorm1, nNorm2, normFVal1, normFVal2;
    static Complex tmp21, tmp22;
    DBL tmpx, tmpy, tmpz, tmpw;

    if (!samePoint)
    {
        Complex tmp1;

        typename QuaternionPwrFractalRules::IterationData *pIterStack =
            reinterpret_cast<typename QuaternionPwrFractalRules::IterationData *>(pIterData);

        nNorm1 = pIterStack[iter].nNorm[0];
        normFVal1 = pIterStack[iter].normFVal[0];

        tmp1 = Complex(x, nNorm1);

        complex_fn::Recip(tmp21, tmp1);

        nNorm2 = pIterStack[iter].nNorm[1];
        normFVal2 = pIterStack[iter].normFVal[1];

        tmp22 = pIterStack[iter].expVal;

        if (nNorm1 != 0.0)
        {
            ny1 = y / nNorm1;
            nz1 = z / nNorm1;
            nw1 = w / nNorm1;
        }

        if (nNorm2 != 0.0)
        {
            ny2 = pIterStack[iter].lg[0].imag() / nNorm2;
            nz2 = pIterStack[iter].lg[1].real() / nNorm2;
            nw2 = pIterStack[iter].lg[1].imag() / nNorm2;
        }
    }

    if (nNorm1 == 0.0)
    {
        AltGenDerivs(rDx, rDy, rDz, rDw, tmp21);
    }
    else
    {
        CalcGenDerivs(rDx, rDy, rDz, rDw, normFVal1, ny1, nz1, nw1, tmp21);
    }

    ComponentComplexMult(rDx, rDy, mExponent);
    if (InfoRulesBase<>::mInfo.funcType.variant == kVar_Left)
    {
        ComponentComplexMult(rDz, rDw, mExponent);
    }
    else
    {
        ComponentComplexMult(rDz, rDw, mExponentConj);
    }

    if (nNorm2 == 0.0)
    {
        AltGenDerivs(rDx, rDy, rDz, rDw, tmp22);
    }
    else
    {
        CalcGenDerivs(rDx, rDy, rDz, rDw, normFVal2, ny2, nz2, nw2, tmp22);
    }
}

template <class Estimator>
bool QuaternionPwrFractalRules<Estimator>::
DiscontinuityCheck(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL& rDist,
                   DBL tx, DBL ty, DBL tz, DBL tw, DBL px, DBL py, DBL pz, DBL pw,
                   int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    Complex tmp, tPt, pPt;
    DBL dist, nc;

    typename QuaternionPwrFractalRules::IterationData *pTIterStack =
        reinterpret_cast<typename QuaternionPwrFractalRules::IterationData *>(pTIterData),
        *pPIterStack = reinterpret_cast<typename QuaternionPwrFractalRules::IterationData *>(pPIterData);

    tPt = Complex(tx, pTIterStack[iter].nNorm[0]);

    pPt = Complex(px, pPIterStack[iter].nNorm[0]);

    if (complex_fn::NegReal_DTest(tmp, dist, tPt, pPt, mExponent))
    {
        rDx = tmp.real();
        if (pTIterStack[iter].nNorm[0] == 0.0)
        {
            rDy = tmp.imag();
            rDz = rDw = 0.0;
        }
        else
        {
            nc = tmp.imag() / pTIterStack[iter].nNorm[0];
            rDy = ty * nc;
            rDz = tz * nc;
            rDw = tw * nc;
        }
        rDist = dist;
        return true;
    }
    else
        return false;
}

void QuaternionDispatchInit() {
    static const MagicRulesDispatch<QuaternionSqrFractalRules,
                                    magic::Estim_OrigSpecial<DistEstimatorSpecialOrig_QuatSqr> >
        QuatSqrDispatch(CreateFuncType(kQuaternion, kFunc_Sqr, kVar_Normal));

    static const MagicRulesDispatch<QuaternionCubeFractalRules,
                                    magic::Estim_OrigSpecial<DistEstimatorSpecialOrig_QuatCube> >
        QuatCubeDispatch(CreateFuncType(kQuaternion, kFunc_Cube, kVar_Normal));

    static const MagicRulesDispatch<QuaternionRecipFractalRules>
        QuatRecipDispatch(CreateFuncType(kQuaternion, kFunc_Reciprocal, kVar_Normal));

    static const MagicRulesDispatch<QuaternionFuncFractalRules>
        QuatFuncDispatch(Func_FuncTypeSet(kQuaternion), -1);

    static const MagicRulesDispatch<QuaternionPwrFractalRules>
        QuatPwrDispatch(CreateFuncTypeSet(2, CreateFuncType(kQuaternion, kFunc_Pwr, kVar_Left),
                                          CreateFuncType(kQuaternion, kFunc_Pwr, kVar_Right)));
}

}
