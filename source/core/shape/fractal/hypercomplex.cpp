//******************************************************************************
///
/// @file core/shape/fractal/hypercomplex.cpp
///
/// This module implements hypercomplex Julia fractals.
///
/// This file was written by Pascal Massimino.
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
#include "core/shape/fractal/hypercomplex.h"

#include "core/math/complexfn.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/dispatch.h"
#include "core/shape/fractal/distestimator.h"
#include "core/shape/fractal/estimmagic.h"
#include "core/shape/fractal/magicimpl.h"
#include "core/shape/fractal/util.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

static inline void ComplexSqrAdd(Complex& rC, const Complex& c, const Complex& a);
static inline void ComplexSqrDeriv(Complex& rC, const Complex& c);
static inline void ComplexCubeAdd(Complex& rC, const Complex& c, const Complex& a);
static inline void ComplexCubeDeriv(Complex& rC, const Complex& c);
static inline void ComplexRecipAdd(Complex& rC, const Complex& c, const Complex& a);
static inline void ComplexRecipDeriv(Complex& rC, const Complex& c);

static inline void ComplexSqrAdd(Complex& rC, const Complex& c, const Complex& a) {
    Complex temp;
    complex_fn::Sqr(temp, c);
    rC = temp + a;
}

static inline void ComplexSqrDeriv(Complex& rC, const Complex& c)
{
    rC = 2.0 * c;
}

static inline void ComplexCubeAdd(Complex& rC, const Complex& c, const Complex& a)
{
    DBL tmpxx, tmpyy;
    tmpxx = Sqr(c.real());
    tmpyy = Sqr(c.imag());
    rC = Complex(c.real() * (tmpxx - 3 * tmpyy), c.imag() * (3 * tmpxx - tmpyy)) + a;
}

static inline void ComplexCubeDeriv(Complex& rC, const Complex& c)
{
    Complex temp;
    complex_fn::Sqr(temp, c);
    rC = 3.0 * temp;
}

static inline void ComplexRecipAdd(Complex& rC, const Complex& c, const Complex& a)
{
    Complex temp;
    complex_fn::Recip(temp, c);
    rC = temp + a;
}

static inline void ComplexRecipDeriv(Complex& rC, const Complex& c)
{
    DBL mod = complex_fn::Norm(c);

    rC = Complex((1.0 - 2.0 * c.real() / mod) / mod, 2.0 * c.real() * c.imag() / Sqr(mod));
}

template <class Estimator>
inline void HypercomplexSqrFractalRules<Estimator>::
IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    ComplexSqrAdd(rC0, rC0, mDuplexJuliaParm[0]);
    ComplexSqrAdd(rC1, rC1, mDuplexJuliaParm[1]);
}

template <class Estimator>
inline void HypercomplexSqrFractalRules<Estimator>::
ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
               int iter, const Fractal *pFractal, void *pIterData) const
{
    Complex tmp0, tmp1;

    ComplexSqrDeriv(tmp0, c0);
    ComplexSqrDeriv(tmp1, c1);

    complex_fn::Mult(rD0, rD0, tmp0);
    complex_fn::Mult(rD1, rD1, tmp1);
}

template <class Estimator>
inline void HypercomplexCubeFractalRules<Estimator>::
IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    ComplexCubeAdd(rC0, rC0, mDuplexJuliaParm[0]);
    ComplexCubeAdd(rC1, rC1, mDuplexJuliaParm[1]);
}

template <class Estimator>
inline void HypercomplexCubeFractalRules<Estimator>::
ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
               int iter, const Fractal *pFractal, void *pIterData) const
{
    Complex tmp0, tmp1;

    ComplexCubeDeriv(tmp0, c0);
    ComplexCubeDeriv(tmp1, c1);

    complex_fn::Mult(rD0, rD0, tmp0);
    complex_fn::Mult(rD1, rD1, tmp1);
}

template <class Estimator>
inline void HypercomplexRecipFractalRules<Estimator>::
IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
                            int iter, const Fractal *pFractal, void *pIterData) const
{
    ComplexRecipAdd(rC0, rC0, mDuplexJuliaParm[0]);
    ComplexRecipAdd(rC1, rC1, mDuplexJuliaParm[1]);
}

template <class Estimator>
inline void HypercomplexRecipFractalRules<Estimator>::
ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
               int iter, const Fractal *pFractal, void *pIterData) const
{
    Complex tmp0, tmp1;

    ComplexRecipDeriv(tmp0, c0);
    ComplexRecipDeriv(tmp1, c1);

    complex_fn::Mult(rD0, rD0, tmp0);
    complex_fn::Mult(rD1, rD1, tmp1);
}

template <class Estimator>
inline void HypercomplexFuncFractalRules<Estimator>::
IterateCalc(Complex& rC0, Complex& rC1, DBL norm,
            int iter, const Fractal *pFractal, void *pIterData) const
{
    (*(mFunc.pFunc))(rC0, rC0, mExponent);
    (*(mFunc.pFunc))(rC1, rC1, mExponent);

    rC0 += mDuplexJuliaParm[0];

    rC1 += mDuplexJuliaParm[1];
}

template <class Estimator>
inline void HypercomplexFuncFractalRules<Estimator>::
ApplyDerivCalc(Complex& rD0, Complex& rD1, const Complex& c0, const Complex& c1,
               int iter, const Fractal *pFractal, void *pIterData) const
{
    Complex tmp0, tmp1;

    (*(mFunc.pDeriv))(tmp0, c0, mExponent);
    (*(mFunc.pDeriv))(tmp1, c1, mExponent);

    complex_fn::Mult(rD0, rD0, tmp0);
    complex_fn::Mult(rD1, rD1, tmp1);
}

template <class Estimator>
bool HypercomplexFuncFractalRules<Estimator>::
DiscontinuityCheck(Complex& rD0, Complex& rD1, DBL& rDist,
                   const Complex& t0, const Complex& t1, const Complex& p0, const Complex& p1,
                   int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const
{
    Complex tmp;
    DBL dist;
    if ((*(mFunc.pDisc))(tmp, dist, t0, p0, mExponent))
    {
        rD0 = tmp;
        rD1 = 0.0;
        rDist = dist;
        return true;
    }
    else if ((*(mFunc.pDisc))(tmp, dist, t1, p1, mExponent))
    {
        rD0 = 0.0;
        rD1 = tmp;
        rDist = dist;
        return true;
    }
    else
        return false;
}

void HypercomplexDispatchInit() {
    static const MagicRulesDispatch<HypercomplexSqrFractalRules>
        HyperSqrDispatch(CreateFuncType(kHypercomplex, kFunc_Sqr, kVar_Normal));

    static const MagicRulesDispatch<HypercomplexCubeFractalRules>
        HyperCubeDispatch(CreateFuncType(kHypercomplex, kFunc_Cube, kVar_Normal));

    static const MagicRulesDispatch<HypercomplexRecipFractalRules>
        HyperRecipDispatch(CreateFuncType(kHypercomplex, kFunc_Reciprocal, kVar_Normal));

    static const MagicRulesDispatch<HypercomplexFuncFractalRules>
        HyperFuncDispatch(Func_FuncTypeSet(kHypercomplex), -1);
}

}
