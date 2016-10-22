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
    DBL tmp;
    tmp = Sqr(c.x) - Sqr(c.y) + a.x;
    rC.y = 2 * c.x * c.y + a.y;
    rC.x = tmp;
}

static inline void ComplexSqrDeriv(Complex& rC, const Complex& c)
{
    rC.x = 2.0 * c.x;
    rC.y = 2.0 * c.y;
}

static inline void ComplexCubeAdd(Complex& rC, const Complex& c, const Complex& a)
{
    DBL tmpxx, tmpyy;
    tmpxx = Sqr(c.x);
    tmpyy = Sqr(c.y);
    rC.x = c.x * (tmpxx - 3 * tmpyy) + a.x;
    rC.y = c.y * (3 * tmpxx - tmpyy) + a.y;
}

static inline void ComplexCubeDeriv(Complex& rC, const Complex& c)
{
    DBL tmp;
    tmp = Sqr(c.x) - Sqr(c.y);
    rC.y = 6 * c.x * c.y;
    rC.x = 3 * tmp;
}

static inline void ComplexRecipAdd(Complex& rC, const Complex& c, const Complex& a)
{
    DBL mod = Sqr(c.x) + Sqr(c.y);
    if (mod == 0.0)
        return;

    rC.x = c.x / mod + a.x;
    rC.y = -c.y / mod + a.y;
}

static inline void ComplexRecipDeriv(Complex& rC, const Complex& c)
{
    DBL mod = Sqr(c.x) + Sqr(c.y);
    if (mod == 0.0)
        return;

    rC.y = 2 * c.x * c.y / Sqr(mod);
    rC.x = (1 - 2 * c.x / mod) / mod;
}

void HypercomplexSqrFractalRules::
IterateCalc(Duplex& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    ComplexSqrAdd(rC[0], rC[0], mDuplexJuliaParm[0]);
    ComplexSqrAdd(rC[1], rC[1], mDuplexJuliaParm[1]);
}

void HypercomplexSqrFractalRules::
DerivCalc(Duplex& rD, const Duplex& c, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Duplex tmp;

    ComplexSqrDeriv(tmp[0], c[0]);
    ComplexSqrDeriv(tmp[1], c[1]);

    complex_fn::Mult(rD[0], rD[0], tmp[0]);
    complex_fn::Mult(rD[1], rD[1], tmp[1]);
}

void HypercomplexCubeFractalRules::
IterateCalc(Duplex& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    ComplexCubeAdd(rC[0], rC[0], mDuplexJuliaParm[0]);
    ComplexCubeAdd(rC[1], rC[1], mDuplexJuliaParm[1]);
}

void HypercomplexCubeFractalRules::
DerivCalc(Duplex& rD, const Duplex& c, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Duplex tmp;

    ComplexCubeDeriv(tmp[0], c[0]);
    ComplexCubeDeriv(tmp[1], c[1]);

    complex_fn::Mult(rD[0], rD[0], tmp[0]);
    complex_fn::Mult(rD[1], rD[1], tmp[1]);
}

void HypercomplexRecipFractalRules::
IterateCalc(Duplex& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    ComplexRecipAdd(rC[0], rC[0], mDuplexJuliaParm[0]);
    ComplexRecipAdd(rC[1], rC[1], mDuplexJuliaParm[1]);
}

void HypercomplexRecipFractalRules::
DerivCalc(Duplex& rD, const Duplex& c, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Duplex tmp;

    ComplexRecipDeriv(tmp[0], c[0]);
    ComplexRecipDeriv(tmp[1], c[1]);

    complex_fn::Mult(rD[0], rD[0], tmp[0]);
    complex_fn::Mult(rD[1], rD[1], tmp[1]);
}

void HypercomplexFuncFractalRules::
IterateCalc(Duplex& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    (*(mFunc.pFunc))(rC[0], rC[0], mExponent);
    (*(mFunc.pFunc))(rC[1], rC[1], mExponent);

    rC[0].x += mDuplexJuliaParm[0].x;
    rC[0].y += mDuplexJuliaParm[0].y;

    rC[1].x += mDuplexJuliaParm[1].x;
    rC[1].y += mDuplexJuliaParm[1].y;
}

void HypercomplexFuncFractalRules::
DerivCalc(Duplex& rD, const Duplex& c, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Duplex tmp;

    (*(mFunc.pDeriv))(tmp[0], c[0], mExponent);
    (*(mFunc.pDeriv))(tmp[1], c[1], mExponent);

    complex_fn::Mult(rD[0], rD[0], tmp[0]);
    complex_fn::Mult(rD[1], rD[1], tmp[1]);
}

bool HypercomplexFuncFractalRules::
DiscontinuityCheck(Duplex& rD, DBL& rDist, const Duplex& t, const Duplex& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    Complex tmp;
    DBL dist;
    if ((*(mFunc.pDisc))(tmp, dist, t[0], p[0], mExponent))
    {
        rD[0] = tmp;
        rD[1].x = rD[1].y = 0.0;
        rDist = dist;
        return true;
    }
    else if ((*(mFunc.pDisc))(tmp, dist, t[1], p[1], mExponent))
    {
        rD[0].x = rD[0].y = 0.0;
        rD[1] = tmp;
        rDist = dist;
        return true;
    }
    else
        return false;
}

void HypercomplexDispatchInit() {
    static const RulesDispatch HyperSqrDispatch(MakeCreatorFunc<HypercomplexSqrFractalRules>,
                                                CreateFuncType(kHypercomplex, kFunc_Sqr, kVar_Normal));

    static const RulesDispatch HyperCubeDispatch(MakeCreatorFunc<HypercomplexCubeFractalRules>,
                                                 CreateFuncType(kHypercomplex, kFunc_Cube, kVar_Normal));

    static const RulesDispatch HyperRecipDispatch(MakeCreatorFunc<HypercomplexRecipFractalRules>,
                                                  CreateFuncType(kHypercomplex, kFunc_Reciprocal, kVar_Normal));

    static const RulesDispatch HyperFuncDispatch(MakeCreatorFunc<HypercomplexFuncFractalRules>,
                                                 Func_FuncTypeSet(kHypercomplex), -1);
}

}
