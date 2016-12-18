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
#include "core/shape/fractal/hypercomplex.h"

#include "core/math/complexfn.h"
#include "core/shape/fractal.h"
#include "core/shape/fractal/dispatch.h"

// this must be the last file included
#include "base/povdebug.h"

#define APPLY_DUPLEX_2(func, a, b)              \
    func(AsComplex(a, 0), AsComplex(b, 0));     \
    func(AsComplex(a, 1), AsComplex(b, 1))

#define APPLY_DUPLEX_2_PLUS_1(func, a, b, c)    \
    func(AsComplex(a, 0), AsComplex(b, 0), c);  \
    func(AsComplex(a, 1), AsComplex(b, 1), c)

#define APPLY_DUPLEX_3(func, a, b, c)                           \
    func(AsComplex(a, 0), AsComplex(b, 0), AsComplex(c, 0));    \
    func(AsComplex(a, 1), AsComplex(b, 1), AsComplex(c, 1))


namespace pov
{

static inline void ComplexSqrAdd(Complex rC, const Complex c, const Complex a);
static inline void ComplexCubeAdd(Complex rC, const Complex c, const Complex a);
static inline void ComplexCubeDeriv(Complex rC, const Complex c);
static inline void ComplexRecipAdd(Complex rC, const Complex c, const Complex a);
static inline void ComplexRecipDeriv(Complex rC, const Complex c);
static inline void ComplexMultConj(Complex rC, const Complex c, const Complex a);

static inline void ComplexSqrAdd(Complex rC, const Complex c, const Complex a) {
    DBL tmp;
    tmp = Sqr(c[X]) - Sqr(c[Y]) + a[X];
    rC[Y] = 2 * c[X] * c[Y] + a[Y];
    rC[X] = tmp;
}

static inline void ComplexCubeAdd(Complex rC, const Complex c, const Complex a)
{
    DBL tmpxx, tmpyy;
    tmpxx = Sqr(c[X]);
    tmpyy = Sqr(c[Y]);
    rC[X] = c[X] * (tmpxx - 3 * tmpyy) + a[X];
    rC[Y] = c[Y] * (3 * tmpxx - tmpyy) + a[Y];
}

static inline void ComplexCubeDeriv(Complex rC, const Complex c)
{
    DBL tmp;
    tmp = Sqr(c[X]) - Sqr(c[Y]);
    rC[Y] = 2.0 * c[X] * c[Y];
    rC[X] = tmp;
}

static inline void ComplexRecipAdd(Complex rC, const Complex c, const Complex a)
{
    DBL mod = Sqr(c[X]) + Sqr(c[Y]);
    if (mod == 0.0)
        return;

    rC[X] = c[X] / mod + a[X];
    rC[Y] = -c[Y] / mod + a[Y];
}

static inline void ComplexRecipDeriv(Complex rC, const Complex c)
{
    DBL mod = Sqr(c[X]) + Sqr(c[Y]);
    if (mod == 0.0)
        return;

    rC[Y] = 2 * c[X] * c[Y] / Sqr(mod);
    rC[X] = (1 - 2 * c[X] / mod) / mod;
}

static inline void ComplexMultConj(Complex rC, const Complex c, const Complex a)
{
    DBL tmp = c[X] * a[X] + c[Y] * a[Y];
    rC[Y] = c[Y] * a[X] - c[X] * a[Y];
    rC[X] = tmp;
}

void HypercomplexSqrFractalRules::
IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    APPLY_DUPLEX_3(ComplexSqrAdd, rC, rC, mDuplexJuliaParm);
}

void HypercomplexSqrFractalRules::
GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const
{
    APPLY_DUPLEX_3(ComplexMultConj, rD, rD, c);

    rMult *= 2.0;
}

void HypercomplexCubeFractalRules::
IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    APPLY_DUPLEX_3(ComplexCubeAdd, rC, rC, mDuplexJuliaParm);
}

void HypercomplexCubeFractalRules::
GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Vector4d tmp;

    APPLY_DUPLEX_2(ComplexCubeDeriv, tmp, c);
    APPLY_DUPLEX_3(ComplexMultConj, rD, rD, tmp);

    rMult *= 3.0;
}

void HypercomplexRecipFractalRules::
IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{
    APPLY_DUPLEX_3(ComplexRecipAdd, rC, rC, mDuplexJuliaParm);
}

void HypercomplexRecipFractalRules::
GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Vector4d tmp;

    APPLY_DUPLEX_2(ComplexRecipDeriv, tmp, c);
    APPLY_DUPLEX_3(ComplexMultConj, rD, rD, tmp);
}

void HypercomplexFuncFractalRules::
IterateCalc(Vector4d& rC, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const
{

    APPLY_DUPLEX_2_PLUS_1((*mFunc.pFunc), rC, rC, *mExponent);

    rC += mDuplexJuliaParm;
}

void HypercomplexFuncFractalRules::
GradientCalc(Vector4d& rD, const Vector4d& c, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const
{
    Vector4d tmp;

    APPLY_DUPLEX_2_PLUS_1((*mFunc.pDeriv), tmp, c, *mExponent);
    APPLY_DUPLEX_3(ComplexMultConj, rD, rD, tmp);
}

bool HypercomplexFuncFractalRules::
DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                   int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const
{
    Complex tmp;
    DBL dist;
    if ((*(mFunc.pDisc))(tmp, dist, AsComplex(t, 0), AsComplex(p, 0), *mExponent))
    {
        AssignComplex(AsComplex(rD, 0), tmp);
        rD[Z] = rD[W] = 0.0;
        rDist = dist;
        return true;
    }
    else if ((*(mFunc.pDisc))(tmp, dist, AsComplex(t, 1), AsComplex(p, 1), *mExponent))
    {
        rD[X] = rD[Y] = 0.0;
        AssignComplex(AsComplex(rD, 1), tmp);
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
