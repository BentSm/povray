//******************************************************************************
///
/// @file core/math/complexfn.h
///
/// This file contains prototypes, etc., for complex math (used for julia_fractal).
///
/// @todo   It might be a good thing to replace much of this with C++03 and
/// TR1 ~= C99 ~= boost/math/complex functions.
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

#ifndef POVRAY_CORE_COMPLEXFN_H
#define POVRAY_CORE_COMPLEXFN_H

#include "core/coretypes.h"

namespace pov
{

typedef DBL Complex[2];

namespace complex_fn
{

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

/*****************************************************************************
* Global typedefs
******************************************************************************/

typedef void FuncType(Complex, const Complex, const Complex);
typedef bool DiscontinuityTestFn(Complex, DBL&, const Complex, const Complex, const Complex);

/*****************************************************************************
* Global variables
******************************************************************************/

/*****************************************************************************
* Global functions
******************************************************************************/

inline void Mult(Complex rTarget, const Complex source1, const Complex source2)
{
    DBL tmpx;
    tmpx = source1[X] * source2[X] - source1[Y] * source2[Y];
    rTarget[Y] = source1[X] * source2[Y] + source1[Y] * source2[X];
    rTarget[X] = tmpx;
}

inline void Div(Complex rTarget, const Complex source1, const Complex source2)
{
    DBL mod, tmpx, yxmod, yymod;
    mod = pov::Sqr(source2[X]) + pov::Sqr(source2[Y]);
    if (mod==0)
        return;
    yxmod = source2[X]/mod;
    yymod = - source2[Y]/mod;
    tmpx = source1[X] * yxmod - source1[Y] * yymod;
    rTarget[Y] = source1[X] * yymod + source1[Y] * yxmod;
    rTarget[X] = tmpx;
}

void Exp(Complex rTarget, const Complex source, const Complex unused = NULL);
void Ln(Complex rTarget, const Complex source, const Complex unused = NULL);
void Sin(Complex rTarget, const Complex source, const Complex unused = NULL);
void ASin(Complex rTarget, const Complex source, const Complex unused = NULL);
void Sinh(Complex rTarget, const Complex source, const Complex unused = NULL);
void ASinh(Complex rTarget, const Complex source, const Complex unused = NULL);
void Cos(Complex rTarget, const Complex source, const Complex unused = NULL);
void ACos(Complex rTarget, const Complex source, const Complex unused = NULL);
void ACos_Alt(Complex rTarget, const Complex source, const Complex unused = NULL);
void Cosh(Complex rTarget, const Complex source, const Complex unused = NULL);
void ACosh(Complex rTarget, const Complex source, const Complex unused = NULL);
void ACosh_Alt(Complex rTarget, const Complex source, const Complex unused = NULL);
void Tan(Complex rTarget, const Complex source, const Complex unused = NULL);
void ATan(Complex rTarget, const Complex source, const Complex unused = NULL);
void Tanh(Complex rTarget, const Complex source, const Complex unused = NULL);
void ATanh(Complex rTarget, const Complex source, const Complex unused = NULL);
void Sqrt(Complex rTarget, const Complex source, const Complex unused = NULL);
void Sqr(Complex rTarget, const Complex source, const Complex unused = NULL);
void Recip(Complex rTarget, const Complex source, const Complex unused = NULL);
void Pwr(Complex rTarget, const Complex source1, const Complex source2);

void ASin_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void ASinh_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void Cos_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void ACos_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void ACosh_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void Tan_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void ATan_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void Tanh_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void ATanh_Deriv(Complex rTarget, const Complex source, const Complex unused = NULL);
void Pwr_Deriv(Complex rTarget, const Complex source1, const Complex source2);

bool False_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused = NULL);
bool NegReal_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused = NULL);
bool ASin_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused = NULL);
bool ASinh_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused = NULL);
bool ACos_Alt_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused = NULL);
bool ACosh_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused = NULL);

}

}

#endif
