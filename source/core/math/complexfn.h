/*******************************************************************************
 * complexfn.h
 *
 * This module contains all defines, typedefs, and prototypes for COMPLEXFN.CPP.
 *
 * ---------------------------------------------------------------------------
 * Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
 * Copyright 1991-2013 Persistence of Vision Raytracer Pty. Ltd.
 *
 * POV-Ray is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * POV-Ray is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------------
 * POV-Ray is based on the popular DKB raytracer version 2.12.
 * DKBTrace was originally written by David K. Buck.
 * DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
 * ---------------------------------------------------------------------------
 * $File: //depot/public/povray/3.x/source/backend/math/hcmplx.h $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#ifndef POVRAY_CORE_COMPLEXFN_H
#define POVRAY_CORE_COMPLEXFN_H

#include "core/coretypes.h"

namespace pov
{

namespace complex_fn
{

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

#define cNULL pov::complex_fn::kCDummy

/*****************************************************************************
* Global typedefs
******************************************************************************/

typedef void FuncType(Complex&, const Complex&, const Complex&);
typedef bool DiscontinuityTestFn(Complex&, DBL&, const Complex&, const Complex&, const Complex&);

/*****************************************************************************
* Global variables
******************************************************************************/

/* It really does not matter at all what this is! */
static const Complex kCDummy = {0.0, 0.0};

/*****************************************************************************
* Global functions
******************************************************************************/

inline void Mult(Complex& rTarget, const Complex& source1, const Complex& source2)
{
    DBL tmpx;
    tmpx = source1.x * source2.x - source1.y * source2.y;
    rTarget.y = source1.x * source2.y + source1.y * source2.x;
    rTarget.x = tmpx;
}

inline void Div(Complex& rTarget, const Complex& source1, const Complex& source2)
{
    DBL mod, tmpx, yxmod, yymod;
    mod = pov::Sqr(source2.x) + pov::Sqr(source2.y);
    if (mod==0)
        return;
    yxmod = source2.x/mod;
    yymod = - source2.y/mod;
    tmpx = source1.x * yxmod - source1.y * yymod;
    rTarget.y = source1.x * yymod + source1.y * yxmod;
    rTarget.x = tmpx;
}

void Exp(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Ln(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Sin(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ASin(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Sinh(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ASinh(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Cos(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ACos(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ACos_Alt(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Cosh(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ACosh(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ACosh_Alt(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Tan(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ATan(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Tanh(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ATanh(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Sqrt(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Sqrt_UpperHalfPlane(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Sqr(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Recip(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Pwr(Complex& rTarget, const Complex& source1, const Complex& source2);

void ASin_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ASinh_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Cos_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ACos_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ACosh_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Tan_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ATan_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Tanh_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void ATanh_Deriv(Complex& rTarget, const Complex& source, const Complex& unused = cNULL);
void Pwr_Deriv(Complex& rTarget, const Complex& source1, const Complex& source2);

bool False_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused = cNULL);
bool NegReal_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused = cNULL);
bool ASin_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused = cNULL);
bool ASinh_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused = cNULL);
bool ACos_Alt_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused = cNULL);

}

}

#endif
