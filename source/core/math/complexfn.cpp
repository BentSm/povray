/*******************************************************************************
 * complexfn.cpp
 *
 * This module implements complex functions, their derivatives, and discontinuity
 * tests.
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
 * $File: //depot/public/povray/3.x/source/backend/math/hcmplx.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// configcore.h must always be the first POV file included in core *.cpp files (pulls in platform config)
#include "core/configcore.h"
#include "core/math/complexfn.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

namespace complex_fn
{

/*****************************************************************************
 *
 * FUNCTION  Complex functions
 *
 * INPUT     pointer to source complex number
 *
 * OUTPUT    fn(input)
 *
 * RETURNS   void
 *
 * AUTHOR
 *
 *   Tim Wegner
 *
 * DESCRIPTION  Calculate common functions on complexes
 *   Since our purpose is fractals, error checking is lax.
 *
 * CHANGES
 *
 ******************************************************************************/

void Sqr(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL tmpx;
    tmpx = pov::Sqr(source.x) - pov::Sqr(source.y);
    rTarget.y = 2.0 * source.x * source.y;
    rTarget.x = tmpx;
}

void Sqrt(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL mag;
    DBL theta;

    if(source.x == 0.0 && source.y == 0.0)
    {
        rTarget.x = rTarget.y = 0.0;
    }
    else
    {
        mag   = sqrt(sqrt(pov::Sqr(source.x) + pov::Sqr(source.y)));
        theta = atan2(source.y, source.x) / 2;
        rTarget.y = mag * sin(theta);
        rTarget.x = mag * cos(theta);
    }
}

void Sqrt_UpperHalfPlane(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL mag;
    DBL theta;

    if(source.x == 0.0 && source.y == 0.0)
    {
        rTarget.x = rTarget.y = 0.0;
    }
    else
    {
        mag   = sqrt(sqrt(pov::Sqr(source.x) + pov::Sqr(source.y)));
        theta = atan2(source.y, source.x) / 2;
        if (theta < 0.0)
        {
            rTarget.y = -mag * sin(theta);
            rTarget.x = -mag * cos(theta);
        }
        else
        {
            rTarget.y = mag * sin(theta);
            rTarget.x = mag * cos(theta);
        }
    }
}

void Ln(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL mod, zx, zy;
    mod = pov::Sqr(source.x) + pov::Sqr(source.y);
    zx = 0.5 * log(mod);
    zy = atan2(source.y,source.x);

    rTarget.x = zx;
    rTarget.y = zy;
}

void Recip(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL mod = pov::Sqr(source.x) + pov::Sqr(source.y);
    if (mod == 0.0)
        return;
    rTarget.x = source.x / mod;
    rTarget.y = -source.y / mod;
}

void Exp(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL expx;
    expx = exp(source.x);
    rTarget.x = expx * cos(source.y);
    rTarget.y = expx * sin(source.y);
}

void Sin(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget.x = sin(source.x) * cosh(source.y);
    rTarget.y = cos(source.x) * sinh(source.y);
}

void Sinh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget.x = sinh(source.x) * cos(source.y);
    rTarget.y = cosh(source.x) * sin(source.y);
}


void Cos(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget.x = cos(source.x) * cosh(source.y);
    rTarget.y = -sin(source.x) * sinh(source.y);
}

void Cos_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget.x = -sin(source.x) * cosh(source.y);
    rTarget.y = -cos(source.x) * sinh(source.y);
}

void Cosh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget.x = cosh(source.x) * cos(source.y);
    rTarget.y = sinh(source.x) * sin(source.y);
}


/* rz=Arcsin(z)=-i*Log{i*z+sqrt(1-z*z)} */
void ASin(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex tempz1, tempz2;

    Sqr(tempz1, source);
    tempz1.x = 1 - tempz1.x; tempz1.y = -tempz1.y;
    Sqrt(tempz1, tempz1);

    tempz2.x = -source.y; tempz2.y = source.x;
    tempz1.x += tempz2.x;  tempz1.y += tempz2.y;

    Ln(tempz1, tempz1);
    rTarget.x = tempz1.y;  rTarget.y = -tempz1.x;
}

void ASin_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp.x = 1 - temp.x; temp.y = -temp.y;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

void ACos(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;

    Sqr(temp, source);
    temp.x -= 1;
    /* This is important to get the standard principal value for arccosine. */
    Sqrt_UpperHalfPlane(temp, temp);

    temp.x += source.x; temp.y += source.y;

    Ln(temp, temp);
    rTarget.x = temp.y;  rTarget.y = -temp.x;
}

/* There are some oddities with this function, due to branch cuts. */
void ACos_Alt(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;

    Sqr(temp, source);
    temp.x -= 1;
    Sqrt(temp, temp);

    temp.x += source.x; temp.y += source.y;

    Ln(temp, temp);
    rTarget.x = temp.y;  rTarget.y = -temp.x;
}

void ACos_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp.x -= 1;
    Sqrt(temp, temp);

    Recip(temp, temp);
    rTarget.x = temp.y; rTarget.y = -temp.x;
}

void ASinh(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;

    Sqr(temp, source);
    temp.x += 1;
    Sqrt(temp, temp);
    temp.x += source.x; temp.y += source.y;
    Ln(rTarget, temp);
}

void ASinh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp.x += 1;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

/* rz=Arccosh(z)=Log(z+sqrt(z*z-1)} */
void ACosh(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex tempz;
    Sqr(tempz, source);
    tempz.x -= 1;
    Sqrt_UpperHalfPlane(tempz, tempz);
    tempz.x += source.x; tempz.y += source.y;
    Ln(rTarget, tempz);
}

/* There are some oddities with this function, due to branch cuts. */
void ACosh_Alt(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex tempz;
    Sqr(tempz, source);
    tempz.x -= 1;
    Sqrt(tempz, tempz);
    tempz.x += source.x; tempz.y += source.y;
    Ln(rTarget, tempz);
}

void ACosh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp.x -= 1;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

/* rz=Arctanh(z)=1/2*Log{(1+z)/(1-z)} */
void ATanh(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp0, temp1, temp2;

    if (source.x == 0.0)
    {
        rTarget.x = 0;
        rTarget.y = atan(source.y);
        return;
    }
    else
    {
        if (fabs(source.x) == 1.0 && source.y == 0.0)
        {
            return;
        }
        else if (fabs(source.x) < 1.0 && source.y == 0.0)
        {
            rTarget.x = log((1+source.x)/(1-source.x))/2;
            rTarget.y = 0;
            return;
        }
        else
        {
            temp0.x = 1 + source.x; temp0.y = source.y;
            temp1.x = 1 - source.x; temp1.y = -source.y;
            Div(temp2, temp0, temp1);
            Ln(temp2, temp2);
            rTarget.x = .5 * temp2.x; rTarget.y = .5 * temp2.y;
            return;
        }
    }
}

void ATanh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp.x = 1 - temp.x; temp.y = -temp.y;
    Recip(rTarget, temp);
}

/* rz=Arctan(z)=i/2*Log{(1-i*z)/(1+i*z)} */
void ATan(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp0, temp1, temp2, temp3;
    if (source.x == 0.0 && source.y == 0.0)
        rTarget.x = rTarget.y = 0;
    else if (source.x != 0.0 && source.y == 0.0){
        rTarget.x = atan(source.x);
        rTarget.y = 0;
    }
    else if (source.x == 0.0 && source.y != 0.0){
        temp0.x = source.y;  temp0.y = 0.0;
        ATanh(temp0, temp0);
        rTarget.x = -temp0.y; rTarget.y = temp0.x;
    }
    else if (source.x != 0.0 && source.y != 0.0)
    {
        temp0.x = -source.y; temp0.y = source.x;
        temp1.x = 1 - temp0.x; temp1.y = -temp0.y;
        temp2.x = 1 + temp0.x; temp2.y = temp0.y;

        Div(temp3, temp1, temp2);
        Ln(temp3, temp3);
        rTarget.x = -temp3.y * .5; rTarget.y = .5 * temp3.x;
    }
}

void ATan_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp.x = 1 + temp.x;
    Recip(rTarget, temp);
}

void Tan(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL x, y, sinx, cosx, sinhy, coshy, denom;
    x = 2 * source.x;
    y = 2 * source.y;
    sinx = sin(x); cosx = cos(x);
    sinhy = sinh(y); coshy = cosh(y);
    denom = cosx + coshy;
    if (denom == 0)
        return;
    rTarget.x = sinx / denom;
    rTarget.y = sinhy / denom;
}

void Tan_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL x, y, sinx, cosx, sinhy, coshy, denom;
    x = 2 * source.x;
    y = 2 * source.y;
    sinx = sin(x); cosx = cos(x);
    sinhy = sinh(y); coshy = cosh(y);
    denom = pov::Sqr(cosx + coshy) / 2;
    if (denom == 0)
        return;
    rTarget.x = (cosx*coshy+1)/denom;
    rTarget.y = sinx*sinhy/denom;
}

void Tanh(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL x, y, siny, cosy, sinhx, coshx, denom;
    x = 2 * source.x;
    y = 2 * source.y;
    siny = sin(y); cosy = cos(y);
    sinhx = sinh(x); coshx = cosh(x);
    denom = coshx + cosy;
    if (denom == 0)
        return;
    rTarget.x = sinhx/denom;
    rTarget.y = siny/denom;
}

void Tanh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    DBL x, y, siny, cosy, sinhx, coshx, denom;
    x = 2 * source.x;
    y = 2 * source.y;
    siny = sin(y); cosy = cos(y);
    sinhx = sinh(x); coshx = cosh(x);
    denom = pov::Sqr(coshx + cosy) / 2;
    if (denom == 0)
        return;
    rTarget.x = (coshx*cosy+1)/denom;
    rTarget.y = -sinhx*siny/denom;
}

void Pwr(Complex& rTarget, const Complex& source1, const Complex& source2)
{
    Complex cLog, t;
    DBL e2x;

    if (source1.x == 0 && source1.y == 0)
    {
        rTarget.x = rTarget.y = 0.0;
        return;
    }

    Ln(cLog, source1);
    Mult(t, cLog, source2);

    if (t.x < -690)
        e2x = 0;
    else
        e2x = exp(t.x);
    rTarget.x = e2x * cos(t.y);
    rTarget.y = e2x * sin(t.y);
}

void Pwr_Deriv(Complex& rTarget, const Complex& source1, const Complex& source2)
{
    Complex cLog, t, aMul, temp;
    DBL e2x;

    if (source1.x == 0 && source1.y == 0)
    {
        rTarget.x = rTarget.y = 0.0;
        return;
    }

    Ln(cLog, source1);
    aMul.x = source2.x - 1; aMul.y = source2.y;
    Mult(t, cLog, aMul);

    if (t.x < -690)
        e2x = 0;
    else
        e2x = exp(t.x);
    temp.x = e2x * cos(t.y);
    temp.y = e2x * sin(t.y);

    Mult(rTarget, temp, source2);
}

bool False_DTest(Complex&, DBL&, const Complex&, const Complex&, const Complex&)
{
    return false;
}

bool NegReal_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if (((cT.y >= 0.0) ^ (cP.y >= 0.0)) && (cT.x < 0.0 || cP.x < 0.0))
    {
        tmp = -cP.y / (cT.y - cP.y);
        if ((cT.x < 0.0 && cP.x < 0.0) || ((cT.x - cP.x) * tmp + cP.x < 0.0))
        {
            rDist = tmp;
            rNormal.x = 0;
            rNormal.y = (cT.y >= 0.0 ? -1.0 : 1.0);
            return true;
        }
    }
    return false;
}

bool ASin_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if (((cT.y >= 0.0) ^ (cP.y >= 0.0)) && (fabs(cT.x) > 1.0 || fabs(cT.x) > 1.0))
    {
        tmp = -cP.y / (cT.y - cP.y);
        if (fabs((cT.x - cP.x) * tmp + cP.x) > 1.0)
        {
            rDist = tmp;
            rNormal.x = 0;
            rNormal.y = (cT.y >= 0.0 ? -1.0 : 1.0);
            return true;
        }
    }
    return false;
}

bool ASinh_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if (((cT.x >= 0.0) ^ (cP.x >= 0.0)) && (fabs(cT.y) > 1.0 || fabs(cT.y) > 1.0))
    {
        tmp = -cP.x / (cT.x - cP.x);
        if (fabs((cT.y - cP.y) * tmp + cP.y) > 1.0)
        {
            rDist = tmp;
            rNormal.x = (cT.x >= 0.0 ? -1.0 : 1.0);
            rNormal.y = 0;
            return true;
        }
    }
    return false;
}

bool ACos_Alt_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if ((cT.x >= 0.0) ^ (cP.x >= 0.0))
    {
        rDist = -cP.x / (cT.x - cP.x);
        rNormal.x = (cT.x >= 0.0 ? -1.0 : 1.0);
        rNormal.y = 0;
    }
    else if (((cT.y >= 0.0) ^ (cP.y >= 0.0)) && (cT.x < 1.0 || cP.x < 1.0))
    {
        tmp = -cP.y / (cT.y - cP.y);
        if ((cT.x - cP.x) * tmp + cP.x < 1.0)
        {
            rDist = tmp;
            rNormal.x = 0;
            rNormal.y = (cT.y >= 0.0 ? -1.0 : 1.0);
            return true;
        }
    }
    return false;
}

}

}
