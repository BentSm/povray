//******************************************************************************
///
/// @file core/math/complexfn.cpp
///
/// This module implements complex functions, their derivatives, and discontinuity
/// tests.
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

void Sqr(Complex rTarget, const Complex source, const Complex)
{
    DBL tmpx;
    tmpx = pov::Sqr(source[X]) - pov::Sqr(source[Y]);
    rTarget[Y] = 2.0 * source[X] * source[Y];
    rTarget[X] = tmpx;
}

void Sqrt(Complex rTarget, const Complex source, const Complex)
{
    DBL mag;
    DBL theta;

    if(source[X] == 0.0 && source[Y] == 0.0)
    {
        rTarget[X] = rTarget[Y] = 0.0;
    }
    else
    {
        mag   = sqrt(sqrt(pov::Sqr(source[X]) + pov::Sqr(source[Y])));
        theta = atan2(source[Y], source[X]) / 2;
        rTarget[Y] = mag * sin(theta);
        rTarget[X] = mag * cos(theta);
    }
}

void Ln(Complex rTarget, const Complex source, const Complex)
{
    DBL mod, zx, zy;
    mod = pov::Sqr(source[X]) + pov::Sqr(source[Y]);
    zx = 0.5 * log(mod);
    zy = atan2(source[Y],source[X]);

    rTarget[X] = zx;
    rTarget[Y] = zy;
}

void Recip(Complex rTarget, const Complex source, const Complex)
{
    DBL mod = pov::Sqr(source[X]) + pov::Sqr(source[Y]);
    if (mod == 0.0)
        return;
    rTarget[X] = source[X] / mod;
    rTarget[Y] = -source[Y] / mod;
}

void Exp(Complex rTarget, const Complex source, const Complex)
{
    DBL expx;
    expx = exp(source[X]);
    rTarget[X] = expx * cos(source[Y]);
    rTarget[Y] = expx * sin(source[Y]);
}

void Sin(Complex rTarget, const Complex source, const Complex)
{
    rTarget[X] = sin(source[X]) * cosh(source[Y]);
    rTarget[Y] = cos(source[X]) * sinh(source[Y]);
}

void Sinh(Complex rTarget, const Complex source, const Complex)
{
    rTarget[X] = sinh(source[X]) * cos(source[Y]);
    rTarget[Y] = cosh(source[X]) * sin(source[Y]);
}


void Cos(Complex rTarget, const Complex source, const Complex)
{
    rTarget[X] = cos(source[X]) * cosh(source[Y]);
    rTarget[Y] = -sin(source[X]) * sinh(source[Y]);
}

void Cos_Deriv(Complex rTarget, const Complex source, const Complex)
{
    rTarget[X] = -sin(source[X]) * cosh(source[Y]);
    rTarget[Y] = -cos(source[X]) * sinh(source[Y]);
}

void Cosh(Complex rTarget, const Complex source, const Complex)
{
    rTarget[X] = cosh(source[X]) * cos(source[Y]);
    rTarget[Y] = sinh(source[X]) * sin(source[Y]);
}


/* rz=Arcsin(z)=-i*Log{i*z+sqrt(1-z*z)} */
void ASin(Complex rTarget, const Complex source, const Complex)
{
    Complex temp1, temp2;

    Sqr(temp1, source);
    temp1[X] = 1 - temp1[X]; temp1[Y] = -temp1[Y];
    Sqrt(temp1, temp1);

    temp2[X] = -source[Y]; temp2[Y] = source[X];
    temp1[X] += temp2[X];  temp1[Y] += temp2[Y];

    Ln(temp1, temp1);
    rTarget[X] = temp1[Y];  rTarget[Y] = -temp1[X];
}

void ASin_Deriv(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;
    Sqr(temp, source);
    temp[X] = 1 - temp[X]; temp[Y] = -temp[Y];
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

void ACos(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;

    Sqr(temp, source);
    temp[X] -= 1;
    Sqrt(temp, temp);

    /* Tweak to get the principal value. */
    /// @todo Use @c signbit -- negative zeros make a difference!
    if ((source[X] < 0) ^ (source[Y] < 0))
    {
        temp[X] *= -1; temp[Y] *= -1;
    }

    temp[X] += source[X]; temp[Y] += source[Y];

    Ln(temp, temp);
    rTarget[X] = temp[Y];  rTarget[Y] = -temp[X];
}

/* There are some oddities with this function, due to branch cuts. */
void ACos_Alt(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;

    Sqr(temp, source);
    temp[X] -= 1;
    Sqrt(temp, temp);

    temp[X] += source[X]; temp[Y] += source[Y];

    Ln(temp, temp);
    rTarget[X] = temp[Y];  rTarget[Y] = -temp[X];
}

void ACos_Deriv(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;
    Sqr(temp, source);
    temp[X] -= 1;
    Sqrt(temp, temp);

    Recip(temp, temp);
    rTarget[X] = temp[Y]; rTarget[Y] = -temp[X];
}

void ASinh(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;

    Sqr(temp, source);
    temp[X] += 1;
    Sqrt(temp, temp);
    temp[X] += source[X]; temp[Y] += source[Y];
    Ln(rTarget, temp);
}

void ASinh_Deriv(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;
    Sqr(temp, source);
    temp[X] += 1;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

/* rz=Arccosh(z)=Log(z+sqrt(z*z-1)} */
void ACosh(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;

    Sqr(temp, source);
    temp[X] -= 1;
    Sqrt(temp, temp);

    /* Tweak to get the principal value. */
    /// @todo Use @c signbit -- negative zeros make a difference!
    if (source[X] < 0)
    {
        temp[X] *= -1; temp[Y] *= -1;
    }

    temp[X] += source[X]; temp[Y] += source[Y];
    Ln(rTarget, temp);
}

/* There are some oddities with this function, due to branch cuts. */
void ACosh_Alt(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;
    Sqr(temp, source);
    temp[X] -= 1;
    Sqrt(temp, temp);
    temp[X] += source[X]; temp[Y] += source[Y];
    Ln(rTarget, temp);
}

void ACosh_Deriv(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;
    Sqr(temp, source);
    temp[X] -= 1;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

/* rz=Arctanh(z)=1/2*Log{(1+z)/(1-z)} */
void ATanh(Complex rTarget, const Complex source, const Complex)
{
    Complex temp0, temp1, temp2;

    if (source[X] == 0.0)
    {
        rTarget[X] = 0;
        rTarget[Y] = atan(source[Y]);
        return;
    }
    else
    {
        if (fabs(source[X]) == 1.0 && source[Y] == 0.0)
        {
            return;
        }
        else if (fabs(source[X]) < 1.0 && source[Y] == 0.0)
        {
            rTarget[X] = log((1 + source[X]) / (1 - source[X])) / 2;
            rTarget[Y] = 0;
            return;
        }
        else
        {
            temp0[X] = 1 + source[X]; temp0[Y] = source[Y];
            temp1[X] = 1 - source[X]; temp1[Y] = -source[Y];
            Div(temp2, temp0, temp1);
            Ln(temp2, temp2);
            rTarget[X] = .5 * temp2[X]; rTarget[Y] = .5 * temp2[Y];
            return;
        }
    }
}

void ATanh_Deriv(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;
    Sqr(temp, source);
    temp[X] = 1 - temp[X]; temp[Y] = -temp[Y];
    Recip(rTarget, temp);
}

/* rz=Arctan(z)=i/2*Log{(1-i*z)/(1+i*z)} */
void ATan(Complex rTarget, const Complex source, const Complex)
{
    Complex temp0, temp1, temp2, temp3;
    if (source[X] == 0.0 && source[Y] == 0.0)
        rTarget[X] = rTarget[Y] = 0;
    else if (source[X] != 0.0 && source[Y] == 0.0)
    {
        rTarget[X] = atan(source[X]);
        rTarget[Y] = 0;
    }
    else if (source[X] == 0.0 && source[Y] != 0.0)
    {
        temp0[X] = source[Y];  temp0[Y] = 0.0;
        ATanh(temp0, temp0);
        rTarget[X] = -temp0[Y]; rTarget[Y] = temp0[X];
    }
    else if (source[X] != 0.0 && source[Y] != 0.0)
    {
        temp0[X] = -source[Y]; temp0[Y] = source[X];
        temp1[X] = 1 - temp0[X]; temp1[Y] = -temp0[Y];
        temp2[X] = 1 + temp0[X]; temp2[Y] = temp0[Y];

        Div(temp3, temp1, temp2);
        Ln(temp3, temp3);
        rTarget[X] = -temp3[Y] * .5; rTarget[Y] = .5 * temp3[X];
    }
}

void ATan_Deriv(Complex rTarget, const Complex source, const Complex)
{
    Complex temp;
    Sqr(temp, source);
    temp[X] = 1 + temp[X];
    Recip(rTarget, temp);
}

void Tan(Complex rTarget, const Complex source, const Complex)
{
    DBL x, y, sinx, cosx, sinhy, coshy, denom;
    x = 2 * source[X];
    y = 2 * source[Y];
    sinx = sin(x); cosx = cos(x);
    sinhy = sinh(y); coshy = cosh(y);
    denom = cosx + coshy;
    if (denom == 0)
        return;
    rTarget[X] = sinx / denom;
    rTarget[Y] = sinhy / denom;
}

void Tan_Deriv(Complex rTarget, const Complex source, const Complex)
{
    DBL x, y, sinx, cosx, sinhy, coshy, denom;
    x = 2 * source[X];
    y = 2 * source[Y];
    sinx = sin(x); cosx = cos(x);
    sinhy = sinh(y); coshy = cosh(y);
    denom = pov::Sqr(cosx + coshy) / 2;
    if (denom == 0)
        return;
    rTarget[X] = (cosx * coshy + 1) / denom;
    rTarget[Y] = sinx * sinhy / denom;
}

void Tanh(Complex rTarget, const Complex source, const Complex)
{
    DBL x, y, siny, cosy, sinhx, coshx, denom;
    x = 2 * source[X];
    y = 2 * source[Y];
    siny = sin(y); cosy = cos(y);
    sinhx = sinh(x); coshx = cosh(x);
    denom = coshx + cosy;
    if (denom == 0)
        return;
    rTarget[X] = sinhx / denom;
    rTarget[Y] = siny / denom;
}

void Tanh_Deriv(Complex rTarget, const Complex source, const Complex)
{
    DBL x, y, siny, cosy, sinhx, coshx, denom;
    x = 2 * source[X];
    y = 2 * source[Y];
    siny = sin(y); cosy = cos(y);
    sinhx = sinh(x); coshx = cosh(x);
    denom = pov::Sqr(coshx + cosy) / 2;
    if (denom == 0)
        return;
    rTarget[X] = (coshx * cosy + 1) / denom;
    rTarget[Y] = -sinhx * siny / denom;
}

void Pwr(Complex rTarget, const Complex source1, const Complex source2)
{
    Complex cLog, t;
    DBL e2x;

    if (source1[X] == 0 && source1[Y] == 0)
    {
        rTarget[X] = rTarget[Y] = 0.0;
        return;
    }

    Ln(cLog, source1);
    Mult(t, cLog, source2);

    if (t[X] < -690)
        e2x = 0;
    else
        e2x = exp(t[X]);
    rTarget[X] = e2x * cos(t[Y]);
    rTarget[Y] = e2x * sin(t[Y]);
}

void Pwr_Deriv(Complex rTarget, const Complex source1, const Complex source2)
{
    Complex cLog, t, aMul, temp;
    DBL e2x;

    if (source1[X] == 0 && source1[Y] == 0)
    {
        rTarget[X] = rTarget[Y] = 0.0;
        return;
    }

    Ln(cLog, source1);
    aMul[X] = source2[X] - 1; aMul[Y] = source2[Y];
    Mult(t, cLog, aMul);

    if (t[X] < -690)
        e2x = 0;
    else
        e2x = exp(t[X]);
    temp[X] = e2x * cos(t[Y]);
    temp[Y] = e2x * sin(t[Y]);

    Mult(rTarget, temp, source2);
}

bool False_DTest(Complex, DBL&, const Complex, const Complex, const Complex)
{
    return false;
}

bool NegReal_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused)
{
    DBL tmp;
    if (((cT[Y] >= 0.0) ^ (cP[Y] >= 0.0)) && (cT[X] < 0.0 || cP[X] < 0.0))
    {
        tmp = -cP[Y] / (cT[Y] - cP[Y]);
        if ((cT[X] < 0.0 && cP[X] < 0.0) || ((cT[X] - cP[X]) * tmp + cP[X] < 0.0))
        {
            rDist = tmp;
            rNormal[X] = 0;
            rNormal[Y] = (cT[Y] >= 0.0 ? -1.0 : 1.0);
            return true;
        }
    }
    return false;
}

bool ASin_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused)
{
    DBL tmp;
    if (((cT[Y] >= 0.0) ^ (cP[Y] >= 0.0)) && (fabs(cT[X]) > 1.0 || fabs(cP[X]) > 1.0))
    {
        tmp = -cP[Y] / (cT[Y] - cP[Y]);
        if (fabs((cT[X] - cP[X]) * tmp + cP[X]) > 1.0)
        {
            rDist = tmp;
            rNormal[X] = 0;
            rNormal[Y] = (cT[Y] >= 0.0 ? -1.0 : 1.0);
            return true;
        }
    }
    return false;
}

bool ASinh_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused)
{
    DBL tmp;
    if (((cT[X] >= 0.0) ^ (cP[X] >= 0.0)) && (fabs(cT[Y]) > 1.0 || fabs(cP[Y]) > 1.0))
    {
        tmp = -cP[X] / (cT[X] - cP[X]);
        if (fabs((cT[Y] - cP[Y]) * tmp + cP[Y]) > 1.0)
        {
            rDist = tmp;
            rNormal[X] = (cT[X] >= 0.0 ? -1.0 : 1.0);
            rNormal[Y] = 0;
            return true;
        }
    }
    return false;
}

bool ACos_Alt_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused)
{
    DBL tmp;
    if ((cT[X] >= 0.0) ^ (cP[X] >= 0.0))
    {
        rDist = -cP[X] / (cT[X] - cP[X]);
        rNormal[X] = (cT[X] >= 0.0 ? -1.0 : 1.0);
        rNormal[Y] = 0;
    }
    else if (((cT[Y] >= 0.0) ^ (cP[Y] >= 0.0)) && (cT[X] < 1.0 || cP[X] < 1.0))
    {
        tmp = -cP[Y] / (cT[Y] - cP[Y]);
        if ((cT[X] - cP[X]) * tmp + cP[X] < 1.0)
        {
            rDist = tmp;
            rNormal[X] = 0;
            rNormal[Y] = (cT[Y] >= 0.0 ? -1.0 : 1.0);
            return true;
        }
    }
    return false;
}

bool ACosh_DTest(Complex rNormal, DBL& rDist, const Complex cT, const Complex cP, const Complex unused)
{
    DBL tmp;
    if (((cT[Y] >= 0.0) ^ (cP[Y] >= 0.0)) && (cT[X] < 1.0 || cP[X] < 1.0))
    {
        tmp = -cP[Y] / (cT[Y] - cP[Y]);
        if ((cT[X] - cP[X]) * tmp + cP[X] < 1.0)
        {
            rDist = tmp;
            rNormal[X] = 0;
            rNormal[Y] = (cT[Y] >= 0.0 ? -1.0 : 1.0);
            return true;
        }
    }
    return false;
}

}

}
