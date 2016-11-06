//******************************************************************************
///
/// @file core/math/complexfn.cpp
///
/// This module implements complex functions, their derivatives, and discontinuity
/// tests.  Most are merely hand-offs to the appropriate STL functions.
///
/// @todo   It might be a good thing to replace much of this with C++03 and
/// TR1 ~= C99 ~= boost/math/complex functions. -- In progress...
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

void Sqr(Complex& rTarget, const Complex& source, const Complex&)
{
#ifdef COMPLEXFN_MANUAL_SQR
    rTarget = Complex(pov::Sqr(source.real()) - pov::Sqr(source.imag()),
                      2.0 * source.real() * source.imag());
#else
    rTarget = source * source;
#endif
}

void Sqrt(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = sqrt(source);
}

void Ln(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = log(source);
}

void Recip(Complex& rTarget, const Complex& source, const Complex&)
{
#ifdef COMPLEXFN_MANUAL_RECIP
    rTarget = Conj(source) / Norm(source);
#else
    rTarget = 1.0 / source;
#endif
}

void Exp(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = exp(source);
}

void Sin(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = sin(source);
}

void Sinh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = sinh(source);
}

void Cos(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = cos(source);
}

void Cos_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = -sin(source);
}

void Cosh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = cosh(source);
}

void ASin(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = asin(source);
}

void ASin_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp = 1.0 - temp;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

void ACos(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = acos(source);
}

/* There are some oddities with this function, due to branch cuts. */
void ACos_Alt(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;

    Sqr(temp, source);
    temp -= 1.0;
    Sqrt(temp, temp);

    temp += source;

    Ln(temp, temp);
    rTarget = -kImag_Unit * temp;
}

void ACos_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp -= 1.0;
    Sqrt(temp, temp);

    Recip(temp, temp);
    rTarget = -kImag_Unit * temp;
}

void ASinh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = asinh(source);
}

void ASinh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp += 1.0;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

void ACosh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = acosh(source);
}

/* There are some oddities with this function, due to branch cuts. */
void ACosh_Alt(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp -= 1.0;
    Sqrt(temp, temp);
    temp += source;
    Ln(rTarget, temp);
}

void ACosh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp -= 1.0;
    Sqrt(temp, temp);

    Recip(rTarget, temp);
}

void ATanh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = atanh(source);
}

void ATanh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp = 1.0 - temp;
    Recip(rTarget, temp);
}

void ATan(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = atan(source);
}

void ATan_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Sqr(temp, source);
    temp += 1.0;
    Recip(rTarget, temp);
}

void Tan(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = tan(source);
}

void Tan_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Cos(temp, source);
    Sqr(temp, temp);
    temp *= -1.0;
    Recip(rTarget, temp);
}

void Tanh(Complex& rTarget, const Complex& source, const Complex&)
{
    rTarget = tanh(source);
}

void Tanh_Deriv(Complex& rTarget, const Complex& source, const Complex&)
{
    Complex temp;
    Cosh(temp, source);
    Sqr(temp, temp);
    Recip(rTarget, temp);
}

void Pwr(Complex& rTarget, const Complex& source1, const Complex& source2)
{
    rTarget = pow(source1, source2);
}

void Pwr_Deriv(Complex& rTarget, const Complex& source1, const Complex& source2)
{
    Complex temp;
    temp = source2 - 1.0;
    Pwr(temp, source1, temp);
    rTarget = source2 * temp;
}

bool False_DTest(Complex&, DBL&, const Complex&, const Complex&, const Complex&)
{
    return false;
}

bool NegReal_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if (((cT.imag() >= 0.0) ^ (cP.imag() >= 0.0)) && ((cT.real() < 0.0) || (cP.real() < 0.0)))
    {
        tmp = -cP.imag() / (cT.imag() - cP.imag());
        if (((cT.real() < 0.0) && (cP.real() < 0.0)) || ((cT.real() - cP.real()) * tmp + cP.real() < 0.0))
        {
            rDist = tmp;
            rNormal = Complex(0.0, ((cT.imag() < 0.0) ? 1.0 : -1.0));
            return true;
        }
    }
    return false;
}

bool ASin_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if (((cT.imag() >= 0.0) ^ (cP.imag() >= 0.0)) && (fabs(cT.real()) > 1.0 || fabs(cP.real()) > 1.0))
    {
        tmp = -cP.imag() / (cT.imag() - cP.imag());
        if (fabs((cT.real() - cP.real()) * tmp + cP.real()) > 1.0)
        {
            rDist = tmp;
            rNormal = Complex(0.0, ((cT.imag() < 0.0) ? 1.0 : -1.0));
            return true;
        }
    }
    return false;
}

bool ASinh_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if (((cT.real() >= 0.0) ^ (cP.real() >= 0.0)) && (fabs(cT.imag()) > 1.0 || fabs(cP.imag()) > 1.0))
    {
        tmp = -cP.real() / (cT.real() - cP.real());
        if (fabs((cT.imag() - cP.imag()) * tmp + cP.imag()) > 1.0)
        {
            rDist = tmp;
            rNormal = Complex(((cT.real() < 0.0) ? 1.0 : -1.0), 0.0);
            return true;
        }
    }
    return false;
}

bool ACos_Alt_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if ((cT.real() >= 0.0) ^ (cP.real() >= 0.0))
    {
        rDist = -cP.real() / (cT.real() - cP.real());
        rNormal = Complex(((cT.real() < 0.0) ? 1.0 : -1.0), 0.0);
    }
    else if (((cT.imag() >= 0.0) ^ (cP.imag() >= 0.0)) && (cT.real() < 1.0 || cP.real() < 1.0))
    {
        tmp = -cP.imag() / (cT.imag() - cP.imag());
        if ((cT.real() - cP.real()) * tmp + cP.real() < 1.0)
        {
            rDist = tmp;
            rNormal = Complex(0.0, ((cT.imag() < 0.0) ? 1.0 : -1.0));
            return true;
        }
    }
    return false;
}

bool ACosh_DTest(Complex& rNormal, DBL& rDist, const Complex& cT, const Complex& cP, const Complex& unused)
{
    DBL tmp;
    if (((cT.imag() >= 0.0) ^ (cP.imag() >= 0.0)) && (cT.real() < 1.0 || cP.real() < 1.0))
    {
        tmp = -cP.imag() / (cT.imag() - cP.imag());
        if ((cT.real() - cP.real()) * tmp + cP.real() < 1.0)
        {
            rDist = tmp;
            rNormal = Complex(0.0, ((cT.imag() < 0.0) ? 1.0 : -1.0));
            return true;
        }
    }
    return false;
}

}

}
