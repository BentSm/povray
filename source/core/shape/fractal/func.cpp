//******************************************************************************
///
/// @file core/shape/fractal/func.cpp
///
/// This module contains the function structs for fractals.
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
#include "core/shape/fractal/func.h"

#include <list>
#include <utility>

#include "base/pov_err.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

using namespace complex_fn;

static const std::map<FractalFuncType, FractalFunc> makeMap(const FractalFunc array[], int len)
{
    std::list<std::pair<FractalFuncType, FractalFunc> > vec;
    int i;
    for (i = 0; i < len; i++)
        vec.push_back(std::make_pair(array[i].type, array[i]));
    return std::map<FractalFuncType, FractalFunc>(vec.begin(), vec.end());
}

static const FractalFunc kFractalFuncArray[] =
{
    { { kComplex,  kFunc_Exp,    kVar_Normal    },  Exp,        Exp,          False_DTest    },
    { { kComplex,  kFunc_Ln,     kVar_Normal    },  Ln,         Recip,        NegReal_DTest  },
    { { kComplex,  kFunc_Sin,    kVar_Normal    },  Sin,        Cos,          False_DTest    },
    { { kComplex,  kFunc_ASin,   kVar_Normal    },  ASin,       ASin_Deriv,   ASin_DTest     },
    { { kComplex,  kFunc_Cos,    kVar_Normal    },  Cos,        Cos_Deriv,    False_DTest    },
    { { kComplex,  kFunc_ACos,   kVar_Normal    },  ACos,       ACos_Deriv,   ASin_DTest     },
    { { kComplex,  kFunc_ACos,   kVar_Alternate },  ACos_Alt,   ACos_Deriv,   ACos_Alt_DTest },
    { { kComplex,  kFunc_Tan,    kVar_Normal    },  Tan,        Tan_Deriv,    False_DTest    }, /* See note */
    { { kComplex,  kFunc_ATan,   kVar_Normal    },  ATan,       ATan_Deriv,   ASinh_DTest    },
    { { kComplex,  kFunc_Sinh,   kVar_Normal    },  Sinh,       Cosh,         False_DTest    },
    { { kComplex,  kFunc_ASinh,  kVar_Normal    },  ASinh,      ASinh_Deriv,  ASinh_DTest    },
    { { kComplex,  kFunc_Cosh,   kVar_Normal    },  Cosh,       Sinh,         False_DTest    },
    { { kComplex,  kFunc_ACosh,  kVar_Normal    },  ACosh,      ACosh_Deriv,  ASin_DTest     },
    { { kComplex,  kFunc_ACosh,  kVar_Alternate },  ACosh_Alt,  ACosh_Deriv,  ACos_Alt_DTest },
    { { kComplex,  kFunc_Tanh,   kVar_Normal    },  Tanh,       Tanh_Deriv,   False_DTest    }, /* See note */
    { { kComplex,  kFunc_ATanh,  kVar_Normal    },  ATanh,      ATanh_Deriv,  ASin_DTest     },
    { { kComplex,  kFunc_Pwr,    kVar_Normal    },  Pwr,        Pwr_Deriv,    NegReal_DTest  }
};

/* Note: Certain of above (namely, tan and tanh) have discontinuities, but still
   use the always-false discontinuity test.  This is intentional.  (See the next
   comment if you are interested in a somewhat-technical explanation.) */

/* Somewhat-technical explanation: The functions only have point discontinuities
   where the modulus of the function value approaches infinity as the test point
   approaches the discontinuity.  Thus the discontinuities cannot be on the
   boundary of the fractal (since the bailout would always be reached "prior to"
   "reaching" a discontinuity).  [More precisely, for a fixed bailout and specific
   discontinuity, there is a disc of positive radius centered at the
   discontinuity such that the modulus of the function exceeds the bailout value
   at every point of the (punctured) disc.] */

const std::map<FractalFuncType, FractalFunc>& AllFractalFuncs() {
    static const std::map<FractalFuncType, FractalFunc> kFractalFuncs = makeMap(kFractalFuncArray, 17);
    return kFractalFuncs;
}

const FractalFunc& FractalFuncForType(const FractalFuncType& funcType) {
    std::map<FractalFuncType, FractalFunc>::const_iterator funcPairIter;
    FractalFuncType cType = funcType;
    cType.algebra = kComplex;
    funcPairIter = AllFractalFuncs().find(cType);
    if (funcPairIter == AllFractalFuncs().end())
        throw POV_EXCEPTION_STRING("Function type/variant combination unknown in fractal.");
    return funcPairIter->second;
}

}
