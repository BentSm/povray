//******************************************************************************
///
/// @file core/shape/fractal.h
///
/// Declarations related to the fractal set geometric primitives.
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

#ifndef POVRAY_CORE_FRACTAL_H
#define POVRAY_CORE_FRACTAL_H

// Module config header file must be the first file included within POV-Ray unit header files
#include "core/configcore.h"

#include "core/math/matrix.h"
#include "core/math/vector.h"
#include "core/scene/object.h"
#include "core/shape/fractal/types.h"

namespace pov
{

//##############################################################################
///
/// @addtogroup PovCoreShape
///
/// @{

/*****************************************************************************
* Global preprocessor defines
******************************************************************************/

/*****************************************************************************
* Global typedefs
******************************************************************************/

class Fractal;

class Fractal : public ObjectBase
{
public:
    Vector3d Center;
    Vector4d Julia_Parm;
    Vector4d Slice;               /* vector perpendicular to slice plane */
    DBL SliceDist;                /* distance from slice plane to origin */
    FractalTransformMethod TransformMethod;
    DBL Bailout;
    DBL Exit_Value;
    int Num_Iterations;           /* number of iterations */
    DBL Precision;                /* Precision value */
    int Discontinuity_Test;
    EstimatorType Distance_Estimator;
    DBL Jump_Max;
    DBL Jump_Max_Lower;
    DBL Jump_Decay;
    DBL Jump_Min;
    FractalFuncType Func_Type;
    Complex exponent;             /* exponent of power function */
    DBL Radius_Squared;           /* For F_Bound(), if needed */
    FractalRulesPtr Rules;
    FractalSpacePtr RulesSpace;

    Fractal();
    virtual ~Fractal();

    virtual ObjectPtr Copy();

    virtual bool All_Intersections(const Ray&, IStack&, TraceThreadData *);
    virtual bool Inside(const Vector3d&, TraceThreadData *) const;
    virtual void Normal(Vector3d&, Intersection *, TraceThreadData *) const;
    virtual void Translate(const Vector3d&, const TRANSFORM *);
    virtual void Rotate(const Vector3d&, const TRANSFORM *);
    virtual void Scale(const Vector3d&, const TRANSFORM *);
    virtual void Transform(const TRANSFORM *);
    virtual void Compute_BBox();

    static const int kNumIterStacks = 3;

    int SetUp_Fractal();
    const FractalDataSizes& IterationDataSizes() const;

};

/// @}
///
//##############################################################################

}

#endif // POVRAY_CORE_FRACTAL_H
