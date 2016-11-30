//******************************************************************************
///
/// @file core/shape/fractal.cpp
///
/// Implementation of the fractal set geometric primitives.
///
/// @author Pascal Massimino
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

// Unit header file must be the first file included within POV-Ray *.cpp files (pulls in config)
#include "core/shape/fractal.h"

#include "base/pov_err.h"

#include "core/bounding/boundingbox.h"
#include "core/math/complexfn.h"
#include "core/math/matrix.h"
#include "core/math/vector.h"
#include "core/render/ray.h"
#include "core/scene/tracethreaddata.h"
#include "core/shape/fractal/dispatch.h"
#include "core/shape/fractal/space.h"
#include "core/shape/fractal/types.h"
#include "core/shape/fractal/util.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
* Local preprocessor defines
******************************************************************************/

/*****************************************************************************
* Local variables
******************************************************************************/

static const Vector4d kDummyVector;

const DBL Fractal_Tolerance = 1e-7;

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

#define SWAP(a,b) tmpStack = a; a = b; b = tmpStack

bool Fractal::All_Intersections(const Ray& ray, IStack& Depth_Stack, TraceThreadData *Thread)
{
    bool Intersection_Found;
    bool LastIsInside = false;
    bool CurrentIsInside, NextIsInside;
    DBL Depth, Depth_Max;
    DBL Dist, Dist_Next, LenSqr, LenInv;

    Vector4d IPoint, Mid_Point, Next_Point;
    Vector3d Real_Normal, F_Normal;
    Vector3d Direction0, IPoint0, Real_Pt;
    Vector4d Direction;
    BasicRay New_Ray;
    int cIter, nIter, lIter, tIter;
    FractalIterData *cStack = &(Thread->Fractal_IterData[0]), *nStack = &(Thread->Fractal_IterData[1]),
        *lStack = &(Thread->Fractal_IterData[2]), *tmpStack, *tStack, *pStack;

    Thread->Stats()[Ray_Fractal_Tests]++;

    if (Test_Flag(this, DEGENERATE_FLAG))
        return false;

    Intersection_Found = false;

    /* Get into Fractal's world. */

    if (Trans != NULL)
    {
        MInvTransDirection(Direction0, ray.Direction, Trans);
        LenSqr = Direction0.lengthSqr();

        if (LenSqr == 0.0)
        {
            return (false);
        }

        if (LenSqr != 1.0)
        {
            LenInv = 1.0 / sqrt(LenSqr);
            Direction0 *= LenInv;
        }
        else
            LenInv = 1.0;

        New_Ray.Direction = Direction0;
        MInvTransPoint(New_Ray.Origin, ray.Origin, Trans);
    }
    else
    {
        Direction0 = ray.Direction;
        New_Ray = ray;
        LenInv = 1.0;
    }

    Direction = RulesSpace->TransformDirTo4D(Direction0);

    /* Bound fractal. */

    if (!RulesSpace->Bound(New_Ray, this, &Depth, &Depth_Max))
    {
        return (false);
    }

    if (Depth_Max < Fractal_Tolerance)
    {
        return (false);
    }

    if (Depth < Fractal_Tolerance)
    {
        Depth = Fractal_Tolerance;
    }

    /* Jump to starting point */

    Next_Point = RulesSpace->TransformTo4D(New_Ray.Origin) + Direction * Depth;

    cIter = Rules->Iterate(Next_Point, this, Direction, &Dist, cStack);

    CurrentIsInside = (cIter == Num_Iterations + 1);

    /* Light ray starting inside ? */

    if (CurrentIsInside)
    {
        Next_Point += (2.0 * Fractal_Tolerance) * Direction;

        Depth += 2.0 * Fractal_Tolerance;

        if (Depth > Depth_Max)
        {
            return (false);
        }

        cIter = Rules->Iterate(Next_Point, this, Direction, &Dist, cStack);

        CurrentIsInside = (cIter == Num_Iterations + 1);
    }

    /* Ok. Trace it */

    while (Depth <= Depth_Max)
    {
        /*
         * Get close to the root: Advance with Next_Point, keeping track of last
         * position in IPoint...
         */

        while (1)
        {
            if (Depth == Depth_Max)
            {
                if (Intersection_Found)
                    Thread->Stats()[Ray_Fractal_Tests_Succeeded]++;
                return (Intersection_Found);
            }

            if (Dist < Precision)
                Dist = Precision;

            /* Make sure the outer edge gets checked. */
            if (Depth + Dist > Depth_Max)
            {
                Dist = Depth_Max - Depth;
                Depth = Depth_Max;
            }
            else
            {
                Depth += Dist;
            }

            IPoint = Next_Point;
            Next_Point += Dist * Direction;

            nIter = Rules->Iterate(Next_Point, this, Direction, &Dist_Next, nStack);

            NextIsInside = (nIter == Num_Iterations + 1);

            if (NextIsInside != CurrentIsInside)
            {
                /* Set surface was crossed... */

                Depth -= Dist;
                break;
            }
            else
            {
                cIter = nIter;
                SWAP(cStack, nStack);
                Dist = Dist_Next; /* not reached */
            }
        }

        /* then, polish the root via bisection method... */

        while (Dist > Fractal_Tolerance)
        {
            Dist *= 0.5;
            Mid_Point = IPoint + Dist * Direction;

            lIter = Rules->Iterate(Mid_Point, this, kDummyVector, NULL, lStack);

            LastIsInside = (lIter == Num_Iterations + 1);

            if (LastIsInside == CurrentIsInside)
            {
                IPoint = Mid_Point;

                Depth += Dist;

                cIter = lIter;
                SWAP(cStack, lStack);

                if (Depth > Depth_Max)
                {
                    if (Intersection_Found)
                        Thread->Stats()[Ray_Fractal_Tests_Succeeded]++;
                    return (Intersection_Found);
                }
            }
            else
            {
                nIter = lIter;
                SWAP(nStack, lStack);
            }
        }

        if (!CurrentIsInside) /* IPoint is outside */
        {
            Depth += Dist;

            tIter = cIter;

            pStack = cStack;
            tStack = nStack;

            if (!LastIsInside) /* Mid_Point == IPoint */
            {
                IPoint += Dist * Direction;

                Rules->Iterate(IPoint, this, kDummyVector, NULL, nStack);
            }
            else
            {
                IPoint = Mid_Point;
            }
        }
        else
        {
            tIter = nIter;

            pStack = nStack;
            tStack = cStack;

            if (!LastIsInside) /* Mid_Point isn't inside the set */
            {
                Rules->Iterate(IPoint, this, kDummyVector, NULL, cStack);
            }
        }

        IPoint0 = New_Ray.Origin + Direction0 * Depth;

        if (Trans != NULL)
        {
            MTransPoint(Real_Pt, IPoint0, Trans);
            Rules->CalcNormal(F_Normal, tIter, this, tStack, pStack);
            MTransNormal(Real_Normal, F_Normal, Trans);
        }
        else
        {
            Real_Pt = IPoint0;
            Rules->CalcNormal(Real_Normal, tIter, this, tStack, pStack);
        }

        if (Clip.empty() || Point_In_Clip(Real_Pt, Clip, Thread))
        {
            Real_Normal.normalize();
            Depth_Stack->push(Intersection(Depth * LenInv, Real_Pt, Real_Normal, this));
            Intersection_Found = true;

            /* If fractal isn't used with CSG we can exit now. */

            if (!(Type & IS_CHILD_OBJECT))
            {
                break;
            }
        }

        /* Start over where work was left */

        IPoint = Next_Point;
        Dist = Dist_Next;
        CurrentIsInside = NextIsInside;

    }

    if (Intersection_Found)
        Thread->Stats()[Ray_Fractal_Tests_Succeeded]++;
    return (Intersection_Found);
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

bool Fractal::Inside(const Vector3d& IPoint, TraceThreadData *Thread) const
{
    bool Result;
    Vector3d New_Point;

    if (Test_Flag(this, DEGENERATE_FLAG))
        return Test_Flag(this, INVERTED_FLAG);

    if (Trans != NULL)
    {
        MInvTransPoint(New_Point, IPoint, Trans);

        Result = (Rules->Iterate(RulesSpace->TransformTo4D(New_Point), this, kDummyVector, NULL,
                                 &(Thread->Fractal_IterData[0])) == Num_Iterations + 1);
    }
    else
    {
        Result = (Rules->Iterate(RulesSpace->TransformTo4D(IPoint), this, kDummyVector, NULL,
                                 &(Thread->Fractal_IterData[0])) == Num_Iterations + 1);
    }

    if (Test_Flag(this, INVERTED_FLAG))
    {
        return (!Result);
    }
    else
    {
        return (Result);
    }
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

void Fractal::Normal(Vector3d& Result, Intersection *Intersect, TraceThreadData *) const
{
    Result = Intersect->INormal;
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

void Fractal::Translate(const Vector3d&, const TRANSFORM *tr)
{
    Transform(tr);
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

void Fractal::Rotate(const Vector3d&, const TRANSFORM *tr)
{
    Transform(tr);
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

void Fractal::Scale(const Vector3d&, const TRANSFORM *tr)
{
    Transform(tr);
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*   Mar 1996 : Moved call to Recompute_BBox to Compute_Fractal_BBox() (TW)
*
******************************************************************************/

void Fractal::Transform(const TRANSFORM *tr)
{
    if(Trans == NULL)
        Trans = Create_Transform();

    Compose_Transforms(Trans, tr);

    Compute_BBox();
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*   Mar 1996 : Added call to recompute_BBox() to bottom (TW)
*
******************************************************************************/

void Fractal::Compute_BBox()
{
    if (!RulesSpace)
        return;

    if (!RulesSpace->Compute_BBox(BBox, this))
    {
        ;// TODO MESSAGE        Warning("Degenerate julia_fractal.");

        // This is basically superfluous, but it doesn't hurt.
        Set_Flag(this, DEGENERATE_FLAG);

        // This should cause any bounding checks to fail quickly.
        // (Is this okay to do?)
        Make_BBox(BBox, BOUND_HUGE, BOUND_HUGE, BOUND_HUGE,
                  -2.0 * BOUND_HUGE, -2.0 * BOUND_HUGE, -2.0 * BOUND_HUGE);
    }
    else
    {
            Recompute_BBox(&BBox, Trans);
    }
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

Fractal::Fractal() : ObjectBase(BASIC_OBJECT)
{
    Trans = NULL;

    Center = Vector3d(0.0, 0.0, 0.0);

    Julia_Parm = Vector4d(1.0, 0.0, 0.0, 0.0);

    Slice = Vector4d(0.0, 0.0, 0.0, 1.0);
    SliceDist = 0.0;
    TransformMethod = kTransformProjection;

    Bailout = 0.0;
    Exit_Value = 0.0;

    Num_Iterations = 20;

    Precision = 1.0 / 20.0;

    Discontinuity_Test = -1;

    Distance_Estimator = kDefaultEstimator;

    Jump_Max = 30.0;
    Jump_Max_Lower = 1.0;
    Jump_Decay = -1.0;
    Jump_Min = 2.0;

    Func_Type.algebra = kQuaternion;
    Func_Type.type = kFunc_Sqr;
    Func_Type.variant = kVar_Normal;

    Rules.reset();
    RulesSpace.reset();

    Radius_Squared = 0.0;
    exponent[X] = 0.0;
    exponent[Y] = 0.0;

    InitDispatch();
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

ObjectPtr Fractal::Copy()
{
    Fractal *New = new Fractal();
    Destroy_Transform(New->Trans);
    *New = *this;
    New->Trans = Copy_Transform(Trans);
    New->Rules = Rules;
    New->RulesSpace = RulesSpace;

    return (New);
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

Fractal::~Fractal()
{}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Pascal Massimino
*
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Dec 1994 : Creation.
*
******************************************************************************/

int Fractal::SetUp_Fractal()
{
    FractalConstructorData ctorData;
    DBL R;

    RulesSpace = static_cast<FractalSpacePtr>(new FractalSpace(TransformMethod, Func_Type.algebra, Slice, SliceDist));

    ctorData.juliaParm = Julia_Parm;
    ctorData.estimatorType = Distance_Estimator;
    ctorData.funcType = Func_Type;
    AssignComplex(ctorData.exponent, exponent);
    ctorData.space = RulesSpace;

    /* ... And this is [one reason] why all that dispatch stuff is nice! */
    Rules = RulesDispatch::CreateNew(ctorData);

    if (Discontinuity_Test < 0)
    {
        Discontinuity_Test = 1;
    }

    if (Discontinuity_Test > 0)
    {
        if (Rules->Info().discontinuitySupport == kDiscontinuityNotImplemented)
        {
            // Currently, this warning gets issued by the parser.
            ;// TODO MESSAGE        Warning("Discontinuity testing not supported for fractal type.");

            Discontinuity_Test = -1;
        }
        else if (Rules->Info().discontinuitySupport == kDiscontinuityUnneeded)
        {
            // Discontinuity testing would have no effect, so we turn it off for speed's sake.
            Discontinuity_Test = 0;
        }
    }

    // The automatic settings for the Newton estimator are rather ad hoc...
    if (Rules->Info().estimatorType == kNewtonEstimator && Jump_Decay < 0.0)
    {
        if (pow(0.75, Num_Iterations) * Jump_Max > Jump_Max_Lower)
        {
            Jump_Decay = pow(Jump_Max / Jump_Max_Lower, -1.0 / Num_Iterations);
        }
        else
        {
            Jump_Decay = 0.75;
        }
    }

    if (Bailout > 0.0)
    {
        R = Bailout;
    }
    else if (Func_Type.algebra == kQuaternion && (Func_Type.type == kFunc_Sqr || Func_Type.type == kFunc_Cube))
    {
        R = 1.0 + Julia_Parm.length();
        if (R > 2.0)
        {
            R = 2.0;
        }

        Bailout = R;
    }
    else
    {
        Bailout = R = 4.0;
    }

    Exit_Value = Sqr(Bailout);

    /* To make sure the outside of the fractal doesn't get cut off. */
    R += 2 * Fractal_Tolerance;

    Radius_Squared = Sqr(R);

    Compute_BBox();

    return Num_Iterations;
}

const FractalDataSizes& Fractal::IterationDataSizes() const
{
    return Rules->Info().sizes;
}

}
