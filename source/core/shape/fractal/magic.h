//******************************************************************************
///
/// @file core/shape/fractal/magic.h
///
/// This module contains prototypes for the generic implementation of
/// FractalRules subclasses in magic.cpp.
///
/// (This was originally named for the templating 'magic' that it used, but
/// could equally well be described as the 'magic' that makes everything work!)
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

#ifndef POVRAY_CORE_FRACTAL_MAGIC_H
#define POVRAY_CORE_FRACTAL_MAGIC_H

#include "core/coretypes.h"
#include "core/math/vector.h"
#include "core/shape/fractal/distestimator.h"
#include "core/shape/fractal/types.h"
#include "core/shape/fractal/util.h"

namespace pov
{

class MagicFractalRulesBase : public FractalRules
{
public:
    typedef NilFractalData FixedData;
    typedef NilFractalData MainIterData;
    typedef NilFractalData AuxIterData;

    MagicFractalRulesBase(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                          const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        mSpace4D(*(data.space)), mEstimator(estimator),
        mInfo(CreateRulesInfo(data.funcType, discontinuitySupport, sizes, estimator.eType)) { }
    virtual const FractalRulesInfo& Info() const { return mInfo; }
    virtual const FractalSpace& Space() const { return mSpace4D; }

    static const DistanceEstimator& GetEstimatorFromType(EstimatorType estimatorType, EstimatorType defaultEstimator = kNewtonEstimator,
                                                         EstimatorType legacyEstimator = kOrigNewtonEstimator,
                                                         const DistanceEstimator& (*ExtraEstimators)(EstimatorType eType) = NULL);

protected:
    const FractalSpace& mSpace4D;
    const FractalRulesInfo mInfo;
    const DistanceEstimator& mEstimator;

};

/// Abstract class providing the standard implementation of the fractal iteration mechanisms.
///
/// Currently, for this class to be correct mathematically, the following must hold:
///   - The underlying vector space must be @f$\mathbb{R}^4@f$ (or a subspace thereof),
///   - The bailout condition must be @f$\|f^n(v)\|_2^2<\mathtt{Exit\_Value}@f$, and
///   - The functions IterateCalc(), GradientCalc(), and DiscontinuityCheck() must perform the functions
///     specified in their documentation.
///   .
class MagicFractalRules : public MagicFractalRulesBase
{
public:
    typedef Vector4d MainIterData;

    MagicFractalRules(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                      const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        MagicFractalRulesBase(data, discontinuitySupport, sizes, estimator) { }

    virtual int Iterate(const Vector4d& iPoint, const Fractal *pFractal, const Vector4d& direction,
                        DBL *pDist, FractalIterData *pIterData) const;
    virtual DBL CalcDirDeriv(const Vector4d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData,
                            FractalIterData *pPIterData) const;

    /// Performs one iteration of the function to be iterated, i.e., @f$\mathtt{rV}\gets f(\mathtt{rV})@f$.
    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    /// Performs one iteration of the gradient calculation, left-multiplying @c rD by [a positive multiple of] the transpose of
    /// the function's Jacobian at @c v.
    ///
    /// Specifically, @f$\mathtt{rD}\gets cJ_f(\mathtt{v})^T\mathtt{rD}@f$ and @f$\mathtt{rMult}\gets\frac{1}{c}\mathtt{rMult}@f$
    /// (with @f$c>0@f$).
    virtual void GradientCalc(Vector4d& rD, const Vector4d& v, int iter, DBL& rMult, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    /// Returns a bool indicating whether a path from @c p to @c t crosses a discontinuity of the function being iterated.
    ///
    /// If it does, it stores the normal of the discontinuity's surface at the intersection (oriented towards @c t)
    /// in @c rD.
    ///
    /// Also, it stores the distance from @c p, relative to the distance to @c t, to the intersection in @c rDist.
    /// **Note, however, that this latter property is unused by the current implementation, and may change or be removed.**
    virtual bool DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;
};

}

#endif
