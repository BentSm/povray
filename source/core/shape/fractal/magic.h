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

class MagicFractalSpace : public FractalSpace
{
public:
    MagicFractalSpace(FractalTransformMethod transformMethod, FractalAlgebra algebra,
                      const Vector4d& slice, DBL sliceDist);
    const Vector4d& transformedX() const { return tX; }
    const Vector4d& transformedY() const { return tY; }
    const Vector4d& transformedZ() const { return tZ; }
    const Vector4d& transformed0() const { return t0; }

    const Vector4d TransformTo4D(const Vector3d& point) const;
    const Vector4d TransformDirTo4D(const Vector3d& point) const;
    bool Bound(const BasicRay& ray, const Fractal *pFractal, DBL *pDepthMin, DBL *pDepthMax) const;
    bool Compute_BBox(BoundingBox& BBox, const Fractal *pFractal) const;

protected:
    Vector4d tX, tY, tZ, t0, sliceNorm;
    Vector3d sliceNorm3d;
    DBL sliceDistNorm;
    FractalTransformMethod mTransformMethod;
    FractalAlgebra mAlgebra;
};

class MagicRulesBase : public FractalRules
{
public:
    typedef NilFractalData FixedData;
    typedef NilFractalData MainIterData;
    typedef NilFractalData AuxIterData;

    MagicRulesBase(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                   const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        mSpace4D(data.transformMethod, data.funcType.algebra, data.slice, data.sliceDist),
        mInfo(CreateRulesInfo(data.funcType, discontinuitySupport, sizes, estimator.eType)),
        mEstimator(estimator) { }
    virtual const FractalRulesInfo& Info() const { return mInfo; }
    virtual const FractalSpace *GetSpace() const { return &mSpace4D; }

protected:
    const MagicFractalSpace mSpace4D;
    const FractalRulesInfo mInfo;
    const DistanceEstimator& mEstimator;

    static const DistanceEstimator& GetEstimatorFromType(EstimatorType estimatorType, EstimatorType defaultEstimator = kNewtonEstimator,
                                                         EstimatorType legacyEstimator = kOrigNewtonEstimator,
                                                         const DistanceEstimator& (*ExtraEstimators)(EstimatorType eType) = NULL);
};

class MagicQuaternionFractalRules : public MagicRulesBase
{
public:
    typedef Vector4d MainIterData;

    MagicQuaternionFractalRules(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                                const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        mJuliaParm(data.juliaParm), MagicRulesBase(data, discontinuitySupport, sizes, estimator) {}

    virtual int Iterate(const Vector4d& iPoint, const Fractal *pFractal, const Vector4d& direction,
                        DBL *pDist, FractalIterData *pIterData) const;
    virtual DBL CalcDirDeriv(const Vector4d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData,
                            FractalIterData *pPIterData) const;

    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual void DirDerivCalc(Vector4d& rD, const Vector4d& v, int iter, bool samePoint, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual bool DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;

protected:
    const Vector4d mJuliaParm;

};

class MagicHypercomplexFractalRules : public MagicRulesBase
{
public:
    typedef Vector4d MainIterData;

    MagicHypercomplexFractalRules(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                                  const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        mDuplexJuliaParm(DuplexFromHypercomplex(data.juliaParm)),
        MagicRulesBase(data, discontinuitySupport, sizes, estimator) {}

    virtual int Iterate(const Vector4d& iPoint, const Fractal *pFractal, const Vector4d& direction,
                        DBL *pDist, FractalIterData *pIterData) const;
    virtual DBL CalcDirDeriv(const Vector4d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const;
    virtual void CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData,
                            FractalIterData *pPIterData) const;

    virtual void IterateCalc(Vector4d& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual void DerivCalc(Vector4d& rD, const Vector4d& v, int iter, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual bool DiscontinuityCheck(Vector4d& rD, DBL& rDist, const Vector4d& t, const Vector4d& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;

protected:
    const Vector4d mDuplexJuliaParm;

};

}

#endif
