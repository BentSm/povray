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

#ifndef POVRAY_CORE_FRACTAL_MAGIC_H
#define POVRAY_CORE_FRACTAL_MAGIC_H

#include "core/coretypes.h"
#include "core/math/vector.h"
#include "core/shape/fractal/types.h"
#include "core/shape/fractal/util.h"

namespace pov
{

class MagicRulesBase : public FractalRules
{
public:
    typedef NilData FixedData;
    typedef NilData MainIterData;
    typedef NilData AuxIterData;

    MagicRulesBase(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                   const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        mInfo(CreateRulesInfo(data.funcType, discontinuitySupport, sizes, estimator.eType)),
        mEstimator(estimator) { }

    virtual bool Bound(const BasicRay& ray, const Fractal *pFractal, DBL *pDepthMin, DBL *pDepthMax) const;
    virtual const FractalRulesInfo& Info() const { return mInfo; }

protected:
    const FractalRulesInfo mInfo;
    const DistanceEstimator& mEstimator;

    static const DistanceEstimator& GetEstimatorFromType(EstimatorType estimatorType, EstimatorType defaultEstimator = kNewtonEstimator,
                                                         EstimatorType legacyEstimator = kOrigNewtonEstimator,
                                                         const DistanceEstimator& (*ExtraEstimators)(EstimatorType eType) = NULL);
};

class MagicQuaternionFractalRules : public MagicRulesBase
{
public:
    typedef VECTOR_4D MainIterData;

    MagicQuaternionFractalRules(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                                const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        INIT_VECTOR_4D(mJuliaParm, data.juliaParm[X], data.juliaParm[Y], data.juliaParm[Z], data.juliaParm[W]),
        MagicRulesBase(data, discontinuitySupport, sizes, estimator) {}

    virtual int Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction,
                        DBL *pDist, FractalIterData *pIterData) const;
    virtual void CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData,
                            FractalIterData *pPIterData) const;

    virtual DBL CalcDirDeriv(const Vector3d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const;

    virtual void IterateCalc(VECTOR_4D& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual void DirDerivCalc(VECTOR_4D& rD, const VECTOR_4D& v, int iter, bool samePoint, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual bool DiscontinuityCheck(VECTOR_4D& rD, DBL& rDist, const VECTOR_4D& t, const VECTOR_4D& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;

protected:
    const IVECTOR_4D mJuliaParm;

};

class MagicHypercomplexFractalRules : public MagicRulesBase
{
public:
    typedef Duplex MainIterData;

    MagicHypercomplexFractalRules(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport,
                                  const FractalDataSizes& sizes, const DistanceEstimator& estimator) :
        INIT_DUPLEX(mDuplexJuliaParm, CreateComplex(data.juliaParm[X] - data.juliaParm[W], data.juliaParm[Y] + data.juliaParm[Z]),
                    CreateComplex(data.juliaParm[X] + data.juliaParm[W], data.juliaParm[Y] - data.juliaParm[Z])),
        MagicRulesBase(data, discontinuitySupport, sizes, estimator) {}

    virtual int Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction,
                        DBL *pDist, FractalIterData *pIterData) const;
    virtual void CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, FractalIterData *pTIterData,
                            FractalIterData *pPIterData) const;

    virtual DBL CalcDirDeriv(const Vector3d& dir, int nMax, const Fractal *pFractal, FractalIterData *pIterData) const
    {
        Vector3d normal;
        CalcNormal(normal, nMax, pFractal, pIterData, NULL);
        return dot(dir, normal);
    }

    virtual void IterateCalc(Duplex& rV, DBL norm, int iter, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual void DerivCalc(Duplex& rD, const Duplex& v, int iter, const Fractal *pFractal, FractalIterData *pIterData) const = 0;
    virtual bool DiscontinuityCheck(Duplex& rD, DBL& rDist, const Duplex& t, const Duplex& p,
                                    int iter, const Fractal *pFractal, FractalIterData *pTIterData, FractalIterData *pPIterData) const;

protected:
    const IDuplex mDuplexJuliaParm;

};

}

#endif
