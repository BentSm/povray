//******************************************************************************
///
/// @file core/shape/fractal/magic.h
///
/// This module contains prototypes for the generic implementation of
/// FractalRules subclasses in magicimpl.h.
///
/// (This was originally named for the templating 'magic' that it uses, but
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
#include "core/shape/fractal/cxx03compat.h"
#include "core/shape/fractal/types.h"

namespace pov
{

template <class BaseRules = FractalRules>
class InfoRulesBase : virtual public BaseRules
{
public:
    InfoRulesBase(const FractalFuncType& funcType, EstimatorType estimatorType,
                  DiscontinuitySupportLevel discontinuitySupport, int iterationDataSize) :
        mInfo(CreateRulesInfo(funcType, estimatorType, discontinuitySupport, iterationDataSize)) {}
    virtual const FractalRulesInfo& Info() const { return mInfo; }

protected:
    const FractalRulesInfo mInfo;
};

template <template <class> class RulesClass, class Estimator, class BaseRules = FractalRules>
class MagicFractalRulesBase : public InfoRulesBase<BaseRules>
{
public:
    MagicFractalRulesBase(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport) :
        InfoRulesBase<BaseRules>(data.funcType, Estimator::eType, discontinuitySupport,
                                 sizeof(typename RulesClass<Estimator>::IterationData)) {}
    virtual void CalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pTIterData,
                            void *pPIterData) const;
    virtual bool Bound(const BasicRay& ray, const Fractal *pFractal, DBL *pDepthMin, DBL *pDepthMax) const;
};

template <template <class> class RulesClass, class Estimator, class BaseRules = FractalRules>
class MagicQuaternionFractalRules : public MagicFractalRulesBase<RulesClass, Estimator, BaseRules> {
public:
    MagicQuaternionFractalRules(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport) :
        INIT_QUATERNION(mJuliaParm, data.juliaParm[X], data.juliaParm[Y], data.juliaParm[Z], data.juliaParm[W]),
        MagicFractalRulesBase<RulesClass, Estimator, BaseRules>(data, discontinuitySupport) {}
    virtual int Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction,
                        DBL *pDist, void *pIterData) const;
    DBL CalcDirDeriv(const Vector3d& dir, int nMax, const Fractal *pFractal, void *pIterData) const;

    inline void NonVirtualCalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pIterData) const
    { TemplatedCalcNormal<false>(rResult, nMax, pFractal, pIterData, NULL); }
    inline void NonVirtualCalcNormalWDisc(Vector3d& rResult, int nMax, const Fractal *pFractal,
                                          void *pTIterData, void *pPIterData) const
    { TemplatedCalcNormal<true>(rResult, nMax, pFractal, pTIterData, pPIterData); }

    bool DiscontinuityCheck(DBL& rDx, DBL& rDy, DBL& rDz, DBL& rDw, DBL& rDist,
                            DBL tx, DBL ty, DBL tz, DBL tw, DBL px, DBL py, DBL pz, DBL pw,
                            int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const;

protected:
    const Quaternion mJuliaParm;

    template <bool disc>
    void TemplatedCalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pTIterData,
                             void *pPIterData = NULL) const;

};

template <template <class> class RulesClass, class Estimator, class BaseRules = FractalRules>
class MagicHypercomplexFractalRules : public MagicFractalRulesBase<RulesClass, Estimator, BaseRules> {
public:
    MagicHypercomplexFractalRules(const FractalConstructorData& data, DiscontinuitySupportLevel discontinuitySupport) :
        INIT_DUPLEX(mDuplexJuliaParm, CreateComplex(data.juliaParm[X] - data.juliaParm[W], data.juliaParm[Y] + data.juliaParm[Z]),
                    CreateComplex(data.juliaParm[X] + data.juliaParm[W], data.juliaParm[Y] - data.juliaParm[Z])),
        MagicFractalRulesBase<RulesClass, Estimator, BaseRules>(data, discontinuitySupport) {}
    virtual int Iterate(const Vector3d& iPoint, const Fractal *pFractal, const Vector3d& direction,
                        DBL *pDist, void *pIterData) const;
    inline DBL CalcDirDeriv(const Vector3d& dir, int nMax, const Fractal *pFractal, void *pIterData) const
    {
        Vector3d normal;
        NonVirtualCalcNormal(normal, nMax, pFractal, pIterData);
        return dot(dir, normal);
    }

    inline void NonVirtualCalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pIterData) const
    { TemplatedCalcNormal<false>(rResult, nMax, pFractal, pIterData, NULL); }
    inline void NonVirtualCalcNormalWDisc(Vector3d& rResult, int nMax, const Fractal *pFractal,
                                          void *pTIterData, void *pPIterData) const
    { TemplatedCalcNormal<true>(rResult, nMax, pFractal, pTIterData, pPIterData); }

    bool DiscontinuityCheck(Complex& rD0, Complex& rD1, DBL& rDist,
                            const Complex& t0, const Complex& t1, const Complex& p0, const Complex& p1,
                            int iter, const Fractal *pFractal, void *pTIterData, void *pPIterData) const;

protected:
    const Duplex mDuplexJuliaParm;

    template <bool disc>
    void TemplatedCalcNormal(Vector3d& rResult, int nMax, const Fractal *pFractal, void *pTIterData,
                             void *pPIterData = NULL) const;
};

}

#endif
