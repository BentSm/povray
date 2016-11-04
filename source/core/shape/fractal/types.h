//******************************************************************************
///
/// @file core/shape/fractal/types.h
///
/// This module contains types for fractals and other general miscellany.
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

#ifndef POVRAY_CORE_FRACTAL_TYPES_H
#define POVRAY_CORE_FRACTAL_TYPES_H

#include "core/coretypes.h"

#ifndef FRACTAL_USE_CXX11
#include <cstdarg>
#endif

namespace pov
{

enum FractalAlgebra
{
    kComplex = -2,
    kInvalidAlgebra = -1,
    kQuaternion,
    kHypercomplex
};

enum FractalFunc_FuncType
{
    kFunc_Invalid = -1,
    kFunc_Sqr,
    kFunc_Cube,
    kFunc_Reciprocal,
    kFunc_Exp,
    kFunc_Ln,
    kFunc_Sin,
    kFunc_ASin,
    kFunc_Cos,
    kFunc_ACos,
    kFunc_Tan,
    kFunc_ATan,
    kFunc_Sinh,
    kFunc_ASinh,
    kFunc_Cosh,
    kFunc_ACosh,
    kFunc_Tanh,
    kFunc_ATanh,
    kFunc_Pwr
};

enum FractalFunc_VariantType
{
    kVar_Invalid = -1,
    kVar_Normal = 0,
    kVar_Left = 0,
    kVar_Alternate = 1,
    kVar_Right = 1
};

enum EstimatorType
{
    kInvalidEstimator = -1,
    kNoEstimator = 0,
    kDefaultEstimator = 1,
    kNewtonEstimator = 2,
    kOrigNewtonEstimator = 3,
    kOrigSpecialEstimators = 7,
    kLegacyEstimator = 15
};

enum DiscontinuitySupportLevel
{
    kDiscontinuityUnneeded = -1,
    kDiscontinuityNotImplemented = 0,
    kDiscontinuitySupported = 1
};

struct FractalFuncType
{
    FractalAlgebra algebra;
    FractalFunc_FuncType type;
    FractalFunc_VariantType variant;
};

inline bool operator<(const FractalFuncType& a, const FractalFuncType& b)
{
    return (a.algebra < b.algebra || (a.algebra == b.algebra && (a.type < b.type || (a.type == b.type && a.variant < b.variant))));
}

class FractalDataStore
{
public:
    FractalDataStore() : mData(NULL), mSize(-1) { }

    FractalDataStore(int dataSize) : mData(NULL) { AllocateNewData(dataSize); }

    ~FractalDataStore() { DeallocateData(); }

    void AllocateNewData(int dataSize)
    {
        DeallocateData();
        if (dataSize >= 0)
            mData = operator new(dataSize);
        mSize = dataSize;
    }

    void ReserveData(int dataSize)
    {
        if (mSize < dataSize)
            AllocateNewData(dataSize);
    }

    void DeallocateData()
    {
        operator delete(mData);
        mData = NULL;
        mSize = -1;
    }

    void *data() const { return mData; }

    int size() const { return mSize; }

protected:
    void *mData;
    int mSize;
};

struct FractalIterData
{
    FractalDataStore fixed, mainIter, auxIter;

    void ReserveStores(const FractalDataSizes& tSizes)
    {
        fixed.ReserveData(tSizes.fixedSize);
        mainIter.ReserveData(tSizes.mainIterSize);
        auxIter.ReserveData(tSizes.auxIterSize);
    }

    void ReserveStoresForIters(const FractalDataSizes& indSizes, int maxIters)
    {
        FractalDataSizes tSizes = { indSizes.fixedSize, indSizes.mainIterSize * (maxIters + 1),
                                    indSizes.auxIterSize * (maxIters + 1) };
        ReserveStores(tSizes);
    }

    const FractalDataSizes sizes()
    {
        FractalDataSizes sizes = { fixed.size(), mainIter.size(), auxIter.size() };
        return sizes;
    }
};

struct FractalRulesInfo
{
    FractalFuncType funcType;
    DiscontinuitySupportLevel discontinuitySupport;
    FractalDataSizes sizes;
    EstimatorType estimatorType;
};

struct NilData { };

typedef Complex Duplex[2];

typedef DBL EstimatorFunc(const FractalRules *pRules, DBL norm, int iters, const Vector3d& direction,
                          const Fractal *pFractal, FractalIterData *pIterData);

struct DistanceEstimator
{
    EstimatorFunc *pEstim;
    EstimatorType eType;
};

struct FractalConstructorData
{
    FractalFuncType funcType;
    VECTOR_4D juliaParm;
    Complex exponent;
    EstimatorType estimatorType;
};

/* If we don't have C++11, we can't initialize const arrays... */
#ifndef FRACTAL_USE_CXX11

template <class T, unsigned N>
class IArray
{
public:
    IArray(const T t0, ...)
    {
        va_list vals;
        int i;
        va_start(vals, t0);
        value[0] = t0;
        for (i = 1; i < N; i++)
            value[i] = va_arg(vals, const T);
        va_end(vals);
    }
    T& operator[](unsigned k) { return value[k]; }
    const T& operator[](unsigned k) const { return value[k]; }

    typedef T Array_Type[N];

    operator Array_Type&() { return value; }
    operator const Array_Type&() const { return value; }

private:
    Array_Type value;
};

typedef IArray<Complex, 2> IDuplex;
typedef IArray<DBL, 4> IVECTOR_4D;

#else

typedef Duplex IDuplex;
typedef VECTOR_4D IVECTOR_4D;

#endif

}

#endif
