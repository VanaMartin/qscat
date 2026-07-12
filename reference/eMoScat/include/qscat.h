/*!
    Master Header fo the QSCAT problem:

    This is just a helping header to provide easy to use
    inclusion on the QSCAT project classes.
*/

#ifndef INCLUDE_QSCAT_H_
#define INCLUDE_QSCAT_H_

#include <algorithm>
#include <string>
#include <cassert>
#include <stdio.h>
#include <complex>				// Complex algebra
#include <cmath>
#include <fstream>
#include <omp.h>				// Parallelization library

#include "blas.h"

#ifndef def_float
    typedef double                      def_float;
#endif

#ifndef dfloat       // TODO: cleanup def float and dfloat
    typedef double                      dfloat;
    typedef std::complex<double>        dcomp;
#endif

#ifndef qfloat       // TODO: cleanup def float and dfloat
    typedef double                      qfloat;
    typedef std::complex<double>        qcomp;
#endif
///// ARRAYS /////

namespace QSCAT
{ // Some Forward definitions for simplicity
template<typename T>
class Buffer;
template<typename T>
class Vector;
template<typename T>
class VectorScalarMultiple;
template<typename T>
class ShallowVector;
template<typename T>
class Matrix;
template<typename T>
class EigenSystem;

} // namespace QSCAT

#include "ScalarMultiple.h"
#include "Arrays/Buffer.h"
#include "Arrays/Vector.h"
#include "Arrays/EigenSystem.h"
#include "Arrays/Matrix.h"
#include "Arrays/RowCompressedMatrix.h"

namespace QSCAT
{ // Most common types

/// \brief floating point buffer
typedef Buffer<dfloat>      dBuffer;

/// \brief complex floating point buffer
typedef Buffer<dcomp>       zBuffer;

/// \brief integral type vector
/// \details some operations cannot be performed effectively, done via for loops
typedef Vector<blas_int>    iVector;

/// \brief floating point vector
typedef Vector<dfloat>      dVector;

/// \brief complex floating point vector
typedef Vector<dcomp>       zVector;

/// \brief integral type matrix
typedef Matrix<blas_int>    iMatrix;

/// \brief floating point matrix
typedef Matrix<dfloat>      dMatrix;

/// \brief complex floating point matrix
typedef Matrix<dcomp>       zMatrix;

/// \brief floating point eigen system
typedef EigenSystem<dfloat>             dEigenSystem;

/// \brief complex floating point eigensystem
typedef EigenSystem<dcomp>              zEigenSystem;

/// \brief floating point matrix in row compressed format (csr)
typedef RowCompressedMatrix<dfloat>     dRCMatrix;

/// \brief complex floating point matrix in row compressed format (csr)
typedef RowCompressedMatrix<dcomp>      zRCMatrix;
} // namspace QSCAT

#include "bessel.h"
#include "common.h"
#include "input.h"

///// FEM_DVR_ECS /////

#include "FemDvrEcs/FemDvrFunctions.h"
#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"
#include "FemDvrEcs/Projector.h"
#include "FemDvrEcs/OperatorDiagonal.h"
#include "FemDvrEcs/OperatorFull.h"
#include "FemDvrEcs/OperatorRowCompressed.h"
#include "FemDvrEcs/KineticEnergy.h"
#include "FemDvrEcs/DiscreteStates.h"
#include "FemDvrEcs.h"

namespace QSCAT
{
    OperatorFull buildFullHamiltonian(FemDvrEcsGrid& g, GridVector& v, const dfloat& mu);
}   // namspace QSCAT


#include "FemDvrEcs2d/FemDvrEcsGrid2d.h"
#include "FemDvrEcs2d/GridVector2d.h"
#include "FemDvrEcs2d/ShallowGridVector2d.h"
//#include "FemDvrEcs2d/DoubleGridVector2D.h"
#include "FemDvrEcs2d/OperatorRowCompressed2d.h"
#include "FemDvrEcs2d/OperatorFull2d.h"
#include "FemDvrEcs2d/ZoomFilter.h"
#include "FemDvrEcs2d/EquidistantProjector2d.h"
//#include "FemDvrEcs2d/DoubleOperator2DRC.h"
#include "FemDvrEcs2d.h"

#include "picojson/pjson.h"
#include "pjinput.h"
#include "interface.h"

#include "Model2d/TestFunctionInterface2d.h"

#include "Model2d/TestFunction2d.h"
#include "Model2d/DiracTestFunction2d.h"
#include "Model2d/FluxTestFunction2d.h"
#include "Model2d/MultiTestFunction2d.h"
#include "Model2d/TimeDependentModel2d.h"
#include "Model2d/TimeIndependentModel2d.h"
#include "Model2d/CoupledModel2d.h"

#include "ModelLCP/SMatrix.h"
#include "ModelLCP/ModelLCP.h"

#endif  // INCLUDE_QSCAT_H_
