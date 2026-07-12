#ifndef INCLUDE_QSCAT_ARRAYS_H_
#define INCLUDE_QSCAT_ARRAYS_H_

#include <algorithm>
#include <string>"
#include <cassert>
#include "stdio.h"
#include "blas.h"

#ifndef def_float       // TODO: cleanup def float and dfloat
    typedef double                      def_float;
#endif
#ifndef dfloat       // TODO: cleanup def float and dfloat
    typedef double                      dfloat;
    typedef std::complex<double>        dcomp;
#endif

/// Simple mathematical array templates
/*!
    The classes are closely connected to fast algebra via BLAS and
    solver methods via LAPACK. All necessary mathematical operations
    should be placed as class methods in here, so the scientific code
    can be clear of purely mathematical (and implementaitonal) problems.
*/
namespace QSCAT
{ // Forward definitions for simplicity
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

/** \addtogroup Interface
* @{ */

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

/** @} */

} // namspace QSCAT

#endif  // INCLUDE_QSCAT_ARRAYS_H_
