#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <complex>

#include "common.h"
#include "Arrays.h"

#include "Arrays/Buffer.hpp"
#include "Arrays/Vector.hpp"
//#include "Arrays/VectorScalarMultiple.hpp"
#include "Arrays/EigenSystem.hpp"
#include "Arrays/Matrix.hpp"
#include "Arrays/RowCompressedMatrix.hpp"

namespace QSCAT
{
    // Specify templates for common classes
    template class Buffer<dfloat>;
    template class Buffer<dcomp>;
    template class Vector<blas_int>;
    template class Vector<dfloat>;
    template class Vector<dcomp>;
    template class EigenSystem<dfloat>;
    template class EigenSystem<dcomp>;
    template class Matrix<blas_int>;
    template class Matrix<dfloat>;
    template class Matrix<dcomp>;
    template class RowCompressedMatrix<dfloat>;
    template class RowCompressedMatrix<dcomp>;

    // Specify functions for common classes
    template ScalarMultiple<dfloat, Vector<dfloat> > operator* (const dfloat& scalar, Vector<dfloat>& object);
    template ScalarMultiple<dfloat, Vector<dfloat> > operator* (Vector<dfloat>& object, const dfloat& scalar);
    template ScalarMultiple<dcomp, Vector<dcomp> > operator* (const dcomp& scalar, Vector<dcomp>& object);
    template ScalarMultiple<dcomp, Vector<dcomp> > operator* (Vector<dcomp>& object, const dcomp& scalar);

    template RowCompressedMatrix<dcomp> TensorSum(const RowCompressedMatrix<dcomp>& A, const RowCompressedMatrix<dcomp>& B);

}   // namespace QSCAT
