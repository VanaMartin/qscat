#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <complex>
#include <cassert>
#include <math.h>
#include <stdlib.h>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "common.h"
#include "Arrays.h"
#include "input.h"
#include "bessel.h"
#include "coulomb.h"

#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"
#include "FemDvrEcs/OperatorDiagonal.h"
#include "FemDvrEcs/OperatorRowCompressed.h"
#include "FemDvrEcs/KineticEnergy.h"

namespace QSCAT
{
 // internals

    void OperatorRowCompressed::initialize(const FemDvrEcsGrid& grid, blas_int num_nonzeros, blas_int length, blas_int shift)    // Shallow initialization : body is set to zero
    {
        rows_ = grid.nb();
        columns_ = rows_;
        if (num_nonzeros) {
            body_ = RowCompressedMatrix<dcomp>(rows_, columns_, num_nonzeros);
        }
        start_ = shift;
        size_ = length;
        grid_ = grid;
        *init_ = true;
    }

 // constructors

    OperatorRowCompressed::OperatorRowCompressed() : Object(), body_(), grid_()
    {
        *init_ = false;
        rows_ = 0;
        columns_ = 0;
        size_ = 0;
        start_ = 0;
        max_size_ = 0;
        body_ = RowCompressedMatrix<dcomp>();
    }
    OperatorRowCompressed::OperatorRowCompressed(const FemDvrEcsGrid& grid) : Object()
    {
        assert(grid.init());
      //
        initialize(grid, 0, grid.nb(), 0);
    }
    OperatorRowCompressed::OperatorRowCompressed(const FemDvrEcsGrid& grid, blas_int length, blas_int shift) : Object()
    {
        assert(grid.init());
        assert(length > 0);
        assert(shift >= 0);
        assert(length + shift <= grid.nb());
      //
        initialize(grid, 0, length, shift);
    }
    OperatorRowCompressed::OperatorRowCompressed(const FemDvrEcsGrid& grid, blas_int num_nonzeros, blas_int length, blas_int shift) : Object()
    {
        assert(grid.init());
        assert(length > 0);
        assert(shift >= 0);
        assert(length + shift <= grid.nb());
        assert(num_nonzeros >= 0);
      //
        initialize(grid, num_nonzeros, length, shift);
    }
    OperatorRowCompressed::OperatorRowCompressed(const OperatorRowCompressed& old) :
        Object(old),
        rows_(old.rows_),
        columns_(old.columns_),
        max_size_(old.max_size_),
        size_(old.size_),
        start_(old.start_),
        grid_(old.grid_),
        body_(old.body_)
    {}
    OperatorRowCompressed::~OperatorRowCompressed()
    {
        decref();
    }

 // accessors

    blas_int OperatorRowCompressed::get_size() const
    {
        assert(init());
      //
        return size_;
    }
    blas_int OperatorRowCompressed::get_shift() const
    {
        assert(init());
      //
        return start_;
    }
    const FemDvrEcsGrid& OperatorRowCompressed::get_grid() const
    {
        assert(init());
      //
        return grid_;
    }
    dcomp OperatorRowCompressed::operator() (int i, int j) const
    {
        assert(init());
        assert(i < grid_.get_size());
        assert(j < grid_.get_size());
      //
        return body_.get_element(i, j);
    }

 // modifiers

    OperatorRowCompressed& OperatorRowCompressed::swap(OperatorRowCompressed& rhs)
    {
        Object::swap(rhs);
        body_.swap(rhs.body_);
        std::swap(max_size_, rhs.max_size_);
        std::swap(rows_, rhs.rows_);
        std::swap(columns_, rhs.columns_);
        std::swap(size_, rhs.size_);
        std::swap(start_, rhs.start_);
        grid_.swap(rhs.grid_);
        return *this;
    }
    OperatorRowCompressed& OperatorRowCompressed::set_kinetic_term(dfloat mu)
    {
        assert(init());
      //
        body_ = generateKineticTermRCM(grid_, mu);
        start_ = 0;
        size_ = max_size_;
        return *this;
    }
    OperatorRowCompressed& OperatorRowCompressed::complex_conjugate()
    {
        assert(init());
      //
        body_.complex_conjugate();
        return *this;
    }
    OperatorRowCompressed& OperatorRowCompressed::conjugate()
    {
        assert(init());
      //
        body_.conjugate();
        return *this;
    }

 // operators

    OperatorRowCompressed OperatorRowCompressed::copy() const
    {
        OperatorRowCompressed out;
        *out.init_ = *init_;
        out.max_size_ = max_size_;
        out.rows_ = rows_;
        out.columns_ = columns_;
        out.start_ = start_;
        out.size_ = size_;
        out.grid_ = grid_.copy();
        out.body_ = body_.copy();
        return out;
    }

    OperatorRowCompressed& OperatorRowCompressed::operator= (OperatorRowCompressed rhs)
    {
        this->swap(rhs);
        return *this;
    }
    OperatorRowCompressed& OperatorRowCompressed::operator+= (const GridVector& rhs)
    {
        assert(init());
        assert(grid_ == rhs.get_grid());
      //
        body_.add_vector_to_diagonal(rhs.function_values());
        return *this;
    }
    OperatorRowCompressed& OperatorRowCompressed::operator+= (const dcomp& alpha)
    {
        assert(init());
      //
        body_.add_to_diagonal(alpha);
        return *this;
    }
    OperatorRowCompressed& OperatorRowCompressed::operator*= (const dcomp& alpha)
    {
        assert(init());
      //
        body_ *= alpha;
        return *this;
    }

 // custom operations

    void OperatorRowCompressed::gemv(const dcomp alpha, const GridVector& x, const dcomp beta, GridVector& y) const
    {
        assert(init());
        assert(x.init());
        assert(y.init());
        assert(grid_ == x.get_grid());
        assert(grid_ == y.get_grid());
      //
        body_.gemv(alpha, x.body(), beta, y.body());
    }
    void OperatorRowCompressed::LU_factorize()
    {
        assert(init());
        assert(rows_ == columns_);
      //
        body_.LU_factorize();
    }
    void OperatorRowCompressed::LU_back_substitution(GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        body_.LU_back_substitution(rhs.body());
    }

    const RowCompressedMatrix<dcomp>& OperatorRowCompressed::body() const
    {
        return body_;
    }

} // namespace QSCAT
