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
#include "FemDvrEcs/KineticEnergy.h"
#include "FemDvrEcs/OperatorFull.h"

namespace QSCAT
{

  // internals
    void OperatorFull::initialize(const FemDvrEcsGrid& grid, blas_int length, blas_int start)
    {
        assert(length > 0);
        assert(start >= 0);
        assert(length + start <= grid.nb());
      //
        grid_ = grid;
        size_ = length;
        start_ = start;
        max_size_ = grid.nb();
        body_ = zMatrix(size_, size_);
        body_.fill(dcomp(0));
        *init_ = true;
        transposed_ = 'N';
        factorized_ = false;
    }

    // constructors
    OperatorFull::OperatorFull() : Object()
    {
        size_ = 0;
        start_ = 0;
        max_size_ = 0;
        *init_ = false;
        transposed_ = 'N';
        factorized_ = false;
    }
    OperatorFull::OperatorFull(const FemDvrEcsGrid& grid, blas_int length, blas_int start) : Object()
    {
        initialize(grid, length, start);
    }
    OperatorFull::OperatorFull(const FemDvrEcsGrid& grid) : Object()
    {
        initialize(grid, grid.nb(), 0);
    }
    OperatorFull::OperatorFull(const OperatorFull& old):
        Object(old),
        size_(old.size_),
        start_(old.start_),
        max_size_(old.max_size_),
        grid_(old.grid_),
        body_(old.body_),
        transposed_(old.transposed_),
        factorized_(old.factorized_)
    {}
    OperatorFull::~OperatorFull()
    {
        decref();
    }

    // TODO check consistency of indices for shifted operators
    // accessors
    const dcomp& OperatorFull::operator[] (blas_int i) const
    {
        assert(init());
        assert(i >= 0);
        assert(i < size_ * size_);
      //
        return body_[i];
    }
    dcomp& OperatorFull::operator[] (blas_int i)
    {
        assert(init());
        assert(i >= 0);
        assert(i < size_*size_);
      //
        return body_[i];
    }
    const FemDvrEcsGrid& OperatorFull::grid() const
    {
        assert(init());
      //
        return grid_;
    }
    blas_int OperatorFull::get_size() const
    {
        assert(init());
      //
        return size_;
    }
    blas_int OperatorFull::get_shift() const
    {
        assert(init());
      //
        return start_;
    }

    // modifiers
    OperatorFull& OperatorFull::transpose()
    {
        assert(init());
      //
        if (transposed_=='N') {
            transposed_='T';
        } else {
            transposed_='N';
        }
        return *this;
    }
    OperatorFull& OperatorFull::set_identity()
    {
        assert(init());
      //
        body_.set_identity(size_);
        transposed_ = 'N';
        factorized_ = false;
        return *this;
    }
    OperatorFull& OperatorFull::swap(OperatorFull& rhs)
    {
        Object::swap(rhs);
        std::swap(size_, rhs.size_);
        std::swap(start_, rhs.start_);
        std::swap(max_size_, rhs.max_size_);
        grid_.swap(rhs.grid_);
        std::swap(transposed_, rhs.transposed_);
        std::swap(factorized_, rhs.factorized_);
        body_.swap(rhs.body_);
        return *this;
    }
    OperatorFull& OperatorFull::outer_product(const GridVector& rhs1, const GridVector& rhs2)
    {
        assert(init());
        assert(grid_ == rhs1.get_grid());
        assert(grid_ == rhs2.get_grid());
      //
        dcomp alpha = dcomp(1.0);
        dcomp beta = dcomp(0.0);

        blas::matrix_matrix('N', 'N', size_, size_, 1, alpha, beta, &rhs1[start_], &rhs2[start_], &body_[0]);
        transposed_ = 'N';
        factorized_ = false;
        return *this;
    }
    OperatorFull& OperatorFull::inverse()
    {
        assert(init());
      //
        zMatrix aux;
        aux.set_identity(size_);
        blas::lapack_solve(size_, size_, &body_[0], &aux[0]);
        body_.swap(aux);
        return *this;
    }
    OperatorFull& OperatorFull::LU_factorize()
    {
        assert(0);    // NIY
        return *this;
    }
    OperatorFull& OperatorFull::add_kinetic_term(const dfloat& mu)
    {
        assert(init());
      //
        body_ += generateKineticTerm( grid_, mu);
        return *this;
    }

    // operators

    OperatorFull OperatorFull::copy() const
    {
        OperatorFull out;
        *out.init_ = *init_;
        out.max_size_ = max_size_;
        out.start_ = start_;
        out.size_ = size_;
        out.transposed_ = transposed_;
        out.factorized_ = factorized_;
        out.grid_ = grid_.copy();
        out.body_ = body_.copy();
        return out;
    }
    OperatorFull& OperatorFull::operator=  (OperatorFull tmp)
    {
        this->swap(tmp);
        return *this;
    }
    OperatorFull& OperatorFull::operator=  (const OperatorDiagonal& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        blas_int length = min(size_, rhs.get_size()) - abs(start_ - rhs.get_shift());
        body_.fill(dcomp(0));
        blas::copy(length, &rhs[max(0ll, rhs.get_shift() - start_)], 1, &body_[max(0ll, start_ - rhs.get_shift())], size_ + 1);
        factorized_ = false;
        transposed_ = 'N';
        return *this;
    }
    OperatorFull& OperatorFull::operator=  (const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
     //
        body_.fill(dcomp(0));
        for (blas_int i=0; i<size_; ++i){
            body_[i*(size_+1)] = rhs.f(i + start_);
        }
        factorized_ = false;
        transposed_ = 'N';
        return *this;
    }
    OperatorFull& OperatorFull::operator=  (const zVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(size_ <= rhs.get_size());
      //
        body_.fill(dcomp(0));
        blas::copy(size_, &rhs[0], 1, &body_[0], size_+1);
        factorized_ = false;
        return *this;
    }
    OperatorFull& OperatorFull::operator=  (const dcomp& alpha)
    {
        assert(init());
      //
        //body_.fill(Z(0));
        //blas::copy(size_, &alpha, 0, &body_[0], size_+1);
        body_.fill(alpha);
        factorized_ = false;
        transposed_ = 'N';
        return *this;
    }
    OperatorFull& OperatorFull::operator+= (const OperatorFull& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        if (size_ == rhs.size_ && start_ == rhs.start_){
            blas::axpy(size_*size_, dcomp(1), &(rhs.body_[0]), 1, &body_[0], 1);
        } else {
            blas_int start = max(start_, rhs.start_);
            blas_int length = min(start_ + size_ - start, rhs.start_ + rhs.size_ - start);
            for (blas_int i=0; i<length; ++i){
                blas::axpy(length, dcomp(1), &(rhs.body_[max(start-rhs.get_shift(),0ll) + i*rhs.size_]), 1, &body_[max(start-start_,0ll) + i*size_], 1);
            }
        }
        return *this;
    }
    OperatorFull& OperatorFull::operator+= (const OperatorDiagonal& rhs)                // Adds the diagonal operator
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        blas_int start = max(start_, rhs.get_shift());
        blas_int length = min(start_ + size_ - start, rhs.get_shift() + rhs.get_size() - start);
        blas::axpy(length, dcomp(1), &rhs[max(start-rhs.get_shift(),0ll)], 1, &body_[max(start-start_,0ll)], size_+1);
        return *this;
    }
    OperatorFull& OperatorFull::operator+= (const zVector& rhs)                         // Adds the Vector to the diagonal
    {
        assert(init());
        assert(rhs.init());
        assert(size_ <= rhs.get_size());
      //
        blas::axpy(size_, dcomp(1), &rhs[0], 1, &body_[0], size_+1);
        return *this;
    }
    OperatorFull& OperatorFull::operator+= (const GridVector& rhs)                      // Adds the functional values to the diagonal
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        //#pragma omp parallel for
        for (blas_int i=0; i<size_; ++i){
            body_[i*(size_+1)] += rhs.f(start_+i);
        }
        return *this;
    }
    OperatorFull& OperatorFull::operator+= (const dcomp& alpha)                         // Adds the Identity Matrix multiplied by alpha
    {
        assert(init());
      //
        blas::axpy(size_, dcomp(1), &alpha, 0, &body_[0], size_+1);
        return *this;
    }
    OperatorFull& OperatorFull::operator-= (const OperatorFull& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        if (size_ == rhs.size_ && start_ == rhs.start_){
            blas::axpy(size_*size_, dcomp(-1), &(rhs.body_[0]), 1, &body_[0], 1);
        } else {
            blas_int start = max(start_, rhs.start_);
            blas_int length = min(start_ + size_ - start, rhs.start_ + rhs.size_ - start);
            for (blas_int i=0; i<length; ++i){
                blas::axpy(length, dcomp(-1), &(rhs.body_[max(start-rhs.get_shift(),0ll) + i*rhs.size_]), 1, &body_[max(start-start_,0ll) + i*size_], 1);
            }
        }
        return *this;
    }
    OperatorFull& OperatorFull::operator-= (const OperatorDiagonal& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        blas_int start = max(start_, rhs.get_shift());
        blas_int length = min(start_ + size_ - start, rhs.get_shift() + rhs.get_size() - start);
        blas::axpy(length, dcomp(-1), &rhs[max(start-rhs.get_shift(),0ll)], 1, &body_[max(start-start_,0ll)], size_+1);
        return *this;
    }
    OperatorFull& OperatorFull::operator-= (const zVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(size_ <= rhs.get_size());
      //
        blas::axpy(size_, dcomp(-1), &rhs[0], 1, &body_[0], size_+1);
        return *this;
    }
    OperatorFull& OperatorFull::operator-= (const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        //#pragma omp parallel for
        for (blas_int i=0; i<size_; ++i){
            body_[i*(size_+1)] -= rhs.f(start_+i);
        }
        return *this;
    }
    OperatorFull& OperatorFull::operator-= (const dcomp& alpha)
    {
        assert(init());
      //
        blas::axpy(size_, dcomp(-1), &alpha, 0, &body_[0], size_+1);
        return *this;
    }
    OperatorFull& OperatorFull::operator*= (const OperatorFull& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
        assert(size_==rhs.size_);
      //
        dcomp alpha = dcomp(1);
        dcomp beta  = dcomp(0);
        zMatrix out(size_, size_);

        blas::matrix_matrix(rhs.transposed_, transposed_, size_, size_, size_, alpha, beta, &(rhs.body_[0]), &body_[0], &out[0]);
        body_.swap(out);
        return *this;
    }
    OperatorFull& OperatorFull::operator*= (const OperatorDiagonal& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid() );
      //
        blas_int start = max(start_,rhs.get_shift());
        blas_int length = min(start_ + size_ - start, rhs.get_shift() + rhs.get_size() - start);

        if (transposed_=='N'){
            if (start >= 0 && length > 0){
                for (blas_int i=0; i<length; ++i){
                    blas::scale(length, &body_[start - start_ + i*size_], rhs[start - rhs.get_shift() + i]);
                }
            }
        } else {
            if (start >= 0 && length > 0){
                for (blas_int i=0; i<length; ++i){
                    blas::scale(length, &body_[start - start_ + i], size_, rhs[start - rhs.get_shift() + i]);
                }
            }
        }
        return *this;
    }
    OperatorFull& OperatorFull::operator*= (const dcomp& alpha)
    {
        assert(init());
      //
        body_ *= alpha;
        return *this;
    }
    GridVector OperatorFull::operator* (const GridVector& rhs) const
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        GridVector out = rhs.copy();
        out.fill(0);
        blas::matrix_vector(transposed_, size_, size_, dcomp(1), dcomp(0), &body_[0], &rhs[start_], &out[start_]);
        return out;
    }

    // custom operations
    EigenSystem<dcomp> OperatorFull::eigen_system()
    {
        assert(init());
      //
        return body_.get_eigen_system();
    }
    GridVector& OperatorFull::back_substitution(GridVector& rhs)
    {
        assert(init());
        assert(grid_ == rhs.get_grid());
      //
        blas::lapack_solve(size_, 1, &body_[0], &rhs[start_]);
        return rhs;
    }
    GridVector& OperatorFull::LU_back_substitution(GridVector& rhs)
    {
        assert(0);   // NIY
        return rhs;
    }
    GridVector& OperatorFull::smart_back_substitution(GridVector& rhs)
    {
        if (factorized_) {
            LU_back_substitution(rhs);
        } else {
            back_substitution(rhs);
        }
        return rhs;
    }
    GridVector& OperatorFull::gemv(const dcomp alpha, const GridVector &x, const dcomp beta, GridVector& y)
    {
        assert(init());
        assert(x.init());
        assert(y.init());
        assert(grid_ == x.get_grid());
        assert(grid_ == y.get_grid());
        assert(size_ == max_size_); // gemv implemented only for full operators (TODO)

        body_.gemv(alpha, x.body(), beta, y.body());
        return y;
    }

    // TODO
    // Operator binary oprators
    OperatorFull operator* (OperatorFull lhs, OperatorDiagonal& rhs)
    {
        return lhs*=rhs;
    }
    OperatorFull operator* (OperatorDiagonal& lhs, OperatorFull rhs)
    {
        return ((rhs.transpose()) *= lhs).transpose();
    }
    OperatorFull operator* (OperatorFull lhs, OperatorFull& rhs)
    {
        lhs*=rhs;
        return lhs;
    }
    GridVector operator* (GridVector v, OperatorFull& O)
    {
        O.transpose();
        v = O*v;
        O.transpose();
        return v;
    }
} // namespace QSCAT
