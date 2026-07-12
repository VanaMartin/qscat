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

namespace QSCAT
{

 // internals

    void OperatorDiagonal::initialize(const FemDvrEcsGrid& grid, blas_int length, blas_int start)
    {
        assert(grid.init());
        assert(length > 0);
        assert(start >= 0);
        assert(length + start <= grid.nb());
    //
        start_ = start;
        size_ = length;
        body_ = zVector(size_);
        grid_ = grid;
        max_size_ = grid.nb();
        *init_ = true;
    }

 // constructors

    OperatorDiagonal::OperatorDiagonal() : Object()
    {
        grid_ = FemDvrEcsGrid();
        start_ = 0;
        size_ = 0;
        max_size_ = 0;
        *init_ = false;
        body_ = zVector();
    }
    OperatorDiagonal::OperatorDiagonal(const FemDvrEcsGrid& grid) : Object()
    {
        assert(grid.init());
      //
        initialize(grid, grid.nb(), 0);
    }
    OperatorDiagonal::OperatorDiagonal(const FemDvrEcsGrid& grid, blas_int length, blas_int start) : Object()
    {
        assert(grid.init());
        assert(length>0);
        assert(start>=0);
      //
        initialize(grid, length, start);
    }
    OperatorDiagonal::OperatorDiagonal(const OperatorDiagonal& old) :
        Object(old),
        max_size_(old.max_size_),
        start_(old.start_),
        size_(old.size_),
        body_(old.body_),
        grid_(old.grid_)
    {}
    OperatorDiagonal::~OperatorDiagonal()
    {
        decref();
    }

 // accessors
    const dcomp& OperatorDiagonal::operator[] (blas_int i) const
    {
        assert(init());
        assert(i - start_ >= 0);
        assert(i - start_ < size_);
      //
        return body_[i - start_];
    }
    dcomp& OperatorDiagonal::operator[] (blas_int i)
    {
        assert(init());
        assert(i-start_ >= 0);
        assert(i-start_ < size_);
      //
        return body_[i-start_];
    }
    blas_int OperatorDiagonal::get_size() const
    {
        assert(init());
      //
        return size_;
    }
    blas_int OperatorDiagonal::get_shift() const
    {
        assert(init());
      //
        return start_;
    }
    const FemDvrEcsGrid& OperatorDiagonal::get_grid() const
    {
        assert(init());
      //
        return grid_;
    }

 // modifiers

    OperatorDiagonal& OperatorDiagonal::identity()
    {
        assert(init());
      //
        body_ = zVector(max_size_);
        size_ = max_size_;
        start_ = 0;
        body_.fill(dcomp(1));
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::swap(OperatorDiagonal& rhs)
    {
        Object::swap(rhs);
        std::swap(size_, rhs.size_);
        std::swap(start_, rhs.start_);
        std::swap(max_size_, rhs.max_size_);
        grid_.swap(rhs.grid_);
        body_.swap(rhs.body_);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::inverse()
    {
        assert(init());
      //
        for (blas_int i=0; i<size_; ++i){
            body_[i] = dcomp(1)/body_[i];
        }
        return *this;
    }

 // operators

    OperatorDiagonal& OperatorDiagonal::operator=  (OperatorDiagonal tmp)
    {
        swap(tmp);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator=  (const GridVector& rhs)
    {
        assert(init());                              // Considering only initialized operators
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());               // Only if grids are equivalent
        assert(start_ + size_ <= rhs.get_size());   // Considering only Vectors covering all of the length
      //
        for (blas_int i=0; i<size_; ++i){
            body_[i] = rhs.f(start_ + i);
        }
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator=  (const zVector& rhs)
    {
        assert(init());                              // Considering only initialized operators
        assert(rhs.init());
        assert(size_ <= rhs.get_size());   // Considering only Vectors covering all of the length
      //
        blas::copy(size_, &rhs[0], 1, &body_[0], 1);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator=  (const dcomp& alpha)
    {
        assert(init());
      //
        body_.fill(alpha);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator+= (const OperatorDiagonal& rhs)
    {
        assert(init());
        assert(rhs.init());
      //
        blas_int length = min(size_,rhs.size_) - abs(start_ - rhs.start_);
        blas::axpy(length, dcomp(1), &(rhs.body_[max(0ll,rhs.start_-start_)]), 1, &body_[max(0ll, start_-rhs.start_)], 1);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator+= (const zVector& rhs)
    {
        assert(init());
        assert(rhs.init());
      //
        blas_int length = min(size_, rhs.get_size());
        blas::axpy(length, dcomp(1), &(rhs[0]), 1, &body_[0], 1);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator+= (const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
        assert(start_ + size_ <= rhs.get_size());
      //
        for (blas_int i=0; i<size_; ++i){
            body_[i] += rhs.f(i+start_);
        }
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator+= (const dcomp& alpha)
    {
        assert(init());
      //
        body_ += alpha;
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator-= (const OperatorDiagonal& rhs)
    {
        assert(init());
        assert(rhs.init());
      //
        blas_int length = min(size_,rhs.size_) - abs(start_ - rhs.start_);

        blas::axpy(length, dcomp(-1), &(rhs.body_[max(0ll,rhs.start_-start_)]), 1, &body_[max(0ll, start_-rhs.start_)], 1);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator-= (const zVector& rhs)
    {
        assert(init());
        assert(rhs.init());
      //
        blas_int length = min(size_, rhs.get_size());
        blas::axpy(length, dcomp(-1), &(rhs[0]), 1, &body_[0], 1);
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator-= (const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
        assert(start_ + size_ <= rhs.get_size());
      //
        for (blas_int i=0; i<size_; ++i){
            body_[i] -= rhs.f(i+start_);
        }
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator-= (const dcomp& alpha)
    {
        assert(init());

        body_ -= alpha;
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator*= (const OperatorDiagonal& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        blas_int length = min(size_, rhs.size_) - abs(start_ - rhs.start_);
        dcomp* aux = new dcomp[length];

        blas::ewxy(length, &body_[max(0ll, start_-rhs.start_)], &(rhs.body_[max(0ll, rhs.start_-start_)]), aux);
        blas::copy(length, aux, &body_[max(0ll, start_-rhs.start_)]);
        delete[] aux;
        return *this;
    }
    OperatorDiagonal& OperatorDiagonal::operator*= (const dcomp& alpha)
    {
        assert(init());
      //
        body_ *= alpha;
        return *this;
    }
    GridVector OperatorDiagonal::operator* (const GridVector& rhs) const
    {
        assert(init());
        assert(grid_ == rhs.get_grid());
        assert(start_+size_ <= rhs.get_size());
        // FIXME
        //for (blas_int i=0; i<size_; ++i){
        //    rhs[i+s] *= body_[i];
        //}
        GridVector out = rhs.copy();

        dcomp zero = dcomp(0);
        if (start_) {   // fill with zeros on uncovered space
            blas::copy(start_, &zero, 0, &out[0], 1);
        }
        dcomp *aux = new dcomp[size_];
        blas::ewxy(size_, &body_[0], &out[start_], aux);
        blas::copy(size_, aux, &out[start_]);
        delete[] aux;
        if (start_ + size_ < max_size_) {
            blas::copy(max_size_ - start_ - size_, &zero, 0, &out[start_+size_], 1);
        }
        return out;
    }
} // namespace QSCAT
