#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <complex>
#include <string>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "bessel.h"
#include "common.h"
#include "blas.h"
#include "Arrays.h"
#include "input.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d/FemDvrEcsGrid2d.h"
#include "FemDvrEcs2d/GridVector2d.h"
#include "FemDvrEcs2d/ShallowGridVector2d.h"
#include "FemDvrEcs2d/DoubleGridVector2d.h"

namespace QSCAT
{

  // constructors
    DoubleGridVector2d::DoubleGridVector2d()
    {
        init_ = false;
    }
    DoubleGridVector2d::DoubleGridVector2d(FemDvrEcsGrid2d& grid) :
        grid1_(grid),
        grid2_(grid),
        body_(2*grid.get_size())
    {
        assert(grid.init());
      //
        init_ = true;
    }
    DoubleGridVector2d::DoubleGridVector2d(FemDvrEcsGrid2d& grid1, FemDvrEcsGrid2d& grid2) :
        grid1_(grid1),
        grid2_(grid2),
        body_(grid1.get_size() + grid2.get_size())
    {
        assert(grid1.init());
        assert(grdi2.init());
      //
        init_ = true;
    }
    DoubleGridVector2d::DoubleGridVector2d(FemDvrEcsGrid2d& grid, GridVector2d& psi, blas_int i=0) :
        grid1_(grid),
        grid2_(grid),
        body_(2*grid.get_size())
    {
        assert(i<2);
        assert(grid.init());
        assert(psi.init());
        assert(psi.get_grid() == grid);
      //
        body_.partial_assign(grid.get_size(), (i==1)? grid.get_size() : 0, 1, psi.body(), 0, 1);
        init_ = true;
    }
    DoubleGridVector2d::DoubleGridVector2d(FemDvrEcsGrid2d& grid1, FemDvrEcsGrid2d& grid2, GridVector2d& psi1, GridVector2d& psi2) :
        grid1_(grid1),
        grid2_(grid2),
        body_(grid1.get_size()+grid2.get_size())
    {
        assert(grid1.init());
        assert(grid2.init());
        assert(psi1.init());
        assert(psi2.init());
        assert(grid1 == psi1.get_grid());
        assert(grid2 == psi2.get_grid());
      //
        body_.partial_assign(grid1_.get_size(), 0, 1, psi1.body(), 0, 1);
        body_.partial_assign(grid2_.get_size(), grid1_.get_size(), 1, psi2.body(), 0, 1);
        init_ = true;
    }
    DoubleGridVector2d::DoubleGridVector2d(const DoubleGridVector2d& old) :
        grid1_(old.grid1_),
        grid2_(old.grid2_),
        body_(old.body_),
        init_(old.init_)
    {}

  // accessors
    bool DoubleGridVector2d::init() const
    {
        return init_;
    }
    dcomp& DoubleGridVector2d::operator[] (blas_int i)
    {
        assert(init_);
        assert(i < get_size());
      //
        return body_[i];
    }
    const dcomp& DoubleGridVector2d::operator[] (blas_int i) const
    {
        assert(init_);
        assert(i<get_size());
      //
        return body_[i];
    }
    blas_int DoubleGridVector2d::get_size1() const
    {
        assert(init_);
        return grid1_.get_size();
    }
    blas_int DoubleGridVector2d::get_size2() const
    {
        assert(init_);
      //
        return grid2_.get_size();
    }
    blas_int DoubleGridVector2d::get_size() const
    {
        assert(init_);
      //
        return get_size1() + get_size2();
    }
    void DoubleGridVector2d::f(const dcomp& val, blas_int i)
    {
        assert(init_);
        assert(i < get_size());
      //
        if (i<get_size1()) {
            body_[i] = val * sqrt(grid1_.wz(i));
        } else {
            body_[i] = val * sqrt( grid2_.wz(i-get_size1()) );
        }
    }
    dcomp DoubleGridVector2d::f(blas_int i) const
    {
        assert(init_);
        assert(i < get_size());
      //
        if (i< get_size1()) {
            return body_[i] / sqrt(grid1_.wz(i));
        } else {
            return body_[i] / sqrt(grid2_.wz(i - get_size1()));
        }
        return dcomp(0); // Security only
    }
    const zVector& DoubleGridVector2d::body() const
    {
        assert(init_);
      //
        return body_;
    }
    zVector& DoubleGridVector2d::body()
    {
        assert(init_);
      //
        return body_;
    }
    const FemDvrEcsGrid2d& DoubleGridVector2d::get_grid1() const
    {
        assert(init_);
      //
        return grid1_;
    }
    const FemDvrEcsGrid2d& DoubleGridVector2d::get_grid2() const
    {
        assert(init_);
      //
        return grid2_;
    }

  // modifiers
    DoubleGridVector2d& DoubleGridVector2d::swap(DoubleGridVector2d &rhs)
    {
        grid1_.swap(rhs.grid1_);
        grid2_.swap(rhs.grid2_);
        body_.swap(rhs.body_);
        std::swap(init_, rhs.init_);
        return *this;
    }

  // operators
    dcomp DoubleGridVector2d::operator* (const DoubleGridVector2d& rhs) const
    {
        assert(init_);
        assert(rhs.init_);
        assert(grid1_ == rhs.grid1_);
        assert(grid2_ == rhs.grid2_);
      //
        return body_*rhs.body_;
    }
    DoubleGridVector2d& DoubleGridVector2d::operator= (DoubleGridVector2d tmp)
    {
        return swap(tmp);
    }
    DoubleGridVector2d& DoubleGridVector2d::operator+= (const DoubleGridVector2d& rhs)
    {
        assert(init_);
        assert(rhs.init_);
        assert(grid1_ == rhs.grid1_);
        assert(grid2_ == rhs.grid2_);
      //
        body_ += rhs.body_;
        return *this;
    }
    DoubleGridVector2d& DoubleGridVector2d::operator-= (const DoubleGridVector2d& rhs)
    {
        assert(init_);
        assert(rhs.init_);
        assert(grid1_ == rhs.grid1_);
        assert(grid2_ == rhs.grid2_);
      //
        body_ -= rhs.body_;
        return *this;
    }
    DoubleGridVector2d& DoubleGridVector2d::operator*= (const dcomp& alpha)
    {
        assert(init_);
      //
        body_*=alpha;
        return *this;
    }

  // custom operators
    DoubleGridVector2d& DoubleGridVector2d::axpy(const dcomp& alpha, const DoubleGridVector2d& x)
    {
        assert(init_);
        assert(x.init_);
        assert(grid1_ == x.grid1_);
        assert(grid2_ == x.grid2_);
      //
        body_.axpy(alpha, x.body_);
        return *this;
    }
    DoubleGridVector2d& DoubleGridVector2d::ax(const dcomp& alpha, const DoubleGridVector2d& x)
    {
        assert(init_);
        assert(x.init_);
        assert(grid1_ == x.grid1_);
        assert(grid2_ == x.grid2_);
      //
        body_.ax(alpha, x.body());
        return *this;
    }

}   // namespace QSCAT
