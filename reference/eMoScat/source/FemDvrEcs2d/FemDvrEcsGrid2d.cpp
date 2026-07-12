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

namespace QSCAT
{
  //Two dimensional grid type. The ordering of the grid is given by
  //repeating the X coordinate (as a column of a matrix) for each point of the Y coordinate. Stored
  //in a matrix, the elements a(i,j) = a_x(i)*a_y(j)

  // internals
    void FemDvrEcsGrid2d::initialize(const FemDvrEcsGrid& gx, const FemDvrEcsGrid& gy)
    {
        assert(gx.init());
        assert(gy.init());
      //
        xgrid_ = gx;
        ygrid_ = gy;
        xbasis_size_ = gx.nb();
        ybasis_size_ = gy.nb();
        basis_size_ = xbasis_size_ * ybasis_size_;
        *init_ = true;
    }
    void FemDvrEcsGrid2d::clean()
    {}
    bool FemDvrEcsGrid2d::save_bin_body(std::ofstream& file) const
    {
        assert(init());
      //
        bool stat = file.is_open();
        if (stat) {
            file.write((char*) &xbasis_size_, sizeof(blas_int));
            file.write((char*) &ybasis_size_, sizeof(blas_int));
            file.write((char*) &basis_size_, sizeof(blas_int));
            stat = xgrid_.save_binary(file);
            if(stat) stat = ygrid_.save_binary(file);
        }
        return stat;
    }
    bool FemDvrEcsGrid2d::read_bin_body(std::ifstream & file)
    {
        bool stat = file.is_open();
        if (stat) {
            file.read((char*) &xbasis_size_, sizeof(blas_int));
            file.read((char*) &ybasis_size_, sizeof(blas_int));
            file.read((char*) &basis_size_,  sizeof(blas_int));
            stat = xgrid_.read_binary(file);
            if(stat) stat = ygrid_.read_binary(file);
            *init_ = true;
        }
        return stat;
    }

  // constructors
    FemDvrEcsGrid2d::FemDvrEcsGrid2d() : Object()
    {
        xbasis_size_ = 0;
        ybasis_size_ = 0;
        basis_size_ = 0;
        *init_ = false;
    }
    FemDvrEcsGrid2d::FemDvrEcsGrid2d(const FemDvrEcsGrid& gx, const FemDvrEcsGrid& gy) : Object()
    {
        assert(gx.init());
        assert(gy.init());
      //
        initialize(gx,gy);
    }
    FemDvrEcsGrid2d::FemDvrEcsGrid2d(const FemDvrEcsGrid2d& old) :
        Object(old),
        xgrid_(old.xgrid_),
        ygrid_(old.ygrid_),
        xbasis_size_(old.xbasis_size_),
        ybasis_size_(old.ybasis_size_),
        basis_size_(old.basis_size_)
    {}
    FemDvrEcsGrid2d::~FemDvrEcsGrid2d()
    {
        decref();
        clean();
    }

    blas_int FemDvrEcsGrid2d::get_xsize() const
    {
        assert(init());
      //
        return xbasis_size_;
    }
    blas_int FemDvrEcsGrid2d::get_ysize() const
    {
        assert(init());
      //
        return ybasis_size_;
    }
    blas_int FemDvrEcsGrid2d::get_size() const
    {
        assert(init());
      //
        return basis_size_;
    }
    const dfloat& FemDvrEcsGrid2d::xr(blas_int i) const
    {
        assert(init());
        assert(i < xbasis_size_);
      //
        return xgrid_.xr(i);
    }
    const dfloat& FemDvrEcsGrid2d::yr(blas_int i) const
    {
        assert(init());
        assert(i < ybasis_size_);
      //
        return ygrid_.xr(i);
    }
    const dcomp& FemDvrEcsGrid2d::xz(blas_int i) const
    {
        assert(init());
        assert(i < xbasis_size_);
      //
        return xgrid_.x(i);
    }
    const dcomp& FemDvrEcsGrid2d::yz(blas_int i) const
    {
        assert(init());
        assert(i < ybasis_size_);
      //
        return ygrid_.x(i);
    }
    const FemDvrEcsGrid& FemDvrEcsGrid2d::get_xgrid() const
    {
        assert(init());
      //
        return xgrid_;
    }
    const FemDvrEcsGrid& FemDvrEcsGrid2d::get_ygrid() const
    {
        assert(init());
      //
        return ygrid_;
    }
    dcomp FemDvrEcsGrid2d::wz(blas_int i, blas_int j) const
    {
        assert(init());
        assert(i < xbasis_size_);
        assert(j < ybasis_size_);
      //
        return xgrid_.w(i) * ygrid_.w(j);
    }
    dcomp FemDvrEcsGrid2d::wz(blas_int i) const
    {
        assert(init());
        assert(i < xbasis_size_ * ybasis_size_);
      //
        blas_int l = i/xbasis_size_;
        blas_int k = i % xbasis_size_;
        return xgrid_.w(k) * ygrid_.w(l);
    }
    const dcomp& FemDvrEcsGrid2d::xwz(blas_int i) const
    {
        assert(init());
        assert(i < xbasis_size_);
      //
        return xgrid_.w(i);
    }
    const dcomp& FemDvrEcsGrid2d::ywz(blas_int i) const
    {
        assert(init());
        assert(i < ybasis_size_);
      //
        return ygrid_.w(i);
    }
    blas_int FemDvrEcsGrid2d::x_element_end(const dfloat& X) const
    {
        assert(init());
      //
        return xgrid_.get_element_end_x(X);
    }
    blas_int FemDvrEcsGrid2d::y_element_end(const dfloat& X) const
    {
        assert(init());
      //
        return ygrid_.get_element_end_x(X);
    }
    blas_int FemDvrEcsGrid2d::get_real_xsize() const
    {
        assert(init());
      //
        return xgrid_.nr();
    }
    blas_int FemDvrEcsGrid2d::get_real_ysize() const
    {
        assert(init());
      //
        return ygrid_.nr();
    }

  // modifiers
    FemDvrEcsGrid2d& FemDvrEcsGrid2d::swap(FemDvrEcsGrid2d& rhs)
    {
        Object::swap(rhs);
        xgrid_.swap(rhs.xgrid_);
        ygrid_.swap(rhs.ygrid_);
        std::swap(xbasis_size_, rhs.xbasis_size_);
        std::swap(ybasis_size_, rhs.ybasis_size_);
        std::swap(basis_size_,  rhs.basis_size_);
        return *this;
    }
    FemDvrEcsGrid2d& FemDvrEcsGrid2d::operator= (FemDvrEcsGrid2d tmp)
    {
        swap(tmp);
        return *this;
    }

    bool FemDvrEcsGrid2d::operator== (const FemDvrEcsGrid2d& rhs) const
    {
        bool stat = (xgrid_ == rhs.xgrid_);
        if (stat) stat = (ygrid_ == rhs.ygrid_);
        return stat;
    }
    bool FemDvrEcsGrid2d::operator!= (const FemDvrEcsGrid2d& rhs) const
    {
        return !((*this)==rhs);
    }
}   // QSCAT
