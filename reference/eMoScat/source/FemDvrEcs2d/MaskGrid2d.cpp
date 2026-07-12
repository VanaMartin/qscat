#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <complex>
#include <string>
#include <limits>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "bessel.h"
#include "common.h"
#include "blas.h"
#include "Arrays.h"
#include "input.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"

namespace QSCAT
{
    // constructors
    MaskGrid2d::MaskGrid2d()
    {
        init_ = false;
    }
    MaskGrid2d::MaskGrid2d(const FemDvrEcsGrid2d& grid, const iMatrix& mask)
    {
        assert( grid.init() );
        assert( mask.columns() == grid.get_xgrid().tnel() );
        assert( mask.rows() == grid.get_ygrid().tnel() );
      //
        grid_ = grid;
        mask_ = mask;

        x_size_ = grid.get_xgrid().nb();
        y_size_ = grid.get_ygrid().nb();
        full_size_ = x_size_ * y_size_;

        x_elements_ = grid_.get_xgrid().tnel();
        y_elements_ = grid_.get_ygrid().tnel();
        x_quadrature_ = grid_.get_xgrid().quadrature();
        y_quadrature_ = grid_.get_ygrid().quadrature();

      // first determine the size of the compression
        compressed_size_ = 0;
        total_elements_ = 0;
        blas_int xp, yp, xpyp; // connection controllers

        for (blas_int i=0; i<y_elements_; ++i) {
            for (blas_int j=0; j<x_elements_; ++j) {
                if (mask(i,j)!=0) {    // if the element is used
                    yp = (i>0)? mask(i-1,j) : 0;                // connected to previous?
                    xp = (j>0)? mask(i,j-1) : 0;                // connected to previous?
                    xpyp = (xp && yp)? mask_(i-1, j-1) : 0;     // corner basis function is present

                  // We treat the element as not connected to the next, if connected the next element includes the basis function
                    compressed_size_ += (x_quadrature_ - 2 + xp) * (y_quadrature_ - 2 + yp) - ((xp && yp)?  1-xpyp : 0);
                    total_elements_++;
                }
            }
        }

      // second determine the column positions and row sizes
        columns_ = Vector<blas_int>(compressed_size_);   // column index in full representation
        row_index_ = Vector<blas_int>(y_size_ + 1);      // how many basis functions in each row

        blas_int idc=0;
        blas_int idr=0;
        row_index_[0] = 0;
        for (blas_int i=0; i<y_elements_; ++i) {
            for (blas_int qy= (i)? 0:1; qy < y_quadrature_-1; ++qy) {    // for all quadrature point except the last (to be solved in next element (ommits first row)
                idr++;
                row_index_[idr] = row_index_[idr-1];
                for (blas_int j=0; j<x_elements_; ++j) {
                    if (mask(i,j)!=0) {                             // if the element is used
                        xp = (j>0)? mask(i,j-1) : 0;                // connected to x-previous?
                        yp = (i>0)? mask(i-1,j) : 0;                // connected to y-previous?
                        xpyp = (xp && yp)? mask_(i-1, j-1) : 0;     // corner point is present

                        if (qy == 0) {
                            if (yp!=0) {
                                for (blas_int qx=1-xpyp; qx<(x_quadrature_-1); ++qx) {
                                    columns_[idc++] = j * (x_quadrature_ - 1) + qx - 1;     // minus one for the very first element
                                    row_index_[idr]++;
                                }
                            }
                        } else {
                            for (blas_int qx=1-xp; qx<x_quadrature_-1; ++qx) {
                                columns_[idc++] = j * (x_quadrature_ - 1) + qx - 1;         // minus one for the very first element
                                row_index_[idr]++;
                            }
                        }
                    }
                }
            }
        }

      //
        assert(idc == compressed_size_);
      //
        cout << "Compressed grid built: compression ratio = " << double(compressed_size_) / double(full_size_) << endl;
        init_=true;
    }
    MaskGrid2d::MaskGrid2d(const MaskGrid2d& old) :
        grid_(old.grid_),
        mask_(old.mask_),
        x_size_(old.x_size_),
        y_size_(old.y_size_),
        full_size_(old.full_size_),
        row_index_(old.row_index_),
        columns_(old.columns_),
        x_elements_(old.x_elements_),
        y_elements_(old.y_elements_),
        x_quadrature_(old.x_quadrature_),
        y_quadrature_(old.y_quadrature_),
        total_elements_(old.total_elements_),
        compressed_size_(old.compressed_size_),
        init_(old.init_)
    {}


    // accessors
    bool MaskGrid2d::init() const
    {
        return init_;
    }
    blas_int MaskGrid2d::get_size() const
    {
        assert(init_);
      //
        return compressed_size_;
    }
    blas_int MaskGrid2d::get_xsize() const
    {
        assert(init_);
      //
        return grid_.get_xsize();
    }
    blas_int MaskGrid2d::get_ysize() const
    {
        assert(init_);
      //
        return grid_.get_ysize();
    }
    const FemDvrEcsGrid2d& MaskGrid2d::full_grid() const
    {
        assert(init_);
      //
        return grid_;
    }
    const iVector& MaskGrid2d::get_columns() const
    {
        assert(init_);
      //
        return columns_;
    }
    const iVector& MaskGrid2d::get_row_index() const
    {
        assert(init_);
      //
        return row_index_;
    }
    const iMatrix& MaskGrid2d::get_mask() const
    {
        assert(init_);
      //
        return mask_;
    }
    dcomp MaskGrid2d::wz(blas_int i) const
    {
        assert( init_ );
      //
        blas_int row = 0;
        while (row_index_[row+1] < i) row++;
        return grid_.wz( row * grid_.get_xsize() + columns_[i] );
    }

    // modifiers
    MaskGrid2d& MaskGrid2d::swap(MaskGrid2d& rhs)
    {
        grid_.swap(rhs.grid_);
        mask_.swap(rhs.mask_);
        std::swap(x_size_,rhs.x_size_);
        std::swap(y_size_,rhs.y_size_);
        std::swap(full_size_,rhs.full_size_);
        row_index_.swap(rhs.row_index_);
        columns_.swap(rhs.columns_);
        std::swap(x_elements_,rhs.x_elements_);
        std::swap(y_elements_,rhs.y_elements_);
        std::swap(x_quadrature_,rhs.x_quadrature_);
        std::swap(y_quadrature_,rhs.y_quadrature_);
        std::swap(total_elements_,rhs.total_elements_);
        std::swap(compressed_size_,rhs.compressed_size_);
        std::swap(init_,rhs.init_);

        return *this;
    }


    // operators
    MaskGrid2d& MaskGrid2d::operator= (MaskGrid2d tmp)
    {
        return this->swap(tmp);
    }
    bool MaskGrid2d::operator== (const MaskGrid2d& rhs) const
    {
        assert( init_ );
        assert( rhs.init_ );
      //
        blas_int err=0;

        if ( x_size_ != rhs.x_size_ ) err++;
        if ( y_size_ != rhs.y_size_ ) err++;
        if ( compressed_size_ != rhs.compressed_size_ ) err++;
        if ( total_elements_ != rhs.total_elements_ ) err++;
        if ( x_quadrature_ != rhs.x_quadrature_ ) err++;
        if ( y_quadrature_ != rhs.y_quadrature_ ) err++;
        if ( x_elements_ != rhs.x_elements_ ) err++;
        if ( y_elements_ != rhs.y_elements_ ) err++;

        // TODO use hash on mask ----- most important part of comparison

        return (err==0)? true : false;
    };

}   // namespace QSCAT
