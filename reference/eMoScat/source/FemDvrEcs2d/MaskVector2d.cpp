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
    // Constructors
    MaskVector2d::MaskVector2d()
    {
        init_ = false;
    }
    MaskVector2d::MaskVector2d(const MaskGrid2d& grid)
    {
        assert( grid.init() );
      //
        grid_ = grid;
        body_ = zVector(grid_.get_size());

        init_ = true;
    }
    MaskVector2d::MaskVector2d(const MaskGrid2d& grid, const GridVector& x, const GridVector& y)
    {
        assert( grid.init() );
        assert( x.init() );
        assert( y.init() );
        assert( x.get_size() == grid.get_xsize() );
        assert( y.get_size() == grid.get_ysize() );
      //
        grid_ = grid;
        body_ = zVector( grid_.get_size() );

        blas_int i=0;
        const iVector& columns = grid_.get_columns();
        const iVector& row_index = grid_.get_row_index();
        for (blas_int j=0; j<grid_.get_size(); ++j) {
            if (j == row_index[i+1] ) ++i;
            body_[j] = x[ columns[j] ] * y[i];
        }

        init_ = true;
    }
    MaskVector2d::MaskVector2d(const MaskVector2d& old) :
        init_(old.init_),
        body_(old.body_),
        grid_(old.grid_)
    {}


    // Accessors
    bool MaskVector2d::init() const
    {
        return init_;
    }
    blas_int MaskVector2d::get_size() const
    {
        return grid_.get_size();
    }
    const MaskGrid2d& MaskVector2d::get_grid() const
    {
        assert( init_ );
      //
        return grid_;
    }
    const zVector& MaskVector2d::body() const
    {
        assert(init_);
      //
        return body_;
    }
    zVector& MaskVector2d::body()
    {
        assert( init_ );
      //
        return body_;
    }
    dcomp MaskVector2d::f(blas_int i) const
    {
        assert(init_);
        assert(i < body_.get_size() );
      //
        return body_[i] / sqrt( grid_.wz(i) );   // TODO move to grid_.wz(i)
    }
    void MaskVector2d::f(const dcomp& val,  blas_int i)
    {
        assert( init_ );
        assert( i < body_.get_size() );
      //
        body_[i] = val * sqrt( grid_.wz(i) );   // TODO move to grid_.wz(i)
    }
    zVector MaskVector2d::function_values() const
    {
        assert(init_);
      //
        zVector out(get_size());
        for (blas_int i=0; i<get_size(); ++i)
            out[i] = f(i);
        return out;
    }

    // Modifiers
    GridVector2d& MaskVector2d::get(GridVector2d& dst) const
    {
        assert(init_);
        assert(dst.init());
        assert(dst.get_grid() == grid_.full_grid());    // consistencty check
      //
        const iVector& cols = grid_.get_columns();
        const iVector& row_index = grid_.get_row_index();

        blas_int row_size = grid_.full_grid().get_xsize();

        for (blas_int i=0; i<row_index.get_size()-1; ++i) {
            for (blas_int j=row_index[i]; j<row_index[i+1]; ++j) {
                dst[i*row_size + cols[j]] = body_[j];
            }
        }

        return dst;
    }
    MaskVector2d& MaskVector2d::set(const GridVector2d& src)
    {
        assert(init_);
        assert(src.init());
        assert(src.get_grid() == grid_.full_grid());
      //

        const iVector& cols = grid_.get_columns();
        const iVector& row_index = grid_.get_row_index();

        blas_int row_size = grid_.full_grid().get_xsize();

        for (blas_int i=0; i<row_index.get_size()-1; ++i) {
            for (blas_int j=row_index[i]; j<row_index[i+1]; ++j) {
                body_[j] = src[i*row_size + cols[j]];
            }
        }

        return *this;
    }
    MaskVector2d& MaskVector2d::swap(MaskVector2d& rhs)
    {
        std::swap(init_, rhs.init_);
        body_.swap(rhs.body_);
        grid_.swap(rhs.grid_);
        return *this;
    }
    MaskVector2d& MaskVector2d::copy(const MaskVector2d& src)
    {
        init_ = src.init_;
        body_ = src.body_;
        grid_ = src.grid_;
        return *this;
    }


    // Operators
    dcomp MaskVector2d::operator* (const MaskVector2d& rhs) const
    {
        assert( init_ );
        assert( rhs.init_ );
        assert( grid_ == rhs.grid_ );
      //
        return body_ * rhs.body_;
    }

    void MaskVector2d::save(const char*filename) const
    {
        assert(init_);
      //
        FILE * file;
        dcomp v;
        dfloat aux = 0.; //std::numeric_limits<def_float>::quiet_NaN();

        fopen_s(&file,filename,"w");
        fprintf(file,"#Function on FEM-DVR-ECS grid of %lld basis functions.\n", get_size());
        fprintf(file,"#Coordinate X   \tCoordinate Y   \tReal part of z     \tImaginary part of z\n");
        const iVector& columns = grid_.get_columns();
        const iVector& row_index = grid_.get_row_index();
        const FemDvrEcsGrid2d& grid = grid_.full_grid();

        blas_int nbx = grid_.get_xsize();
        blas_int nby = grid_.get_ysize();
        blas_int size = grid_.get_size();
        blas_int pos=0;
        blas_int y = 0;
        for (blas_int i=0; i<nby; ++i) {
            for (blas_int j=0; j<nbx; ++j) {
                if (pos<size && j==columns[pos]) {
                    v = this->f(pos++);
                    fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\n", grid.yr(i), grid.xr(j), real(v), imag(v)); //
                } else {
                    //fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\n", grid.yr(i), grid.xr(j), aux, aux); //
                    fprintf(file, "%.12E\t%.12E\t%s\t%s\n", grid.yr(i), grid.xr(j), "nan", "nan"); //
                }
            }
            fprintf(file, "\n");
        }

        //for (int i=0; i<row_index.get_size()-1; ++i) {
        //    for (int j=row_index[i]; j<row_index[i+1]; ++j) {
        //        v = this->f(j);
        //        fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\n", grid.yr(i), grid.xr(columns[j]), real(v), imag(v)); //
        //    }
        //    fprintf(file, "\n");
        //}

        fclose(file);
    }
//    template<>
//    void GridVector2d<def_float,def_comp>::save_equidistant(const char *name, int x_points, const double& x_min, const double& x_max, int y_points, const double& y_min, const double& y_max) const
//    {
//        assert(init_);
//      //
//        FILE * file;
//        def_comp v;
//        fopen_s(&file,name,"w");
//
//    // Auxiliary variables declarations
//        Vector<double> Yr(y_points+1,y_min,y_max, true);
//        Vector<double> Xr(x_points+1,x_min,x_max, true);
//
//
//        int x_start, x_end, y_start, y_end;
//        (x_min <= grid_.xr(0))? x_start = 1: x_start = 0;
//        (x_max >= grid_.xr(get_xsize()-1))? x_end = x_points-1: x_end = x_points;
//        (y_min <= grid_.yr(0))? y_start = 1: y_start = 0;
//        (y_max >= grid_.yr(get_ysize()-1))? y_end = y_points-1: y_end = y_points;
//        //fprintf(file,"#Function on 2d range (%f,%f)x(%f,%f) enumerated on %dx%d points.\n", x_min, x_max, y_min, y_max, x_end - x_start, y_end - y_start);
//        //fprintf(file,"#Coordinate X   \tCoordinate Y   \tReal part of z     \tImaginary part of z\n");
//        for (int i=0;i<=y_points;++i){
//            if (i<y_start || i > y_end){
//                for (int j=0;j<=x_points;++j){
//                    fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\n", Yr[i], Xr[j], 0.0, 0.0, 0.0, 0.0);
//                }
//            } else {
//                for (int j=0;j<=x_points;++j){
//                    if (j<x_start || j > x_end){
//                        fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\n", Yr[i], Xr[j], 0.0, 0.0, 0.0, 0.0);
//                    } else {
//                        v = evaluate(Xr[j],Yr[i]);
//                        fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\n",Yr[i],Xr[j], real(v), imag(v), pow(abs(v),2), arg(v));
//                    }
//                }
//
//            }
//            fprintf(file, "\n");
//        }
//        fclose(file);
//    }
} // namespace QSCAT
