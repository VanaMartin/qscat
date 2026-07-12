#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <complex>
#include <string>
#include <cmath>

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
    EquidistantProjector2d::EquidistantProjector2d(const FemDvrEcsGrid2d& grid, size_t x_sampling, size_t y_sampling, dfloat xa, dfloat xb, dfloat ya, dfloat yb) : Object()
    {
     //
        //assert(grid.get_xgrid().xr(0) <= xa);
        assert(grid.get_xgrid().xr(grid.get_xsize()-1) >= xb);
        //assert(grid.get_ygrid().xr(0) <= ya);
        assert(grid.get_ygrid().xr(grid.get_ysize()-1) >= yb);
        assert(x_sampling > 0);
        assert(y_sampling > 0);
      //
        grid_ = grid;
        x_samples_ = x_sampling;
        y_samples_ = y_sampling;
        values_ = zVector(x_samples_*y_samples_);
        values_.fill(0);
        x_start_ = xa;
        x_end_ = xb;
        y_start_ = ya;
        y_end_ = yb;

        x_coordinate_ = dVector(x_samples_, xa, xb, true);
        y_coordinate_ = dVector(y_samples_, ya, yb, true);

        const FemDvrEcsGrid& px = grid_.get_xgrid();   // Auxiliary pointer;
        const FemDvrEcsGrid& py = grid_.get_ygrid();   // Auxiliary pointer;

        size_t nnz = px.quadrature() * py.quadrature() * x_samples_ * y_samples_;  // each point gets whole element contributions

        size_t rows = x_samples_ * y_samples_;  // result size
        size_t columns = grid_.get_size();      // input size

        body_ = RowCompressedMatrix<dcomp>(rows, columns, nnz);

        size_t position = 0;
        blas_int nqx = px.quadrature();
        blas_int nqy = py.quadrature();

        zVector y_hlp(nqy);
        zVector x_hlp(nqx);

        // TODO Transpose output XxY

        for (size_t pos_y=0; pos_y<y_samples_; ++pos_y) {
            dfloat y = y_coordinate_[pos_y];
            blas_int el_y = py.get_element_index(y);
            blas_int el_ys = py.get_element_start(el_y);
            blas_int el_ye = py.get_element_end(el_y);
            for (blas_int i=el_ys; i<=el_ye; ++i){
                y_hlp[i - el_ys] = py.basis_function_value(i,y,el_ys,el_ye);
            }
            for (size_t pos_x=0; pos_x<x_samples_; ++pos_x) {
                body_.row_index(pos_y*x_samples_ + pos_x) = position;    // starting index of given row

                dfloat x = x_coordinate_[pos_x];
                blas_int el_x = px.get_element_index(x);
                blas_int el_xs = px.get_element_start(el_x);
                blas_int el_xe = px.get_element_end(el_x);
                for (blas_int i=el_xs; i<=el_xe; ++i){
                    x_hlp[i - el_xs] = px.basis_function_value(i,x,el_xs,el_xe);
                }
                //for (blas_int i=el_ys; i<=el_ye; ++i){
                //    for (blas_int j=el_xs; j<=el_xe; ++j){
                for (blas_int i=0; i<nqy; ++i){
                    for (blas_int j=0; j<nqx; ++j){
                        //if (i < get_ysize() && j < get_xsize())
                        body_.nonzeros(position) = y_hlp[i] * x_hlp[j];
                        body_.columns(position) = (el_ys+i) * grid_.get_xsize() + el_xs + j;
                        position++;
                    }
                }
            }
        }
        body_.row_index(y_samples_ * x_samples_) = position;

        *init_ = true;
        incref();
    }
    EquidistantProjector2d::~EquidistantProjector2d()
    {
        *init_ = false;
        decref();
    }
    void EquidistantProjector2d::operator<< (const GridVector2d& state) const
    {
        assert(init());
        assert(state.get_grid() == grid_);
      //
        body_.gemv(1.0, state.body(), 0.0, values_);
    }

    void EquidistantProjector2d::export_state_hsv(GridVector2d& state, const char *filename, dfloat magnitude, dfloat phase) const
    {
        *this << state;
        export_state_hsv(filename, magnitude, phase);
    }
    void EquidistantProjector2d::export_state(GridVector2d& state, const char* filename) const
    {
        *this << state;
        export_state(filename);
    }

    void EquidistantProjector2d::set_zoom_filter(dfloat x_start, dfloat x_end, dfloat y_start, dfloat y_end, dfloat magnitude)
    {
        blas_int xa = 0, xb = x_samples_-1, ya = 0, yb = y_samples_-1;
        while (x_start > x_coordinate_[xa] && xa < x_samples_) ++xa;
        while (x_end < x_coordinate_[xb-1] && xb > 0) --xb;
        while (y_start > y_coordinate_[ya] && ya < y_samples_) ++ya;
        ya = y_samples_ - 1 - ya;
        while (y_end < y_coordinate_[yb] && yb > 0) --yb;
        yb = y_samples_ - 2 - yb;

        zoom_ = ZoomFilter(xa, xb, ya, yb, magnitude);
    }

    void EquidistantProjector2d::export_state(const char* filename) const
    {
        FILE * file;
        fopen_s(&file, filename,"w");

        //fprintf(file,"#Function on 2d range (%f,%f)x(%f,%f) enumerated on %dx%d points.\n", x_min, x_max, y_min, y_max, x_end - x_start, y_end - y_start);
        //fprintf(file,"#Coordinate X   \tCoordinate Y   \tReal part of z     \tImaginary part of z\n");
        for (blas_int i=0; i<y_samples_; ++i){
            for (blas_int j=0; j<x_samples_; ++j){
                dcomp v = values_[i*x_samples_ + j];
                fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\n", y_coordinate_[i], x_coordinate_[j], real(v), imag(v), pow(abs(v),2), arg(v));
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }

    void EquidistantProjector2d::export_state_hsv( const char *filename, dfloat magnitude, dfloat phase) const
    {

      // TODO clean this method

        std::ofstream file;
        file.open(filename, std::ios::out | std::ios::binary );

        // write header
        char aux[4] = {'B', 'M', 0, 0};
        file.write((char*) aux, 2);

        blas_int row_size = ( ( 24 * x_samples_ + 31) / 32 ) * 4;
        blas_int h_size = 14 + 12;
        blas_int b_size = h_size + row_size * y_samples_;

        file.write((char*) &b_size, 4);
        file.write((char*) aux, 4);
        file.write((char*) &h_size, 4);
        // DIB header
        h_size = 12;
        short pix_size = x_samples_;
        file.write((char*) &h_size, 4);
        file.write((char*) &pix_size, 2);
        pix_size = y_samples_;
        file.write((char*) &pix_size, 2);
        pix_size = 1;
        file.write((char*) &pix_size, 2);
        pix_size = 24;
        file.write((char*) &pix_size, 2);

        char *row = new char[row_size];
        memset(row, 0, row_size);

        double z, h, s, v, f;
        blas_int ii;
        char p, q, t, vc;
        for (blas_int i=0; i<y_samples_; ++i){
            char* rp = row;
            for (blas_int j=0; j<x_samples_; ++j){
                def_comp val = values_[i*x_samples_ + j] * std::exp(imu * phase);
                double mag;
                if (zoom_.init()) {
                   mag = magnitude;
                } else {
                    if ( i >= zoom_.ya && i <= zoom_.yb && j >= zoom_.xa && j <= zoom_.xb)
                        mag = zoom_.mag;
                    else
                        mag = magnitude;
                }
                double rv = abs(values_[i*x_samples_ + j + ((j<x_samples_-1)?1:0)])/mag;
                double lv = abs(values_[i*x_samples_ + j - ((j>0)?1:0)]) / mag;
                double tv = abs(values_[(i+((i<y_samples_-1)?1:0))*x_samples_ + j])/mag;
                double bv = abs(values_[(i-((i>0)?1:0))*x_samples_ + j]) / mag;
                double xv = abs(val) / mag;

                h = (arg(val) + pi) / (2*pi);                                         // hue
                //s = 1.0 - std::atan(max(abs(val) / mag - 1.0, 0.0)) / pi * 2.;  // saturation
                //v = min(abs(val) / mag, 1.0);                                   // value

                z = 1 - 1 / ( std::pow(abs(val)/mag,2) + 1);
                s = min( 2. * (1. - z), 1.0 );
                v = min( 2. * z, 1.0 );

                double n=10.;
                int levels = min(abs( int(xv*n) - int(lv*n) ) + abs( int(xv*n) - int(rv*n) ) + abs(int(xv*n) - int(tv*n)) + abs(int(xv*n)-int(bv*n)), 1);
                double correction = 0.25 * double(levels);

                s = max( s - correction, 0.0);

                if (s == 0.0) {
                    for (ii=0; ii<3; ++ii)
                        *rp++ = 255;
                } else {
                    ii = h * 6;
                    f = h*6 - ii;
                    p = int(255 * v * (1. - s));
                    q = int(255 * v * (1. - s*f));
                    t = int(255 * v * (1. - s*(1.-f)));
                    vc = int(255 * v);
                    ii %= 6;
                    switch (ii) {
                        case 1:     // v, t, p
                            *rp++ = vc;
                            *rp++ = t;
                            *rp++ = p;
                            break;
                        case 2:     // q, v, p
                            *rp++ = q;
                            *rp++ = vc;
                            *rp++ = p;
                            break;
                        case 3:     // p, v, t
                            *rp++ = p;
                            *rp++ = vc;
                            *rp++ = t;
                            break;
                        case 4:     // p, q, v
                            *rp++ = p;
                            *rp++ = q;
                            *rp++ = vc;
                            break;
                        case 5:     // t, p, v
                            *rp++ = t;
                            *rp++ = p;
                            *rp++ = vc;
                            break;
                        default:    // v, p, q
                            *rp++ = vc;
                            *rp++ = p;
                            *rp++ = q;

                    }
                }
            }
            file.write(row, row_size);
        }
        file.close();
    }
} // namespace QSCAT
