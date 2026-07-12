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
    MaskOperator2d::MaskOperator2d()
    {
        init_ = false;
    }
    MaskOperator2d::MaskOperator2d(const MaskGrid2d& grid)
    {
        assert(grid.init());
      //
        grid_ = grid;

        init_ = true;
    }
    // Accessors
    bool MaskOperator2d::init() const
    {
        return init_;
    }
    const MaskGrid2d& MaskOperator2d::get_grid() const
    {
        assert(init_);
      //
        return grid_;
    }
    // Modifiers
    MaskOperator2d& MaskOperator2d::set_kinetic_term(dfloat mu_x, dfloat mu_y)
    {
        assert(init_);
      //
        const FemDvrEcsGrid& xgrid = grid_.full_grid().get_xgrid();
        const FemDvrEcsGrid& ygrid = grid_.full_grid().get_ygrid();
        RowCompressedMatrix<dcomp> Tx = generateKineticTermRCM(xgrid, mu_x);
        RowCompressedMatrix<dcomp> Ty = generateKineticTermRCM(ygrid, mu_y);

        RowCompressedMatrix<dcomp> ref = TensorSum(Tx, Ty);

        blas_int tnx = xgrid.tnel();
        blas_int tny = ygrid.tnel();

        blas_int nbx = xgrid.nb();
        blas_int nby = ygrid.nb();

        blas_int nqx = xgrid.quadrature();
        blas_int nqy = ygrid.quadrature();

      // for gathering compression statistics
        blas_int max = Tx.num_nonzeros() * nby + Ty.num_nonzeros() * nbx - nbx * nby;

      // first determine the total of nonzero values
        blas_int px, py; // previous
        blas_int nx, ny; // next

        blas_int nnz = 0;
        blas_int size = 0;

        const iVector& cols = grid_.get_columns();
        const iVector& row_index = grid_.get_row_index();
        const iMatrix& mask = grid_.get_mask();

        size = grid_.get_size();

      // build body

        blas_int count = 0;
        blas_int y=0;
        for (blas_int b=0; b<size; ++b) {
            // actuall indices of position in grid: x,y
            if (row_index[y+1] == b) ++y;
            blas_int x = cols[b];

            blas_int idx = x / (nqx - 1);
            blas_int idy = y / (nqy - 1);

            blas_int bx = ( (x+1) % (nqx-1) == 0 )? 1 : 0;   // x-border
            blas_int by = ( (y+1) % (nqy-1) == 0 )? 1 : 0;   // y-border

            blas_int ysize = Ty.row_index(y+1) - Ty.row_index(y);
            blas_int xsize = Tx.row_index(x+1) - Tx.row_index(x);
            //int ysize = Ty.row_index(y+1) - Ty.row_index(y) + ((idy==0)? 1:0);
            //int xsize = Tx.row_index(x+1) - Tx.row_index(x) + ((idx==0)? 1:0);

            /*
                cross border: cE resp. rE - marks the columns
                resp. rows that surely exists (otherwise the
                basis function would not be defined at all)

                NOTE: if the border points exist, both elements
                have to be present! :)
            */
            if (bx && by) {
                /*     cross border:

                          cE     cE
                              +       : idy-1
                              |
                          idx | idx+1 : idy     rE
                              |
                  idx-1 +-----+-----+  idx+2
                              |
                          idx | idx+1 : idy+1   rE
                              |
                              +       : idy+2
                */
                if (idy > 0 && ( mask(idy-1, idx)==0 || mask(idy-1,idx+1)==0 ) ) ysize--;
                if (idx > 0 && ( mask(idy, idx-1)==0 || mask(idy+1, idx-1)==0 ) ) xsize--;
                if (idy < tny-2 && ( mask(idy+2, idx)==0 || mask(idy+2,idx+1)==0 ) ) ysize--;
                if (idx < tnx-2 && ( mask(idy, idx+2)==0 || mask(idy+1,idx+2)==0 ) ) xsize--;
            } else if (bx) {
                /*      x-only border

                          cE     cE
                          idx + idx+1       : idy-1
                              |
                  idx-1 +.....|.....+ idx+2 : idy   rE
                              |
                          idx + idx+1       : idy+1
                */
                if (idy > 0 && ( mask(idy-1, idx)==0 || mask(idy-1,idx+1)==0 ) ) ysize--;
                if (idx > 0 && mask(idy, idx-1)==0 ) xsize--;
                if (idy < tny-1 && ( mask(idy+1, idx)==0 || mask(idy+1, idx+1)==0 ) ) ysize--;
                if (idx < tnx-2 && mask(idy, idx+2)==0 ) xsize--;
            } else if (by) {
                /*      y-only border

                          cE
                          idx
                           +         : idy-1
                           .
                           .         : idy      rE
                           .
                  idx-1 +-----+ idx+1
                           .
                           .         : idy+1    rE
                           .
                           +         : idy+2
                */
                if (idy > 0 && mask(idy-1, idx)==0 ) ysize--;
                if (idx > 0 && ( mask(idy, idx-1)==0 || mask(idy+1, idx-1)==0 ) ) xsize--;
                if (idy < tny-2 && mask(idy+2, idx)==0 ) ysize--;
                if (idx < tnx-1 && ( mask(idy, idx+1)==0 || mask(idy+1, idx)==0 ) ) xsize--;
            } else {
                /*      internal

                          cE
                          idx
                           +            : idy-1
                           .
                  idx-1 +.....+ idx+1   : idy   rE
                           .
                           +            : idy+1
                */
                if (idy > 0 && mask(idy-1, idx)==0 ) ysize--;
                if (idx > 0 && mask(idy, idx-1)==0 ) xsize--;
                if (idy < tny-1 && mask(idy+1, idx)==0 ) ysize--;
                if (idx < tnx-1 && mask(idy, idx+1)==0 ) xsize--;
            }
            count += xsize + ysize - 1;
        }

        cout << "Kinetic term compression ratio " << count << "/" << max << " = " << float(count) / float(max) << endl;

        body_ = RowCompressedMatrix<dcomp>(size, size, count);

        assert( body_.num_nonzeros() == count);

        blas_int pos = 0;
        blas_int col = 0;

        blas_int offset = 0; // points to the beginning of current element-row in compressed basis

        body_.row_index(0) = 0;
        y = 0;
        for (blas_int b=0; b<size; ++b) {
            // actuall indices of position in grid: x,y
            if (row_index[y+1] == b) ++y;
            blas_int x = cols[b];

            blas_int idx = x / (nqx - 1);
            blas_int idy = y / (nqy - 1);

            blas_int bx = ( (x+1) % (nqx-1) == 0 )? 1 : 0;   // x-border
            blas_int by = ( (y+1) % (nqy-1) == 0 )? 1 : 0;   // y-border

            blas_int xstart = Tx.row_index(x);
            blas_int xend = Tx.row_index(x+1);
            //int xstart = (idx==0)? -1 : idx*(nqx-1)-1;       // start on the full x-grid
            //int xend = min(nbx-1, (idx+1+bx)*(nqx-1)-1);    // end on the full x-grid

            blas_int ystart = Ty.row_index(y);
            blas_int yend = Ty.row_index(y+1);
            //int ystart = (idy==0)? -1 : idy*(nqy-1)-1;       // start on the full y-grid
            //int yend = min(nby-1, (idy+1+by)*(nqy-1)-1);    // end on the full y-grid

            blas_int ysize = Ty.row_index(y+1) - Ty.row_index(y);
            blas_int xsize = Tx.row_index(x+1) - Tx.row_index(x);

            if (bx && by) {
                /*     cross border:

                          cE     cE
                              +       : idy-1
                              |
                          idx | idx+1 : idy     rE
                              |
                  idx-1 +-----+-----+  idx+2
                              |
                          idx | idx+1 : idy+1   rE
                              |
                              +       : idy+2
                */
                if (idy > 0 && ( mask(idy-1, idx)==0 || mask(idy-1,idx+1)==0 ) ) ystart++;
                if (idx > 0 && ( mask(idy, idx-1)==0 || mask(idy+1, idx-1)==0 ) ) xstart++;
                if (idy < tny-2 && ( mask(idy+2, idx)==0 || mask(idy+2,idx+1)==0 ) ) yend--;
                if (idx < tnx-2 && ( mask(idy, idx+2)==0 || mask(idy+1,idx+2)==0 ) ) xend--;

                if (idy > 0 && ( mask(idy-1, idx)==0 || mask(idy-1,idx+1)==0 ) ) ysize--;
                if (idx > 0 && ( mask(idy, idx-1)==0 || mask(idy+1, idx-1)==0 ) ) xsize--;
                if (idy < tny-2 && ( mask(idy+2, idx)==0 || mask(idy+2,idx+1)==0 ) ) ysize--;
                if (idx < tnx-2 && ( mask(idy, idx+2)==0 || mask(idy+1,idx+2)==0 ) ) xsize--;
            } else if (bx) {
                /*      x-only border

                          cE     cE
                          idx + idx+1       : idy-1
                              |
                  idx-1 +.....|.....+ idx+2 : idy   rE
                              |
                          idx + idx+1       : idy+1
                */
                if (idy > 0 && ( mask(idy-1, idx)==0 || mask(idy-1,idx+1)==0 ) ) ystart++;
                if (idx > 0 && mask(idy, idx-1)==0 ) xstart++;
                if (idy < tny-1 && ( mask(idy+1, idx)==0 || mask(idy+1, idx+1)==0 ) ) yend--;
                if (idx < tnx-2 && mask(idy, idx+2)==0 ) xend--;

                if (idy > 0 && ( mask(idy-1, idx)==0 || mask(idy-1,idx+1)==0 ) ) ysize--;
                if (idx > 0 && mask(idy, idx-1)==0 ) xsize--;
                if (idy < tny-1 && ( mask(idy+1, idx)==0 || mask(idy+1, idx+1)==0 ) ) ysize--;
                if (idx < tnx-2 && mask(idy, idx+2)==0 ) xsize--;
            } else if (by) {
                /*      y-only border

                          cE
                          idx
                           +         : idy-1
                           .
                           .         : idy      rE
                           .
                  idx-1 +-----+ idx+1
                           .
                           .         : idy+1    rE
                           .
                           +         : idy+2
                */
                if (idy > 0 && mask(idy-1, idx)==0 ) ystart++;
                if (idx > 0 && ( mask(idy, idx-1)==0 || mask(idy+1, idx-1)==0 ) ) xstart++;
                if (idy < tny-2 && mask(idy+2, idx)==0 ) yend--;
                if (idx < tnx-1 && ( mask(idy, idx+1)==0 || mask(idy+1, idx)==0 ) ) xend--;

                if (idy > 0 && mask(idy-1, idx)==0 ) ysize--;
                if (idx > 0 && ( mask(idy, idx-1)==0 || mask(idy+1, idx-1)==0 ) ) xsize--;
                if (idy < tny-2 && mask(idy+2, idx)==0 ) ysize--;
                if (idx < tnx-1 && ( mask(idy, idx+1)==0 || mask(idy+1, idx)==0 ) ) xsize--;
            } else {
                /*      internal

                          cE
                          idx
                           +            : idy-1
                           .
                  idx-1 +.....+ idx+1   : idy   rE
                           .
                           +            : idy+1
                */
                if (idy > 0 && mask(idy-1, idx)==0 ) ystart++;
                if (idx > 0 && mask(idy, idx-1)==0 ) xstart++;
                if (idy < tny-1 && mask(idy+1, idx)==0 ) yend--;
                if (idx < tnx-1 && mask(idy, idx+1)==0 ) xend--;

                if (idy > 0 && mask(idy-1, idx)==0 ) ysize--;
                if (idx > 0 && mask(idy, idx-1)==0 ) xsize--;
                if (idy < tny-1 && mask(idy+1, idx)==0 ) ysize--;
                if (idx < tnx-1 && mask(idy, idx+1)==0 ) xsize--;
            }

            assert(xsize == xend - xstart);
            assert(ysize == yend - ystart);

    //        if (by) {
    //            if (idx>0 && ( mask(idy,idx-1)==0 || mask(idy+1,idx-1)==0 ) ) xstart++;    // shift the x-start point
    //            if (idx==tnx-1-bx || mask(idy,idx+1+bx)==0 || mask(idy+1,idx+1+bx)==0 ) xend--;    // shift the x-end point
    //        } else {
    //            if (idx>0 && mask(idy, idx-1)==0) xstart++;             // shift the x-start point
    //            if (idx==tnx-1-bx || mask(idy, idx+1+bx)==0) xend--;    // shift the x-end point
    //        }
    //
    //        if (bx) {
    //            if (idy>0 && ( mask(idy-1,idx)==0 || mask(idy-1,idx+1)==0 ) ) ystart++;    // shift the y-start point
    //            if (idy==tny-1-by || mask(idy+1+by,idx)==0 || mask(idy+1+by,idx+1)==0 ) yend--;    // shift the y-end point
    //        } else {
    //            if (idy>0 && mask(idy-1, idx)==0) ystart++;             // shift the y-start point
    //            if (idy==tny-1-by || mask(idy+1+by, idx)==0) yend--;    // shift the y-end point
    //        }

    //        if (bx && by) {
    //            if (idx>0 && ( mask(idy,idx-1)==0 || mask(idy+1,idx-1)==0 ) ) xstart++;    // shift the x-start point
    //            if (idy>0 && ( mask(idy-1,idx)==0 || mask(idy-1,idx+1)==0 ) ) ystart++;    // shift the y-start point
    //            if (idx==tnx-2 || mask(idy,idx+2)==0 || mask(idy+1,idx+2)==0 ) xend--;    // shift the x-end point
    //            if (idy==tny-2 || mask(idy+2,idx)==0 || mask(idy+2,idx+1)==0 ) yend--;    // shift the y-end point
    //        } else {
    //            if (idx>0 && mask(idy, idx-1)==0) xstart++;             // shift the x-start point
    //            if (idy>0 && mask(idy-1, idx)==0) ystart++;             // shift the y-start point
    //            if (idx==tnx-1-bx || mask(idy, idx+1+bx)==0) xend--;    // shift the x-end point
    //            if (idy==tny-1-by || mask(idy+1+by, idx)==0) yend--;    // shift the y-end point
    //        }

            for (blas_int i=ystart; i<yend; ++ i) {
                // get column
                if (Ty.columns(i)==y) { // diagonal
                    for (col=row_index[y]; cols[col]<Tx.columns(xstart); ++col) {} // TODO: get rid of this (too slow)
                    assert(col<row_index[y+1]);

                    for (blas_int j=xstart; j<xend; ++j) {
                        body_.nonzeros(pos) = Tx.nonzeros(j);
                        if (Tx.columns(j)==x)
                            body_.nonzeros(pos) += Ty.nonzeros(i);
                        body_.columns(pos++) = col++;
                    }
                } else {
                    for (col=row_index[Ty.columns(i)]; cols[col]<x; ++col) {} // TODO: get rid of this (too slow)
                    assert(col<row_index[Ty.columns(i)+1]);

                    body_.nonzeros(pos) = Ty.nonzeros(i);
                    body_.columns(pos++) = col;
                }
            }

            body_.row_index(b+1) = pos;
        }

    //    int l = 32;
    //    int s = 16;
    //
    //    for (int i=0; i<10; ++i) {
    //        int q = body_.row_index(row_index[i] + s);
    //        for (int j=q; j<q+l; ++j) {
    //            cout << body_.columns(j) << " ";
    //        }
    //        cout << endl;
    //        for (int j=q; j<q+l; ++j) {
    //            cout << ref.columns(j) << " ";
    //        }
    //        cout << endl;
    //        cout << endl;
    //    }

        return *this;
    }
    void MaskOperator2d::LU_factorize()
    {
        assert( init_ );
      //
        body_.LU_factorize();
    }
    // Operators
    MaskOperator2d& MaskOperator2d::operator*= (const dcomp& alpha)
    {
        assert(init_);
      //
        body_ *= alpha;
        return *this;
    }
    MaskOperator2d& MaskOperator2d::operator+= (const dcomp& alpha)
    {
        assert(init_);
      //
        body_.add_to_diagonal( alpha );
        return *this;
    }
    MaskOperator2d& MaskOperator2d::operator+= (const MaskVector2d& rhs)
    {
        assert(init_);
        assert( grid_ == rhs.get_grid() );
      //
        body_.add_vector_to_diagonal(rhs.function_values());

        return *this;
    }

    MaskVector2d MaskOperator2d::operator* (const MaskVector2d& rhs) const
    {
        assert(init_);
        assert(grid_ == rhs.get_grid() );
      //
        MaskVector2d out(grid_);
        this->gemv(1.0, rhs, 0.0, out);
        return out;
    }
    void MaskOperator2d::gemv(const dcomp& alpha, const MaskVector2d& x, const dcomp& beta, MaskVector2d& y) const
    {
        assert(init_);
        assert( grid_ == x.get_grid() );
        assert( grid_ == y.get_grid() );
      //
        body_.gemv(alpha, x.body(), beta, y.body());
    }
    void MaskOperator2d::LU_back_substitution(MaskVector2d& x)
    {
        assert( init_ );
        assert( x.init() );
        assert( grid_ == x.get_grid() );
      //
        body_.LU_back_substitution(x.body());
    }

}   // namespace QSCAT
