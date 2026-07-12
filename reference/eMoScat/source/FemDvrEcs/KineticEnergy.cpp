#include <complex>
#include "common.h"
#include "Arrays.h"
#include "input.h"
#include "bessel.h"
#include "coulomb.h"

#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"
#include "FemDvrEcs/KineticEnergy.h"

namespace QSCAT
{
zMatrix generateKineticTerm(const FemDvrEcsGrid& grid, dfloat mu)
{
    // TODO remove commented body assignement lines : the competence to choose matrix ordering is now fully on matrix class

    zMatrix body(grid.nb(), grid.nb());

    dfloat coef;                     // Coefficient for auxiliary use
    blas_int nq, nel, nb;            // shortcuts for clarity

    nq  = grid.quadrature();    // number of quadrature points
    nel = grid.tnel();          // number of elements
    nb  = grid.nb();            // number of basis functions

    blas_int offset; //, rowoffset, coloffset;   // Auxiliary integers
    blas_int start, end;
    blas_int i, k1, j, ii, jj, l;

    dcomp hz, hlp, old_corner;      // Auxiliary complex numbers
    zVector wze(nq);          // scaled weights in one element (!= wz at the endpoints!!!)
    zMatrix dBF(nq, nq);      // derivatives of basis functions in one element

    coef = 0.5/mu;

    for (i=1;i<nb*nb;++i){
        body[i] = 0.0;          // zeros everywhere
    }

    offset = 0;             // position of the first point of the element in KE%A
    old_corner = zzero;     // contribution of the previous element to the overlapping point

    for (blas_int k=0; k<nel; ++k) {
        start = (k==0)? 1 : 0;                  // in the first element skip the first point (boundary condition psi(0) = 0)
        end = (k==nel-1)? nq-1:nq;              // in the last element skip the last point (boundary condition psi(x_max) = 0)
        hz = def_float(0.5) * grid.aaz(k+1);    // a half of the element interval
        for (i=0; i<nq; ++i) {
            wze[i] = hz * grid.dvr().w(i);
        }
        for (i=0; i<nq*nq; ++i){
            dBF[i] = grid.dlp()[i] / hz;
        }
        k1 = k * (nq - 1) - 1;            // index of the first point of the k-th element in xz(:)
        if (k1 < 0) {k1 = 0;}               // correction for the first element
        for (i=start; i<end; ++i) {
            hlp = def_float(1.0) / sqrt(grid.w(k1 + i - start));
            for (j=0; j<nq; ++j){
                dBF(i,j) = dBF(i,j) * hlp;      // normalization factor of basis functions
            }
        }
        k1 = (k) * (nq - 1) - 1;            // index of the first point of the k-th element in xz(:)
        if (k1 < 0) {k1 = 0;}               // correction for the first element
        for (i=start; i<end; ++i){
            ii = offset + i - start;        // current row in the matrix
            for (j=start; j<end; ++j){
                jj = offset + j - start;    // current column in the matrix
                hlp = 0.0;
                for (l=0; l<nq; ++l) {
                    hlp += wze[l] * dBF(i,l) * dBF(j,l);
                    // add the conjg for logical reasons?
                }
                //body[jj*nb + ii] = coef * hlp;
                //body[ii*nb + jj] = body[jj*nb + ii];
                body(ii, jj) = coef * hlp;
                body(jj, ii) = body(ii, jj);
            }
        }
        //body[offset*nb + offset] = body[offset*nb + offset] + old_corner;
        body(offset, offset) = body(offset, offset) + old_corner;
        offset = offset + end - start - 1;
        // old_corner = body[offset*nb + offset];
        old_corner = body(offset, offset);
    }
    return body;
}

RowCompressedMatrix<dcomp> generateKineticTermRCM(const FemDvrEcsGrid& grid, dfloat mu)    // Generate One dimensional kinetic term
{
// code
    blas_int basis_size = grid.nb();                             // total basis functions (points on grid)
    blas_int quadrature = grid.quadrature();                     // quadrature order (functions on each element)
    blas_int total_elements = grid.tnel();
    blas_int num_nonzeros = quadrature*quadrature*total_elements - 4*quadrature + 2 - total_elements + 1;
    RowCompressedMatrix<dcomp> out(basis_size, basis_size, num_nonzeros);

    dfloat coef;                                                 // Coefficient for auxiliary use
    blas_int offset, rowoffset, coloffset;                       // Auxiliary integers
    blas_int start, end;
    blas_int i, k1, j, nep, pos, l;

    coef = 0.5/mu;

    dcomp hz, hlp, old_corner;                                  // Auxiliary complex numbers
    zVector wze(quadrature);                              // scaled weights in one element (!= wz at the endpoints!!!)
    zMatrix dBF(quadrature, quadrature);                  // derivatives of basis functions in one element

    out.row_index(0)=0;

    offset = 0;                                             // position of the first point in the element in KE_SLU(:)
    coloffset = 0;                                          // (row index of the first point in the element) - 1
    rowoffset = 1;                                          // index in the array KE%colIndex (i.e. column index + 1)
    old_corner = zzero;                                     // contribution of the previous element to the overlapping point


    for (blas_int k=0; k<total_elements; ++k){
        start = (k==0)? 1 : 0;                                  // in the first element skip the first point (boundary condition psi(0) = 0)
        end = (k==total_elements-1)? quadrature-1:quadrature;   // in the last element skip the last point (boundary condition psi(x_max) = 0)
        hz = dfloat(0.5) * (grid.aaz(k+1));                          // a half of the element interval
        for (i=0; i<quadrature; ++i) {
            wze[i] = hz * (grid.dvr().w(i));
        }
        for (i=0; i<quadrature*quadrature; ++i){
            dBF[i] = grid.dlp()[i] / hz;
        }
        k1 = (k) * (quadrature - 1) - 1;                    // index of the first point of the k-th element in xz(:)
        if (k1 < 0) {k1 = 0;}                               // correction for the first element
        for (i=start; i<end; ++i) {
            hlp = dfloat(1) / sqrt(grid.w(k1 + i - start));
            for (j=0; j<quadrature; ++j){
                dBF(i,j) = dBF(i,j) * hlp;                  // normalization factor of basis functions
            }
        }
        nep = end - start;                                  // number of points in the k-th element
        pos = offset;                                       // position in nze(.)
        for (i=start; i<end; ++i){
            for (j=start; j<end; ++j){
                if (j>=i){
                    hlp = 0.0;
                    for (l=0; l<quadrature; ++l) {
                        hlp += coef * (wze[l]*dBF(i,l)*dBF(j,l));
                    }
                    out.nonzeros(pos) = hlp;                // Assigning value
                } else {
                    out.nonzeros(pos) = out.nonzeros(offset + (j - start) * nep + i - start);
                }
                out.columns(pos) = coloffset + j - start;
                ++pos;                                      // Shifting position to next one
            }
            out.row_index(rowoffset) = pos;
            ++rowoffset;
        }
        --rowoffset;                                                // Returning back one row
        out.nonzeros(offset) = out.nonzeros(offset) + old_corner;   // overlapping point
        offset = pos - 1;                                           // shift to next element
        coloffset = coloffset + nep - 1;                            // setting the column offset
        old_corner = out.nonzeros(offset);                          // Storing the overlaping value for next element
    }

    return out;
}
} // namespace QSCAT
