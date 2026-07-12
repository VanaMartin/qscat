#include <fstream>
#include <string>
#include <cassert>
#include <complex>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "common.h"

#include "mkl_types.h"
#include "mkl_blas.h"
#include "mkl_dfti.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"

#include "Arrays/Matrix.h"

namespace QSCAT
{
    template<>
    void Matrix<def_float>::save(const char *name) const
    {
        assert(init());
        assert(rows_);
        assert(columns_);
        assert(array_);
      //
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Matrix of %lld rows and %lld columns\n", rows_, columns_);
        for (blas_int i=0; i<rows_; ++i){
            for (blas_int j=0; j<columns_; ++j) {
                fprintf(file, "%.12E", (*this)(i,j));
                fprintf(file, (j!=columns_-1)? "\t" : "\n");
            }
        }
        fclose(file);
        return;
    }

    template<>
    void Matrix<def_comp>::save(const char *name) const
    {
        assert(init());
        assert(rows_);
        assert(columns_);
        assert(array_);
      //
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Complex Matrix of %lld rows and %lld columns\n", rows_, columns_);
        def_comp val;
        for (blas_int i=0; i<rows_; ++i){
            for (blas_int j=0; j<columns_; ++j) {
                val = (*this)(i,j);
                fprintf(file, "%.12E\t%.12E", real(val), imag(val));
                fprintf(file, (j!=columns_-1)? "\t" : "\n");
            }
        }
        fclose(file);
        return;
    }

    template<>
    void Matrix<def_float>::save(const Vector<def_float>& range, const char *name) const
    {
        assert(init());
        assert(rows_);
        assert(columns_);
        assert(array_);
        assert(range.get_size() >= rows_);
      //
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Matrix of %lld rows and %lld columns, with xrange in first column\n", rows_, columns_+1);
        for (blas_int i=0; i<rows_; ++i){
            fprintf(file, "%.12E\t", range[i]);
            for (blas_int j=0; j<columns_; ++j) {
                fprintf(file, "%.12E", (*this)(i,j));
                fprintf(file, (j!=columns_-1)? "\t" : "\n");
            }
        }
        fclose(file);
        return;
    }

    template<>
    void Matrix<def_comp>::save(const Vector<def_float>& range, const char *name) const
    {
        assert(init());
        assert(rows_);
        assert(columns_);
        assert(array_);
      //
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Complex Matrix of %lld rows and %lld columns, with xrange in first column\n", rows_, columns_+1);
        def_comp val;
        for (blas_int i=0; i<rows_; ++i){
            fprintf(file, "%.12E\t", range[i]);
            for (blas_int j=0; j<columns_; ++j) {
                val = (*this)(i,j);
                fprintf(file, "%.12E\t%.12E", real(val), imag(val));
                fprintf(file, (j!=columns_-1)? "\t" : "\n");
            }
        }
        fclose(file);
        return;
    }
} // namespace QSCAT
