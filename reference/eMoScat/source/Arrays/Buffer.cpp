#include <cassert>
#include <string>
#include <fstream>
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

#include "Arrays/Buffer.h"

// TODO clean unnecessary includes

namespace QSCAT
{
    template<>
    void Buffer<def_float>::save (const char* filename) const 
    {
        FILE * file;
        fopen_s(&file,filename,"w");
        fprintf(file,"#Vector of %lld values\n", num_values_);
        fprintf(file,"#Index\tValue \n");
        for (blas_int i=0; i<num_values_; ++i) {
            fprintf(file, "%lld\t%.12E\n", i, begin_[i]);
        }
        fclose(file);
    }
    template<>
    void Buffer<def_comp>::save (const char* filename) const 
    {
        FILE * file;
        fopen_s(&file,filename,"w");
        fprintf(file,"#Vector of %lld values\n", num_values_);
        fprintf(file,"#Index\tReal part of z     \tImaginary part of z\n");
        for (blas_int i=0; i<num_values_; ++i) {
            fprintf(file, "%lld\t%.12E\t%.12E\n", i, real(begin_[i]), imag(begin_[i]));
        }
        fclose(file);
    }

    template<>
    void Buffer<def_float>::save_range (def_float x, def_float y, const char* filename) const 
    {
        FILE * file;
        def_float r = (y - x)/(num_values_-1);
        fopen_s(&file,filename,"w");
        fprintf(file,"#Vector of %lld values\n", num_values_);
        fprintf(file,"#Index\tValue \n");
        for (blas_int i=0; i<num_values_; ++i){
            fprintf(file, "%.12E\t%.12E\n", x+i*r, *(begin_+i));
        }
        fclose(file);
    }
    template<>
    void Buffer<def_comp>::save_range (def_float x, def_float y, const char* filename) const 
    {
        FILE * file;
        def_float r = (y-x)/(num_values_-1);
        fopen_s(&file,filename,"w");
        fprintf(file,"#Vector of %lld values\n", num_values_);
        fprintf(file,"#Index\tReal part of z \tImaginary part of z\n");
        for (blas_int i=0; i<num_values_; ++i) {
            fprintf(file, "%.12E\t%.12E\t%.12E\n", x+i*r, real(begin_[i]), imag(begin_[i]));
        }
        fclose(file);
    }
}; // namespace QSCAT
