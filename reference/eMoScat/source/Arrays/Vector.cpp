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

#include "Arrays/Vector.h"

namespace QSCAT
{
    template<>
    void Vector<blas_int>::save(const char*name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Vector of %lld values\n", size_);
        fprintf(file,"#Index\tvalue\n");
        for (blas_int i=0;i<size_;++i){
            fprintf(file, "%lld\t%lld\n", i, array_[i]);
        }
        fclose(file);
        return;
    }
    template<>
    void Vector<def_comp>::save(const char*name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Vector of %lld values\n", size_);
        fprintf(file,"#Index\tReal part of z     \tImaginary part of z\n");
        for (blas_int i=0;i<size_;++i){
            fprintf(file, "%lld\t%.12E\t%.12E\n", i, real(array_[i]), imag(array_[i]));
        }
        fclose(file);
        return;
    }
    template<>
    void Vector<def_float>::save(const char*name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Vector of %lld values\n", size_);
        fprintf(file,"#Index\tReal part of z     \tImaginary part of z\n");
        for (blas_int i=0;i<size_;++i){
            fprintf(file, "%lld\t%.12E\n", i, array_[i]);
        }
        fclose(file);
        return;
    }
    template<>
    void Vector<def_comp>::save(Vector<double> & X, const char * name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Vector of %lld values\n", size_);
        fprintf(file,"#Index\tReal part\tImaginary part\n");
        for (blas_int i=0;i<size_;++i){
            fprintf(file, "%.12E\t%.12E\t%.12E\n", X[i], real(array_[i]), imag(array_[i]));
        }
        fclose(file);
        return;
    }
    template<>
    void Vector<def_float>::save(Vector<double> & X, const char * name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Vector of %lld values\n", size_);
        fprintf(file,"#Index\tReal part\tImaginary part\n");
        for (blas_int i=0;i<size_;++i){
            fprintf(file, "%.12E\t%.12E\n", X[i], array_[i]);
        }
        fclose(file);
        return;
    }

    template<>
    void Vector<def_comp>::save(Vector<dComplex> & X, const char * name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Vector of %lld values\n", size_);
        fprintf(file,"#Index\tReal part\tImaginary part\n");
        for (blas_int i=0;i<size_;++i){
            fprintf(file, "%.12E\t%.12E\t%.12E\n", abs(X[i]), real(array_[i]), imag(array_[i]));
        }
        fclose(file);
        return;
    }
    template<>
    void Vector<def_float>::save(Vector<dComplex> & X, const char * name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Vector of %lld values\n", size_);
        fprintf(file,"#Index\tReal part\tImaginary part\n");
        for (blas_int i=0;i<size_;++i){
            fprintf(file, "%.12E\t%.12E\n", abs(X[i]), array_[i]);
        }
        fclose(file);
        return;
    }

    template<>
    void Vector<float>::print() const
    {
        assert(init());
      //
        printf("[ ");
        for (int i=0; i<size_ - 1; ++i)
            printf("%f, ", array_[i]);
        printf("%f ]\n", array_[size_-1]);
    }

    template<>
    void Vector<def_float>::print() const
    {
        assert(init());
      //
        printf("[ ");
        for (int i=0; i<size_ - 1; ++i)
            printf("%f, ", array_[i]);
        printf("%f ]\n", array_[size_-1]);
    }

    template<>
    void Vector<def_comp>::print() const
    {
        assert(init());
      //
        printf("[ ");
        for (int i=0; i<size_ - 1; ++i)
            printf("(%f,%f), ", real(array_[i]), imag(array_[i]));
        printf("(%f,%f) ]\n", real(array_[size_-1]), imag(array_[size_-1]));
    }

    template<>
    void SaveMultipleVectors(int N, int M, Vector<dComplex> ** X, const char * name)
    {
    // N stands for number of Vectors, M denotes the number of values
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Multiple Vectors stored in one file.");
        fprintf(file,"#%d Vectors\n", N);
        fprintf(file,"#%d values in each\n", M);
        for (blas_int i=0; i<M; ++i){
            for (blas_int j=0; j<N; ++j){
                if (j < N - 1) {
                // Not last Vector
                    fprintf(file, "%.12E\t%.12E\t", real((*(X[j]))[i]), imag((*(X[j]))[i]));
                } else {
                // last Vector
                    fprintf(file, "%.12E\t%.12E", real((*(X[j]))[i]), imag((*(X[j]))[i]));
                }
            }
            fprintf(file,"\n");
        }
        fclose(file);
        return;
    }
    void SaveMultipleVectorsAbs2(const int & N, const int & M, Vector<dComplex> ** X, const char * name)
    {
    // N stands for number of Vectors, M denotes the number of values
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Multiple Vectors stored in one file.");
        fprintf(file,"#%d Vectors\n", N);
        fprintf(file,"#%d values in each\n", M);
        for (blas_int i=0; i<M; ++i){
            for (blas_int j=0; j<N; ++j){
                if (j < N - 1) {
                // Not last Vector
                    fprintf(file, "%.12E\t", std::pow(std::abs((*(X[j]))[i]),2));
                } else {
                // last Vector
                    fprintf(file, "%.12E", std::pow(std::abs((*(X[j]))[i]),2));
                }
            }
            fprintf(file,"\n");
        }
        fclose(file);
        return;
    }
} // namespace QSCAT
