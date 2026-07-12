#include <iostream>
#include <fstream>
#include <stdio.h>
#include <complex>
#include <cassert>
#include <string>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "bessel.h"
#include "common.h"
#include "blas.h"
#include "Arrays.h"
#include "input.h"
#include "fem_dvr_ecs.h"

namespace QSCAT
{

    template<>
    void poly_grid_vector<double,std::complex<double> >::Save(const char * filename)
    {
        std::complex<double> val;
        FILE * file;
        fopen_s(&file,filename,"w");
        fprintf(file,"#%d-Poly-vector of %d values\n", n, base);
        fprintf(file,"#Coordinate");
        for (int i=0;i<n;++i){
            fprintf(file,"#\tReal Part\tImaginary Part");
        }
        fprintf(file,"#\n");
        for (int j=0;j<base;++j){
            fprintf(file,"%.12E", grid->Xr(j));
            for (int i=0;i<n;++i){
                val = F(i,j);
                fprintf(file, "\t%.12E\t%.12E", real(val), imag(val) );
            }
            fprintf(file,"#\n");
        }
        fclose(file);
    }
    template<>
    void poly_grid_vector<double,double>::Save(const char * filename)
    {
        double val;
        FILE * file;
        fopen_s(&file,filename,"w");
        fprintf(file,"#%d-Poly-vector of %d values\n", n, base);
        fprintf(file,"#Coordinate");
        for (int i=0;i<n;++i){
            fprintf(file,"#\tReal Part");
        }
        fprintf(file,"#\n");
        for (int j=0;j<base;++j){
            fprintf(file,"%.12E", grid->Xr(j));
            for (int i=0;i<n;++i){
                val = F(j,i);
                fprintf(file, "\t%.12E", val);
            }
            fprintf(file,"#\n");
        }
        fclose(file);
    }
    template<>
    void poly_grid_vector<def_float, def_comp>::SaveTransposed(const Vector<def_comp>& X, const char* filename)
    {
        std::complex<double> val;
        FILE * file;
        fopen_s(&file,filename,"w");
        fprintf(file,"#%d values of transposed %d-Poly-vector\n", base, n);
        fprintf(file,"#Coordinate");
        for (int i=0;i<n;++i){
            fprintf(file,"#\tReal Part\tImaginary Part");
        }
        fprintf(file,"#\n");
        for (int j=0;j<n;++j){
            fprintf(file,"%.12E\t%.12E", real(X[j]), imag(X[j]));
            for (int i=0;i<base;++i){
                val = F(j,i);
                fprintf(file, "\t%.12E\t%.12E", real(val), imag(val) );
            }
            fprintf(file,"#\n");
        }
        fclose(file);
    }
    template<>
    void poly_grid_vector<def_float, def_float>::SaveTransposed(const Vector<def_float>& X, const char* filename)
    {
        double val;
        FILE * file;
        fopen_s(&file,filename,"w");
        fprintf(file,"#%d values of transposed %d-Poly-vector\n", base, n);
        fprintf(file,"#Coordinate");
        for (int i=0;i<n;++i){
            fprintf(file,"#\tReal Part\tImaginary Part");
        }
        fprintf(file,"#\n");
        for (int j=0;j<n;++j){
            fprintf(file,"%.12E", X[j]);
            for (int i=0;i<base;++i){
                val = F(j,i);
                fprintf(file, "\t%.12E", val);
            }
            fprintf(file,"#\n");
        }
        fclose(file);
    }
} // namespace QSCAT
