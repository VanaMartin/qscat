#include <stdio.h>
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
//#include "solver.h"           // Customized version of mkl_solver.h

#include "Arrays/RowCompressedMatrix.h"

namespace QSCAT
{
    template<>
    void RowCompressedMatrix<dComplex>::save(const char * name) const
    {
        FILE * file;
        fopen_s(&file,name,"w");
        fprintf(file,"#Row Compressed matrix of %lld x %lld size.\n", num_rows_, num_columns_);
        fprintf(file,"#Total number of %lld nonzero elements. \n", num_nonzeros_);
        fprintf(file,"#Index\tValue\n");
        for (blas_int i=0;i<num_nonzeros_;++i){
            fprintf(file, "%.12E\t%.12E\t%lld\n", real(nonzeros_[i]), imag(nonzeros_[i]), columns_[i]);
        }
        fprintf(file, "\n");
        for (blas_int i=0;i<num_rows_+1;++i){
            fprintf(file, "%lld\n", row_index_[i]);
        }
        fclose(file);
        return;
    }
} // namespace QSCAT
