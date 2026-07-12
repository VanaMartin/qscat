#ifndef _MKL_H_
#define _MKL_H_

#if !defined( __MKL_SOLVER_H )

#define __MKL_SOLVER_H

#include "mkl_dss.h"
//#include "mkl_pardiso.h"
#include "mkl_rci.h"

#endif

#include "mkl_vml.h"
#include "mkl_vsl.h"
#include "mkl_service.h"
#include "mkl_spblas.h"

#endif /* _MKL_H_ */

typedef std::complex<double> doublecomplex;

#if !defined( __MKL_DSS_PARDISO_H )

#define __MKL_DSS_PARDISO_H

#include "mkl_dss.h"
#include "mkl_types.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define SOLVER_HANDLE_DIMENSION                      64
#define SOLVER_HANDLE_INIT_LOCATION                   0
#define SOLVER_HANDLE_S_BASIC_LOCATION                1
#define SOLVER_HANDLE_BACKLINK_LOCATION               2
#define SOLVER_HANDLE_FACT_ADR_LOCATION               3
#define SOLVER_HANDLE_STRUC_FI_LOCATION               7

#define AUXINFO_DIMENSION                            64
#define S_BASIC_COMI_MATRIX_TYPE_LOCATION            11
#define S_BASIC_COMI_TMPSIZ_LOCATION                 36
#define S_BASIC_COMI_NSUP_LOCATION                   30
#define S_BASIC_COMI_SYM_LOCATION                    18
#define S_BASIC_COMI_GSSUB_LOCATION                  35
#define S_BASIC_COMR_TIME_ADJ_LOCATION               0
#define S_BASIC_COMR_TIME_REORD_LOCATION             1
#define S_BASIC_COMR_TIME_SYMFCT_LOCATION            2
#define S_BASIC_COMR_TIME_SCAT_A_LOCATION            4
#define S_BASIC_COMR_TIME_NUMFCT_LOCATION            5
#define S_BASIC_COMR_TIME_SOLVE_LOCATION             6
#define S_BASIC_COMR_NUMFLOP_LOCATION                14
#define S_BASIC_DISP_XLNZ_LOCATION                   63
#define S_BASIC_DISP_XSUP_LOCATION                   60
#define S_BASIC_DISP_PIVOT_LOCATION                  67
#define S_BASIC_DISP_XUNZ_LOCATION                   64
#define S_BASIC_DISP_XLINDX_LOCATION                 61
#define S_BASIC_DISP_LINDX_LOCATION                  62
#define S_BASIC_DISP_FIN_INT_LOCATION                73
#define S_BASIC_DISP_FIN_NMOD_LOCATION               65
#define S_BASIC_DISP_FIN_SNODE_LOCATION              68

#define L_COMI  100
#define L_COMR  20
#define L_DISP  100
typedef struct
{
	_LONG_t comi[L_COMI];
	double comr[L_COMR];
	_LONG_t disp[L_DISP];
} s_basic;

typedef struct _sparseDataStruct {
        MKL_INT          facilityHandle;
        _INTEGER_t   nRows;
        _INTEGER_t   nCols;
        _INTEGER_t   nNonZeros;
        _INTEGER_t   origNNonZeros;
        _INTEGER_t * majorIndex;
        _INTEGER_t * minorIndex;
        _INTEGER_t * permVector;
        _INTEGER_t * origMajorIndex;
        _INTEGER_t * origMinorIndex;
        _INTEGER_t * fill;

        MKL_INT          matrixType;
        MKL_INT          matrixStructure;
        MKL_INT          valueStructure;
        MKL_INT          dataType;
        void       * matrixValues;
        void       * origMatrixValues;
        _INTEGER_t   nRhs;
        MKL_INT          maxFacStore;
        MKL_INT          matrixNumber;
        void       * rhsValues;
        void       * computedSolution;
        void       * solverHandle[SOLVER_HANDLE_DIMENSION];
        MKL_INT          nCpus;
        _INTEGER_t   phase;
        _INTEGER_t   auxInfoArray[AUXINFO_DIMENSION];
        MKL_INT          solverMessageLevel;
        MKL_INT          messageLevel;
        MKL_INT          terminationLevel;
        MKL_INT          option1;		/* for debug only */
        MKL_INT          solverState;
        MKL_INT          detSign;
        double       detPow;
        double       detRealBase;
        double       detImgBase;
        MKL_INT          nPosEig;
        MKL_INT          nNegEig;
        MKL_INT          nZeroEig;
} typeDssHandle;


#if defined( _WIN32 ) || defined( _WIN64 )
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif

/* PARDISO prototype. */

extern "C" void pardisoinit (void   *, MKL_INT *,   MKL_INT *);
extern "C" void pardiso_64     (void   *, long long int *,   long long int *, long long int *,    long long int *, long long int *,
                  doublecomplex *, long long int *,   long long int *, long long int *,   long long int *, long long int *,
                  long long int *, doublecomplex *, doublecomplex *, long long int *);
extern "C" void pardiso     (void   *, int *,   int *, int *,    int *, int *,
                  doublecomplex *, int *,   int *, int *,   int *, int *,
                  int *, doublecomplex *, doublecomplex *, int *);
}
#endif /* __cplusplus */
