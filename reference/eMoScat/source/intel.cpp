#include <iostream>
#include <fstream>
#include <complex>
#include <string.h>

#include "common.h"

#include "intel.h"
#include <cassert>

using namespace std;

namespace QSCAT
{
    /*
        This  namespace contains the classes  and functions necessary to access
        the sparse solver  routines of Intel Math Kernel Library.  The contains
        is adapted  to the electron molecule  scattering problems, but could be
        possibly used for other codes.
    */
    char CSRMV[6] = {'G', 'U', 'N', 'C'};       // The constant Matrix descriptor for sparse matrix vector product

    /*
        Definitions of FORTRAN-like variables for Intel MKL, added for security reasons
    */
    #ifndef MKL_INT
    #define MKL_INT int
    #endif

    #ifndef lapack_int
    #define lapack_int MKL_INT
    #endif

    #ifndef lapack_logical
    #define lapack_logical lapack_int
    #endif

    #ifndef lapack_complex_float
    #define lapack_complex_float   MKL_Complex8
    #endif

    #ifndef lapack_complex_double
    #define lapack_complex_double   MKL_Complex16
    #endif

    //const lapack_int LAPACK_COL_MAJOR = 102;


    /*
        BLAS  functions:  General  set  of Parallelized  Basic  Linear Algebra
        functions provided by Intel MKL Libraries. The further generalizations
        are expected. Namely implementation of PBLAS methods for multi-machine
        computations.
    */
// Array copy method

    void blas::copy(blas_int N,const blas_int * X, blas_int * Y)
    {
        assert(X);
        assert(Y);
        assert(N>0);
    //
        for (blas_int i=0;i<N;++i){
            Y[i] = X[i];
        }
        return;
    }
    void blas::copy(blas_int N,const double * X, double * Y)
    {
        assert(X);
        assert(Y);
        assert(N>0);
    //
        MKL_INT n = N;
        cblas_dcopy (n, X, 1, Y, 1);
        return;
    }
    void blas::copy(blas_int N, const dComplex * X, dComplex * Y)
    {
        assert(X);
        assert(Y);
        assert(N>0);
    //
        MKL_INT n = N;
        cblas_zcopy (n, X, 1, Y, 1);
        return;
    }

    void blas::copy(blas_int N,const blas_int * X, blas_int incx, blas_int * Y, blas_int incy)
    {
        assert(X);
        assert(Y);
        assert(N>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        for (blas_int i=0;i<N;++i){
            Y[i*incy] = X[i*incx];
        }
        return;
    }

    void blas::copy(blas_int N, const double* X, blas_int incx, double* Y, blas_int incy)
    {
        assert(X);
        assert(Y);
        assert(N>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        MKL_INT n = N;
        cblas_dcopy (n, X, incx, Y, incy);
        return;
    }
    void blas::copy(blas_int N, const dComplex* X, blas_int incx, dComplex * Y, blas_int incy)
    {
        assert(X);
        assert(Y);
        assert(N>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        MKL_INT n = N;
        cblas_zcopy (n, X, incx, Y, incy);
        return;
    }

// TODO investigate strange invalid reads in intel library for larger arrays
// Vector operation sub(y) := sub(y) + a*sub(x)
    void blas::axpy (blas_int n, const blas_int a, const blas_int * x, blas_int * y)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
    //
        for (blas_int k=0; k<n; ++k) {
            y[k] += a * x[k];
        }
        return;
    }
    void blas::axpy (blas_int n, const dComplex a, const dComplex * x, dComplex * y)
    {
        assert(x);
        assert(y);
        assert(n>0);
    //
        cblas_zaxpy(n,&a,x,1,y,1);
        return;
    }
    void blas::axpy (blas_int n, const double a, const dComplex * x, dComplex * y)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
    //
        dComplex an = a;
        cblas_zaxpy(n,&an,x,1,y,1);
        return;
    }
    void blas::axpy (blas_int n, const double a, const double * x, double * y)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
    //
        cblas_daxpy(n,a,x,1,y,1);
        return;
    }
  //
    void blas::axpy (blas_int n, const blas_int a, const blas_int * x, blas_int incx, blas_int * y, blas_int incy)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        for (blas_int k=0; k<n; ++k) {
            y[k*incy] += a * x[k*incx];
        }
        return;
    }
    void blas::axpy (blas_int n, const dComplex a, const dComplex * x, blas_int incx, dComplex * y, blas_int incy)
    {
        assert(x);
        assert(y);
        assert(n>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        cblas_zaxpy(n,&a,x,incx,y,incy);
        return;
    }
    void blas::axpy (blas_int n, const double a, const dComplex * x, blas_int incx, dComplex * y, blas_int incy)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        dComplex an = a;
        cblas_zaxpy(n,&an,x,incx,y,incy);
        return;
    }
    void blas::axpy (blas_int n, const double a, const double * x, blas_int incx, double * y, blas_int incy)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        cblas_daxpy(n,a,x,incx,y,incy);
        return;
    }

// Vector dot product
    dComplex blas::dotproduct(blas_int n,const dComplex * x, const dComplex * y)
    {
        assert(n>0);
        assert(x);
        assert(y);
    //
        dComplex out;
        cblas_zdotc_sub(n,x,1,y,1,&out);
        return out;
    }
    double blas::dotproduct(blas_int n,const double * x, const double * y)
    {
        assert(n>0);
        assert(x);
        assert(y);
    //
        double out;
        out = cblas_ddot(n,x,1,y,1);
        return out;
    }
// Vector reduction
    dComplex blas::reduct(blas_int n,const dComplex * x, const dComplex * y)
    {
        assert(n>0);
        assert(x);
        assert(y);
    //
        dComplex out;
        cblas_zdotu_sub(n,x,1,y,1,&out);
        return out;

    }
    double blas::reduct(blas_int n,const double * x, const double * y)
    {
        return blas::dotproduct(n, x, y);
    }
// vector x vector - element wise multiplication
    void blas::ewxy(blas_int n, const double * x, const double * y, double * o)
    {
        assert(n>0);
        assert(x);
        assert(y);
        assert(o);
    //
        vdMul( n, x, y, o);
    }
    void blas::ewxy(blas_int n, const dComplex * x, const dComplex * y, dComplex * o)
    {
        assert(n>0);
        assert(x);
        assert(y);
        assert(o);
    //
        vzMul( n, (MKL_Complex16*) x, (MKL_Complex16*) y, o);
    }
    void blas::subewxy(blas_int n, blas_int i, blas_int j, double * x, double * y, double * o)
    {
        assert(n>0);
        assert(i>=0);
        assert(j>=0);
        assert(x);
        assert(y);
        assert(o);
    //
        vdMul( n,  &(x[i]), &(y[j]), o);
    }
    void blas::subewxy(blas_int n, blas_int i, blas_int j, dComplex * x, dComplex * y, dComplex * o)
    {
        assert(n>0);
        assert(i>=0);
        assert(j>=0);
        assert(x);
        assert(y);
        assert(o);
    //
        vzMul( n, (MKL_Complex16*) &(x[i]), (MKL_Complex16*) &(y[j]), o);
    }
// Partial dot product for projections
    dComplex blas::partial_dotproduct(blas_int n,const dComplex * x, blas_int incX, const dComplex * y, blas_int incY)
    {
        assert(n>0);
        assert(incX>=0);
        assert(incY>=0);
        assert(x);
        assert(y);
    //
        dComplex out;
        cblas_zdotc_sub(n,x,incX,y,incY,&out);
        return out;
    }
    double blas::partial_dotproduct(blas_int n,const double * x, blas_int incX, const double * y, blas_int incY)
    {
        assert(n>0);
        assert(incX>=0);
        assert(incY>=0);
        assert(x);
        assert(y);
    //
        double out;
        out = cblas_ddot(n,x,incX,y,incY);
        return out;
    }
// Vector simple scaling
    void blas::scale(blas_int N, dComplex * a, dComplex alpha)
    {
        assert(N>0);
        assert(a);
    ////
        const void * Alpha = &alpha;
        cblas_zscal(N, Alpha, a, 1);
    }
    void blas::scale(blas_int N, dComplex * a, double alpha)
    {
        assert(N>0);
        assert(a);
    ////
        cblas_zdscal(N, alpha, a, 1);
    }
    void blas::scale(blas_int N, double * a, double alpha)
    {
        assert(N>0);
        assert(a);
    ////
        cblas_dscal(N, alpha, a, 1);
    }
    void blas::scale(blas_int N, dComplex * a, blas_int inc, dComplex alpha)
    {
        assert(N>0);
        assert(a);
        assert(inc>=0);
    ////
        cblas_zscal(N, &alpha, a, inc);
    }
    void blas::scale(blas_int N, dComplex * a, blas_int inc, double alpha)
    {
        assert(N>0);
        assert(a);
        assert(inc>=0);
    ////
        cblas_zdscal(N, alpha, a, inc);
    }
    void blas::scale(blas_int N, double * a, blas_int inc, double alpha)
    {
        assert(N>0);
        assert(a);
        assert(inc>=0);
    ////
        cblas_dscal(N, alpha, a, inc);
    }

    void blas::conj(blas_int N, double *a)
    {
        // real : nothing to do here, but operation is legal
    }
    void blas::conj(blas_int N, dComplex *a)
    {
        // first blas_interpret the poblas_inter as double
        double *x = (double*) a;
        // shift the poblas_inter by one (imaginary part)
        x++;
        // finally scale every second value (imaginary part) with -1.0
        blas::scale(N, x, 2, -1.0);
    }

// Matrix vector multiplication
    void blas::matrix_vector(const char & trans, blas_int m, blas_int n, double alpha, double beta, const double * a, const double * x, double * y)
    {
        assert(trans=='N'||trans=='T'||trans=='C');
        assert(m>0);
        assert(n>0);
        assert(a);
        assert(x);
        assert(y);
    // code
        CBLAS_TRANSPOSE C;
        switch (trans){
            case 'N': C = CblasNoTrans;
                break;
            case 'T': C = CblasTrans;
                break;
            case 'C': C = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                C = CblasNoTrans;
        }

        double Alpha = alpha;
        double Beta = beta;
    // ColumnMajor version for matrix definition conventions
        cblas_dgemv(CblasColMajor,C,m,n,Alpha,a,m,x,1,Beta,y,1);
        return;
    }
    void blas::matrix_vector(const char & trans, blas_int m, blas_int n, dComplex alpha, dComplex beta, const dComplex * a, const dComplex * x, dComplex * y)
    {
        assert(trans=='N'||trans=='T'||trans=='C');
        assert(m>0);
        assert(n>0);
        assert(a);
        assert(x);
        assert(y);
    // code
        CBLAS_TRANSPOSE C;
        switch (trans){
            case 'N': C = CblasNoTrans;
                break;
            case 'T': C = CblasTrans;
                break;
            case 'C': C = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                C = CblasNoTrans;
        }

        dComplex Alpha = alpha;
        dComplex Beta = beta;

        cblas_zgemv(CblasColMajor,C,m,n,&Alpha,a,m,x,1,&Beta,y,1);
        return;
    }

    // Sub matrix case
    void blas::sub_matrix_vector(const char & trans, blas_int m, blas_int n, blas_int lda, double alpha, double beta, double * a, double * x, double * y)
    {
        assert(trans=='N'||trans=='T'||trans=='C');
        assert(m>0);
        assert(n>0);
        assert(lda>0);
        assert(a);
        assert(x);
        assert(y);
    // code
        CBLAS_TRANSPOSE C;
        switch (trans){
            case 'N': C = CblasNoTrans;
                break;
            case 'T': C = CblasTrans;
                break;
            case 'C': C = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                C = CblasNoTrans;
        }

        double Alpha = alpha;
        double Beta = beta;
    // ColumnMajor version for matrix definition conventions
        cblas_dgemv(CblasColMajor,C,m,n,Alpha,a,lda,x,1,Beta,y,1);
        return;
    }
    void blas::sub_matrix_vector(const char & trans, blas_int m, blas_int n, blas_int lda, dComplex alpha, dComplex beta, dComplex * a, dComplex * x, dComplex * y)
    {
        assert(trans=='N'||trans=='T'||trans=='C');
        assert(m>0);
        assert(n>0);
        assert(lda>0);
        assert(a);
        assert(x);
        assert(y);
    // code
        CBLAS_TRANSPOSE C;
        switch (trans){
            case 'N': C = CblasNoTrans;
                break;
            case 'T': C = CblasTrans;
                break;
            case 'C': C = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                C = CblasNoTrans;
        }

        dComplex Alpha = alpha;
        dComplex Beta = beta;

        cblas_zgemv(CblasColMajor,C,m,n,&Alpha,a,lda,x,1,&Beta,y,1);
        return;
    }

// Sparse matrix vector multiplication
    void blas::RCmatrix_vector(char trans, blas_int m, double alpha, const double *a, const blas_int *ia, const blas_int *ja, const double *x, double beta, double *y)
    {
    // Assertions
        assert(m>0);
        assert(a);
        assert(ia);
        assert(ja);
        assert(x);
        assert(y);
    // Code starts here

        mkl_dcsrmv(&trans, &m, &m, &alpha, CSRMV, (double*) a, (blas_int*) ja, (blas_int*) &(ia[0]), (blas_int*) &(ia[1]), (double*) x, &beta, y);
        //mkl_dcsrmv(&trans, &m, &m, alpha, CSRMV, a, ja, &(ia[0]), &(ia[1]), x, beta, y);
    }
    void blas::RCmatrix_vector(char trans, blas_int m, dComplex alpha, const dComplex *a, const blas_int *ia, const blas_int *ja, const dComplex *x, dComplex beta, dComplex *y)
    {
    // Assertions
        assert(m>0);
        assert(a);
        assert(ia);
        assert(ja);
        assert(x);
        assert(y);
    // Code starts here

        mkl_zcsrmv(&trans, &m, &m, &alpha, CSRMV, (dComplex*) a, (blas_int*) ja, (blas_int*) &(ia[0]), (blas_int*) &(ia[1]), (dComplex*) x, &beta, y);
        //mkl_zcsrmv(&trans, &m, &m, alpha, CSRMV, a, ja, &(ia[0]), &(ia[1]), x, beta, y);
    }
// Full Matrix-Matrix operations
    void blas::matrix_matrix(const char & transa, const char & transb, blas_int m, blas_int n, const blas_int k, double alpha, double beta, const double*A, const double*B, double*C)
    {
    // Assertions
        assert(transa=='N'||transa=='T'||transa=='C');
        assert(transb=='N'||transb=='T'||transb=='C');
        assert(m>0);
        assert(n>0);
        assert(k>0);
        assert(A);
        assert(B);
        assert(C);
    // Code starts here
        CBLAS_TRANSPOSE TA, TB;
        switch (transa){
            case 'N': TA = CblasNoTrans;
                break;
            case 'T': TA = CblasTrans;
                break;
            case 'C': TA = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                TA = CblasNoTrans;
        }
        switch (transb){
            case 'N': TB = CblasNoTrans;
                break;
            case 'T': TB = CblasTrans;
                break;
            case 'C': TB = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                TB = CblasNoTrans;
        }
        cblas_dgemm(CblasColMajor,TA, TB, m,n,k,alpha,A,(transa=='N')? m:k,B,(transb=='N')? k:n, beta, C,m);
    }
    void blas::matrix_matrix(const char & transa, const char & transb, blas_int m, blas_int n, const blas_int k, dComplex alpha, dComplex beta, const dComplex*A, const dComplex*B, dComplex*C)
    {
    // Assertions
        assert(transa=='N'||transa=='T'||transa=='C');
        assert(transb=='N'||transb=='T'||transb=='C');
        assert(m>0);
        assert(n>0);
        assert(k>0);
        assert(A);
        assert(B);
        assert(C);
    // Code starts here
        CBLAS_TRANSPOSE TA, TB;
        switch (transa){
            case 'N': TA = CblasNoTrans;
                break;
            case 'T': TA = CblasTrans;
                break;
            case 'C': TA = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                TA = CblasNoTrans;
        }
        switch (transb){
            case 'N': TB = CblasNoTrans;
                break;
            case 'T': TB = CblasTrans;
                break;
            case 'C': TB = CblasConjTrans;
                break;
            default: std::cout << "No transposition identified. Setting to 'N'";
                TB = CblasNoTrans;
        }
        cblas_zgemm(CblasColMajor,TA, TB, m,n,k,&alpha,A,(transa=='N')? m:k,B,(transb=='N')? k:n, &beta, C,m);
    }
    // Eigenvalue and eigen vector functions
    //
    //    Complex<double> case: calling Intel ZGEEV procedure. See the MKL manual
    //    for detailed information.
    //
    void blas::eigen(dComplex * a, dComplex * W, blas_int n)
    {
    // Assertions
        assert(a);
        assert(W);
        assert(n>0);
    // Code starts here
        char jobvl;
        char jobvr;
        lapack_int N;
        lapack_int lda;
        lapack_int ldvl;
        lapack_int ldvr;
        //MKL_Complex16 * w = new MKL_Complex16[n];     // The array of eigen values
        dComplex * w = new dComplex[n];
        MKL_Complex16 * vr = new MKL_Complex16[n*n];    // Right hand eigenvector space pointer
        MKL_Complex16 * vl = new MKL_Complex16;         // Left hand eigenvector space pointer
        MKL_Complex16 * work = new MKL_Complex16[1];    // The function workspace
        MKL_INT lwork;
        double * rwork = new double[2*n];
        MKL_INT info;

        jobvl = 'N';                        // Left hand eigenvectors disabled
        jobvr = 'V';                        // Right hand eigenvectors enabled
        N = n;                              // The order of the matrix elements stored as a vector
        lda = n;                            // Array dimension
        ldvl = 1;                           // The leading dimension for left eigenvectors to be computed
        ldvr = n;                           // The leading dimension for right eigenvectors to be computed
        lwork = -1;

    // The first run returns the appropriate space of the workspace variable
        zgeev(&jobvl, &jobvr, &N, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

        lwork = blas_int(abs(work[0]));              // Proper workspace dimension
        delete[] work;                      // Cleaning old workspace
        work = new MKL_Complex16[lwork+1];  // Allocating new workspace

    // The second call invokes the solver to compute the eigenvalues
        zgeev(&jobvl, &jobvr, &N, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

    // Overwriting the matrix a with the eigenvalues ordered by real value of energy
        blas_int * order = new blas_int[n];
        blas::order_values(w, order, n);
        // TODO copy()
        for (blas_int i=0;i<n;++i){
            W[i] = w[i];
            for (blas_int j=0;j<n;++j){
                a[i*n + j] = vr[order[i]*n + j];
            }
        }

        delete[] work;
        delete[] rwork;
        delete[] w ;
        delete[] vr;
        delete vl;
        return;
    }

// Lapack linear equations solver
    void blas::lapack_solve(blas_int n, blas_int nrhs, double * a, double * b)
    {
    // Assertions
        assert(n>0);
        assert(nrhs>0);
        assert(a);
        assert(b);
    // Code starts here
        MKL_INT * ipiv = new MKL_INT[n];
        MKL_INT info;
        dgesv( (MKL_INT*) &n, (MKL_INT*) &nrhs, a, (MKL_INT*) &n, ipiv, b, (MKL_INT*) &n, &info );
        delete[] ipiv;
    }
    void blas::lapack_solve(blas_int n, blas_int nrhs, dComplex * a, dComplex * b)
    {
    // Assertions
        assert(n>0);
        assert(nrhs>0);
        assert(a);
        assert(b);
    // Code starts here
        MKL_INT * ipiv = new MKL_INT[n];
        MKL_INT info;
        zgesv( (MKL_INT*) &n, (MKL_INT*) &nrhs, (MKL_Complex16*) a, (MKL_INT*) &n, ipiv, (MKL_Complex16*) b, (MKL_INT*) &n, &info );
        delete[] ipiv;
    }
// LU factorizations, Full matrix case
    void blas::lu_factorize(blas_int n, double * a, blas_int * ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(ipiv);
    // Code starts here
        lapack_int info;
    // LAPACK - retrieves the LU factorization of matrix A
        dgetrf_( &n, &n, a, &n, ipiv, &info );
    }
    void blas::lu_factorize(blas_int n, dComplex * a,  blas_int * ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(ipiv);
    // Code starts here
        lapack_int info;
    // LAPACK - retrieves the LU factorization of matrix A
        zgetrf_( &n, &n, a, &n, ipiv, &info );
    }
// LU back substitutions, Full matrix case
    void blas::lu_back_subst(char & trans, blas_int n, double * a, double * b, blas_int * ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(b);
        assert(ipiv);
    // Code starts here
        lapack_int nrhs = 1;
        lapack_int info;

    // LAPACK - retrieves the back substitution for the factorized matrix A and one right-hand side: L*U*x = b, x->b;
        dgetrs( &trans, &n, &nrhs, a, &n, ipiv, b, &n, &info );

    }
    void blas::lu_back_subst(char & trans, blas_int n, dComplex * a, dComplex * b, blas_int * ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(b);
        assert(ipiv);
    // Code starts here

        lapack_int nrhs = 1;
        lapack_int info;
    // LAPACK - retrieves the back substitution for the factorized matrix A and one right-hand side: L*U*x = b, x->b;
        zgetrs( &trans, &n, &nrhs, a, &n, ipiv, b, &n, &info );
    }



// LU factorization, Row Compressed Matrix case
    void blas::lu_factorize_RCM(_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, dComplex * a, blas_int * ia, blas_int * ja)
    {
    // Assertions
        assert(pt);
        assert(iparm);
        assert(n>0);
        assert(a);
        assert(ia);
        assert(ja);
    // Code starts here
        memset(iparm, 0, 64*sizeof(blas_int));
        memset(pt, 0, 64*sizeof(_MKL_DSS_HANDLE_t));

        MKL_INT error  = 0;
        MKL_INT  mtype = 13;
        iparm[0]  = 0;

        pardisoinit (pt,  &mtype, iparm);

        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT phase;
        MKL_INT nrhs = 0;
        MKL_INT msglvl = 0;
        dComplex ddum;
        MKL_INT idum = 0;

        phase = 11;

        //pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
        pardiso_64 (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

         if (error != 0) {
            printf("\nERROR during symbolic factorization: %d", int(error));
            exit(1);
        }
        //printf("\nReordering completed ... ");
        //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
        //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

        phase = 22;

        //pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
        pardiso_64 (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    }
// LU Back Substitution, Row Compressed Matrix case TODO optimize this method (allocations)
    void blas::lu_back_subst_RCM (_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, dComplex * a, blas_int * ia, blas_int * ja, dComplex * b)
    {

        MKL_INT error  = 0;
        MKL_INT mtype = 13;
        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT phase = 33;
        MKL_INT nrhs = 1;
        MKL_INT msglvl = 0;
        dComplex * x = new dComplex[n];
        MKL_INT idum;

        iparm[1] = 0; /* Fill-in reordering from METIS */
        /* Numbers of processors, value of OMP_NUM_THREADS */
        iparm[2] = 3;
        iparm[7] = 1; // Maximum # of iterative refinements

        //pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error, dparm);
        pardiso_64 (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);


        for (blas_int i=0;i<n;++i){
            b[i] = x[i];
        }

        delete[] x;
        return;
    }

// Releasing internal memmory of factorized matrix
    void blas::Pardiso_clean (_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, dComplex * a, blas_int * ia, blas_int * ja)
    {
    // Assertions
        assert(pt);
        assert(iparm);
        assert(n>0);
        assert(a);
        assert(ia);
        assert(ja);
    // Code starts here
        MKL_INT error  = 0;
        MKL_INT  mtype = 13;
        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT phase = -1;     // Release internal memory
        MKL_INT nrhs = 1;
        MKL_INT msglvl = 0;
        dComplex ddum;
        MKL_INT idum;

        pardiso_64 (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        return;
    }
} // namespace QSCAT
