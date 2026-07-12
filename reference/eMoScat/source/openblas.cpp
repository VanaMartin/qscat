#include <iostream>
#include <fstream>
#include <complex>

#include "common.h"

#include "openblas.h"
#include "cblas.h"
#include "lapacke.h"
#include <cassert>

using namespace std;
namespace QSCAT
{
namespace blas
{
// Simplified call
    template<>
    void copy(int N, const int *X, int *Y)
    {
        assert(X);
        assert(Y);
        assert(N>0);
    //
        for (int i=0;i<N;++i){
            Y[i] = X[i];
        }
        return;
    }
    template<>
    void copy(int N,const double *X, double *Y)
    {
        assert(X);
        assert(Y);
        assert(N>0);
    //
        cblas_dcopy (N, X, 1, Y, 1);
        return;
    }
    template<>
    void copy(int N, const dComplex *X, dComplex *Y)
    {
        assert(X);
        assert(Y);
        assert(N>0);
    //
        cblas_zcopy (N, (double*) X, 1, (double*) Y, 1);
        return;
    }
// Full Call
    template<>
    void copy(int N,const int *X, int incx, int *Y, int incy)
    {
        assert(X);
        assert(Y);
        assert(N>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        for (int i=0;i<N;++i){
            Y[i*incy] = X[i*incx];
        }
        return;
    }
    template<>
    void copy(int N, const double *X, int incx, double *Y, int incy)
    {
        assert(X);
        assert(Y);
        assert(N>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        cblas_dcopy (N, X, incx, Y, incy);
        return;
    }
    template<>
    void copy(int N, const dComplex *X, int incx, dComplex *Y, int incy)
    {
        assert(X);
        assert(Y);
        assert(N>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        cblas_zcopy (N, (double*) X, incx, (double*) Y, incy);
        return;
    }
// Vector operation sub(y) := sub(y) + a*sub(x)
    template<>
    void axpy (int n, const dComplex a, const dComplex *x, dComplex *y)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
    //
        cblas_zaxpy(n,(double*) &a, (double*) x,1, (double*) y,1);
        return;
    }
    template<>
    void axpy (int n, const double a, const dComplex *x, dComplex *y)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
    //
        dComplex an = a;
        cblas_zaxpy(n,(double*) &an, (double*) x,1, (double*) y,1);
        return;
    }
    template<>
    void axpy (int n, const double a, const double *x, double *y)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
    //
        cblas_daxpy(n,a,x,1,y,1);
        return;
    }
    template<>
    void axpy (int n, const dComplex a, const dComplex *x, int incx, dComplex *y, int incy)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        cblas_zaxpy(n, (double*) &a, (double*) x, incx, (double*) y, incy);
        return;
    }
    template<>
    void axpy (int n, const double a, const dComplex *x, int incx, dComplex *y, int incy)
    {
        assert(a);
        assert(x);
        assert(y);
        assert(n>0);
        assert(incx>=0);
        assert(incy>=0);
    //
        dComplex an = a;
        cblas_zaxpy(n, (double*) &an, (double*) x, incx, (double*) y, incy);
        return;
    }
    template<>
    void axpy (int n, const double a, const double *x, int incx, double *y, int incy)
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

    template<>
    void swap(int n, int *x, int incx, int *y, int incy)
    {
        assert(x);
        assert(y);
        assert(incx>=0);
        assert(incy>=0);

        if (n) {
            //cblas_iswap(n, x, incx, y, incy);
            int aux;
            for (int i=0; i<n; ++i){
                aux = x[i*incx];
                x[i*incx] = y[i*incy];
                y[i*incy] = aux;
            }
        }
    }
    template<>
    void swap(int n, double *x, int incx, double *y, int incy)
    {
        assert(x);
        assert(y);
        assert(incx>=0);
        assert(incy>=0);

        if (n) {
            cblas_dswap(n, x, incx, y, incy);
        }
    }
    template<>
    void swap(int n, dComplex *x, int incx, dComplex *y, int incy)
    {
        assert(x);
        assert(y);
        assert(incx>=0);
        assert(incy>=0);

        if (n) {
            cblas_zswap(n, (double*) x, incx, (double*) y, incy);
        }
    }


// Vector dot product
    template<>
    dComplex dotproduct(int n, const dComplex *x, const dComplex *y)
    {
        assert(n>0);
        assert(x);
        assert(y);
    //
        dComplex out;
        cblas_zdotc_sub(n, (double*) x, 1, (double*) y, 1, (openblas_complex_double*) &out);
        return out;
    }
    template<>
    double dotproduct(int n, const double *x, const double *y)
    {
        assert(n>0);
        assert(x);
        assert(y);
    //
        double out;
        out = cblas_ddot(n,x,1,y,1);
        return out;
    }
// vector x vector - element wise multiplication
    template<>
    void ewxy(int n, double *x, double *y, double *o)
    {
        assert(n>0);
        assert(x);
        assert(y);
        assert(o);
    //
        // FIXME vdMul( n, x, y, o);
        #pragma omp parallel for
        for (int i=0; i<n; ++i){
            o[i] = x[i]*y[i];
        }
    }
    template<>
    void ewxy(int n, dComplex *x, dComplex *y, dComplex *o)
    {
        assert(n>0);
        assert(x);
        assert(y);
        assert(o);
    //
        // FIXME vzMul( n, (MKL_Complex16*) x, (MKL_Complex16*) y, o);
        #pragma omp parallel for
        for (int i=0; i<n; ++i){
            o[i] = x[i]*y[i];
        }
    }
    template<>
    void subewxy(int n, int i, int j, double *x, double *y, double *o)
    {
        assert(n>0);
        assert(i>=0);
        assert(j>=0);
        assert(x);
        assert(y);
        assert(o);
    //
        // FIXME vdMul( n,  &(x[i]), &(y[j]), o);
        #pragma omp parallel for
        for (int i=0; i<n; ++i){
            o[i] = x[i]*y[i];
        }
    }
    template<>
    void subewxy(int n, int i, int j, dComplex *x, dComplex *y, dComplex *o)
    {
        assert(n>0);
        assert(i>=0);
        assert(j>=0);
        assert(x);
        assert(y);
        assert(o);
    //
        // FIXME vzMul( n, (MKL_Complex16*) &(x[i]), (MKL_Complex16*) &(y[j]), o);
        #pragma omp parallel for
        for (int i=0; i<n; ++i){
            o[i] = x[i]*y[i];
        }
    }
// Partial dot product for projections
    template<>
    dComplex partial_dotproduct(int n,const dComplex *x, int incX, const dComplex *y, int incY)
    {
        assert(n>0);
        assert(incX>=0);
        assert(incY>=0);
        assert(x);
        assert(y);
    //
        dComplex out;
        cblas_zdotc_sub(n, (double*) x, incX, (double*) y, incY, (openblas_complex_double*) &out);
        return out;
    }
    template<>
    double partial_dotproduct(int n,const double *x, int incX, const double *y, int incY)
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
    template<>
    void scale(int N, dComplex *a, dComplex alpha)
    {
        assert(N>0);
        assert(a);
    ////
        cblas_zscal(N, (const double*) &alpha, (double*) a, 1);
    }
    template<>
    void scale(int N, dComplex *a, double alpha)
    {
        assert(N>0);
        assert(a);
    ////
        cblas_zdscal(N, alpha, (double*) a, 1);
    }
    template<>
    void scale(int N, double *a, double alpha)
    {
        assert(N>0);
        assert(a);
    ////
        cblas_dscal(N, alpha, a, 1);
    }
    template<>
    void scale(int N, dComplex *a, const int inc, dComplex alpha)
    {
        assert(N>0);
        assert(a);
        assert(inc>=0);
    ////
        cblas_zscal(N, (const double*) &alpha, (double*) a, inc);
    }
    template<>
    void scale(int N, dComplex *a, const int inc, double alpha)
    {
        assert(N>0);
        assert(a);
        assert(inc>=0);
    ////
        cblas_zdscal(N, alpha, (double*) a, inc);
    }
    template<>
    void scale(int N, double *a, const int inc, double alpha)
    {
        assert(N>0);
        assert(a);
        assert(inc>=0);
    ////
        cblas_dscal(N, alpha, a, inc);
    }

// Matrix vector multiplication
    template<>
    void matrix_vector(const char & trans, int m, int n, double alpha, double beta, const double *a, const double *x, double *y)
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

    // ColumnMajor version for matrix definition conventions
        cblas_dgemv(CblasColMajor, C, m, n, alpha, a, m, x, 1, beta, y, 1);
        return;
    }
    template<>
    void matrix_vector(constzo char & trans, int m, int n, dComplex alpha, dComplex beta, const dComplex *a, const dComplex *x, dComplex *y)
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

        cblas_zgemv(CblasColMajor, C, m, n, (double*) &alpha, (double*) a, m, (double*) x, 1, (double*) &beta, (double*) y, 1);
        return;
    }

    // Sub matrix case
    template<>
    void sub_matrix_vector(const char & trans, int m, int n, int lda, double alpha, double beta, double *a, double *x, double *y)
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

    // ColumnMajor version for matrix definition conventions
        cblas_dgemv(CblasColMajor, C, m, n, alpha, a, lda, x, 1, beta, y, 1);
        return;
    }
    template<>
    void sub_matrix_vector(const char & trans, int m, int n, int lda, dComplex alpha, dComplex beta, dComplex *a, dComplex *x, dComplex *y)
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

        cblas_zgemv(CblasColMajor, C, m, n, (double*) &alpha, (double*) a, lda, (double*) x, 1, (double*) &beta, (double*) y, 1);
        return;
    }

// Sparse matrix vector multiplication
    template<>
    void RCmatrix_vector(char & trans, int & m, double *alpha, double *a,int *ia, int *ja, double *x, double *beta, double *y)
    {
    // Assertions
//      assert(m>0);
//      assert(alpha);
//      assert(a);
//      assert(ia);
//      assert(ja);
//      assert(x);
//      assert(beta);
//      assert(y);
//  // Code starts here
//      mkl_dcsrmv(&trans, &m, &m, alpha, CSRMV, a, ja, &(ia[0]), &(ia[1]), x, beta, y);
        cout << "Error, Sparse BLAS not implemented yet!" << endl;
        exit(-1);
    }
    template<>
    void RCmatrix_vector(char & trans, int & m, dComplex *alpha, dComplex *a,int *ia, int *ja, dComplex *x, dComplex *beta, dComplex *y)
    {
    // Assertions
        assert(m>0);
//      assert(alpha);
//      assert(a);
//      assert(ia);
//      assert(ja);
//      assert(x);
//      assert(beta);
//      assert(y);
//  // Code starts here
//      mkl_zcsrmv(&trans, &m, &m, alpha, CSRMV, a, ja, &(ia[0]), &(ia[1]), x, beta, y);
        cout << "Error, Sparse BLAS not implemented yet!" << endl;
        exit(-1);
    };
// Full Matrix-Matrix operations
    template<>
    void matrix_matrix(const char & transa, const char & transb, int m, int n, const int k, double alpha, double beta, const double*A, const double*B, double*C)
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
        cblas_dgemm(CblasColMajor,TA, TB, m, n, k, alpha, A, (transa=='N')? m:k, B, (transb=='N')? k:n, beta, C, m);
    }
    template<>
    void matrix_matrix(const char & transa, const char & transb, int m, int n, const int k, dComplex alpha, dComplex beta, const dComplex *A, const dComplex *B, dComplex *C)
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
        cblas_zgemm(CblasColMajor, TA, TB, m, n, k, (double*) &alpha, (double*) A, (transa=='N')? m:k, (double*) B, (transb=='N')? k:n, (double*) &beta, (double*) C, m);
    }


// Eigenvalue and eigen vector functions
    /*
        Auxiliary function  for ordering values  and storing the permutation as
        a  vector.  The ordering  is performed  through  the real  value of the
        source.
    */
    //template<typename T>
    //void order_values(T * W,int *P,int n){
    // Preparing permutation vector
    //  for (int i=0;i<n;++i){
    //      P[i] = i;
    //  }

    //  int perm;
    //  double val;
    //  val = real(W[0]);       // starting value
    //  for(int i=0;i<n;++i){   // getting the order
    //      perm = i;
    //      for (int j=i;j<n;++j) {
    //          if (std::real(W[perm])>std::real(W[j])) {
    //              perm = j;
    //          }
    //      }
    //      std::swap(W[i],W[perm]);
    //      std::swap(P[i],P[perm]);
    //  }
    //  return;
    //}
    /*
        Default  case  for  the  LAPACK  generalized  eigenvalue  problem  as a
        template function.  The  default matrix type  is deprecated. Only cases
        such as complex<double> or double are allowed.
    */
    /*
        Complex<float> case: calling Intel ZGEEV procedure. See the MKL manual
        for detailed information.
    */
/*  void eigen(comp0 * a, comp0 * W, int n) {
        char jobvl;
        char jobvr;
        lapack_int N;
        lapack_int lda;
        lapack_int ldvl;
        lapack_int ldvr;
        MKL_Complex8 * w = new MKL_Complex8[n];     // The array of eigen values
        MKL_Complex8 * vr = new MKL_Complex8[n*n];  // Right hand eigenvector space pointer
        MKL_Complex8 * vl = new MKL_Complex8;           // Left hand eigenvector space pointer
        MKL_Complex8 * work = new MKL_Complex8[1];  // The function workspace
        MKL_INT lwork;
        float * rwork = new float[2*n];
        MKL_INT info;

        jobvl = 'N';                        // Left hand eigenvectors disabled
        jobvr = 'V';                        // Right hand eigenvectors enabled
        N = n;                              // The order of the matrix elements stored as a vector
        lda = n;                            // Array dimension
        ldvl = 1;                           // The leading dimension for left eigenvectors to be computed
        ldvr = n;                           // The leading dimension for right eigenvectors to be computed
        lwork = -1;

    // The first run returns the appropriate space of the workspace variable
        cgeev(&jobvl, &jobvr, &N, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

        lwork = abs(work[0]);               // Proper workspace dimension
        delete[] work;                      // Cleaning old workspace
        work = new MKL_Complex8[lwork+1];   // Allocating new workspace

    // The second call invokes the solver to compute the eigenvalues
        cgeev(&jobvl, &jobvr, &N, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

    // Overwriting the matrix a with the eigenvalues ordered by real value of energy
        int *order = new int[n];
        order_values(w, order, n);

        for (int i=0;i<n;++i){
            W[i] = w[i];
            for (int j=0;j<n;++j){
                a[i*n + j] = vr[order[i]*n + j];
            }
        }

        delete[] work;
        delete[] rwork;
        delete[] w ;
        delete[] vr;
        delete[] vl;
        return;
    }
*/
    /*
        Complex<double> case: calling Intel ZGEEV procedure. See the MKL manual
        for detailed information.
    */
    template<>
    void eigen(dComplex *a, dComplex *W, int n)
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
        dComplex *w = new dComplex[n];
        lapack_complex_double * vr = new lapack_complex_double[n*n];    // Right hand eigenvector space pointer
        lapack_complex_double * vl = NULL; // new lapack_complex_double;            // Left hand eigenvector space pointer
        //lapack_complex_double * work = new lapack_complex_double[1];  // The function workspace
        //int lwork;
        //double *rwork = new double[2*n];
        int info;

        jobvl = 'N';                        // Left hand eigenvectors disabled
        jobvr = 'V';                        // Right hand eigenvectors enabled
        N = n;                              // The order of the matrix elements stored as a vector
        lda = n;                            // Array dimension
        ldvl = 1;                           // The leading dimension for left eigenvectors to be computed
        ldvr = n;                           // The leading dimension for right eigenvectors to be computed
        //lwork = -1;

    // The first run returns the appropriate space of the workspace variable
        info = LAPACKE_zgeev(LAPACK_COL_MAJOR, jobvl, jobvr, N, (lapack_complex_double*) a, lda, (lapack_complex_double*) w, vl, ldvl, vr, ldvr);

        if (info) {
            cout << "Zgeev failed with argument " << info << endl;
            exit(info);
        }
//      lwork = std::abs((std::complex<double>) work[0]);               // Proper workspace dimension
//      delete[] work;                      // Cleaning old workspace
//      work = new lapack_complex_double[lwork+1];  // Allocating new workspace

    // The second call invokes the solver to compute the eigenvalues
//      info = LAPACKE_zgeev(LAPACK_COL_MAJOR, jobvl, jobvr, N, (lapack_complex_double*) a, lda, (lapack_complex_double*) w, vl, ldvl, vr, ldvr, work, lwork, rwork);

    // Overwriting the matrix a with the eigenvalues ordered by real value of energy
        int *order = new int[n];
        order_values(w, order, n);
        // TODO copy()
        for (int i=0;i<n;++i){
            W[i] = w[i];
            for (int j=0;j<n;++j){
                a[i*n + j] = vr[order[i]*n + j];
            }
        }

        //delete[] work;
        //delete[] rwork;
        delete[] w ;
        delete[] vr;
        delete vl;
        return;
    }

// void    ZGESV( MKL_INT *n, MKL_INT *nrhs, MKL_Complex16 *a, MKL_INT *lda, MKL_INT *ipiv, MKL_Complex16 *b, MKL_INT *ldb, MKL_INT *info );

// Lapack linear equations solver
    template<>
    void lapack_solve(int n, int nrhs, double *a, double *b)
    {
    // Assertions
        assert(n>0);
        assert(nrhs>0);
        assert(a);
        assert(b);
    // Code starts here
        int * ipiv = new int[n];
        int info;
        info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, a, n, ipiv, b, n);
        if (info) {
            cout << "lapack solver failed with argument " << info << endl;
            exit(info);
        }
        delete[] ipiv;
    }
    template<>
    void lapack_solve(int n, int nrhs, dComplex *a, dComplex *b)
    {
    // Assertions
        assert(n>0);
        assert(nrhs>0);
        assert(a);
        assert(b);
    // Code starts here
        int * ipiv = new int[n];
        int info;
        info = LAPACKE_zgesv( LAPACK_COL_MAJOR, n, nrhs, (lapack_complex_double*) a, n, ipiv, (lapack_complex_double*) b, n );
        if (info) {
            cout << "lapack solver failed with argument " << info << endl;
            exit(info);
        }
        delete[] ipiv;
    }
// LU factorizations, Full matrix case
    template<>
    void lu_factorize(int n, double *a, int *ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(ipiv);
    // Code starts here
        lapack_int N = n;
        lapack_int info;
    // LAPACK - retrieves the LU factorization of matrix A
        info = LAPACKE_dgetrf( LAPACK_COL_MAJOR, N, N, a, N, ipiv );
        if (info) {
            cout << "lapack solver failed with argument " << info << endl;
            exit(info);
        }
    }
    template<>
    void lu_factorize(int n, dComplex *a,  int *ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(ipiv);
    // Code starts here
        lapack_int N = n;
        lapack_int info;
    // LAPACK - retrieves the LU factorization of matrix A
        info = LAPACKE_zgetrf( LAPACK_COL_MAJOR, N, N, (lapack_complex_double*) a, N, ipiv );
        if (info) {
            cout << "lapack solver failed with argument " << info << endl;
            exit(info);
        }
    }
// LU back substitutions, Full matrix case
    template<>
    void lu_back_subst(char & trans, int n, double *a, double *b, int *ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(b);
        assert(ipiv);
    // Code starts here
        lapack_int N = n;
        lapack_int nrhs = 1;
        lapack_int info;

    // LAPACK - retrieves the back substitution for the factorized matrix A and one right-hand side: L*U*x = b, x->b;
        info = LAPACKE_dgetrs( LAPACK_COL_MAJOR, trans, N, nrhs, a, N, ipiv, b, N );
        if (info) {
            cout << "lapack solver failed with argument " << info << endl;
            exit(info);
        }

    }
    template<>
    void lu_back_subst(char & trans, int n, dComplex *a, dComplex *b, int *ipiv)
    {
    // Assertions
        assert(n>0);
        assert(a);
        assert(b);
        assert(ipiv);
    // Code starts here

        lapack_int N = n;
        lapack_int nrhs = 1;
        lapack_int info;
    // LAPACK - retrieves the back substitution for the factorized matrix A and one right-hand side: L*U*x = b, x->b;
        info = LAPACKE_zgetrs( LAPACK_COL_MAJOR, trans, N, nrhs, (lapack_complex_double*) a, N, ipiv, (lapack_complex_double*) b, N );
        if (info) {
            cout << "lapack solver failed with argument " << info << endl;
            exit(info);
        }
    }



// LU factorization, Row Compressed Matrix case
    template<>
    void lu_factorize_RCM(_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, dComplex *a, int *ia, int *ja)
    {
    // Assertions
//      assert(pt);
//      assert(iparm);
//      assert(n>0);
//      assert(a);
//      assert(ia);
//      assert(ja);
//  // Code starts here
//      MKL_INT error  = 0;
//      MKL_INT  mtype = 13;
//      iparm[0]  = 0;
//
//      pardisoinit (pt,  &mtype, iparm);
//
//      MKL_INT maxfct = 1;
//      MKL_INT mnum = 1;
//      MKL_INT phase;
//      MKL_INT nrhs = 0;
//      MKL_INT msglvl = 0;
//      dComplex ddum;
//      int idum = 0;
//
//      phase = 11;
//
//      //pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
//      pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//
//       if (error != 0) {
//          printf("\nERROR during symbolic factorization: %d", error);
//          exit(1);
//      }
//      //printf("\nReordering completed ... ");
//      //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
//      //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
//
//      phase = 22;
//
//      //pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
//      pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        cout << "Error, Sparse BLAS not implemented yet!" << endl;
        exit(-1);
    }
// LU Back Substitution, Row Compressed Matrix case
    template<>
    void lu_back_subst_RCM (_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, dComplex *a, int *ia, int *ja, dComplex *b)
    {

//      MKL_INT error  = 0;
//      MKL_INT  mtype = 13;
//      MKL_INT solver = 0;
//      MKL_INT maxfct = 1;
//      MKL_INT mnum = 1;
//      MKL_INT phase = 33;
//      MKL_INT nrhs = 1;
//      MKL_INT msglvl = 0;
//      dComplex *x = new comp[n];
//      int idum;
//
//      iparm[1] = 2; /* Fill-in reordering from METIS */
//      /* Numbers of processors, value of OMP_NUM_THREADS */
//      iparm[2] = 2;
//      iparm[7] = 1; // Maximum # of iterative refinements
//
//      //pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error, dparm);
//      pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
//
//
//      for (int i=0;i<n;++i){
//          b[i] = x[i];
//      }
//
//      delete[] x;
//      return;
        cout << "Error, Sparse BLAS not implemented yet!" << endl;
        exit(-1);
    }

// Releasing internal memmory of factorized matrix
    template<>
    void Pardiso_clean (_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, dComplex *a, int *ia, int *ja)
    {
//  // Assertions
//      assert(pt);
//      assert(iparm);
//      assert(n>0);
//      assert(a);
//      assert(ia);
//      assert(ja);
//  // Code starts here
//      MKL_INT error  = 0;
//      MKL_INT  mtype = 13;
//      MKL_INT solver = 0;
//      MKL_INT maxfct = 1;
//      MKL_INT mnum = 1;
//      MKL_INT phase = -1;     // Release internal memory
//      MKL_INT nrhs = 1;
//      MKL_INT msglvl = 0;
//      dComplex ddum;
//      int idum;
//
//      pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//      return;
        cout << "Error, Sparse BLAS not implemented yet!" << endl;
        exit(-1);
    }
} // namespace blas
} // namespace QSCAT
