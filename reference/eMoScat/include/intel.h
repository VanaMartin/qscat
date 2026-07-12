#if !defined(dComplex)
typedef std::complex<double> dComplex;
typedef std::complex<float> isComp;
#endif

#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>

#include "mkl_types.h"
#include "mkl_blas.h"
#include "mkl_dfti.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#include "solver.h"         // Pardiso wrapper

typedef MKL_INT blas_int;

using namespace std;

namespace QSCAT
{
namespace blas {
    /*
        This  namespace contains the classes  and functions necessary to access
        the sparse solver  routines of Intel Math Kernel Library.  The contains
        is adapted  to the electron molecule  scattering problems, but could be
        possibly used for other codes.
    */
    //char CSRMV[6] = {'G', 'U', 'N', 'C'};     // The constant Matrix descriptor for sparse matrix vector product

    /*
        Definitions of FORTRAN-like variables for Intel MKL, added for security reasons
    */
    //#ifndef MKL_INT
    //#define MKL_INT int
    //#endif


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

    const lapack_int LAPACK_COL_MAJOR = 102;


    /*
        BLAS  functions:  General  set  of Parallelized  Basic  Linear Algebra
        functions provided by Intel MKL Libraries. The further generalizations
        are expected. Namely implementation of PBLAS methods for multi-machine
        computations.
    */
// Array copy method
    template<typename T>
    void copy(blas_int N, const T * X, T * Y)
    {
        Y[0:N] = X[0:N];
        //std::cout << "The copy operation is not supported for this type of variable!" << std::endl;
        //exit(1024);
    }
    template<>
    void copy(blas_int N, const blas_int * X, blas_int * Y);
    template<>
    void copy(blas_int N, const double * X, double * Y);
    template<>
    void copy(blas_int N, const dComplex * X, dComplex * Y);
    template<typename T>
    void copy(blas_int N, const T * X, blas_int incx, T * Y, blas_int incy)
    {
        for (blas_int i=0; i < N; ++i)
            Y[i*incy] = X[i*incx];
        //std::cout << "The copy operation is not supported for this type of variable!" << std::endl;
        //exit(1024);
    }
    template<>
    void copy(blas_int N, const blas_int * X, blas_int incx, blas_int * Y, blas_int incy);
    template<>
    void copy(blas_int N, const double* X, blas_int incx, double* Y, blas_int incy);
    template<>
    void copy(blas_int N, const dComplex* X, blas_int incx, dComplex * Y, blas_int incy);
// Vector operation sub(y) := sub(y) + a*sub(x)
    template<typename T, typename V>
    void axpy(blas_int n, const T a, const V * x, V * y) // Unsupported case
    {
        std::cout << "Error: The vector operation 'ax+y' is not supported for given types'";
        exit(114);
    }
    template<>  // blas_integer case
    void axpy (blas_int n, const blas_int a, const blas_int * x, blas_int * y);
    template<>  // Complex vs. complex case
    void axpy (blas_int n, const dComplex a, const dComplex * x, dComplex * y);
    template<>
    void axpy (blas_int n, const double a, const dComplex * x, dComplex * y);
    template<>
    void axpy (blas_int n, const double a, const double * x, double * y);
  //
    template<typename T, typename V>
    void axpy(blas_int n, const T a, const V * x, blas_int incx, V * y, blas_int incy)  // Unsupported case
    {
        std::cout << "Error: The vector operation 'ax+y' with general increment is not supported for given types'";
        exit(114);
    }
    template<>  // Complex vs. complex case
    void axpy (blas_int n, const blas_int a, const blas_int * x, blas_int incx, blas_int * y, blas_int incy);
    template<>  // Complex vs. complex case
    void axpy (blas_int n, const dComplex a, const dComplex * x, blas_int incx, dComplex * y, blas_int incy);
    template<>
    void axpy (blas_int n, const double a, const dComplex * x, blas_int incx, dComplex * y, blas_int incy);
    template<>
    void axpy (blas_int n, const double a, const double * x, blas_int incx, double * y, blas_int incy);

// Vector dot product
    template<typename T>
    T dotproduct(blas_int n,const T * x, const T * y)
    {
        std::cout << "Error: The dot product is not yet implemented for this type of vector!" << std::endl;
        exit(115);
    }
    template<>
    dComplex dotproduct(blas_int n,const dComplex * x, const dComplex * y);
    template<>
    double dotproduct(blas_int n,const double * x, const double * y);
// Vector reduction
    template<typename T>
    T reduct(blas_int n,const T * x, const T * y)
    {
        std::cout << "Error: The dot product is not yet implemented for this type of vector!" << std::endl;
        exit(115);
    }
    template<>
    dComplex reduct(blas_int n,const dComplex * x, const dComplex * y);
    template<>
    double reduct(blas_int n,const double * x, const double * y);
// vector x vector - element wise multiplication
    template<typename T>
    void ewxy(blas_int n, const T * x, const T * y, T * o)
    {
        std::cout << "Error: The element-wise multiplication is not yet implemented for this type of vector!" << std::endl;
        exit(115);
    }
    template<>
    void ewxy(blas_int n, const double * x, const double * y, double * o);
    template<>
    void ewxy(blas_int n, const dComplex * x, const dComplex * y, dComplex * o);
    template<typename T>
    void subewxy(blas_int N, blas_int i, blas_int j, T * x, T * y, T * o)
    {
        std::cout << "Error: The element wise sub vector multiplication was not yet implemented for this type of vector." << std::endl;
        exit(116);
    }
    template<>
    void subewxy(blas_int n, blas_int i, blas_int j, double * x, double * y, double * o);
    template<>
    void subewxy(blas_int n, blas_int i, blas_int j, dComplex * x, dComplex * y, dComplex * o);
// Partial dot product for projections
    template<typename T>
    T partial_dotproduct(blas_int n,const T * x, blas_int incX, const T * y, blas_int incY)
    {
        std::cout << "Error: The dot product is not yet implemented for this type of vector!" << std::endl;
        exit(115);
    }
    template<>
    dComplex partial_dotproduct(blas_int n,const dComplex * x, blas_int incX, const dComplex * y, blas_int incY);
    template<>
    double partial_dotproduct(blas_int n,const double * x, blas_int incX, const double * y, blas_int incY);
// Vector simple scaling
    template<typename T, typename Z>
    void scale(blas_int N, Z * a, T alpha)
    {
        std::cout << "Error: The scaling is not yet implemented for this type of vector!" << std::endl;
        exit(116);
    }
    template<>
    void scale(blas_int N, dComplex * a, dComplex alpha);
    template<>
    void scale(blas_int N, dComplex * a, double alpha);
    template<>
    void scale(blas_int N, double * a, double alpha);
    template<typename T, typename Z>
    void scale(blas_int N, Z * a, blas_int inc, T alpha)
    {
        std::cout << "Error: The scaling is not yet implemented for this type of vector!" << std::endl;
        exit(116);
    }
    template<>
    void scale(blas_int N, dComplex * a, blas_int inc, dComplex alpha);
    template<>
    void scale(blas_int N, dComplex * a, blas_int inc,  double alpha);
    template<>
    void scale(blas_int N, double * a, blas_int inc, double alpha);

    template<typename T>
    void conj(blas_int N, T *a)
    {
        std::cout << "Error: The conjugation is not yet implemented for this type of vector!" << std::endl;
        exit(116);
    }
    template<>
    void conj(blas_int N, double *a);
    template<>
    void conj(blas_int N, dComplex *a);

// Vector Swap procedure
    template<typename T>
    void swap(blas_int n, T *a, T *b)
    {
        T aux;
        for (blas_int i=0; i<n; ++i){
            aux=a[i];
            a[i]=b[i];
            b[i]=aux;
        }
    }
    template<>
    void swap(blas_int n, double *a, double *b)
    {
        cblas_dswap(n, a, 1, b, 1);
    }
    template<>
    void swap(blas_int n, dComplex *a, dComplex *b)
    {
        cblas_zswap(n, a, 1, b, 1);
    }
    template<typename T>
    void swap(blas_int n, T *a, blas_int inca, T *b, blas_int incb)
    {
        T aux;
        for (blas_int i=0; i<n; ++i){
            aux=a[i*inca];
            a[i*inca]=b[i*incb];
            b[i*incb]=aux;
        }
    }
    template<>
    void swap(blas_int n, double *a, blas_int inca, double *b, blas_int incb)
    {
        cblas_dswap(n, a, inca, b, incb);
    }
    template<>
    void swap(blas_int n, dComplex *a, blas_int inca, dComplex *b, blas_int incb)
    {
        cblas_zswap(n, a, inca, b, incb);
    }

// Matrix vector multiplication
    template<typename T>
    void matrix_vector(const char & trans, blas_int m, blas_int n, T alpha, T beta, const T * a, const T * x, T * y)
    {
        std::cout << "Error: The Matrix-vector product is not yet implemented for this type of vector or matrix!" << std::endl;
        exit(215);
    }
    template<>
    void matrix_vector(const char & trans, blas_int m, blas_int n, double alpha, double beta, const double * a, const double * x, double * y);
    template<>
    void matrix_vector(const char & trans, blas_int m, blas_int n, dComplex alpha, dComplex beta, const dComplex * a, const dComplex * x, dComplex * y);

    // Sub matrix case
    template<typename T>
    void sub_matrix_vector(const char & trans, blas_int m, blas_int n, blas_int lda, T alpha, T beta, T * a, T * x, T * y)
    {
        std::cout << "Error: The Matrix-vector product is not yet implemented for this type of vector or matrix!" << std::endl;
        exit(215);
    }
    template<>
    void sub_matrix_vector(const char & trans, blas_int m, blas_int n, blas_int lda, double alpha, double beta, double * a, double * x, double * y);
    template<>
    void sub_matrix_vector(const char & trans, blas_int m, blas_int n, blas_int lda, dComplex alpha, dComplex beta, dComplex * a, dComplex * x, dComplex * y);

// Sparse matrix vector multiplication
    template<typename Z>
    void RCmatrix_vector(char trans, blas_int m, Z alpha, const Z *a, const blas_int *ia, const blas_int *ja, const Z *x, Z beta, Z *y)
    {
        std::cout << "Error: The Sparse Matrix-vector product is not yet implemented for this type of vector or matrix!" << std::endl;
        exit(315);
    }
    template<>
    void RCmatrix_vector(char trans, blas_int m, double alpha, const double *a, const blas_int *ia, const blas_int *ja, const double *x, double beta, double *y);
    template<>
    void RCmatrix_vector(char trans, blas_int m, dComplex alpha, const dComplex *a, const blas_int *ia, const blas_int *ja, const dComplex *x, dComplex beta, dComplex *y);
// Full Matrix-Matrix operations
    template<typename Z>
    void matrix_matrix(const char & transa, const char & transb, blas_int m, blas_int n, const blas_int k, Z alpha, Z beta, const Z*A, const Z*B, Z*C)
    {
        std::cout << "Error: The matrix-matrix operation is not yet defined for this type of matrices:"  << std::endl;
        exit(315);
    }
    template<>
    void matrix_matrix(const char & transa, const char & transb, blas_int m, blas_int n, const blas_int k, double alpha, double beta, const double*A, const double*B, double*C);
    template<>
    void matrix_matrix(const char & transa, const char & transb, blas_int m, blas_int n, const blas_int k, dComplex alpha, dComplex beta, const dComplex*A, const dComplex*B, dComplex*C);
// Eigenvalue and eigen vector functions
    /*
        Auxiliary function  for ordering values  and storing the permutation as
        a  vector.  The ordering  is performed  through  the real  value of the
        source.
    */
    template<typename T>
    void order_values(T * W,blas_int * P,blas_int n)
    {
    // Preparing permutation vector
        for (blas_int i=0;i<n;++i){
            P[i] = i;
        }

        blas_int perm;
        for(blas_int i=0;i<n;++i){   // getting the order
            perm = i;
            for (blas_int j=i;j<n;++j) {
                if (std::real(W[perm])>std::real(W[j])) {
                    perm = j;
                }
            }
            std::swap(W[i],W[perm]);
            std::swap(P[i],P[perm]);
        }
        return;
    }
    /*
        Default  case  for  the  LAPACK  generalized  eigenvalue  problem  as a
        template function.  The  default matrix type  is deprecated. Only cases
        such as complex<double> or double are allowed.
    */
    template<typename T>
    void eigen(T * a, T * w, blas_int n)
    {
        std::cout << "The eigenvalue problem was not defined for this type of matrix!" << std::endl;
    }
    /*
        Complex<float> case: calling Intel ZGEEV procedure. See the MKL manual
        for detailed information.
    */
//  template<>
//  void eigen(isComp * a, isComp * W, int n);
    /*
        Complex<double> case: calling Intel ZGEEV procedure. See the MKL manual
        for detailed information.
    */
    template<>
    void eigen(dComplex * a, dComplex * W, blas_int n);

// void    ZGESV( MKL_INT *n, MKL_INT *nrhs, MKL_Complex16 *a, MKL_INT *lda, MKL_INT *ipiv, MKL_Complex16 *b, MKL_INT *ldb, MKL_INT *info );

// Lapack linear equations solver
    template<typename T>
    void lapack_solve(blas_int n, blas_int nrhs, T * a, T * b)
    {
        std::cout << "Error. The LAPACK SOLVER is not defined for this type of matrix!" << std::endl;
        exit(260);
    }
    template<>
    void lapack_solve(blas_int n, blas_int nrhs, double * a, double * b);
    template<>
    void lapack_solve(blas_int n, blas_int nrhs, dComplex * a, dComplex * b);
// LU factorizations, Full matrix case
    template<typename T>
    void lu_factorize(blas_int n, T * a, blas_int * ipiv)
    {
        std::cout << "Error. The LU decomposition is not defined for this type of matrix!" << std::endl;
        exit(260);
    }
    template<>
    void lu_factorize(blas_int n, double * a, blas_int * ipiv);
    template<>
    void lu_factorize(blas_int n, dComplex * a,  blas_int * ipiv);
// LU back substitutions, Full matrix case
    template<typename T>
    void lu_back_subst(char & trans, blas_int n, T * a, T * b, blas_int * ipiv)
    {
        std::cout << "Error. The LU back substitution is not defined for this type of matrix!" << std::endl;
        exit(261);
    }
    template<>
    void lu_back_subst(char & trans, blas_int n, double * a, double * b, blas_int * ipiv);
    template<>
    void lu_back_subst(char & trans, blas_int n, dComplex * a, dComplex * b, blas_int * ipiv);

// LU factorization, Row Compressed Matrix case
    template<typename Z>
    void lu_factorize_RCM(_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, Z * a, blas_int * ia, blas_int * ja)
    {
        std::cout << "Error! The sparse solver has no interface for given variable type!" << std::endl;
        exit(2337);
    }
    template<>
    void lu_factorize_RCM(_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, dComplex * a, blas_int * ia, blas_int * ja);
// LU Back Substitution, Row Compressed Matrix case
    template<typename Z>
    void lu_back_subst_RCM (_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, Z * a, blas_int * ia, blas_int * ja, Z * b)
    {
    }
    template<>
    void lu_back_subst_RCM (_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, dComplex * a, blas_int * ia, blas_int * ja, dComplex * b);

// Releasing internal memmory of factorized matrix
    template<typename Z>
    void Pardiso_clean (_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, Z * a, blas_int * ia, blas_int * ja)
    {
    }
    template<>
    void Pardiso_clean (_MKL_DSS_HANDLE_t * pt, blas_int * iparm, blas_int & n, dComplex * a, blas_int * ia, blas_int * ja);

}
} // namespace QSCAT
