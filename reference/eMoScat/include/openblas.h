#include <stdlib.h>
#include <iostream>
#include <complex>

#define _MKL_DSS_HANDLE_t void*
#ifndef dComplex
typedef std::complex<double> dComplex;
#endif


using namespace std;

namespace blas {
    /*  
        This  namespace contains the classes  and functions necessary to access
        the sparse solver  routines of Intel Math Kernel Library.  The contains 
        is adapted  to the electron molecule  scattering problems, but could be
        possibly used for other codes.
    */
    

    /*  
        BLAS  functions:  General  set  of Parallelized  Basic  Linear Algebra 
        functions provided by Intel MKL Libraries. The further generalizations
        are expected. Namely implementation of PBLAS methods for multi-machine
        computations.
    */
// Array copy method
    template<typename T>
    void copy(int N, const T * X, T * Y){
        std::cout << "The copy operation is not supported for this type of variable!" << std::endl;
        exit(1024);
    }
    template<>
    void copy(int N,const int *X, int *Y);
    template<>
    void copy(int N,const double *X, double *Y);
    template<>
    void copy(int N, const dComplex *X, dComplex *Y);
    template<typename T>
    void copy(int N, const T * X, int incx, T * Y, int incy){
        std::cout << "The copy operation is not supported for this type of variable!" << std::endl;
        exit(1024);
    }
    template<>
    void copy(int N,const int *X, int incx, int *Y, int incy);
    template<>
    void copy(int N, const double* X, int incx, double* Y, int incy);
    template<>
    void copy(int N, const dComplex* X, int incx, dComplex *Y, int incy);
// Vector operation sub(y) := sub(y) + a*sub(x)
    template<typename T, typename V>
    void axpy(int n, const T a, const V * x, V * y) { // Unsupported case 
        std::cout << "Error: The vector operation 'axpy' is not supported for given types'";
        exit(114);
    }
    template<>  // Complex vs. complex case
    void axpy (int n, const dComplex a, const dComplex *x, dComplex *y);
    template<>
    void axpy (int n, const double a, const dComplex *x, dComplex *y);
    template<>
    void axpy (int n, const double a, const double *x, double *y);
  //    
    template<typename T, typename V>
    void axpy(int n, const T a, const V * x, int incx, V * y, int incy) { // Unsupported case 
        std::cout << "Error: The vector operation 'axpy' with general increment is not supported for given types'";
        exit(114);
    }
    template<>  // Complex vs. complex case
    void axpy (int n, const dComplex a, const dComplex *x, int incx, dComplex *y, int incy);
    template<>
    void axpy (int n, const double a, const dComplex *x, int incx, dComplex *y, int incy);
    template<>
    void axpy (int n, const double a, const double *x, int incx, double *y, int incy);

    template<typename T>
    void swap(int n, T *x, int incx, T *y, int incy){
        std::cout << "Error: The vector operation is not supported for given types'";
        exit(1166);
    }
    template<>
    void swap(int n, int *x, int incx, int *y, int incy);
    template<>
    void swap(int n, double *x, int incx, double *y, int incy);
    template<>
    void swap(int n, dComplex *x, int incx, dComplex *y, int incy);

// Vector dot product
    template<typename T> 
    T dotproduct(int n, const T * x, const T * y){
        std::cout << "Error: The dot product is not yet implemented for this type of vector!" << std::endl;
        exit(115);
    }
    template<> 
    dComplex dotproduct(int n, const dComplex *x, const dComplex *y);
    template<> 
    double dotproduct(int n,const double *x, const double *y);
// vector x vector - element wise multiplication
    template<typename T>
    void ewxy(int n, T * x, T * y, T * o){
        std::cout << "Error: The element-wise multiplication is not yet implemented for this type of vector!" << std::endl;
        exit(115);
    }
    template<>
    void ewxy(int n, double *x, double *y, double *o);
    template<>
    void ewxy(int n, dComplex *x, dComplex *y, dComplex *o);
    template<typename T>
    void subewxy(int N, int i, int j, T * x, T * y, T * o){
        std::cout << "Error: The element wise sub vector multiplication was not yet implemented for this type of vector." << std::endl;
        exit(116);
    }
    template<> 
    void subewxy(int n, int i, int j, double *x, double *y, double *o);
    template<>
    void subewxy(int n, int i, int j, dComplex *x, dComplex *y, dComplex *o);
// Partial dot product for projections
    template<typename T> 
    T partial_dotproduct(int n,const T * x, int incX, const T * y, int incY){
        std::cout << "Error: The dot product is not yet implemented for this type of vector!" << std::endl;
        exit(115);
    }
    template<> 
    dComplex partial_dotproduct(int n,const dComplex *x, int incX, const dComplex *y, int incY);
    template<> 
    double partial_dotproduct(int n,const double *x, int incX, const double *y, int incY);
// Vector simple scaling
    template<typename T, typename Z>
    void scale(int N, Z * a,T alpha){
        std::cout << "Error: The scaling is not yet implemented for this type of vector!" << std::endl;
        exit(116);
    }
    template<>
    void scale(int N, dComplex *a, dComplex alpha);
    template<>
    void scale(int N, dComplex *a, double alpha);
    template<>
    void scale(int N, double *a, double alpha);
    template<typename T, typename Z>
    void scale(int N, Z * a, const int inc, T alpha){
        std::cout << "Error: The scaling is not yet implemented for this type of vector!" << std::endl;
        exit(116);
    }
    template<>
    void scale(int N, dComplex *a, const int inc, dComplex alpha);
    template<>
    void scale(int N, dComplex *a, const int inc,  double alpha);
    template<>
    void scale(int N, double *a, const int inc, double alpha);


// Matrix vector multiplication
    template<typename T>
    void matrix_vector(const char & trans, int m, int n, T alpha, T beta, const T * a, const T * x, T * y){
        std::cout << "Error: The Matrix-vector product is not yet implemented for this type of vector or matrix!" << std::endl;
        exit(215);
    }
    template<>
    void matrix_vector(const char & trans, int m, int n, double alpha, double beta, const double *a, const double *x, double *y);
    template<>
    void matrix_vector(const char & trans, int m, int n, dComplex alpha, dComplex beta, const dComplex *a, const dComplex *x, dComplex *y);

    // Sub matrix case
    template<typename T>
    void sub_matrix_vector(const char & trans, int m, int n, int lda, T alpha, T beta, T * a, T * x, T * y){
        std::cout << "Error: The Matrix-vector product is not yet implemented for this type of vector or matrix!" << std::endl;
        exit(215);
    }
    template<>
    void sub_matrix_vector(const char & trans, int m, int n, int lda, double alpha, double beta, double *a, double *x, double *y);
    template<>
    void sub_matrix_vector(const char & trans, int m, int n, int lda, dComplex alpha, dComplex beta, dComplex *a, dComplex *x, dComplex *y);

// Sparse matrix vector multiplication
    template<typename Z>
    void RCmatrix_vector(char & trans, int & m, Z* alpha, Z*a, int *ia, int *ja, Z * x, Z * beta, Z * y) {
        std::cout << "Error: The Sparse Matrix-vector product is not yet implemented for this type of vector or matrix!" << std::endl;
        exit(315);
    }
    template<>
    void RCmatrix_vector(char & trans, int & m, double *alpha, double *a,int *ia, int *ja, double *x, double *beta, double *y);
    template<>
    void RCmatrix_vector(char & trans, int & m, dComplex *alpha, dComplex *a,int *ia, int *ja, dComplex *x, dComplex *beta, dComplex *y);
// Full Matrix-Matrix operations
    template<typename Z>
    void matrix_matrix(const char & transa, const char & transb, int m, int n, const int k, Z alpha, Z beta, const Z*A, const Z*B, Z*C){
        std::cout << "Error: The matrix-matrix operation is not yet defined for this type of matrices:"  << std::endl;
        exit(315);
    };
    template<>
    void matrix_matrix(const char & transa, const char & transb, int m, int n, const int k, double alpha, double beta, const double*A, const double*B, double*C);
    template<>
    void matrix_matrix(const char & transa, const char & transb, int m, int n, const int k, dComplex alpha, dComplex beta, const dComplex *A, const dComplex *B, dComplex *C);
// Eigenvalue and eigen vector functions
    /*
        Auxiliary function  for ordering values  and storing the permutation as 
        a  vector.  The ordering  is performed  through  the real  value of the 
        source.
    */
    template<typename T>
    void order_values(T * W,int *P,int n){
    // Preparing permutation vector
        for (int i=0;i<n;++i){
            P[i] = i;
        }

        int perm;
        //double val;
        //val = real(W[0]);     // starting value
        for(int i=0;i<n;++i){   // getting the order
            perm = i;
            for (int j=i;j<n;++j) {
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
    void eigen(T * a, T * w, int n) {
        std::cout << "The eigenvalue problem was not defined for this type of matrix!" << std::endl;
    }
    /*  
        Complex<float> case: calling Intel ZGEEV procedure. See the MKL manual
        for detailed information.
    */
    template<>
    void eigen(dComplex *a, dComplex *W, int n);

// void    ZGESV( MKL_INT *n, MKL_INT *nrhs, MKL_Complex16 *a, MKL_INT *lda, MKL_INT *ipiv, MKL_Complex16 *b, MKL_INT *ldb, MKL_INT *info );

// Lapack linear equations solver
    template<typename T>
    void lapack_solve(int n, int nrhs, T * a, T * b){
        std::cout << "Error. The LAPACK SOLVER is not defined for this type of matrix!" << std::endl;
        exit(260);
    }
    template<>
    void lapack_solve(int n, int nrhs, double *a, double *b);
    template<>
    void lapack_solve(int n, int nrhs, dComplex *a, dComplex *b);
// LU factorizations, Full matrix case  
    template<typename T>
    void lu_factorize(int n, T * a, int *ipiv){
        std::cout << "Error. The LU decomposition is not defined for this type of matrix!" << std::endl;
        exit(260);
    }
    template<>
    void lu_factorize(int n, double *a, int *ipiv);
    template<>
    void lu_factorize(int n, dComplex *a,  int *ipiv);
// LU back substitutions, Full matrix case  
    template<typename T>
    void lu_back_subst(char & trans, int n, T * a, T * b, int *ipiv){
        std::cout << "Error. The LU back substitution is not defined for this type of matrix!" << std::endl;
        exit(261);
    }
    template<>
    void lu_back_subst(char & trans, int n, double *a, double *b, int *ipiv);
    template<>
    void lu_back_subst(char & trans, int n, dComplex *a, dComplex *b, int *ipiv);
    
// LU factorization, Row Compressed Matrix case 
    template<typename Z>
    void lu_factorize_RCM(_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, Z * a, int *ia, int *ja){
        std::cout << "Error! The sparse solver has no interface for given variable type!" << std::endl;
        exit(2337);
    }
    template<>
    void lu_factorize_RCM(_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, dComplex *a, int *ia, int *ja);
// LU Back Substitution, Row Compressed Matrix case     
    template<typename Z>
    void lu_back_subst_RCM (_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, Z * a, int *ia, int *ja, Z * b) {
    }
    template<>
    void lu_back_subst_RCM (_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, dComplex *a, int *ia, int *ja, dComplex *b);

// Releasing internal memmory of factorized matrix
    template<typename Z>
    void Pardiso_clean (_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, Z * a, int *ia, int *ja) {
    }
    template<>
    void Pardiso_clean (_MKL_DSS_HANDLE_t * pt, int *iparm, int & n, dComplex *a, int *ia, int *ja);

}
