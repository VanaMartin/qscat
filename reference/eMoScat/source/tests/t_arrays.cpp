#include <iostream>             // Input/Output library
#include <cassert>

#include <fstream>
#include <stdio.h>
#include <string>
#include <complex>              // Complex algebra

#include "common.h"
#include "bessel.h"     // Solutions to the coulomb potential Schrödinger equation (sph. Bessel functions)

#include "Arrays.h"     // Library to be tested
#include "interface.h"

#include "tests/templates.hpp"

using namespace QSCAT;

int main(){

    int N = 65535;
    cout << endl << "Testing vectors of size " << N << endl;

    {
        cout << endl;
        iVector X(N);
        iVector Y(N);
        for (blas_int i=0; i<N; ++i) {
            X[i] = i;
            Y[i] = i + N;
        }

        TESTS::SECTION = "Integer Vector";
        TEST(copy, X);
        TEST(fill, X, 0);
        TEST(swap, X, Y);
        TEST(inplace_add, X, Y);
        TEST(inplace_sub, X, Y);
        //TEST("Integer vector", TESTS::inplace_scaling, X, 4); TODO fix this one
    }

    {
        cout << endl;
        dVector X(N);
        dVector Y(N);
        for (blas_int i=0; i<N; ++i) {
            X[i] = i;
            Y[i] = i + N;
        }

        TESTS::SECTION = "Double Vector";
        TEST(copy, X);
        TEST(fill, X, 0.0);
        TEST(swap, X, Y);
        TEST(inplace_add, X, Y);
        TEST(inplace_sub, X, Y);
        TEST(inplace_scaling, X, 4.0);
    }

    {
        cout << endl;
        zVector X(N);
        zVector Y(N);
        for (blas_int i=0; i<N; ++i) {
            X[i] = i;
            Y[i] = i + N;
        }

        TESTS::SECTION = "Complex Double Vector";
        TEST(copy, X);
        TEST(fill, X, 0.0);
        TEST(swap, X, Y);
        TEST(inplace_add, X, Y);
        TEST(inplace_sub, X, Y);
        TEST(inplace_scaling, X, 4.0);
    }

    blas_int M = 1024;
    N = 2048;
    cout << endl << "Testing matrices of size " << M << "x" << N << endl;

    {
        cout << endl;
        dMatrix A(M,N);
        dMatrix B(M,N);

        for (blas_int i=0; i<A.rows(); ++i) {
            for (blas_int j=0; j<A.columns(); ++j){
                A(i,j) = i * 0.03 - j * 0.01;
                B(i,j) = i * 0.07 - j * 0.09;
            }
        }

        TESTS::SECTION = "Double Matrix";
        TEST(copy, A);
        TEST(fill, A, 0.0);
        TEST(swap, A, B);
        TEST(inplace_add, A, B);
        TEST(inplace_sub, A, B);
        TEST(inplace_scaling, A, 4.0);
    }

    {
        cout << endl;
        zMatrix A(M,N);
        zMatrix B(M,N);

        for (blas_int i=0; i<A.rows(); ++i) {
            for (blas_int j=0; j<A.columns(); ++j){
                A(i,j) = i * 0.03 - j * 0.01;
                B(i,j) = i * 0.07 - j * 0.09;
            }
        }

        TESTS::SECTION = "Complex Double Matrix";
        TEST(copy, A);
        TEST(fill, A, 0.0);
        TEST(swap, A, B);
        TEST(inplace_add, A, B);
        TEST(inplace_sub, A, B);
        TEST(inplace_scaling, A, 4.0);
    }

    //M = 1024;
    //N = 2048;

    M = 1024;
    N = 2048;
    cout << endl << "Testing row compressed matrices of size " << M << "x" << N << endl;

    blas_int NNZ = 16 * M;

    {
        cout << endl;

        zRCMatrix A(M,N,NNZ);
        blas_int pos=0;
        for (blas_int i=0; i<M; ++i) {
            A.row_index(i) = pos;
            for (blas_int n=0; n<16; ++n) {
                A.nonzeros(pos) = 1.0 * (i%50) / (n+1) + 1.0 * (i/8);
                A.columns(pos) = (i/8) * 16 + n;
                pos++;
            }
        }
        if (pos != NNZ)
            cout << "ERR" << endl;

        A.row_index(M) = pos;
     //
        zVector x(N);
        for (blas_int i=0; i<N; ++i) {
            x[i] = 1.0;
        }
     //
        zVector y(M);
        for (blas_int i=0; i<M; ++i) {
            y[i] = 0.5;
        }
     //
        dcomp alpha = {1.5,0.0}, beta = {0.0, 1.0};
        TEST(general_matrix_vector, alpha, A, x, beta, y);
    }

    cout << endl << "TOTAL " << TEST_COUNT << " TESTS PASSED, " << ERROR_COUNT << " ERRORS DETECTED" << endl << endl;

    return 0;
}
