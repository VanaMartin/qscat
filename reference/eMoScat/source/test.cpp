#include <iostream>				// Input/Output library
#include <complex>				// Complex algebra
#include <math.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <omp.h>				// Parallelization library

//#ifdef linux
//	#include <sec_stream.h>
//#endif

#include "common.h"		// Common functions, parameters and classes
#include "bessel.h"		// Solutions to the coulomb potential Schrödinger equation (sph. Bessel functions)
#include "coulomb.h"
#include "blas.h"		// Intel BLAS, SPARSE-BLAS LAPACK and PARDISO wrapper

#include "arrays.h"
#include "input.h"
#include "potentials.h"
#include "fem_dvr_ecs.h"
#include "fem_dvr_ecs_2D.h"
#include "interface.h"
#include "module_LCP.h"
//#include "module_NRM.h"
#include "model_2D.h"
#include "fem_dvr_ecs_3D.h"

#include "conjugate_gradients.hpp"

using std::cout;
using std::endl;
using std::exp;
using std::abs;
using std::sqrt;
using std::pow;
using std::arg;

typedef FEM_DVR_ECS_3D::FemDvrEcsGrid3D<def_float, def_comp> femGrid3D;
typedef FEM_DVR_ECS_3D::GridVector3D<def_float, def_comp> gVector3D;
typedef FEM_DVR_ECS_3D::Operator3D<def_float, def_comp> zOperator3D;

//typedef FEM_DVR_ECS_2D::PreconditionerOperator2D<def_float,def_comp> Preconditioner2D;

typedef DiagonalPreconditioner<def_comp, femGrid, zOperatorC, gVector> DP1D;
typedef DiagonalPreconditioner<def_comp, femGrid2D, zOperator2D, gVector2D> DP2D;
typedef DiagonalPreconditioner<def_comp, femGrid3D, zOperator3D, gVector3D> DP3D;

class P1D {
    femGrid g_;
    zOperatorC M_;
    zOperatorC iM_;
 public:
    P1D(femGrid& g, double mass) : g_(g), M_(g), iM_(g)
    {
        M_.set_kinetic_term(mass);
        iM_.set_kinetic_term(mass);
        iM_.LU_factorize();
    }
    void gemv(def_comp alpha, const gVector& x, def_comp beta, gVector& y)
    {
        M_.gemv(alpha, x, beta, y);
    }
    void igemv(def_comp alpha, const gVector& x, def_comp beta, gVector& y)
    {
        if (beta == 0.0) {
            y = x;
            iM_.LU_back_substitution(y);
            if (alpha != 1.0)
                y *= alpha;
        } else {
            gVector a(y);
            if (beta != 1.0)
                a *= beta;
            y = x;
            iM_.LU_back_substitution(y);
            if (alpha != 1.0)
                y *= alpha;
            y += a;
        }
    }
};

class P2D {
    femGrid2D g_;
    zMatrix *P_;
    zMatrix *iP_;
 public:
    P2D(const femGrid2D& g, const zOperator2D& H) : g_(g), P_(NULL), iP_(NULL)
    {
        int ny = g.get_ysize();
        int nx = g.get_xsize();
        P_ = new zMatrix[ny];
        iP_ = new zMatrix[ny]; 
        for (int i=0; i<ny; ++i) {
            P_[i] = zMatrix(nx, nx);
            P_[i].fill(0);
            for (int j=0; j<nx; ++j) {
                for (int k=H.body().row_index(i*nx + j); k<H.body().row_index(i*nx + j +1); ++k){
                    int l = H.body().columns(k) - i*nx;
                    if (l < nx && l >= 0) {
                        P_[i](j,l) = H.body().nonzeros(k);
                    }
                }           
            }
            iP_[i] = P_[i]; 
            iP_[i].inverse();
        }
    }
    void gemv(def_comp alpha, const gVector2D& x, def_comp beta, gVector2D& y)
    {
        blas_int shift;
        for (blas_int i=0; i<g_.get_ysize(); ++i) {
            shift = i * g_.get_xsize();
            const ShallowVector<def_comp> p1(g_.get_xsize(), const_cast<def_comp*>(&x[shift]) );
            ShallowVector<def_comp> p2(g_.get_xsize(), &y[shift]);
            P_[i].gemv(alpha, p1, beta, p2);
        }
    }
    void igemv(def_comp alpha, const gVector2D& x, def_comp beta, gVector2D& y)
    {
        blas_int shift;
        for (blas_int i=0; i<g_.get_ysize(); ++i){
            shift = i * g_.get_xsize();
            const ShallowVector<def_comp> p1(g_.get_xsize(), const_cast<def_comp*>(&x[shift]) );
            ShallowVector<def_comp> p2(g_.get_xsize(), &y[shift]);
            iP_[i].gemv(alpha, p1, beta, p2);
        }
    }
};

class P3D {
    femGrid3D g_;
    zMatrix *P_;
    zMatrix *iP_;
 public:
    P3D(const femGrid3D& g, const zOperator3D& H) : g_(g), P_(NULL), iP_(NULL)
    {
        blas_int nz = g.get_zsize();
//        int ny = g.get_ysize();
//        int nx = g.get_xsize();
//        P_ = new zMatrix[ny];
//        iP_ = new zMatrix[ny]; 
//        for (int i=0; i<ny; ++i) {
//            P_[i] = zMatrix(nx, nx);
//            P_[i].fill(0);
//            for (int j=0; j<nx; ++j) {
//                for (int k=H.body().row_index(i*nx + j); k<H.body().row_index(i*nx + j +1); ++k){
//                    int l = H.body().columns(k) - i*nx;
//                    if (l < nx && l >= 0) {
//                        P_[i](j,l) = H.body().nonzeros(k);
//                    }
//                }           
//            }
//            iP_[i] = P_[i]; 
//            iP_[i].inverse();
//        }
    }
    void gemv(def_comp alpha, const gVector2D& x, def_comp beta, gVector2D& y)
    {
//        blas_int shift;
//        for (blas_int i=0; i<g_.get_ysize(); ++i) {
//            shift = i * g_.get_xsize();
//            const ShallowVector<def_comp> p1(g_.get_xsize(), const_cast<def_comp*>(&x[shift]) );
//            ShallowVector<def_comp> p2(g_.get_xsize(), &y[shift]);
//            P_[i].gemv(alpha, p1, beta, p2);
//        }
    }
    void igemv(def_comp alpha, const gVector2D& x, def_comp beta, gVector2D& y)
    {
//        blas_int shift;
//        for (blas_int i=0; i<g_.get_ysize(); ++i){
//            shift = i * g_.get_xsize();
//            const ShallowVector<def_comp> p1(g_.get_xsize(), const_cast<def_comp*>(&x[shift]) );
//            ShallowVector<def_comp> p2(g_.get_xsize(), &y[shift]);
//            iP_[i].gemv(alpha, p1, beta, p2);
//        }
    }
};


int main(int argc, char** argv) {

    if (false) {
        int N = 2000;
        zMatrix A(N,N);
        zVector x(N);

        A.fill(1.0 + imu);
        double dh = 1.0 / (N-1);
        for (int i=0; i<N; ++i)
            A[N*i + i] = -2.0 + ( 2 * exp(- 2 * i*dh) -  exp( - i * dh ));

        //for (int i=0; i<N; ++i) {
        //    for (int j=0; j<N; ++j) {
        //        cout << A(i,j) << ", ";
        //    }
        //    cout << endl;
        //}

        x.fill(1.0);
        x[N/2] = 2.0;

        zVector y1 = compute_biconjugate_gradients<def_float,def_comp,zMatrix,zVector>(A, x, 1e-15);

        //zVector y2 = COCG<def_float, def_comp, zMatrix, zVector>(A, x, 1e-15);
        //for (int i=0; i<N; ++i)
        //    cout << y1[i] << "  " << y2[i] << endl;

        exit(0);
    }

 // energy problem parametrisation:
    parameters2D m2dp("input/N2/2D_model.txt");
    parametersEvolution& ep(m2dp.evol_par);
    dVector erg(ep.steps, ep.e_min, ep.e_max, false);

    def_float dt = ep.dt;
    int innerSpan = ep.loop;
    int CNorder = ep.pade;

 // grids:
    parametersMultiGrid gp("input/N2/grids.txt");
    femGrid gx(gp.gp[0]); // electronic
    femGrid gy(gp.gp[1]); // nuclear
    //femGrid gm(gp.gp[2]); // electronic extended real region for Moeller evolution


    def_comp e = 0.08;

    // 1D case

    if (false) {
        zOperatorC H(gy);
        H.set_kinetic_term(1.0);

        gVector pot(gy);

        for (int i=0; i<gy.nb(); ++i){
            //if (i<gy.nr())
                pot.f( exp(- 2.0 * ( gy.x(i) - 1.0 ) - 2.0 * exp( -(gy.x(i) - 1.0) ) ), i );
            //else
            //    pot.f(0.0, i);
        }

        H += pot;

        gVector psi0(gy);

        // v_int |psi0>
        for  (int i=0; i<gy.nb(); ++i) {
            //if (i<=gy.nr()-1)
                psi0.f( bessel::s_jEn(gy.x(i), sqrt(2*e), 1.0, 0.0) * pot.f(i), i);
            //else
            //    psi0.f(0.0,i);
        }

        H *= -1.0;

        H += e;


        //gVector res = compute_biconjugate_gradients<def_float,def_comp,zOperatorC,gVector>(H, psi0, 1e-15, psi0);

        //gVector res = COCG<def_float,def_comp,zOperatorC,gVector>(H, psi0, 1e-15);

        //P1D P(gy, 1.0);
        //gVector res = PCOCG<def_float,def_comp,zOperatorC,P1D, gVector>(H, P, psi0, 1e-15);

        DP1D PP(gy, H);
        gVector res = PCOCG<def_float,def_comp,zOperatorC,DP1D, gVector>(H, PP, psi0, 1e-15);

//        for (int s=0; s<100; s++) {
//            def_comp step = 0.00001;
//
//            e += step;
//            for  (int i=0; i<gy.nb(); ++i) {
//                if (i<=gy.nr())
//                    psi0.f( bessel::s_jEn(gy.x(i), sqrt(2*e), 1.0, 0.0) * pot.f(i), i);
//                else
//                    psi0.f(0.0,i);
//            }
//            H += step;
//
//            cout << "e = " << e << endl;
//            res = compute_biconjugate_gradients<def_float,def_comp,zOperatorC,gVector>(H, psi0, 1e-15, res);
//        }
        res.save("res.asc");

        zOperatorC H2(H);

        gVector psc(psi0);
        gVector psi(psi0);

        H2.LU_back_substitution(psc);

        psc.save("test.asc");
        //psi -= psc;

        //cout << psc * psi0 << endl;


        //cout << res * psi0 << endl;
        exit(0);
    }

    // 2D case
    femGrid2D g(gx, gy);

    if (false) {
      // potential
        gVector2D pot(g);
        gVector2D vint(g);
        for (int i=0; i<gy.nb(); ++i){
            for (int j=0; j<gx.nb(); ++j){
                pot.f(potentials::potential_2D<def_float,def_comp>(gx.x(j),gy.x(i),m2dp), i*gx.nb() + j);
                vint.f(potentials::V_int<def_float,def_comp>(gx.x(j),gy.x(i),m2dp), i*gx.nb() + j);
            }
        }
      // Operator
        zOperator2D H(g);
        H.set_kinetic_term(1.0, m2dp.mu);
        H += pot;

      // initial state
        // make 1D potential
        gVector aux1(gy);
        aux1.function_fill(potentials::V_zero<def_float,def_comp>,m2dp);
        // build operator and its eigensystem
        zOperatorF H1 = makeHamiltonianF(gy, aux1, m2dp.mu);
        zEigenSystem eSys = H1.eigen_system();
        def_comp ierg = real(eSys.eigen_value(0));
        // get eigenvector
        gVector Y(gy);
        eSys.eigen_vector(Y.body(), 0);
        gVector X(gx);
        for (int i=0; i<gx.nb(); ++i){
            //if (i<=gx.nr())
                X.f( bessel::s_jEn(gx.x(i), sqrt(2*e), 1.0, 1.0), i);
            //else
            //    X.f(0.0,i);
        }

      // H -> E - H
        H *= -1.0;
        H += e + ierg;
        cout << e + ierg << endl;
        // make 2D state
        gVector2D psi0(g,X,Y);
        psi0.element_wise_multiplication(vint);

        //Preconditioner2D P(g,1.0);
        //gVector2D res = compute_preconditioned_biconjugate_gradients<def_float,def_comp,zOperator2D,Preconditioner2D,gVector2D>(H,P,psi0, 1e-6);

        //gVector2D res = compute_biconjugate_gradients<def_float,def_comp,zOperator2D,gVector2D>(H,psi0, 1e-6);

        //gVector2D res = COCG<def_float,def_comp,zOperator2D,gVector2D>(H,psi0, 1e-15);
        //res.save("test2d_cocg.txt");

        //gVector2D res = PCOCG<def_float,def_comp,zOperator2D,Preconditioner2D,gVector2D>(H,P,psi0, 2e-10);

        DP2D P(g, H);
        gVector2D res = PCOCG<def_float,def_comp,zOperator2D,DP2D,gVector2D>(H,P,psi0, 1e-15);

        //P2D PP(g, pot, m2dp.mu, real(e+ierg));
        //P2D PP(g, H);
        //gVector2D res = PCOCG<def_float,def_comp,zOperator2D,P2D,gVector2D>(H,PP,psi0, 1e-15);

        res.save("test2d.txt");

        zOperator2D H2(H);

        res = psi0;
        H2.LU_back_substitution(res);
        res.save("ref2d.txt");

        exit(0);
    }

    // 3D case
    if (true) {
      // grid
        FEM_DVR_ECS_3D::FemDvrEcsGrid3D<def_float,def_comp> g3(gx, gy, gy);     // TODO FIX order of grids
      // potential
        FEM_DVR_ECS_3D::GridVector3D<def_float, def_comp> pot(g3);
        FEM_DVR_ECS_3D::GridVector3D<def_float, def_comp> vint(g3);
        FEM_DVR_ECS_3D::GridVector3D<def_float, def_comp> psi0(g3);
        
        // 1D potentials
        gVector xPot(gy);
        xPot.function_fill(potentials::V_zero<def_float,def_comp>,m2dp);
        // build operator and its eigensystem
        zOperatorF Hc = makeHamiltonianF(gy, xPot, m2dp.mu);
        zEigenSystem eSys = Hc.eigen_system();
        def_comp ierg = real(eSys.eigen_value(0));
        // get eigenvector
        gVector Y(gy);
        eSys.eigen_vector(Y.body(), 0);
        gVector X(gx);
        for (int i=0; i<gx.nb(); ++i) {
            if (i < gx.nr())
                X.f( bessel::s_jEn(gx.x(i), sqrt(2*e), 1.0, 1.0), i);
            else
                X.f(0,i);
        } 

        cout << "States: x=" << X*X << ", y=" << Y*Y << endl;

        def_comp r, v, vi;
        for (int i=0; i<gy.nb(); ++i){
            for (int j=0; j<gy.nb(); ++j){
                r = sqrt(pow( abs(gy.x(i)),2 ) + pow( abs(gy.x(j)),2) ) * exp(imu * (arg(gy.x(i))+arg(gy.x(j))));
                v = potentials::V_zero(gy.x(i), m2dp) + potentials::V_zero(gy.x(j), m2dp);
                for (int k=0; k<gx.nb(); ++k) {
                    vi = potentials::V_int<def_float,def_comp>(gx.x(k), gy.x(i), m2dp);
                    vint.f( vi, i, j, k );
                    //pot.f(potentials::potential_2D<def_float,def_comp>(gx.x(k), r, m2dp), (i*gy.nb() + j)*gy.nb() + k );
                    pot.f( vi + v, (i*gy.nb() + j)*gx.nb() + k);
                    psi0.f( Y.f(i) * Y.f(j) * X.f(k), i, j, k );
                }
            }
        }
        def_comp sum = 0;
        for (int i=0; i<gy.nb(); ++i){
            for (int j=0; j<gy.nb(); ++j){
                for (int k=0; k<gx.nb(); ++k) {
                    //sum += std::pow(abs(psi0(i,j,k)),2);
                    sum += std::pow(abs(psi0[(i*gy.nb() + j)*gx.nb() + k]),2);
                }
            }
        }
        cout << "done: " << psi0*psi0 << " : " << sum << endl;
        FEM_DVR_ECS_3D::Operator3D<def_float,def_comp> H(g3);
        
        cout << "Operator" << endl;

        H.set_kinetic_term(1.0, m2dp.mu, m2dp.mu);
        H += pot;
        H *= -1;
        H += e + 2*ierg;

        cout << "Solver" << endl;
        psi0.element_wise_multiplication(vint);


        DP3D P(g3, H);
        gVector3D res = PCOCG<def_float, def_comp, zOperator3D, DP3D, gVector3D>(H,P,psi0, 1e-9);
    }
}
