#include <iostream>				// Input/Output library
#include <cassert>
#include <complex>				// Complex algebra

#include "cuda.h"		// nVidia CUDA BLAS and SPARSE-BLAS routines

#include "common.h"		// Common functions, parameters and classes
#include "bessel.h"		// Solutions to the coulomb potential Schrödinger equation (sph. Bessel functions)

#include "interface.h"


using namespace QSCAT;

femGrid make_grid() 
{
	int nq = 8;
	int tnel = 20; 
	iVector nel(3);
	nel[0]=0; nel[1]=20; nel[2]=0;
	dVector aa(21);
	double theta = 35;
	aa[0] = -1.0;
	for (int i=1; i<21; ++i){
		aa[i] = 0.1;
	}

	return femGrid(nq, tnel, nel, aa, theta);
}


int doubleGridVector2dTest() 
{
//    cout << "Double grid vector test" << endl;
//    // Make grid
//    fGrid grid;
//    cout << "phase 1 ..." << endl;
//    try { grid = make_grid(); }
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 1; }
//    // Make 2D grid
//    cout << "phase 2 ..." << endl;
//    f2Grid g;
//    try { g = f2Grid(grid,grid); }     
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 2; }
//    // make 2D grid vector
//    cout << "phase 3 ..." << endl;
//    g2Vector psi;
//    try { psi = g2Vector(&g); }     
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 3; }
//    // make double 2D grid vector
//    cout << "phase 4 ..." << endl;
//    dg2Vector dpsi;
//    try { dpsi = dg2Vector(g,g); }     
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 4; }
//    // make double 2D operator
//    cout << "phase 5 ..." << endl;
//    dOperator2d dOp;
//    try { dOp = dOperator2d(g,g); }     
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 5; }
//    // Add kinetic term to the operator
//    cout << "phase 6 ..." << endl;
//    try { dOp.AddKineticTerm(1.0, 1.0); }     
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 6; }
//    // Check result size 
//    cout << "Phase 7 ..." << endl;
//    int err = 0; 
//    RCMatrix<comp> kex = FEM_DVR_ECS::generateKineticTermRCM(grid,1.0);
//    RCMatrix<comp> ke = FEM_DVR_ECS_2D::joinKineticTerms<double,comp>(kex,kex);
//    if (dOp.body.NNZ() != 2*ke.NNZ()) 
//        cout << "Error: The actual length of fields differ:" << dOp.body.NNZ() << ":" << 2*ke.NNZ() << endl;
//    // Check result values
//    cout << "Phase 8 ..." << endl;
//    for (int i=0; i<ke.NNZ(); ++i) {
//        if (dOp.body.NZE(i)!=ke.NZE(i)) ++err;
//    }
//    cout << "Phase 9 ..." << endl;
//    for (int i=ke.NNZ(); i<2*ke.NNZ(); ++i) {
//        if (dOp.body.NZE(i) != ke.NZE(i-ke.NNZ())) ++err;
//    }
//    if (err) { std::cout << err << " errors encountered in phase 9" << std::endl; } 
//    // Add potentials 
//    cout << "Phase 10 ..." << endl;
//    g2Vector pot1 = psi;
//    g2Vector pot2 = psi;
//    double sx0 = 0.25;
//    double sy0 = 0.1; 
//    double x;
//    double y;
//    for (int i=0; i<g.NbY(); ++i){
//        for (int j=0; j<g.NbX(); ++j){
//            x = g.Xr(j);
//            y = g.Yr(i);
//            //pot1.F(exp(-x*x/(2*sx0*sx0)) * exp(-y*y/(2*sy0*sy0)), g.NbX()*i + j);
//            //pot2.F(exp(-x*x/(2*sy0*sy0)) * exp(-y*y/(2*sx0*sx0)), g.NbX()*i + j);
//            //pot1.F( 16.*x*x +  4.*y*y, g.NbX()*i + j);
//            //pot2.F(  4.*x*x + 16.*y*y, g.NbX()*i + j);
//            pot1.F( 0.0, g.NbX()*i + j);
//            pot2.F( 0.0, g.NbX()*i + j);
//        }
//    }
//    try { dOp.AddPotential(pot1,pot2); }
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 10; }
//    // Check diagonal values
//    err = 0;
//    cout << "Phase 11 ..." << endl;
//    RCMatrix<comp>& M = dOp.body;
//    for (int i=0; i<M.M(); ++i){    // All rows
//        if (i<g.NB()){              // first potential
//            for (int j=M.RI(i); j<M.RI(i+1); ++j){
//                if (M.C(j)==i) {
//                    if (M.NZE(j) != ke.NZE(j) + pot1.F(i)) 
//                        err++; 
//                }
//            }
//        } else {                    // second potential
//            for (int j=M.RI(i); j<M.RI(i+1); ++j){
//                if (M.C(j)==i) {
//                    if (M.NZE(j) != ke.NZE(j-ke.NNZ()) + pot2.F(i-g.NB())) 
//                        err++; 
//                }
//            }
//        }
//    }
//    if (err) { std::cout << err << " errors encountered in phase 11" << std::endl; } 
//
//    cout << "Phase 12 ..." << endl;
//    for (int i=0; i<g.NbY(); ++i){
//        for (int j=0; j<g.NbX(); ++j){
//            x = g.Xr(j);
//            y = g.Yr(i);
//            //pot1.F(exp(-x*x/(2*sy0*sy0)) * exp(-y*y/(2*sy0*sy0)) * (-imu), g.NbX()*i + j);
//            //pot2.F(exp(-x*x/(2*sy0*sy0)) * exp(-y*y/(2*sy0*sy0)) * (+imu), g.NbX()*i + j);
//            pot1.F(0.0, g.NbX()*i + j);
//            pot2.F(0.0, g.NbX()*i + j);
//        }
//    }
//    try { dOp.AddCoupling(pot1,pot2); }
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 10; }
//    // Check size 
//    cout << "Phase 13 ..." << endl;
//    if (dOp.body.NNZ()!= 2*ke.NNZ() + 2*g.NB()) 
//        cout << "Error: Size of the result differs from expected value" << endl;
//    // Check positions
//    cout << "Phase 14 ..." << endl;
//    err = 0;
//    for (int i=0; i<g.NB(); ++i) {
//        if ( dOp.body.C( dOp.body.RI(i+1)-1 ) != g.NB()+i )
//            ++err;
//    }
//    for (int i=0; i<g.NB(); ++i) {
//        if ( dOp.body.C( dOp.body.RI(i+g.NB()) ) != i )
//            ++err;
//    }
//    if (err) { std::cout << err << " errors encountered in phase 14" << std::endl; } 
//    // Check values
//    cout << "Phase 15 ..." << endl;
//    err = 0;
//    for (int i=0; i<g.NB(); ++i) {
//        if ( dOp.body.NZE( dOp.body.RI(i+1)-1 ) != pot1.F(i) )
//            ++err;
//    }
//    for (int i=0; i<g.NB(); ++i) {
//        if ( dOp.body.NZE( dOp.body.RI(i+g.NB()) ) != pot2.F(i) )
//            ++err;
//    }
//    if (err) { std::cout << err << " errors encountered in phase 15" << std::endl; } 
//     
//    // Check values
//    cout << "Phase 15 ..." << endl;
//    double x0 = 0.5;
//    double y0 = 0.0;
//    double px = 0.25;
//    double py = 0.1;
//    for (int i=0; i<g.NbY(); ++i){
//        for (int j=0; j<g.NbX(); ++j){
//            x = g.Xr(j);
//            y = g.Yr(i);
//            comp val = functions::Gaussian<double,comp>(x, x0, 0.1, px) * functions::Gaussian<double,comp>(y, y0, 0.1, py);
//            dpsi.F( val, i*g.NbX() + j);
//        }
//    }
//    for (int i=0; i<g.NbY(); ++i){
//        for (int j=0; j<g.NbX(); ++j){
//            dpsi.F( 0.0, g.NB() + i*g.NbX() + j);
//        }
//    }
//
//    double dt = 1.0;
//    comp test;
//    try {
//        for (int s=0; s<1; ++s) {
//            comp alpha = 1.0;
//            test = 1.0;
//            dg2Vector aux1 = dpsi;
//            dg2Vector aux2 = dpsi;
//            for (int i=0; i<20; ++i){
//                alpha *= -imu * dt / double(i+1);
//                test += alpha;
//                dOp.Gemv( alpha, aux1, 0.0, aux2);
//                dpsi += aux2;
//                cout << aux2*aux2 << " " <<  dpsi*dpsi << " " << test*conj(test) << endl;
//                aux1.Swap(aux2);
//            }
//        }
//    }
//    catch (exception& e) { cout << "Exception: " << e.what() << endl; return 15; }
//
    return 0; 

}

int main(){
	
//// Testing input procedure
//	parameters::multi_grid<def_float> gp("input/N2/grids.txt");
//// Testing the FEM_DVR_ECS_grid 
//	
//	int nq = 8;
//	int tnel = 20; 
//	vector<int> nel(3);
//	nel[0]=0; nel[1]=20; nel[2]=0;
//	vector<double> aa(21);
//	double theta = 35;
//	aa[0] = 0.0;
//	for (int i=1; i<21; ++i){
//		aa[i] = 1.0;
//	}
//
//	FEM_DVR_ECS::fem_dvr_ecs_grid<double,std::complex<double> > grid_e(nq, tnel, nel, aa, theta);
//
//	for (int i=0; i<grid_e.NB(); ++i){
//		std::cout << grid_e.Xr(i) << ", "; 
//	}
//	std::cout << std::endl;
//
//	nq = 16;
//	tnel = 40; 
//	nel[0]=0; nel[1]=40; nel[2]=0;
//	aa.Set(41);
//	theta = 35;
//	aa[0] = 0.0;
//	for (int i=1; i<41; ++i){
//		aa[i] = 1.0;
//	}
//
//	FEM_DVR_ECS::fem_dvr_ecs_grid<double, comp > grid_n(nq, tnel, nel, aa, theta);
//
//	for (int i=0; i<grid_n.NB(); ++i){
//		std::cout << grid_n.Xr(i) << ", "; 
//	}
//	std::cout << std::endl;
//
//	FEM_DVR_ECS_2D::fem_dvr_ecs_grid_2D<double, comp> G2D(grid_e, grid_n);
    
//    int error = doubleGridVector2dTest();
//    if (error)
//        std::cout << "Error on double Grid vector test phase " << error << std::endl;

}
