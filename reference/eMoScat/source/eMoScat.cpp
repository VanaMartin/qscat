// Main project file of electron-Molecule-Scattering project
// The project is under construction
// It is expected to contain following modules
// 2D model
// LCP
// NRM

#include <iostream>				// Input/Output library
#include <complex>				// Complex algebra
#include <math.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <omp.h>				// Parallelization library

#ifdef linux
	#include <sec_stream.h>
#endif

#include "common.h"		// Common functions, parameters and classes
#include "bessel.h"		// Solutions to the coulomb potential Schrödinger equation (sph. Bessel functions)
#include "blas.h"		// Intel BLAS, SPARSE-BLAS LAPACK and PARDISO wrapper

#include "Arrays.h"
#include "input.h"
#include "potentials.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"
#include "interface.h"
#include "ModelLCP/ModelLCP.h"
//#include "module_NRM.h"
#include "Model2D.h"
#include "pjinput.h"

/*
	TODO LIST:

	-- Check input integrity and correct grid definitions (common problem ending with segfaults)
*/


using namespace std;
using namespace QSCAT;

int main(int argc, char **argv) {

	std::string input_folder;

	bool isLCP = true;
	bool is2D  = true;

	if (argc>1){
		input_folder += "input/";
		input_folder += argv[1];
		input_folder += "/";
	} else {
		input_folder += "input/N2/";
	}

	if (argc>2){
		std::string opt(argv[2]);
		if (opt == "LCP"){
			is2D = false;
		}
		else if (opt == "2D"){
			isLCP = false;
		}
	}

//// CUDA testing area
//	RCMatrix<comp> K(4, 4, 9);
//	// | 1.0     2.0 3.0 |
//	// |     4.0         |
//	// | 5.0     6.0 7.0 |
//	// |     8.0     9.0 |

//	K.wNZE(0) = 1.0;	K.wC(0) = 0;
//	K.wNZE(1) = 2.0;	K.wC(1) = 2;
//	K.wNZE(2) = 3.0;	K.wC(2) = 3;
//	K.wNZE(3) = 4.0;	K.wC(3) = 1;
//	K.wNZE(4) = 5.0;	K.wC(4) = 0;
//	K.wNZE(5) = 6.0;	K.wC(5) = 2;
//	K.wNZE(6) = 7.0;	K.wC(6) = 3;
//	K.wNZE(7) = 8.0;	K.wC(7) = 1;
//	K.wNZE(8) = 9.0;	K.wC(8) = 3;
//
//	K.wRI(0) = 0;
//	K.wRI(1) = 3;
//	K.wRI(2) = 4;
//	K.wRI(3) = 7;
//	K.wRI(4) = 9;

//	int q = 4;

//	comp* X = new comp[q];
//	comp* Y = new comp[q];
//	for (int i=0; i<q; ++i){
//		X[i] = pow(10.0, i-2); //*(1+i);
//		Y[i] = 0.0;
//	}

//	K.CUDA_set();

//	printf("Array Y before matrix application: (");
//	for (int i=0; i<q; ++i){
//		printf(" %f", real(X[i]));
//	}
//	printf(")\n");

//	K.CUDA_host_matmul( zone, X, zzero, Y);
//	K.CUDA_host_matmul( zone, X, zone, Y);

//	printf("Array Y after matrix application: (");
//	for (int i=0; i<q; ++i){
//		printf(" %f", real(Y[i]));
//	}
//	printf(")\n");

//	vector<comp> U(q);
//	U.fill(1.0);
//	vector<comp> V(q);
//	V.fill(0.0);

//	U.CUDA_set();
//	V.CUDA_set();

//	K.MatMulEx(zone,U,zzero,V);
//	V.CUDA_ret();

//	printf("Vector V after matrix application: (");
//	for (int i=0; i<q; ++i){
//		printf(" %f", real(V(i)));
//	}
//	printf(")\n");

//	std::cout << "Press any key to continue ... ";
//	std::cin.ignore();

//	exit(0);

// Parameters input
	parameters::multi_grid<def_float> gp( (input_folder + "grids.txt").c_str() );
	parameters::model_2D<def_float> m2dp( (input_folder + "2D_model.txt").c_str() );
	parameters::LCP<def_float> LCPp( (input_folder + "LCP.txt").c_str() );
	parameters::NRM<def_float> NRMp( (input_folder + "NRM.txt").c_str() );

// 1D declarations

//	FEM_DVR_ECS::FemDvrEcs_grid<def_float,def_comp> fem_x(gp.gp[0]);
//	FEM_DVR_ECS::FemDvrEcs_grid<def_float,def_comp> fem_y(gp.gp[1]);

//	def_float mu = 1.0;
//	FEM_DVR_ECS::Kinetic_Energy_RCM<def_float,def_comp> KE_x(fem_x,mu);
//	FEM_DVR_ECS::Kinetic_Energy_RCM<def_float,def_comp> KE_y(fem_y,mu);

//	arrays::vector<def_comp> fX(fem_x.NB());
//	for (int i=0; i<fem_x.NB();++i){
//		fX(i) = functions::Gaussian<def_float,def_comp>(fem_x.Xr(i),20.0, 2.0, 2.0);
//	}
//	arrays::vector<def_comp> fY(fem_y.NB());
//	for (int i=0; i<fem_y.NB();++i){
//		fY(i) = functions::Sine(fem_y.Xr(i),10.0,1);
//	}

// LCP Approximation declarations

	//LCP::Model_LCP<def_float, def_comp> LCP(m2dp,gp,LCPp);

	//for (int i=0; i<1000; ++i){
	//	LCP.Multistep();
	//}

	//FEM_DVR_ECS::FemDvrEcs_grid<def_float,def_comp> fem_x(gp.gp[0]);
	//FEM_DVR_ECS::FemDvrEcs_grid<def_float,def_comp> fem_y(gp.gp[1]);
    //FEM_DVR_ECS_2D::FemDvrEcs_grid_2D<def_float,def_comp> fem_2D(fem_x,fem_y);
	//FEM_DVR_ECS_2D::grid_vector_2D<def_float,def_comp> X(&fem_2D);

	//fem_x.save_binary("x.grid.bin");
	//fem_y.save_binary("y.grid.bin");

	//LCP.MakePhiD(X);
	//NRM::MakePhiD_const(X,m2dp,fem_2D);

	//X.Save_equidistant( (m2dp.folder + "LCP/PhiD.dat").c_str() , 200, 0.0, 30.0, 400, 0.0, 6.0);
	//X.Save((m2dp.folder + "LCP/PhiD.dat").c_str());

	//std::cout << "Press any key to continue ... ";
	//std::cin.ignore();

	//exit(0);

// NRM Declarations

//	NRM::ModelNRM NRM(NRMp,m2dp,gp);
//	NRM.TimeIndependentSolution();
//	exit(0);

// 2D declarations

	LCP::ModelLCP LCP(m2dp,gp,LCPp);

    gVector2D PhiD = LCP.MakePhiD();
	TimeDependentModel2D M2D(m2dp,gp.gp[0],gp.gp[1], &PhiD);

	if (M2D.ReadBinary((m2dp.folder + "frame.M2D").c_str())){
		std::cout << "The previous run was successfully loaded." << std::endl;
	}


//	for (int i=0; i<10; ++i){
//		NRM.multiStep();
//	}

//	exit(0);

	for (int i=0;i<m2dp.evol_par.tcutoff/(m2dp.evol_par.dt*m2dp.evol_par.loop);++i){
		if (is2D) M2D.multistep();
		if (isLCP) LCP.Multistep();
		M2D.save_binary((m2dp.folder + "frame.M2D").c_str());
	}

	//FEM_DVR_ECS_2D::FemDvrEcs_grid_2D<def_float,def_comp> fem_2D(fem_x,fem_y);
	//FEM_DVR_ECS_2D::Kinetic_Energy_2D_RCM<def_float,def_comp> KE_2D(&fem_2D,KE_x,KE_y);

	//FEM_DVR_ECS_2D::grid_vector_2D<def_float,def_comp> X(&fem_2D,fX,fY);

	//FEM_DVR_ECS::Kinetic_Energy_Full<def_float,def_comp> KEF_y(fem_y,mu);
	//eigen_system<def_comp> Eig = KEF_y.eigen_sys();
	//model_2D::testfunction<def_float,def_comp> Tx(m2dp.evol_par,m2dp.test_par[0],&fem_2D,0,Eig,0.0,m2dp.mu,m2dp.l);

	//KE_2D.LU_Factorize();



// 1D Testing

//	arrays::vector<def_comp> X1(fem_x.NB());
//	for (int i=0; i<fem_x.NB();++i){
//		X1(i) = functions::Gaussian<def_float,def_comp>(fem_x.Xr(i),20.0, 4.0, 1.0);
//		//X(i) = sqrt(def_float(2.0)/def_float(8.0))*sin(1*pi/def_float(8.0)*fem.Xz(i));
//	}

//	FEM_DVR_ECS::grid_vector<def_float,def_comp> Y1(fem_x,X1);
//	FEM_DVR_ECS::grid_vector<def_float,def_comp> Z1(fem_x,X1);
//	FEM_DVR_ECS::grid_vector<def_float,def_comp> W1(fem_x,X1);

//	Z1 = Y1;
//	FEM_DVR_ECS::grid_vector<def_float,def_comp> Q1(Y1);
//
//	for (int i=0;i<fem_x.NB();++i){
//		X1(i) = fem_x.Wz(i);
//		X1.save("wz.dat");
//	}

//	def_comp norm;
//	norm = Y1*Y1;

//	cout << "Norm of the vector: " << norm << endl;

//	//def_float mu;
//	mu = def_float(1);
//	Q1.fill(0.0);
//	for (int i=0;i<Q1.size();++i){
//		if (i>30*7-1 && i <40*7-1) {
//			Q1(i) = 0.0;
//		}
//	}

//	FEM_DVR_ECS::Kinetic_Energy_Full<def_float,def_comp> KE1(fem_x,mu);
//	FEM_DVR_ECS::Kinetic_Energy_RCM<def_float,def_comp> KE1S(fem_x,mu);
//	FEM_DVR_ECS::Kinetic_Energy_RCM<def_float,def_comp> KE1T(fem_x,mu);
//

//	arrays::eigen_system<def_comp> E1(fem_x.NB());
//	E1 = KE1.eigen_sys();
//
//
//
//	//for (int i=0;i<fem.NB();++i){
//	//	X(i) = E.energies[i];
//	//	X.save("erg.dat");
//	//}

//	// Testing Sparse systems

//	double t;
//	double x_0;
//	double sigma, s_0, p_0;
//	def_comp phi;
//	comp val;

//	int order = 3;
//	def_float dt = 0.01;
//	FEM_DVR_ECS::Chebyshev_Full<def_float,def_comp> Ham(fem_x,Q1,mu, dt, 200);
//	FEM_DVR_ECS::Chebyshev_RCM<def_float,def_comp> H(fem_x,Q1,mu,dt,200);
//	FEM_DVR_ECS::Crank_Nicolson_Full<def_float,def_comp> CN(fem_x,Q1,mu,order,dt,KE1);
//	FEM_DVR_ECS::Crank_Nicolson_RCM<def_float,def_comp> CNS(fem_x,Q1,mu,dt,order);

//	norm = W1*W1;
//	cout << norm << endl;

//	def_comp qt;
//	for (int i=0;i<10000;++i){
//		FEM_DVR_ECS::chebyshev_one_step<def_float,def_comp,grid_vector<def_float,def_comp>,Chebyshev_Full<def_float,def_comp>>(Y1,Ham,dt,30); // 0.1, 50
//		FEM_DVR_ECS::chebyshev_one_step<def_float,def_comp,grid_vector<def_float,def_comp>,Chebyshev_RCM<def_float,def_comp>>(Z1,H,dt,30); // 0.1, 50
//		CNS.One_Step(W1);
//		norm = W1*W1;
//		cout << norm << endl;
//		if (i%100 == 0) {
//			t = dt*(i+1);
//			s_0 = 4.0;
//			p_0 = 1.0;
//			x_0 = 20.0 + p_0*t;
//			sigma = s_0*s_0 + t*t/(4.0*s_0*s_0);
//			for (int l=0;l<Q1.size();++l){
//				phi = p_0*(fem_x.Xr(l)-x_0) + pow(fem_x.Xr(l)-x_0, 2)/(8.0*sigma*s_0*s_0)*t + p_0*p_0*t/2.0 - arg(zone/sqrt(1.0 + imu*t/(2*pow(sigma,2))));
//				val = 1.0/sqrt(sqrt(2*pi*sigma))*exp(-pow(fem_x.Xr(l)-x_0,2)/(4.0*sigma))*exp(imu*phi);
//				Q1.F(val,l);
//			}
//			norm = Q1*Q1;
//			cout << norm << endl;
//			Q1.save("output/ref.dat");
//			Y1.save("output/test.dat");
//			Z1.save("output/test2.dat");
//			W1.save("output/test3.dat");
//			std::cin.ignore();
//		}
//	}
//
//	Y1.save("output/test.dat");
//	Z1.save("output/test2.dat");

//	int nb = fem_x.NB();
//	comp c;

	//cout << "Press a key to continue..." << endl;
	//std::cin.ignore();
	return 0;

}
