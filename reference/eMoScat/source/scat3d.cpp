
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

#include "arrays.h"		
#include "input.h"
#include "potentials.h"
#include "fem_dvr_ecs.h"
#include "fem_dvr_ecs_2D.h"
#include "fem_dvr_ecs_3D.h"
#include "interface.h"
#include "module_LCP.h"
#include "module_NRM.h"
#include "module_2D.h"

using namespace std;

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
 
	parameters::multi_grid<def_float> gp( (input_folder + "grids.txt").c_str() );
	//parameters::multi_grid<def_float> gp( (input_folder + "grids.txt").c_str() );
	parameters::model_2D<def_float> m2dp( (input_folder + "2D_model.txt").c_str() );

    
	FEM_DVR_ECS::fem_dvr_ecs_grid<def_float,def_comp> gx(gp.gp[0]);

    FEM_DVR_ECS_3D::fem_dvr_ecs_grid_3D<def_float, def_comp> g(gx, gx, gx);
    FEM_DVR_ECS_3D::Operator3D<def_float, def_comp> H(g);
    H.addKineticTerm(1.0);

    cin.ignore();
}
