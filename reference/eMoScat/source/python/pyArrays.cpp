#include <iostream>		// Input/Output library
#include <complex>		// Complex algebra
#include <math.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <stdexcept>
#include <omp.h>		// Parallelization library

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL emoscat_ARRAY_API

//#include <numpy/numpyconfig.h>
#include <numpy/arrayobject.h>
//#include <numpy/npy_math.h>

//#include "common.h"		// Common functions, parameters and classes
//#include "bessel.h"		// Solutions to the coulomb potential Schrödinger equation (sph. Bessel functions)
//#include "blas.h"		// Intel BLAS, SPARSE-BLAS LAPACK and PARDISO wrapper

//#include "Arrays.h"
//#include "pjinput.h"
#include "potentials.h"
//#include "FemDvrEcs.h"
//#include "FemDvrEcs2d.h"
//#include "interface.h"
//#include "ModelLCP/ModelLCP.h"
//#include "module_NRM.h"
//#include "Model2d.h"

#include "qscat.h"

// GPU enhancement

//#include "libXcuda.hpp"
//#include "Xcuda.h"


#ifndef dComplex
#define dComplex std::complex<double>
#endif

using namespace boost::python;
using namespace QSCAT;

/*
	Before the module definition, the wrapper classes are defined.
*/

#include "python/Arrays.hpp"
#include "python/input.hpp"
#include "python/FemDvrEcs.hpp"
#include "python/FemDvrEcs2d.hpp"
#include "python/functions.hpp"

// libXcuda Wrapper

//#include "python/libXcuda.hpp"


/*
class pyGrid : public FEM_DVR_ECS::FemDvrEcs_grid<double, dComplex>, public wrapper<FEM_DVR_ECS::FemDvrEcs_grid<double, dComplex> > {
public:
	pyGrid() : FEM_DVR_ECS::FemDvrEcs_grid<double, dComplex>() {}
};

class pyGridVector : public FEM_DVR_ECS::grid_vector<double, dComplex>, public wrapper<FEM_DVR_ECS::grid_vector<double, dComplex> > {

};

class pyGrid2D : public FEM_DVR_ECS_2D::FemDvrEcs_grid_2D<double, dComplex>, public wrapper<FEM_DVR_ECS_2D::FemDvrEcs_grid_2D<double, dComplex> > {

};

class pyGridVector2D2D : public FEM_DVR_ECS_2D::grid_vector_2D<double, dComplex>, public wrapper<FEM_DVR_ECS_2D::grid_vector_2D<double, dComplex> > {
};
*/

/*
class pyModel2D : public model_2D::Model, public wrapper<model_2D::Model> {
public:
	//pyModel2D() : model_2D::Model<double, dComplex>() {}
//	Model(parameters::model_2D<T> & m2dp, parameters::grid<T> & gpx, parameters::grid<T> & gpy, grid_vector_2D<T,Z>* phid = NULL);
//	~Model();
//	Model & Set(parameters::model_2D<T> & p, parameters::grid<T> & gpx, parameters::grid<T> & gpy);
//	void DiscreteProjection();
//	void Multistep();
//	bool SaveBinary(const char * name);
//	bool SaveBinary(std::ofstream & file);
//	bool ReadBinary(const char * name);
//	bool ReadBinary(std::ifstream & file);
};
*/

BOOST_PYTHON_MODULE(pyArrays)
{
    _import_array();

	#include "python/pyArrays.hpp"
	#include "python/pyInput.hpp"
	#include "python/pyFemDvrEcs.hpp"
    #include "python/pyFemDvrEcs2d.hpp"

    #include "python/pyFunctions.hpp"

    def("bessel_sje", &sphBesselJEn );
    def("coulomb_sff", &coulomb::sF_en);

    // libXcuda wrapper

    //#include "python/pylibXcuda.hpp"
};



