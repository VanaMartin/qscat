#include <iostream>             // Input/Output library
#include <cassert>
#include <complex>              // Complex algebra

#include "common.h"     // Common functions, parameters and classes
//#include "bessel.h"     // Solutions to the coulomb potential Schrödinger equation (sph. Bessel functions)

//#include "arrays.h"
//#include "input.h"
//#include "fem_dvr_ecs.h"        // Library to be tested
#include "interface.h"

using namespace QSCAT;

int main(){

// Testing input procedure
    parameters::multi_grid<def_float> gp("input/N2/grids.txt");
// Testing the FEM_DVR_ECS_grid

    int nq = 8;
    int tnel = 20;
    iVector nel(3);
    nel[0]=0; nel[1]=20; nel[2]=0;
    dVector aa(21);
    double theta = 35;
    aa[0] = 0.0;
    for (int i=1; i<21; ++i){
        aa[i] = 1.0;
    }

    femGrid g(nq, tnel, nel, aa, theta);

    femGrid gx(g);
    femGrid gy = gx;

    gVector x(g);

    gVector y(x);
    gVector z = x;

    zOperatorF H(g);
    H = 0.0;
    zOperatorF G = H;
    G += H;
    H.add_kinetic_term(10.0);

    for (int i=0; i<g.nb(); ++i) {
        x.f( 1.0 * i, i);
        y.f( 1.0, i);
    }

    H += x;

    z = H * y;

}
