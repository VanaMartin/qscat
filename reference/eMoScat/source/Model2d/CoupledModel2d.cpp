#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <complex>
#include <string>

#ifdef linux
	#include <sec_stream.h>
#endif

#include "bessel.h"
#include "common.h"
#include "interface.h"
#include "potentials.h"
#include "Model2d.h"



def_comp zero(const def_comp& x, void* str) { return 0.0; }

def_comp QSCAT::p_pi_potential(const def_comp& x, const def_comp& y)
{
    return exp( - x*x/3.0 ) / x * 0.48 * sqrt(y) * ( exp( - pow( y - 3.35, 2) / 6.5) );
}

def_comp QSCAT::p_sigma_potential(const def_comp& x, const def_comp& y)
{
    using std::pow;
    using std::exp;
    using std::tanh;

    def_float a1, a2, a3, a4;
    a1 =  1.6435;
    a2 = 6.2;
    a3 = 0.0125;
    a4 = 1.15;

    def_comp Q = (a2 - y - a3 * pow(y, 4.0))/7.0;

    def_comp R = pow( tanh(y / a4), 4);

    def_comp E = exp( - pow(x, 2.0) / 3.0 ) / x;

    return - a1 * (1.0 - tanh(Q)) * R * E;
    //return  0.0/pow(x,2); //(tanh(y-6.2)+1.0) / 2.0 *
}

// Coupled model
/*
MODEL_2D::CoupledModel2D::CoupledModel2D()
{

    // Initialize simple "Half strip" parameters
    parametersMultiGrid gp( "input/test/grd_half_strip.txt" );

    // Manually set rest
    parametersEvolution ep;
    ep.steps = 1000; ep.loop = 100; ep.Q = 2; ep.pade = 5; ep.e_min = 0; ep.e_max = 4;

    parametersTestfunction tp;
    tp.axis = 'x'; tp.method = 'a';  tp.position = 15.0; tp.sigma = 1.0; tp.impulse = 2.0; tp.channels = 1; tp.io = 1;

    parametersInitState ip;
    ip.axis = 'x'; ip.position = 10.0; ip.sigma = 1.0; ip.impulse = -2.0; ip.channel = 0;

    // Set discretization grids
    gx = femGrid(gp.gp[0]);
    gy = femGrid(gp.gp[1]);
    g = femGrid2D(gx, gy);

    // Set evolution parameters
    time = 0;
    dt = 0.1;
    loop = ep.loop;

    // Set dynamics variables
    gVector2D pot1(g), pot2(g), cpl(g);


    for (int i=0; i<gy.nb(); ++i) {
        for (int j=0; j<gx.nb(); ++j) {
            if (gx.xr(i) < 1.0) {
                pot1.f(-2.0, i, j);
                pot2.f(0.0, i, j);
                cpl.f(-0.5, i, j);
            } else if (gx.xr(i) == 1.0) {
                pot1.f(-1.0, i, j);
                pot2.f(1.0, i, j);
                cpl.f(-0.25, i, j);
            } else {
                pot1.f(0.0, i, j);
                pot2.f(2.0, i, j);
                cpl.f(0.0, i, j);
            }
        }
    }

    H = doubleOperator2D(g);
    H.add_kinetic_term(1.0, 1.0);
    H.add_potential(pot1, pot2);
    H.add_coupling(cpl,cpl);

    CN.set(5, dt, H);

    // Set initial state
    psi = doubleGVector2D(g);
    gVector2D hPsi(g);

    for (int i=0; i<gy.nb(); ++i) {
        for (int j=0; j<gx.nb(); ++j) {
            hPsi.f( zGaussian(gx.xr(i), ip.position,  ip.sigma, ip.impulse) * dSine(gy.xr(j), 1.0, 1), i, j);
        }
    }
    def_float init_erg = std::pow(pi/1.0, 2) / 2;

    for (int i=0; i<g.get_size(); ++i) {
        psi[i] = hPsi[i];
        psi[i+g.get_size()] = 0.0;
    }

    dVector erg(ep.steps, ep.e_min, ep.e_max, false);

    std::cout << "Initial energy: " << init_erg << std::endl;

    // Set Scattering variables
    tf = new MultiTestfunction[2];   // Two channels
    tf[0].set(  g, 0.0, 1.0, 1.0, 0,
                ip, ep, init_erg, 0,
                tp, zero, NULL,
                erg, "output/tests/T1/");

    tf[1].set(  g, 0.0, 1.0, 1.0, 0,
                ip, ep, init_erg, -1,
                tp, zero, NULL,
                erg, "output/tests/T2/");

    tf[0].step_buffer(hPsi,loop);

//            multi_testfunction& set(    femGrid2D& Grid, const def_float& MuX, const def_float& MuY, const int& l,
//                                        const parametersInitState& iPar, const parametersEvolution& ePar, const def_comp& init_erg, const int& init_ch,
//                                        const parametersTestfunction& Par, const Z (*pot)(const Z&, void*), void* p_struct,
//                                        const dVector& erg_vector, const std::string& Folder);

    folder = "output/tests";

    for (int j=0; j<100; ++j){
        // Test evolution

        for (int i=0; i<ep.loop; ++i){
            CN.one_step(psi);

            std::cout << psi*psi << std::endl;
            sgVector2D aux1(g, &psi[0]);
            sgVector2D aux2(g, &psi[g.get_size()]);
            char tname[10];
            sprintf(tname, "%0.1f", time + i*dt);
            //aux1.Save( ( folder + "/test1_" + tname + ".dat" ).c_str());
            //aux2.Save( ( folder + "/test2_" + tname + ".dat" ).c_str());
            aux1.save( ( folder + "/test1.dat" ).c_str());
            aux2.save( ( folder + "/test2.dat" ).c_str());

            tf[0].step_buffer(aux1, i);
            tf[1].step_buffer(aux2, i);
        }

        // Test scattering
        tf[0].close_multistep(time, dt);
        tf[1].close_multistep(time, dt);

        time += dt*loop;

        tf[0].cross_sections(time);
        tf[1].cross_sections(time);
    }
}
*/
