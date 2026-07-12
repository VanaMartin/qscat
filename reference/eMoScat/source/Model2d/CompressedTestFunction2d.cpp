#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cassert>
#include <complex>
#include <string>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "bessel.h"
#include "coulomb.h"
#include "common.h"
#include "blas.h"
#include "Arrays.h"
#include "input.h"
#include "potentials.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"
#include "Model2d.h"

namespace QSCAT
{
using std::abs;
using std::sqrt;
using std::cout;
using std::endl;

void CompressedTestfunction::initialize(const parametersEvolution& ep, const parametersTestfunction& p, const mGrid2D& grid, int ch, zEigenSystem& Eig)
{
    // TODO - add more assertions 
    assert(p.axis == 'x' || p.axis == 'y');
    assert(grid.init());
    assert(Eig.init());
    assert(ch >= 0);
  // 
    axis_ = p.axis;      
    in_out_ = p.io;
    channel_ = ch;
    const femGrid2D& g = grid.full_grid();
    size_ = g.get_size();
    xsize_ = g.get_xsize();
    ysize_ = g.get_ysize();
    char name[50];
// Test function builder
    if (axis_ == 'x') {                                  // Wave packet in electronic (x) coordinate, bound state in the nuclear (y) coordinate
        wave_packet_ = gVector(g.get_xgrid());
        for (int i=0; i<xsize_; ++i){
            wave_packet_.f(functions::Gaussian<def_comp>(g.xr(i),p.position,p.sigma,p.impulse),i); 
        }
        bound_state_ = gVector(g.get_ygrid());
        Eig.eigen_vector(bound_state_.body(),ch);
        body_ = mVector2D(grid, wave_packet_, bound_state_);
    } else if (axis_ == 'y'){                            // Wave packet in nuclear (y) coordinate, bound state in the electronic (x) coordinate
        wave_packet_ = gVector(g.get_ygrid());
        for (int i=0; i<ysize_; ++i){
            wave_packet_.f(functions::Gaussian<def_comp>(g.yr(i),p.position,p.sigma,p.impulse),i); 
        }
        bound_state_ = gVector(g.get_xgrid());
        Eig.eigen_vector(bound_state_.body(),ch);
        body_ = mVector2D(grid, bound_state_, wave_packet_);
    }
    
    energy_ =  Eig.eigen_value(ch);
    energy_shift_ = real(Eig.eigen_value(ch) - initial_energy_); // Channel energy shift  

    cout << "Norm of the Testfunction: \t" << body_ * body_ << endl;
    cout << "Energy of the state: " << energy_  << ", Threshold: " << energy_shift_ << endl;

    fourier_coefficients_ = zVector(ep.steps);
// Computing the fourier coefficients and storing them into the vector variable
    def_comp impulse;
    if (axis_ == 'x'){
        gVector aux(g.get_xgrid()); 
        for (int i=0; i<ep.steps; ++i){ // energy distribution
            if (energies_[i] - energy_shift_ > 0){ // Only the values above the threshold are assumed 
                impulse = 1.0 * sqrt(2*mu_x_*(energies_[i] - energy_shift_));
                for (int j=0; j<xsize_; ++j){
                    if (abs(charge_) == 0) { // no coulombic term
                        aux.f(bessel::s_h1En(g.xz(j), impulse, mu_x_, impulse_momentum_)/2.0, j);
                    } else {
                        aux.f(coulomb::sH1_en(g.xz(j), impulse, charge_, mu_x_, impulse_momentum_)/2.0, j); 
                    }
                } // xsize_
                fourier_coefficients_[i] = aux * wave_packet_; // integral ( conj(radial wave) * wave packet ) dx
            } else {
                fourier_coefficients_[i] = 0;   // The zero value can be used to control the utilisation of the other values. 
            }
        }
    } else { // y - axis
        gVector aux(g.get_ygrid());
        for (int i=0; i<ep.steps; ++i){ // energy distribution
            if (energies_[i] - energy_shift_ > 0) {
                impulse = 1.0 * sqrt(2*mu_y_*(energies_[i] - energy_shift_));
                def_float alpha = sqrt(mu_y_/(2*real(impulse)*pi));
                for (int j=0; j<ysize_; ++j){
                    aux.f(alpha*std::exp(imu*impulse*g.yz(j)),j);
                    //aux.F(bessel::s_h1En(g->Yz(j),k,mu_y,0)/2.0, j);
                }
                fourier_coefficients_[i] = aux * wave_packet_;
                } else {
                fourier_coefficients_[i] = 0;                   // The zero value can be used to control the utilisation of the other values. 
            }
        }
    } // x/y axis
    init_ = true;
    
    //sprintf_s(name, "output/STATS/TW%d_body.dat", channel_);
    //body_.save(name);
    //sprintf_s(name, "output/STATS/TWFC%d.dat", channel_);
    //fourier_coefficients_.save(energies_ , name); 
}
void CompressedTestfunction::clean()
{}
CompressedTestfunction::CompressedTestfunction()
{
    init_ = false;
}
CompressedTestfunction::CompressedTestfunction(const parametersEvolution& ep, const parametersTestfunction& p, const mGrid2D& grid, int ch, zEigenSystem& Eig, const def_comp& IE, const def_comp zQ, const def_float& M, const int L, const dVector& En)
{
    // TODO parameters init check
    assert(ch >= 0);
    assert(Eig.init());
    assert(M != 0);
    assert(L >= 0);
    assert(grid.init());
    assert(En.init());
  // 
    energies_ = En;
    mu_x_ = 1.0;
    mu_y_ = M;
    impulse_momentum_ = L;
    charge_ = zQ;
    initial_energy_ = IE;
    initialize(ep,p,grid,ch,Eig);
}
CompressedTestfunction::~CompressedTestfunction()
{
}
void CompressedTestfunction::set(const parametersEvolution& ep, const parametersTestfunction& p, const mGrid2D& grid, const int& ch, zEigenSystem& Eig, const def_comp& IE, const def_comp zQ, const def_float& M, const int L, const dVector& En)
{
    // TODO parameters init check
    assert(ch >= 0);
    assert(Eig.init());
    assert(M != 0);
    assert(L >= 0);
    assert(grid.init());
    assert(En.init());
  // 
    energies_ = En;
    initial_energy_ = IE;
    mu_x_ = 1.0;
    mu_y_ = M;
    impulse_momentum_ = L;
    charge_ = zQ;
    initialize(ep,p,grid,ch,Eig);
}
} // namespcae QSCAT

