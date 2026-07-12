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

void FluxTestFunction2d::clean()
{}

bool FluxTestFunction2d::save_bin_body(std::ofstream & file) const
{
    assert(init_);
  //
    bool stat = file.is_open();
    if (stat && init_) {
        file.write((char*) &axis_, sizeof(char));
        file.write((char*) &channel_, sizeof(int));
        file.write((char*) &element_index_, sizeof(int));
        file.write((char*) &index_start_, sizeof(int));
        file.write((char*) &index_end_, sizeof(int));
        file.write((char*) &position_, sizeof(int));
        file.write((char*) &quadrature_, sizeof(int));
        file.write((char*) &in_out_, sizeof(int));
        file.write((char*) &size_, sizeof(int));
        file.write((char*) &xsize_, sizeof(int));
        file.write((char*) &ysize_, sizeof(int));
        file.write((char*) &mu_x_, sizeof(def_float));
        file.write((char*) &mu_y_, sizeof(def_float));
        file.write((char*) &impulse_momentum_, sizeof(int));
        file.write((char*) &energy_, sizeof(def_comp));
        file.write((char*) &initial_energy_, sizeof(def_comp));
        file.write((char*) &energy_shift_, sizeof(def_float));
        stat = grid_.save_binary(file);
        if (stat) stat = bound_state_.save_binary(file);
        if (stat) stat = projection_.save_binary(file);
        if (stat) stat = phi_out_.save_binary(file);
        if (stat) stat = dphi_out_.save_binary(file);
    }
    return stat;
}

bool FluxTestFunction2d::read_bin_body(std::ifstream & file)
{
    bool stat = file.is_open();
    if (stat) {
        file.read((char*) &axis_, sizeof(char));
        file.read((char*) &channel_, sizeof(int));
        file.read((char*) &element_index_, sizeof(int));
        file.read((char*) &index_start_, sizeof(int));
        file.read((char*) &index_end_, sizeof(int));
        file.read((char*) &position_, sizeof(int));
        file.read((char*) &quadrature_, sizeof(int));
        file.read((char*) &in_out_, sizeof(int));
        file.read((char*) &size_, sizeof(int));
        file.read((char*) &xsize_, sizeof(int));
        file.read((char*) &ysize_, sizeof(int));
        file.read((char*) &mu_x_, sizeof(def_float));
        file.read((char*) &mu_y_, sizeof(def_float));
        file.read((char*) &impulse_momentum_, sizeof(int));
        file.read((char*) &energy_, sizeof(def_comp));
        file.read((char*) &initial_energy_, sizeof(def_comp));
        file.read((char*) &energy_shift_, sizeof(def_float));
        stat = grid_.read_binary(file);
        if (axis_ == 'x'){
            if (stat) stat = bound_state_.read_binary(file);
            if (stat) stat = bound_state_.get_grid() == grid_.get_ygrid();
            if (stat) stat = projection_.read_binary(file);
            if (stat) stat = projection_.get_grid() == grid_.get_xgrid();
        } else {
            if (stat) stat = bound_state_.read_binary(file);
            if (stat) stat = bound_state_.get_grid() == grid_.get_xgrid();
            if (stat) stat = projection_.read_binary(file);
            if (stat) stat = projection_.get_grid() == grid_.get_ygrid();
        }
        if (stat) stat = phi_out_.read_binary(file);
        if (stat) stat = dphi_out_.read_binary(file);
        init_ = true;
    }
    return stat;
}

FluxTestFunction2d::FluxTestFunction2d()
{
    init_ = false;
}


FluxTestFunction2d::FluxTestFunction2d( const pjvalue& params, const femGrid2D& grid, int channel, zEigenSystem& Eig, def_comp initial_energy, def_comp charge, def_float mass, int impulse_momentum, const dVector& energy)
{
    assert(params.isMember("coordinate"));
    assert(params.isMember("direction"));
    assert(params.isMember("wavepacket"));

    assert(params["coordinate"].asString() == "electronic" || params["coordinate"].asString() == "nuclear");
    assert(params["direction"].asString() == "in" || params["direction"].asString() == "out");

    // TODO - add more assertions
    assert(channel >= 0);
    assert(Eig.init());
    assert(mass != 0);
    assert(impulse_momentum >= 0);
    assert(grid.init());
    assert(energy.init());
  //
    grid_ = grid_;
    opened_ = false;
    energies_ = energy;
    mu_x_ = 1.0;
    mu_y_ = mass;
    impulse_momentum_ = impulse_momentum;
    charge_ = charge;
    initial_energy_ = initial_energy;

  //
    if (params["coordinate"].asString() == "electronic") {
        axis_ = 'x';
    } else {
        axis_ = 'y';
    }
    if (params["direction"].asString() == "in") {
        in_out_ = -1;
    } else {
        in_out_ = 1;
    }
    channel_ = channel;
    size_ = grid.get_size();
    xsize_ = grid.get_xsize();
    ysize_ = grid.get_ysize();
    char name[50];
    const femGrid* pg;

// Test function builder
    if (axis_ == 'x') {                                     // Wave packet in electronic (x) coordinate, bound state in the nuclear (y) coordinate
        pg = &grid.get_xgrid();
        element_index_ = pg->get_element_index(params["wavepacket"]["position"].asDouble()); // The surface is shifted  either to the starting or ending element border
        index_start_ = pg->get_element_start(element_index_);
        index_end_ = pg->get_element_end(element_index_);
        position_ = (element_index_ != pg->tnel())? index_end_ : index_start_;

        bound_state_ = gVector(grid.get_ygrid());
        Eig.eigen_vector(bound_state_.body(),channel_);
        projection_ = gVector(grid.get_xgrid());
    } else if (axis_ == 'y'){                               // Wave packet in nuclear (y) coordinate, bound state in the electronic (x) coordinate
        pg = &grid.get_ygrid();
        element_index_ = pg->get_element_index(params["wavepacket"]["position"].asDouble()); // The surface is shifted  either to the starting or ending element border
        index_start_ = pg->get_element_start(element_index_); // The surface is shifted to the starting element border (due to the danger of complex scaling)
        index_end_ = pg->get_element_end(element_index_);
        position_ = (element_index_ != pg->tnel())? index_end_ : index_start_;
        bound_state_ = gVector(grid.get_xgrid());
        Eig.eigen_vector(bound_state_.body(),channel_);
        projection_ = gVector(grid.get_ygrid());
    }

    energy_ =  Eig.eigen_value(channel);
    energy_shift_ = real(energy_ - initial_energy_); // Channel energy shift

    cout << "Energy of the state: " << energy_  << ", Threshold: " << energy_shift_ << endl;

    int esize = energies_.get_size();
    fourier_coefficients_ = zVector(esize);
// Computing the fourier coefficients and storing them into the vector variable
    phi_out_ = zVector(esize);
    dphi_out_ = zVector(esize);
    quadrature_ = index_end_ - index_start_ + 1;
    def_comp impulse;

    if (axis_ == 'x'){
        gVector aux(*pg);
        for (int i=0; i<esize; ++i){
            if (energies_[i] - energy_shift_ > 0){
                impulse = sqrt(2.0 * mu_x_ * (energies_[i] - energy_shift_));
                if (abs(charge_) == 0.0) {
                    phi_out_[i] = sphHankel1En(grid.xz(position_), impulse, mu_x_, impulse_momentum_)/2.0;
                    for (int j=index_start_; j<=index_end_; ++j){
                        aux.f(sphHankel1En(grid.xz(j), impulse, mu_x_, impulse_momentum_)/2.0, j);
                    }
                    dphi_out_[i] = aux.derivative(position_);
                } else {
                    phi_out_[i] = coulomb::sH1_en(grid.xz(position_), impulse, charge_, mu_x_, impulse_momentum_)/2.0;
                    for (int j=index_start_; j<=index_end_; ++j) {
                        aux.f(coulomb::sH1_en(grid.xz(j), impulse, charge_, mu_x_, impulse_momentum_)/2.0, j);
                    }
                    dphi_out_[i] = aux.derivative(position_);
                }
            } else {
                phi_out_[i] = 0;
                dphi_out_[i] = 0;
            }
        }
    } else {
        gVector aux(*pg);
        for (int i=0; i<esize; ++i){
            if (energies_[i] - energy_shift_ > 0) {
                impulse = sqrt(2 * mu_y_ * (energies_[i] - energy_shift_));
                phi_out_[i] = sphHankel1En(grid.yz(position_), impulse, mu_y_, 0)/2.0;
                for (int j=index_start_; j<=index_end_; ++j) {
                    aux.f(sphHankel1En(grid.yz(j), impulse, mu_y_, 0)/2.0, j);
                }
                dphi_out_[i] = aux.derivative(position_);
            } else {
                phi_out_[i] = 0;
                dphi_out_[i] = 0;
            }
        }
    }
    init_ = true;
    sprintf_s(name,"output/STATS/TFC_%d.dat", channel_);
    phi_out_.save(energies_, name);
    sprintf_s(name,"output/STATS/TFDC_%d.dat", channel_);
    dphi_out_.save(energies_, name);

    //buffer_ << 0.0;
    //initialize(ep,p,g,ch,Eig);
}

FluxTestFunction2d::FluxTestFunction2d(const FluxTestFunction2d& old) :
    TestFunctionInterface2d(old),

    element_index_(old.element_index_),
    index_start_(old.index_start_),
    index_end_(old.index_end_),
    position_(old.position_),
    quadrature_(old.quadrature_),

    projection_(old.projection_),
    phi_out_(old.phi_out_),
    dphi_out_(old.dphi_out_)
    //derivative_buffer_(old.derivative_buffer_)
{}

FluxTestFunction2d& FluxTestFunction2d::swap(FluxTestFunction2d& rhs)
{
    TestFunctionInterface2d::swap(rhs);
    std::swap(element_index_,rhs.element_index_);
    std::swap(index_start_,rhs.index_start_);
    std::swap(index_end_,rhs.index_end_);
    std::swap(position_,rhs.position_);
    std::swap(quadrature_,rhs.quadrature_);

    projection_.swap(rhs.projection_);
    phi_out_.swap(rhs.phi_out_);
    dphi_out_.swap(rhs.dphi_out_);
    //derivative_buffer_.swap(rhs.derivative_buffer_);
    return *this;
}

FluxTestFunction2d& FluxTestFunction2d::operator=(FluxTestFunction2d tmp)
{
    return this->swap(tmp);
}


FluxTestFunction2d::~FluxTestFunction2d()
{}

void FluxTestFunction2d::operator<< (const gVector2D& psi)
{
    assert(init_);
    assert(psi.init());
  //
    if (axis_ == 'x'){
        buffer_ << psi.line_projection(bound_state_,'y',position_);   // NOTE: The axis is perpendicular to channel axis
        for (int i=index_start_; i<=index_end_; ++i){
            projection_.f( psi.line_projection(bound_state_,'y',i), i);
        }
        derivative_buffer_ << projection_.derivative(position_);
    } else if (axis_ == 'y'){
        buffer_ << psi.line_projection(bound_state_,'x',position_);   // NOTE: The axis is perpendicular to channel axis
        for (int i=index_start_; i<=index_end_; ++i){
            projection_.f( psi.line_projection(bound_state_,'x',i), i);
        }
        derivative_buffer_ << projection_.derivative(position_);
    }

}

void FluxTestFunction2d::contribution(zMatrix& S, int idx, def_float t, def_float dt, const zVector& ifc)
{
    assert( init_ );
    assert( buffer_.get_size() > 0 );
    assert( ifc.init() );
    assert( S.init() );
  //
    // Starting coefficient build
    if ( coefficients_.get_size() != buffer_.get_size() ) {  // not compatible sizes
        // get best quad order:
        quad_order_=0;
        int base = 1;
        while( (buffer_.get_size()-1) % (base) == 0 ) { quad_order_++; base*=2; }
        coefficients_ = equidistant_quadrature(quad_order_, buffer_.get_size()-1); // assuming one value is always from the previous loop
    }

    assert( coefficients_.get_size() == buffer_.get_size() );   // if fiels are not correct, we are in big
  //
    def_comp Et;
    def_comp core;
    for (int i=0;i<energies_.get_size();++i){
        if (energies_[i] - energy_shift_ > 0){
            for (int j=0; j<buffer_.get_size(); ++j){
                Et = (energies_[i] + initial_energy_) * (t + 1.0*j*dt);
                core = (conj(phi_out_[i])*derivative_buffer_[j] - buffer_[j]*conj(dphi_out_[i])) * exp(imu*Et);
                S(i,idx) += (-imu/(2.0 * reduced_mass() * ifc[i])) * core * dt * coefficients_[j];  //(/2.0) *
            }

        }
    }
    def_comp aux1 = buffer_();
    def_comp aux2 = derivative_buffer_();
    buffer_.clear();
    derivative_buffer_.clear();
    buffer_ << aux1;
    derivative_buffer_ << aux2;
}

void FluxTestFunction2d::print_coefficients(const std::string& path)
{
    char name[50];
    sprintf_s(name, "Phi_%d(R).dat", channel_);
    phi_out_.save((path + name).c_str());
    sprintf_s(name, "dPhi_%d(R).dat", channel_);
    dphi_out_.save(energies_ , (path+name).c_str());
}

} // namespace QSCAT
