#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <complex>

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


/// TODO remove or move elsewhere
void fill(gVector& v, const pjvalue& p)
{
    assert(v.init());
  //
    if ( !p.isMember("function") )
        throw std::invalid_argument( "Invalid function json node provided!" );
    if (p["function"].asString() == "GaussianWavePacket") {
        const femGrid g = v.get_grid();
        def_float x0 = p["position"].asDouble();
        def_float s = p["sigma"].asDouble();
        def_float p0 = p["impulse"].asDouble();
        for (int i=0; i<v.get_size(); ++i)
           v.f( functions::Gaussian<def_comp>(g.x(i),x0,s,p0), i );
    } else {
        throw std::invalid_argument( "Provided function type not recognized!" );
    }
}

void TestFunction2d::clean()
{}

bool TestFunction2d::save_bin_body(std::ofstream& file) const
{
    assert(init_);
  //
    bool stat = file.is_open();
    if (stat && init_) {
        file.write((char*) &axis_, sizeof(char));
        file.write((char*) &channel_, sizeof(int));
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
        if (stat) stat = body_.save_binary(file);
        if (stat) stat = wave_packet_.save_binary(file);
        if (stat) stat = bound_state_.save_binary(file);
        if (stat) stat = fourier_coefficients_.save_binary(file);
        if (stat) stat = energies_.save_binary(file);
    }
    return stat;
}

bool TestFunction2d::read_bin_body(std::ifstream& file)
{
    bool stat = file.is_open();
    if (stat) {
        file.read((char*) &axis_, sizeof(char));
        file.read((char*) &channel_, sizeof(int));
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
        if (stat) stat = body_.read_binary(file);
        if (stat) stat = body_.get_grid() == grid_;
        if (axis_ == 'x'){
            if (stat) stat = wave_packet_.read_binary(file);
            if (stat) stat = wave_packet_.get_grid() == grid_.get_xgrid();
            if (stat) stat = bound_state_.read_binary(file);
            if (stat) stat = bound_state_.get_grid() == grid_.get_ygrid();
        } else {
            if (stat) stat = wave_packet_.read_binary(file);
            if (stat) stat = wave_packet_.get_grid() == grid_.get_ygrid();
            if (stat) stat = bound_state_.read_binary(file);
            if (stat) stat = bound_state_.get_grid() == grid_.get_xgrid();
        }
        if (stat) fourier_coefficients_.read_binary(file);
        if (stat) stat = energies_.read_binary(file);
        init_ = true;
    }
    return stat;
}

TestFunction2d::TestFunction2d()
{
    init_ = false;
    opened_ = false;
}

TestFunction2d::TestFunction2d( const pjvalue& params, const femGrid2D& grid,
                                int channel, zEigenSystem& Eig,
                                def_comp initial_energy, def_comp charge,
                                def_float mass, int impulse_momentum,
                                const dVector& energy)
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
    grid_ = grid;
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

// Test function builder
    if (axis_ == 'x') {                                  // Wave packet in electronic (x) coordinate, bound state in the nuclear (y) coordinate
        wave_packet_ = gVector(grid.get_xgrid());
        fill(wave_packet_, params["wavepacket"]);

        bound_state_ = gVector(grid.get_ygrid());
        Eig.eigen_vector(bound_state_.body(),channel);
        body_ = gVector2D(grid, wave_packet_, bound_state_);
    } else if (axis_ == 'y'){                            // Wave packet in nuclear (y) coordinate, bound state in the electronic (x) coordinate
        wave_packet_ = gVector(grid.get_ygrid());
        fill(wave_packet_, params["wavepacket"]);

        bound_state_ = gVector(grid.get_xgrid());
        Eig.eigen_vector(bound_state_.body(),channel);
        body_ = gVector2D(grid, bound_state_, wave_packet_);
    }

    energy_ =  Eig.eigen_value(channel);
    energy_shift_ = real(energy_ - initial_energy_); // Channel energy shift

    cout << "Norm of the TestFunction2d: \t" << body_ * body_ << endl;
    cout << "Energy of the state: " << energy_  << ", Threshold: " << energy_shift_ << endl;

    int esize = energies_.get_size();
    fourier_coefficients_ = zVector(esize);
// Computing the fourier coefficients and storing them into the vector variable
    def_comp impulse;
    if (axis_ == 'x'){
        gVector aux(grid.get_xgrid());
        for (int i=0; i<esize; ++i){ // energy distribution
            if (energies_[i] - energy_shift_ > 0){ // Only the values above the threshold are assumed
                impulse = 1.0 * sqrt(2*mu_x_*(energies_[i] - energy_shift_));
                for (int j=0; j<xsize_; ++j){
                    if (abs(charge_) == 0) { // no coulombic term
                        aux.f(sphHankel1En(grid.xz(j), impulse, mu_x_, impulse_momentum_)/2.0, j);
                    } else {
                        aux.f(coulomb::sH1_en(grid.xz(j), impulse, charge_, mu_x_, impulse_momentum_)/2.0, j);
                    }
                } // xsize_
                fourier_coefficients_[i] = aux * wave_packet_; // integral ( conj(radial wave) * wave packet ) dx
            } else {
                fourier_coefficients_[i] = 0;   // The zero value can be used to control the utilisation of the other values.
            }
        }
    } else { // y - axis
        gVector aux(grid.get_ygrid());
        for (int i=0; i<esize; ++i){ // energy distribution
            if (energies_[i] - energy_shift_ > 0) {
                impulse = 1.0 * sqrt(2*mu_y_*(energies_[i] - energy_shift_));
                def_float alpha = sqrt(mu_y_/(2*real(impulse)*pi));     // == sqrt(2 m / pi k) * 1/2
                for (int j=0; j<ysize_; ++j){
                    aux.f(alpha*std::exp(imu*impulse*grid.yz(j)),j);
                    //aux.F(bessel::s_h1En(g->Yz(j),k,mu_y,0)/2.0, j);
                }
                fourier_coefficients_[i] = aux * wave_packet_;
                } else {
                fourier_coefficients_[i] = 0;                   // The zero value can be used to control the utilisation of the other values.
            }
        }
    } // x/y axis
    init_ = true;
}

TestFunction2d::TestFunction2d(const TestFunction2d& old) :
    TestFunctionInterface2d(old),
    body_(old.body_),
    wave_packet_(old.wave_packet_)
{}

TestFunction2d& TestFunction2d::swap(TestFunction2d& rhs)
{
    TestFunctionInterface2d::swap(rhs);
    body_.swap(rhs.body_);
    wave_packet_.swap(rhs.wave_packet_);
    return *this;
}

TestFunction2d& TestFunction2d::operator=(TestFunction2d tmp)
{
    return this->swap(tmp);
}

TestFunction2d::~TestFunction2d()
{
}

void TestFunction2d::operator<< (const gVector2D& psi)
{
    assert(init_);
    assert(psi.init());
    assert(size_ == psi.get_size());
  //
    buffer_ <<  body_ * psi;
}

void TestFunction2d::contribution(zMatrix& S, int idx, def_float t, def_float dt, const zVector& ifc)
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

    assert( coefficients_.get_size() == buffer_.get_size() );   // if fields are not correct, we are in big

    def_comp Et;
    def_comp Q;
    def_comp R;
    def_comp val;

    if( this->opened_ ) {
        outfile_->write((char*) &buffer_[1], (buffer_.get_size()-1)*sizeof(def_comp));
        //for (int i=1; i<buffer_.get_size(); ++i)        // skip the first value (from last loop)
            //fprintf(this->outfile_, "%.12e\t%.12e\t%.12e\n", (t + 1.0*i*dt), real(buffer_[i]), imag(buffer_[i]));
        //fflush(this->outfile_);
    }

    for (int i=0; i<energies_.get_size(); ++i){  // energy distribution
        if (energies_[i] - energy_shift_ > 0) {
            Q =  1.0/(2.0*pi* conj(fourier_coefficients_[i])*ifc[i]);
            for (int j=0; j<buffer_.get_size(); ++j){
                Et = (energies_[i] + initial_energy_) * (t + 1.0*j*dt);     // This is the actual energy of the icoming and outcoming state (ECL)
                R = exp(imu*Et) *dt *coefficients_[j];
                S(i,idx) += Q*buffer_[j] * R;
            }
        }
    }
    def_comp aux = buffer_();
    buffer_.clear();
    buffer_ << aux;
}

void TestFunction2d::print_coefficients(const std::string& path)
{
    char name[50];
    sprintf_s(name, "TF%d_body.dat", channel_);
    body_.save((path + name).c_str());
    sprintf_s(name, "TF%d_FC.dat", channel_);
    fourier_coefficients_.save(energies_ , (path+name).c_str());
}
} // namespace QSCAT
