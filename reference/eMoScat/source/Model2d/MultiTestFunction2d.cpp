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
#include "blas.h"
#include "Arrays.h"
#include "input.h"
#include "potentials.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"
#include "Model2d.h"

namespace QSCAT
{
using std::string;
using std::cout;
using std::endl;
using std::sqrt;
using std::abs;

// FIXME
bool MultiTestFunction2d::save_bin_body(std::ofstream& file) const
{
    assert(init_);
  //
    bool stat = file.is_open();
    if(init_ && stat){
        file.write((char*) &channels_, sizeof(int));
        file.write((char*) &initial_channel_, sizeof(int));
        file.write((char*) &axis_, sizeof(char));
        file.write((char*) &method_, sizeof(char));
        file.write((char*) &mu_, sizeof(char));
        if (stat) stat = energies_.save_binary(file);
        if (tp_["method"].asString() == "tannor" || tp_["method"].asString() == "all" ){
            for(int i=0; i<channels_; ++i){
                if (stat) stat = testfunctions_[i].save_binary(file);
            }
            if (stat) stat = tannor_s_matrix_.save_binary(file);
        }
        if (tp_["method"].asString() == "dirac" || tp_["method"].asString() == "all" ){
            for(int i=0; i<channels_; ++i){
                if (stat) stat = dirac_testfunctions_[i].save_binary(file);
            }
            if (stat) stat = dirac_s_matrix_.save_binary(file);
        }
        if (tp_["method"].asString() == "flux" || tp_["method"].asString() == "all" ){
            for(int i=0; i<channels_; ++i){
                if (stat) stat = flux_testfunctions_[i].save_binary(file);
            }
            if (stat) stat = flux_s_matrix_.save_binary(file);
        }
        if (stat) stat = fourier_coefficients_.save_binary(file);
        int nf = folder_.size();
        char * name = new char[nf+1];
        strcpy_s(name, nf+1, folder_.c_str());
        if (stat) file.write((char*) &nf, sizeof(int));
        if (stat) file.write((char*) name, nf*sizeof(char));
        delete[] name;
    }
    return stat;
}

// FIXME
bool MultiTestFunction2d::read_bin_body(std::ifstream & file)
{
    bool stat = file.is_open();
    if(stat){
        if (init_) {
            if (testfunctions_) {
                delete[] testfunctions_;
            }
            if (dirac_testfunctions_) {
                delete[] dirac_testfunctions_;
            }
            if (flux_testfunctions_) {
                delete[] flux_testfunctions_;
            }
        }
        file.read((char*) &channels_, sizeof(int));
        file.read((char*) &initial_channel_, sizeof(int));
        file.read((char*) &axis_, sizeof(char));
        file.read((char*) &method_, sizeof(char));
        file.read((char*) &mu_, sizeof(char));

        if (stat) stat = energies_.read_binary(file);
        if (tp_["method"].asString() == "tannor" || tp_["method"].asString() == "all" ){
            testfunctions_ = new TestFunction2d[channels_];
            for(int i=0; i<channels_; ++i){
                if (stat) stat = testfunctions_[i].read_binary(file);
            }
            if (stat) stat = tannor_s_matrix_.read_binary(file);
        }
        if (tp_["method"].asString() == "dirac" || tp_["method"].asString() == "all" ){
            dirac_testfunctions_ = new DiracTestFunction2d[channels_];
            for(int i=0; i<channels_; ++i){
                if (stat) stat = dirac_testfunctions_[i].read_binary(file);
            }
            if (stat) stat = dirac_s_matrix_.read_binary(file);
        }
        if (tp_["method"].asString() == "flux" || tp_["method"].asString() == "all" ){
            flux_testfunctions_ = new FluxTestFunction2d[channels_];
            for(int i=0; i<channels_; ++i){
                if (stat) stat = flux_testfunctions_[i].read_binary(file);
            }
            if (stat) stat = flux_s_matrix_.read_binary(file);
        }
        if (stat) stat = fourier_coefficients_.read_binary(file);
        int nf;
        if (stat) file.read((char*) &nf, sizeof(int));
        char * name = new char[nf];
        if (stat) file.read((char*) name, nf*sizeof(char));
        folder_.assign(name,nf);
        delete[] name;
        init_ = true;
    }
    return stat;
}

MultiTestFunction2d::MultiTestFunction2d()
    : p_(pjvalue()), tp_(pjvalue())
{
    init_ = false;
    testfunctions_ = 0;
    dirac_testfunctions_ = 0;
    flux_testfunctions_ = 0;
}


MultiTestFunction2d::MultiTestFunction2d(femGrid2D *g, const pjvalue& params, const pjvalue& tfparams, const def_comp& init_erg) :
    p_(params), tp_(tfparams)
{
    folder_ = p_["model"]["folder"].asString();
    channels_ = tp_["channels"].asInt();        // Setting the channel number

    def_comp Q = 0;
    if (p_["model"].isMember("charge"))
        Q = p_["model"]["charge"].asDouble();

    if (tp_["coordinate"].asString() == p_["initial_state"]["coordinate"].asString()){
        initial_channel_ = p_["initial_state"]["channel"].asInt();
    } else {
        initial_channel_ = -1;
    }
    gVector aux;                              // Auxiliary grid_vector pointer
    femGrid grid;                             // Auxiliary grid pointer
// Determination of the transversal Hamiltonian (and its eigenstates)
    if (tp_["coordinate"].asString() == "electronic") {
        mu_ = 1.0;
        grid = femGrid(g->get_ygrid());
        aux = gVector(grid);
        fill_grid_vector(aux, p_["model"]["potential"], MorsePotential);
    } else if (tp_["coordinate"].asString() == "nuclear") {
        mu_ = p_["model"]["reduced_mass"].asDouble();
        grid = femGrid(g->get_xgrid());
        aux = gVector(grid);
        def_float y = g->yr(g->get_ysize()-1);
        fill_grid_vector_xaxis(aux, y, p_["model"]["potential"], Neutral2dPotential);
    }
    zOperatorF H = buildFullHamiltonian(grid, aux, (tp_["coordinate"].asString()=="electronic")? p_["model"]["reduced_mass"].asDouble() : 1.0);
    zEigenSystem eSys = H.eigen_system();

// Declarations of the energy range
    def_float ea = p_["cross_sections"]["range"][0].asDouble(), eb = p_["cross_sections"]["range"][1].asDouble();
    def_float de = p_["cross_sections"]["dE"].asDouble();
    energies_ = dVector(int((eb-ea)/de), ea, eb, false);
// Preparing auxiliary initial state (for initial state Fourier coefitients determination)
    fourier_coefficients_ = zVector(energies_.get_size());
    def_comp k;
    if (p_["initial_state"]["coordinate"].asString() == "electronic") {
        grid = femGrid(g->get_xgrid());
        aux = gVector(grid);
        fill_grid_vector(aux, p_["initial_state"]["wavepacket"], Gaussian);
        gVector bess(grid);
        for (int i=0;i<energies_.get_size();++i){
            k = 1.0 * sqrt(2.0 * energies_[i]);
            for (int j=0;j<grid.nb();++j){
                bess.f(sphBesselJEn(grid.x(j),k,1.0,p_["model"]["potential"]["impulsemomentum"].asInt()), j);
            }
            fourier_coefficients_[i] = bess*aux;
        }
    } else if (p_["initial_state"]["coordinate"].asString() == "nuclear"){
        grid = femGrid(g->get_ygrid());
        aux = gVector(grid);
        fill_grid_vector(aux, p_["initial_state"]["wavepacket"], Gaussian);
        gVector bess(grid);
        def_float mu = p_["model"]["reduced_mass"].asDouble();
        for (int i=0; i<energies_.get_size(); ++i){
            k = 1.0 * sqrt(2 * mu * energies_[i]);

            for (int j=0;j<grid.nb();++j){
                bess.f(sphBesselJEn(grid.x(j),k,mu,0), j);
            }
            fourier_coefficients_[i] = bess * aux;
        }
    }

    if (tp_["method"].asString() == "tannor" || tp_["method"].asString() == "all" ){
        testfunctions_ = new TestFunction2d[channels_];                // Allocating the test function structure for all channels
        char num[30];
        for (int i=0; i<channels_; ++i){
            cout << "Initializing Tannor & Weeks test function in coordinate: " << tp_["coordinate"].asString() << ", channel: " << i << ",..." << endl;
            testfunctions_[i] = TestFunction2d(tp_, *g, i, eSys, init_erg, Q, p_["model"]["reduced_mass"].asDouble(), p_["model"]["potential"]["impulsemomentum"].asInt(), energies_);

            sprintf(num, "/TW/cf_%s%d.txt", (tp_["coordinate"].asString() == "electronic")? "VE" : "DA" ,i);
            testfunctions_[i].set_output(folder_ + num);
            testfunctions_[i].print_coefficients(folder_ + "/TW/");
            cout << "...done." << endl;
        }
        tannor_cross_sections_ = dMatrix(energies_.get_size(), channels_);
        tannor_s_matrix_ = zMatrix(energies_.get_size(), channels_);
        tannor_s_matrix_.fill(0.0);
    }
    if (tp_["method"].asString() == "dirac" || tp_["method"].asString() == "all" ){
        dirac_testfunctions_ = new DiracTestFunction2d[channels_];          // Allocating the test function structure for all channels
        char num[30];
        for (int i=0; i<channels_; ++i){
            cout << "Initializing T&W with Dirac Delta test function in coordinate: " << axis_ << ", channel: " << i << ",..." << endl;
            dirac_testfunctions_[i] = DiracTestFunction2d(tp_, *g, i, eSys, init_erg, Q, p_["model"]["reduced_mass"].asDouble(), p_["model"]["potential"]["impulsemomentum"].asInt(), energies_);
            //sprintf(num, "/TWD/cf_%d.txt",i);
            sprintf(num, "/TWD/cf_%s%d.txt", (tp_["coordinate"].asString() == "electronic")? "VE" : "DA" ,i);
            dirac_testfunctions_[i].set_output(folder_ + num);
            dirac_testfunctions_[i].print_coefficients(folder_ + "/TWD/");
            cout << "...done." << endl;
        }
        dirac_cross_sections_ = dMatrix(energies_.get_size(), channels_);
        dirac_s_matrix_ = zMatrix(energies_.get_size(), channels_);
        dirac_s_matrix_.fill(0.0);
    }
    if (tp_["method"].asString() == "flux" || tp_["method"].asString() == "all" ){
        flux_testfunctions_ = new FluxTestFunction2d[channels_];           // Allocating the test function structure for all channels
        char num[30];
        for (int i=0; i<channels_; ++i){
            cout << "Initializing probability flux test function in coordinate: " << axis_ << ", channel: " << i << ",..." << endl;
            // FIXME
            //flux_testfunctions_[i].set(p.evol_par, p.test_par[AX], g, i, eSys, init_erg, zQ, p.mu, p.l, energies_);
            flux_testfunctions_[i] = FluxTestFunction2d(tp_, *g, i, eSys, init_erg, Q, p_["model"]["reduced_mass"].asDouble(), p_["model"]["potential"]["impulsemomentum"].asInt(), energies_);
            //sprintf(num, "/TF/cf_%d.txt",i);
            sprintf(num, "/TF/cf_%s%d.txt", (tp_["coordinate"].asString() == "electronic")? "VE" : "DA" ,i);
            flux_testfunctions_[i].set_output(folder_ + num);
            flux_testfunctions_[i].print_coefficients(folder_ + "/TF/");
            cout << "...done." << std::endl;
        }
        flux_cross_sections_ = dMatrix(energies_.get_size(), channels_);
        flux_s_matrix_ = zMatrix(energies_.get_size(), channels_);
        flux_s_matrix_.fill(0.0);
    }

    //s_matrix_.set(p.evol_par.steps, p.evol_par.loop, p.evol_par.Q, channels_, method_);
    init_ = true;
    fourier_coefficients_.save(energies_, "output/STATS/IFC.dat");
}

MultiTestFunction2d::~MultiTestFunction2d()
{
    if (init_) {
        if (testfunctions_){
            delete[] testfunctions_;
        }
        if (dirac_testfunctions_) {
            delete[] dirac_testfunctions_;
        }
        if (flux_testfunctions_){
            delete[] flux_testfunctions_;
        }
    }
}

void MultiTestFunction2d::step_buffer(gVector2D& Psi, const int& step)   // Storing multistep correlation functions or surface integrals
{
    if (tp_["method"].asString() == "tannor" || tp_["method"].asString() == "all" ){
        for (int i=0; i<channels_; ++i){
            testfunctions_[i] << Psi;
        }
    }
    if (tp_["method"].asString() == "dirac" || tp_["method"].asString() == "all" ){
        for (int i=0; i<channels_; ++i){
            dirac_testfunctions_[i] << Psi;
        }
    }
    if (tp_["method"].asString() == "flux" || tp_["method"].asString() == "all" ){
        for (int i=0; i<channels_; ++i){
            flux_testfunctions_[i] << Psi;
        }
    }
}

void MultiTestFunction2d::close_multistep(const def_float& t, const def_float& dt)   // Derives the S-matrix elements and stores them into appropriate vectors and files
{
    if (tp_["method"].asString() == "tannor" || tp_["method"].asString() == "all" ){
        for (int i=0; i<channels_; ++i){
            testfunctions_[i].contribution(tannor_s_matrix_, i, t, dt, fourier_coefficients_);
        }
        tannor_s_matrix_.save(energies_, (folder_ + "TW/" + ((tp_["coordinate"].asString()=="electronic")? "S_VE" : "S_DA") + ".dat").c_str());
    }
    if (tp_["method"].asString() == "dirac" || tp_["method"].asString() == "all" ){
        for (int i=0; i<channels_; ++i){
            dirac_testfunctions_[i].contribution(dirac_s_matrix_, i, t, dt, fourier_coefficients_);
        }
        dirac_s_matrix_.save(energies_, (folder_ + "TWD/" + ((tp_["coordinate"].asString()=="electronic")? "S_VE" : "S_DA") + ".dat").c_str());
    }
    if (tp_["method"].asString() == "flux" || tp_["method"].asString() == "all" ){
        for (int i=0; i<channels_; ++i){
            flux_testfunctions_[i].contribution(flux_s_matrix_, i, t, dt, fourier_coefficients_);
        }
        flux_s_matrix_.save(energies_, (folder_ + "TF/" + ((tp_["coordinate"].asString()=="electronic")? "S_VE" : "S_DA") + ".dat").c_str());
    }
}

void MultiTestFunction2d::cross_sections(const def_float& time)        // Derives the cross sections and stores them into a file with appropriate time
{
    std::string chname = (tp_["coordinate"].asString()=="electronic")? "VE" : "DA";

    char tname[10];
    def_float val;
    def_float erg;
    def_float mass = (axis_=='x')? mu_ : 1.0;

    sprintf_s(tname, "%.0f", time);
    // compute Tannor cross sections
    if (tp_["method"].asString() == "tannor" || tp_["method"].asString() == "all" ){
        for (int i=0; i<tannor_s_matrix_.rows(); ++i){
            for (int j=0; j<channels_; ++j) {
                if (j == initial_channel_) {
                    tannor_cross_sections_(i,j) = pow(abs(tannor_s_matrix_(i,j)-1.0),2)*pi/(2.0*mass*energies_[i]);
                } else {
                    tannor_cross_sections_(i,j) = pow(abs(tannor_s_matrix_(i,j)),2)*pi/(2.0*mass*energies_[i]);
                }
            }
        }
        tannor_cross_sections_.save(energies_,  (folder_ + "TW/CS" + chname + "_t_" + tname + ".dat").c_str());
        tannor_cross_sections_.save(energies_,  (folder_ + "TW/CS" + chname + ".dat").c_str());
    }

    // compute Dirac cross sections
    if (tp_["method"].asString() == "dirac" || tp_["method"].asString() == "all" ){
        for (int i=0; i<dirac_s_matrix_.rows(); ++i){
            for (int j=0; j<channels_; ++j) {
                if (j == initial_channel_) {
                    dirac_cross_sections_(i,j) = pow(abs(dirac_s_matrix_(i,j)-1.0),2)*pi/(2.0*mass*energies_[i]);
                } else {
                    dirac_cross_sections_(i,j) = pow(abs(dirac_s_matrix_(i,j)),2)*pi/(2.0*mass*energies_[i]);
                }
            }
        }
        dirac_cross_sections_.save(energies_,  (folder_ + "TWD/CS" + chname + "_t_" + tname + ".dat").c_str());
        dirac_cross_sections_.save(energies_,  (folder_ + "TWD/CS" + chname + ".dat").c_str());
    }

    // compute Flux cross sections
    if (tp_["method"].asString() == "flux" || tp_["method"].asString() == "all" ){
        for (int i=0; i<flux_s_matrix_.rows(); ++i){
            for (int j=0; j<channels_; ++j) {
                if (j == initial_channel_) {
                    flux_cross_sections_(i,j) = pow(abs(flux_s_matrix_(i,j)-1.0),2)*pi/(2.0*mass*energies_[i]);
                } else {
                    flux_cross_sections_(i,j) = pow(abs(flux_s_matrix_(i,j)),2)*pi/(2.0*mass*energies_[i]);
                }
            }
        }
        flux_cross_sections_.save(energies_,  (folder_ + "TF/CS" + chname + "_t_" + tname + ".dat").c_str());
        flux_cross_sections_.save(energies_,  (folder_ + "TF/CS" + chname + ".dat").c_str());
    }
}

} // namespace QSCAT
