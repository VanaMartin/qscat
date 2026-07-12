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

using namespace std;

// FIXME remove
pjvalue dummyPar;
parameters2D dummyOldPar;

const bool HSV = false;
int vib_levels = 15;

//The main class for the two dimensional model of collisions
/*
void TimeDependentModel2D::initialize(parametersGrid& gpx, parametersGrid& gpy)
{
    xgrid_ = femGrid(gpx);
    ygrid_ = femGrid(gpy);
    grid_ = femGrid2D(xgrid_,ygrid_);
    folder_ = params_.folder;
  // Building the initial state
    if (params_.init_par.axis == 'x'){    // The incoming wave packet in the electronic coordinate
      // make potential
        gVector aux1(ygrid_);
        aux1.function_fill(potentials::V_zero<def_float,def_comp>,params_);
        aux1.save((folder_+"pot_x.dat").c_str());
      // build operator and its eigensystem
        zOperatorF H = buildFullHamiltonian(ygrid_, aux1, params_.mu);
        zEigenSystem eSys = H.eigen_system();
        initial_energy_ = eSys.eigen_value(params_.init_par.channel);
      // Save the energies
        dVector outErgs(ygrid_.nb());
        for (int k=0; k<ygrid_.nb(); ++k){
            outErgs[k] = real(eSys.eigen_value(k));
        }
        outErgs.save((folder_ + "V_0_ergs.dat").c_str());
      // Get resonant potential vibrations
        // Get v_res if possible
        gVector vRes(ygrid_);
        if (vRes.read_binary((folder_ + "LCP/bin/vres.bin").c_str()) && vRes.get_grid() == ygrid_){
            for (int i = ygrid_.nr(); i<ygrid_.nb(); ++i){
                vRes.f(vRes.f(ygrid_.nr()-1), i);
            }
            for (int i=0; i<ygrid_.nb(); ++i)
                vRes.f(real(vRes.f(i)),i);
            zOperatorF H2 = buildFullHamiltonian(ygrid_, vRes, params_.mu);
            zEigenSystem eSys2 = H2.eigen_system();
            for (int k=0; k<ygrid_.nb(); ++k){
                outErgs[k] = real(eSys2.eigen_value(k));
            }
            outErgs.save((folder_ + "V_res_ergs.dat").c_str());
            vibrations_ = new gVector[vib_levels];
            for (int k=0; k<vib_levels; ++k) {
                vibrations_[k] = gVector(ygrid_);
                eSys2.eigen_vector(vibrations_[k].body(),k);
            }
            population_file_.open((folder_ + "vibrational_populations.txt").c_str());
        }

      // get eigenvector
        gVector Y(ygrid_);
        eSys.eigen_vector(Y.body(), params_.init_par.channel);
      // build gaussian packet in free coordinate
        gVector X(xgrid_);
        for (int i=0; i<xgrid_.nb(); ++i){
            X.f(zGaussian(xgrid_.xr(i), params_.init_par.position, params_.init_par.sigma, params_.init_par.impulse), i);
        }
      // compose to two dimensional wave function
        psi_ = gVector2D(grid_,X,Y);
    } else if (params_.init_par.axis == 'y'){     // The incoming wave packet in the nuclear coordinate
      // make potential
        gVector aux1(xgrid_);
        aux1.function_xy_fill_x(potentials::V_int<def_float,def_comp>,ygrid_.xr(ygrid_.nb()),params_);
      // build operator and its eigensystem
        zOperatorF H = buildFullHamiltonian(xgrid_, aux1, 1.0);
        zEigenSystem eSys = H.eigen_system();
        initial_energy_ = eSys.eigen_value(params_.init_par.channel);
      // get eigenvector
        gVector X(xgrid_);
        eSys.eigen_vector(X.body(), params_.init_par.channel);
        gVector Y(ygrid_);
        for (int i=0; i<ygrid_.nb(); ++i){
            Y.f(zGaussian(ygrid_.xr(i), params_.init_par.position, params_.init_par.sigma, params_.init_par.impulse), i);
        }
      // compose to two dimensional wave function
        psi_ = gVector2D(grid_,X,Y);
    }
  // Store the initial wave function
    psi_.save((folder_+"Psi.dat").c_str());
    //psi_.save_equidistant((folder_+"Psi0Eq.dat").c_str(), params_.xSamples, params_.xRange[0], params_.xRange[1],
    //                                                     params_.ySamples, params_.yRange[0], params_.yRange[1]);
    std::cout << "Initial state norm:\t" << psi_.norm() << "\t saving into a file." << std::endl;

    initialize_potential();
    potential_.save((folder_+"Potential.dat").c_str());
    potentials::savePotential2D((folder_ + "Veff.dat").c_str(), potentials::potential_2D<def_float,def_comp>, params_, 400, 0.0, 15.0, 400, 0.0, 10.0);
    //potential->Save_equidistant((folder + "Veff.dat").c_str(), 40, 0.0, 2.0, 60, 0.0, 3);

  // Initialization of the test function
    initialize_testfunctions();

  // Set evolution operators
    chebyshev_use_ = false;
    crank_nicolson_use_ = false;

    if (params_.evol_par.evolution_c == 'c') {
        chebyshev_use_ = true;
        zOperator2D H(grid_);
        H.set_kinetic_term(1.0, params_.mu);
        H += potential_;
        chebyshev_ = new Chebyshev2D(params_.evol_par.cheb, params_.evol_par.dt, H);
    } else if (params_.evol_par.evolution_c == 'n') {
        crank_nicolson_use_ = true;
        zOperator2D H(grid_);
        H.set_kinetic_term(1.0, params_.mu);
        H += potential_;
        crank_nicolson_ = new CrankNicolson2D(params_.evol_par.pade, params_.evol_par.dt, H);
    }
    time_ = 0.0;
    dt_ = params_.evol_par.dt;
    loop_size_ = params_.evol_par.loop;

  // Initial contribution to S-matrix
    if (xtestfunction_use_){
        xtestfunction_->step_buffer(psi_,loop_size_);
    }
    if (ytestfunction_use_){
        ytestfunction_->step_buffer(psi_,loop_size_);
    }
    if (xtestfunction_use_){
        xtestfunction_->cross_sections(time_);
    }
    if (ytestfunction_use_){
        ytestfunction_->cross_sections(time_);
    }

    ep_ = new EquidistantProjector2d(grid_, params_.xSamples, params_.ySamples,
                                                     params_.xRange[0], params_.xRange[1],
                                                     params_.yRange[0], params_.yRange[1]);
    hsv_ = new EquidistantProjector2d(grid_, 1200, 600,
                                                      params_.xRange[0], params_.xRange[1],
                                                      params_.yRange[0], params_.yRange[1]);
    // equidistant projection
    *ep_ << psi_;
    ep_->export_state((folder_ + "WF/PsiEq_t_0.dat").c_str());

    // HSV projection
    *hsv_ << psi_;
    char name[50];
    sprintf(name, "bmp/%d.bmp", int(time_));
    if (HSV) hsv_->export_state_hsv((folder_ + name).c_str(), 5e-1, real(initial_energy_) * time_);

  // ready to go!
    init_ = true;
}
*/

void TimeDependentModel2D::initialize_testfunctions()
{
    if (params_.test_par[0].usage=='y'){
        xtestfunction_use_ = true;
        //xtestfunction_ = new MultiTestFunction2d(&grid_, params_, 0, initial_energy_);
    } else {
        xtestfunction_use_ = false;
    }
    if (params_.test_par[1].usage=='y'){
        ytestfunction_use_ = true;
        //ytestfunction_ = new MultiTestFunction2d(&grid_, params_, 1, initial_energy_);
    } else {
        ytestfunction_use_ = false;
    }
}

void TimeDependentModel2D::initialize_potential()
{
    int i, j;
    potential_ = gVector2D(grid_);
    for (i=0; i<ygrid_.nb(); ++i){
        for (j=0; j<xgrid_.nb(); ++j){
            potential_.f(potentials::potential_2D<def_float,def_comp>(xgrid_.xr(j),ygrid_.xr(i),params_), i*xgrid_.nb() + j);
        }
    }
}

void TimeDependentModel2D::clean()
{
    if (init_){
        if (vibrations_)
            delete[] vibrations_;
        if (xtestfunction_use_){
            delete xtestfunction_;
        }
        if (ytestfunction_use_){
            delete ytestfunction_;
        }
        if (chebyshev_use_){
            delete chebyshev_;
        }
        if (crank_nicolson_use_){
            delete crank_nicolson_;
        }
    }
    if (population_file_.is_open())
        population_file_.close();
}

bool TimeDependentModel2D::save_bin_body(std::ofstream& file) const
{
    bool stat = file.is_open();
    if (init_ && stat){
        stat = xgrid_.save_binary(file);
        if (stat) stat = ygrid_.save_binary(file);
        if (stat) stat = grid_.save_binary(file);
        if (stat) stat = potential_.save_binary(file);
        file.write((char*) &xtestfunction_use_, sizeof(bool));
        if (stat && xtestfunction_use_) {
            stat = xtestfunction_->save_binary(file);
        }
        file.write((char*) &ytestfunction_use_, sizeof(bool));
        if (stat && ytestfunction_use_) {
            stat = ytestfunction_->save_binary(file);
        }
        if (stat) stat = psi_.save_binary(file);
        file.write((char*) &initial_energy_, sizeof(def_comp));
        file.write((char*) &chebyshev_use_, sizeof(bool));
        file.write((char*) &crank_nicolson_use_, sizeof(bool));
        // CHEB and CN operators -- to be decided later
        file.write((char*) &time_, sizeof(def_float));
        file.write((char*) &dt_, sizeof(def_float));
        file.write((char*) &loop_size_, sizeof(int));
        int nf = folder_.size()+1;
        char * name = new char[nf];
        strcpy_s(name, nf, folder_.c_str());
        --nf;
        if (stat) file.write((char*) &(nf), sizeof(int));
        if (stat) file.write((char*) name, nf*sizeof(char));
        delete[] name;
    }
    return stat;
}

bool TimeDependentModel2D::read_bin_body(std::ifstream& file)
{
    bool stat = file.is_open();
    if (stat) {
        if (init_) {
            if (xtestfunction_use_){
                delete xtestfunction_;
            }
            if (ytestfunction_use_){
                delete ytestfunction_;
            }
        }
        stat = xgrid_.read_binary(file);
        if (stat) stat = ygrid_.read_binary(file);
        if (stat) stat = grid_.read_binary(file);
        if (stat) stat = potential_.read_binary(file);
        if (stat) stat = potential_.get_grid() == grid_;
        file.read((char*) &xtestfunction_use_, sizeof(bool));
        if (stat && xtestfunction_use_) {
            xtestfunction_ = new MultiTestFunction2d;
            stat = xtestfunction_->read_binary(file);
        }
        file.read((char*) &ytestfunction_use_, sizeof(bool));
        if (stat && ytestfunction_use_) {
            ytestfunction_ = new MultiTestFunction2d;
            stat = ytestfunction_->read_binary(file);
        }
        if (stat) stat = psi_.read_binary(file);
        if (stat) stat = psi_.get_grid() == grid_;
        file.read((char*) &initial_energy_, sizeof(def_comp));
        file.read((char*) &chebyshev_use_, sizeof(bool));
        file.read((char*) &crank_nicolson_use_, sizeof(bool));
        // CHEB and CN operators -- to be decided later
        file.read((char*) &time_, sizeof(def_float));
        file.read((char*) &dt_, sizeof(def_float));
        file.read((char*) &loop_size_, sizeof(int));
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
/*
TimeDependentModel2D::TimeDependentModel2D(parameters2D& m2dp, parametersGrid& gpx, parametersGrid& gpy, gVector2D *phid): params_(m2dp), par_(dummyPar)
{
    vibrations_ = 0;
    initialize(gpx, gpy);

    if (phid) discrete_state_ = *phid;
    if (discrete_state_.init()) {
        radius_operator_ = zOperatorD(ygrid_, ygrid_.nb(), 0);
        for (int i=0; i<ygrid_.nb(); ++i){
            radius_operator_[i] = ygrid_.xr(i);
        }
        discrete_projection();
        normalization_[0] = abs(psi_discrete_*psi_discrete_);
        radius_[0] = real(psi_discrete_*(radius_operator_*psi_discrete_))/normalization_();
    }
}
*/
/// Constructor
TimeDependentModel2D::TimeDependentModel2D(const pjvalue& parameters, gVector2D* phid) : par_(parameters), params_(dummyOldPar)
{
    vibrations_ = 0;

    xgrid_ = grid_from_parameters(par_["grids"][ par_["model"]["electronic_grid"].asString() ]);
    ygrid_ = grid_from_parameters(par_["grids"][ par_["model"]["nuclear_grid"].asString() ]);

    grid_ = femGrid2D(xgrid_,ygrid_);
    folder_ = par_["model"]["folder"].asString();
  // Building the initial state
    if (par_["initial_state"]["coordinate"].asString() == "electronic"){    // The incoming wave packet in the electronic coordinate
      // make potential
        gVector aux1(ygrid_);
        fill_grid_vector(aux1, par_["model"]["potential"], MorsePotential);
        aux1.save((folder_+"pot_x.dat").c_str());
      // build operator and its eigensystem
        zOperatorF H = buildFullHamiltonian(ygrid_, aux1, par_["model"]["reduced_mass"].asDouble());
        zEigenSystem eSys = H.eigen_system();
        initial_energy_ = eSys.eigen_value(par_["initial_state"]["channel"].asInt() );
      // Save the energies
        dVector outErgs(ygrid_.nb());
        for (int k=0; k<ygrid_.nb(); ++k){
            outErgs[k] = real(eSys.eigen_value(k));
        }
        outErgs.save((folder_ + "V_0_ergs.dat").c_str());
      // Get resonant potential vibrations
        // Get v_res if possible
        gVector vRes;
        if (vRes.read_binary((folder_ + "LCP/bin/vres.bin").c_str()) && vRes.get_grid() == ygrid_){
            for (int i = ygrid_.nr(); i<ygrid_.nb(); ++i){
                vRes.f(vRes.f(ygrid_.nr()-1), i);
            }
            for (int i=0; i<ygrid_.nb(); ++i)
                vRes.f(real(vRes.f(i)),i);
            zOperatorF H2 = buildFullHamiltonian(ygrid_, vRes, par_["model"]["reduced_mass"].asDouble());
            zEigenSystem eSys2 = H2.eigen_system();
            for (int k=0; k<ygrid_.nb(); ++k){
                outErgs[k] = real(eSys2.eigen_value(k));
            }
            outErgs.save((folder_ + "V_res_ergs.dat").c_str());
            vibrations_ = new gVector[vib_levels];
            for (int k=0; k<vib_levels; ++k) {
                vibrations_[k] = gVector(ygrid_);
                eSys2.eigen_vector(vibrations_[k].body(),k);
            }
            population_file_.open((folder_ + "vibrational_populations.txt").c_str());
        }

      // get eigenvector
        gVector Y(ygrid_);
        eSys.eigen_vector(Y.body(), par_["initial_state"]["channel"].asInt());
      // build gaussian packet in free coordinate
        gVector X(xgrid_);
        // TODO add more states
        fill_grid_vector(X, par_["initial_state"]["wavepacket"], Gaussian);
      // compose to two dimensional wave function
        psi_ = gVector2D(grid_,X,Y);
    } else if (par_["initial_state"]["coordinate"].asString() == "nuclear"){     // The incoming wave packet in the nuclear coordinate
      // make potential
        gVector aux1(xgrid_);
        fill_grid_vector_xaxis(aux1, ygrid_.xr(ygrid_.nb()), par_["model"]["potential"], LambdaInteraction);
      // build operator and its eigensystem
        zOperatorF H = buildFullHamiltonian(xgrid_, aux1, 1.0);
        zEigenSystem eSys = H.eigen_system();
        initial_energy_ = eSys.eigen_value(par_["initial_state"]["channel"].asInt());
      // get eigenvector
        gVector X(xgrid_);
        eSys.eigen_vector(X.body(), par_["initial_state"]["channel"].asInt());
        gVector Y(ygrid_);
        fill_grid_vector(Y, par_["initial_state"]["wavepacket"], Gaussian);
      // compose to two dimensional wave function
        psi_ = gVector2D(grid_,X,Y);
    }
  // Store the initial wave function
    psi_.save((folder_+"Psi.dat").c_str());
    psi_.save_binary(folder_+"Psi.qbin");

    std::cout << "Initial state norm:\t" << psi_.norm() << "\t saving into a file." << std::endl;

    potential_ = gVector2D(grid_);
    fill_grid_vector_2d(potential_, par_["model"]["potential"], Neutral2dPotential);
    potential_.save((folder_+"Potential.dat").c_str());
    potential_.save_binary(folder_+"Potential.qbin");

  // Initialization of the test function
    {
        if (par_["cross_sections"]["testfunctions"].isMember("dissociative_attachment")){
            xtestfunction_use_ = true;
            // FIXME
//            xtestfunction_ = new MultiTestFunction2d(&grid_, params_, 0, initial_energy_);
            xtestfunction_ = new MultiTestFunction2d(&grid_, par_, par_["cross_sections"]["testfunctions"]["dissociative_attachment"], initial_energy_);
        } else {
            xtestfunction_use_ = false;
        }
        if (par_["cross_sections"]["testfunctions"].isMember("vibrational_excitation")){
            ytestfunction_use_ = true;
            // femGrid2D *g, const pjvalue& params, const pjvalue& tfparams, const def_comp& init_erg
            ytestfunction_ = new MultiTestFunction2d(&grid_, par_, par_["cross_sections"]["testfunctions"]["vibrational_excitation"], initial_energy_);
        } else {
            ytestfunction_use_ = false;
        }
    }

  // Set evolution operators
    chebyshev_use_ = false;
    crank_nicolson_use_ = false;

    if (par_["evolution"]["approximation"].asString() == "Chebyshev") {
        chebyshev_use_ = true;
        zOperator2D H(grid_);
        H.set_kinetic_term(1.0, par_["model"]["reduced_mass"].asDouble());
        H += potential_;
        chebyshev_ = new Chebyshev2D(par_["evolution"]["order"].asInt(), par_["evolution"]["time_step"].asDouble(), H);
    } else if (par_["evolution"]["approximation"].asString() == "Crank-Nicolson") {
        crank_nicolson_use_ = true;
        zOperator2D H(grid_);
        H.set_kinetic_term(1.0, par_["model"]["reduced_mass"].asDouble());
        H += potential_;
        crank_nicolson_ = new CrankNicolson2D(par_["evolution"]["order"].asInt(), par_["evolution"]["time_step"].asDouble(), H);
    }
    time_ = 0.0;
    dt_ = par_["evolution"]["time_step"].asDouble();
    loop_size_ = par_["evolution"]["loop_steps"].asInt();

  // Initial contribution to S-matrix
    if (xtestfunction_use_){
        xtestfunction_->step_buffer(psi_,loop_size_);
    }
    if (ytestfunction_use_){
        ytestfunction_->step_buffer(psi_,loop_size_);
    }
    if (xtestfunction_use_){
        xtestfunction_->cross_sections(time_);
    }
    if (ytestfunction_use_){
        ytestfunction_->cross_sections(time_);
    }

    {   // Equidistant projector initializations
        const pjvalue& elc = par_["storage"]["equidistant_sampling"]["electronic"];
        const pjvalue& nuc = par_["storage"]["equidistant_sampling"]["nuclear"];

        ep_ = new EquidistantProjector2d(grid_, elc["samples"].asInt(), nuc["samples"].asInt(),
                                                elc["range"][0].asDouble(), elc["range"][1].asDouble(),
                                                nuc["range"][0].asDouble(), nuc["range"][1].asDouble() );

        hsv_ = new EquidistantProjector2d(grid_, 1200, 600,
                                                 elc["range"][0].asDouble(), elc["range"][1].asDouble(),
                                                 nuc["range"][0].asDouble(), nuc["range"][1].asDouble() );
    }
    // equidistant projection
    *ep_ << psi_;
    ep_->export_state((folder_ + "WF/PsiEq_t_0.dat").c_str());
    psi_.save_binary(folder_ + "WF/Psi_t_0.0.qbin");

    // HSV projection
    *hsv_ << psi_;
    char name[50];
    sprintf(name, "bmp/%d.bmp", int(time_));
    if (HSV) hsv_->export_state_hsv((folder_ + name).c_str(), 5e-1, real(initial_energy_) * time_);

  // ready to go!
    init_ = true;

    if (phid) discrete_state_ = *phid;
    if (discrete_state_.init()) {
        cout << "Using discrete state with norm =" << (*phid) * (*phid) << endl;
        discrete_state_.save_binary(folder_ + "PhiD.qbin");
        radius_operator_ = zOperatorD(ygrid_, ygrid_.nb(), 0);
        for (int i=0; i<ygrid_.nb(); ++i){
            radius_operator_[i] = ygrid_.xr(i);
        }
        discrete_projection();
        normalization_[0] = abs(psi_discrete_*psi_discrete_);
        radius_[0] = real(psi_discrete_*(radius_operator_*psi_discrete_))/normalization_();
    }
}

TimeDependentModel2D::~TimeDependentModel2D()
{
    clean();
}

TimeDependentModel2D& TimeDependentModel2D::set(parameters2D& p, parametersGrid& gpx, parametersGrid& gpy)
{
    assert(0);
    params_ = p;
    //initialize(gpx, gpy);
    return *this;
}

void TimeDependentModel2D::discrete_projection()
{
    psi_discrete_ = psi_.contraction(discrete_state_, 'x');
}

void TimeDependentModel2D::multistep()
{
    std::cout << "2D " << time_ << " ";
    char t_name[50];
    for (int i=0; i<loop_size_; ++i){
        if (chebyshev_use_){
            chebyshev_->one_step(psi_);
        } else if (crank_nicolson_use_){
            crank_nicolson_->one_step(psi_);
            //std::cout << " " << (*Psi) * (*Psi) << " ";
        }
        if (xtestfunction_use_){
            xtestfunction_->step_buffer(psi_,i);
        }
        if (ytestfunction_use_){
            ytestfunction_->step_buffer(psi_,i);
        }
        if (discrete_state_.init()) {
            discrete_projection();
            sprintf(t_name, "%d", int(time_ + (i+1)*dt_) );
            psi_discrete_.save_binary(folder_ + "WF/psiD_t_" + t_name + ".qbin");
            normalization_ << real(psi_discrete_*psi_discrete_);
            radius_ << real(psi_discrete_*(radius_operator_*psi_discrete_))/normalization_() ;
            if (population_file_.is_open()){
                population_file_ << time_ + (1.0+i)*dt_;
                for (int k=0; k<vib_levels;++k){
                    def_comp val = psi_discrete_ * vibrations_[k];
                    //population_file_ << "\t" << real(val) << "\t" << imag(val) << "\t" << abs(val*val);
                    population_file_ << "\t" << abs(val*val);
                }
                population_file_ << std::endl;
            }
        } else {
            normalization_ << psi_.norm();
        }
        std::cout<<'.';
    }
    if (xtestfunction_use_){
        xtestfunction_->close_multistep(time_,dt_);
    }
    if (ytestfunction_use_){
        ytestfunction_->close_multistep(time_,dt_);
    }
    time_ += loop_size_*dt_;
    if (xtestfunction_use_){
        xtestfunction_->cross_sections(time_);
    }
    if (ytestfunction_use_){
        ytestfunction_->cross_sections(time_);
    }

    sprintf(t_name, "%d", int(time_));
    psi_.save((folder_ + "Psi.dat").c_str());
    psi_.save_binary(folder_ + "WF/Psi_t_" + t_name + ".qbin");
    //psi_.save_equidistant((folder_+"WF/PsiEq_t_"+t_name+".dat").c_str(), params_.xSamples, params_.xRange[0], params_.xRange[1],
    //                                                                     params_.ySamples, params_.yRange[0], params_.yRange[1]);

    // EQ
    *ep_ << psi_;
    ep_->export_state((folder_ + "WF/PsiEq_t_" + t_name + ".dat").c_str());

    // HSV
    if (HSV) {
        *hsv_ << psi_;
        char name[50];
        double amplitude = (time_<1100)? 5e-1: 1e-1;
        amplitude = (time_<2100)? amplitude:1e-2;
        amplitude = (time_<3100)? amplitude:2e-3;
        amplitude = (time_<4100)? amplitude:1e-3;
        //amplitude = (time_<5100)? amplitude:2e-4;

        sprintf(name, "bmp/%d.bmp", int(time_));
        hsv_->export_state_hsv((folder_ + name).c_str(), amplitude, real(initial_energy_) * time_);
    }

    normalization_.save_range(0.0, time_, (folder_+"normalisation.dat").c_str());
    if (discrete_state_.init()) {
        psi_discrete_.save_equidistant((folder_+"WF/PsiD="+t_name+".dat").c_str(), params_.yRange[0], params_.yRange[1], params_.ySamples );
        radius_.save_range(0.0, time_, (folder_+"internuclear_distance.dat").c_str());
    }
    std::cout << "Evolved state norm:\t" << normalization_() << " at time t=" << time_ << std::endl;
}
} // namespace QSCAT
