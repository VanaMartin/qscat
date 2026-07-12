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
#include "interface.h"
#include "module_LCP.h"

using namespace std;
using namespace QSCAT;
using namespace QSCAT::LCP;

// S matrix methods
void S_Matrix::Initialize(parametersLCP& p, def_float & MU, int & iChannel)
{
    if (init) {
        Clean();        // Flush previous data
    }
    mu = MU;
    init_channel = iChannel;
    order = p.order;
    channels = p.ve_channels + p.da_channels;
    ve_channels = p.ve_channels;
    da_channels = p.da_channels;
    steps = p.steps;
    size = p.e_steps;
    Make_coefficients();
    e = dVector(size, p.e_min, p.e_max, false);
    buffer = new zVector[channels];
    s = new zVector[channels];
    for (int i=0; i<channels; ++i){
        buffer[i] = zVector(steps+1);
        buffer[i].fill(0.0);
        s[i] = zVector(size);
        s[i].fill(0.0);
    }
    init = true;
}

void S_Matrix::Clean()
{
    if (init) {
        delete[] buffer;
        delete[] s;
    }
}

void S_Matrix::Make_coefficients()
{
    switch (order){
        case 0:
            coefficients = dVector(steps);
            coefficients.fill(1.0);
            break;
        case 1:
            coefficients = dVector(steps+1);
            coefficients.fill(1.0);
            coefficients[steps-1] = 0.5;
            coefficients[steps] = 0.5;      // The last element has to be added again if the integration continues another time
            break;
        case 2:
            coefficients = dVector(steps+1);
            coefficients[0] = 4.0/3.0;
            for (int i=1;i<steps-1;++i){
                if (i%2==1) {
                    coefficients[i]= 4.0/3.0;
                } else {
                    coefficients[i]= 2.0/3.0;
                }
            }
            coefficients[steps-1] = 1.0/3.0;
            coefficients[steps] = 1.0/3.0;      // The last element has to be added again if the integration continues another time
    }
}

S_Matrix::S_Matrix()
{
    init=false;
}

S_Matrix::S_Matrix(parametersLCP& LCPp, def_float & MU, int & iChannel)
{
    Initialize(LCPp, MU, iChannel);
}

S_Matrix::~S_Matrix()
{
    Clean();
}

S_Matrix& S_Matrix::Set(parametersLCP& LCPp, def_float & MU, int & iChannel)
{
    Initialize(LCPp, MU, iChannel);
    return *this;
}

void S_Matrix::Contribution(gVector& psi, gVector *states, int & i)
{
    for (int j=0; j<ve_channels; ++j){
        // Assuming the states to already contain the factor sqrt(Gamma(R)/2pi)
        buffer[j][i] = states[j]*psi;
    }
    for (int j=ve_channels; j<channels; ++j){
        buffer[j][i] = psi.f(psi.get_grid().nr());
    }
}

void S_Matrix::Close_multistep(def_float & time, def_float & dt, dVector& ergs, const def_float & ierg, const def_float & X)
{
    for (int i=0; i<ve_channels; ++i){      // Contributions to vibrational excitations
        def_float shift = ierg - ergs[i];
        for (int k=0; k<size; ++k){
            if (e[k] + shift>0) {
                if (order > 0) { // Contribution of the previous loop end
                    s[i][k] += (1.0/imu)*exp(imu*(e[k] + ierg)*(time))*buffer[i][steps]*dt*coefficients[steps];
                }
                for (int j=0; j<steps; ++j){
                    s[i][k] += (1.0/imu)*exp(imu*(e[k] + ierg)*(time + (j+1)*dt))*buffer[i][j]*dt*coefficients[j];
                }
            }
        }
        buffer[i][steps] = buffer[i][0];
    }
    for (int i=ve_channels; i<channels; ++i){   // Contributions to dissociative attachment
        def_float shift = ierg - ergs[i];
        for (int k=0; k<size; ++k){
            if (e[k] + shift>0) {
                def_float K = sqrt(2.0*mu*(e[k]+shift));
                if (order > 0) { // Contribution of the previous loop end
                    s[i][k] += sqrt(K/(2.0*pi*mu)) * exp(-imu*K*X) *
                        exp( imu*(e[k] + ierg )*(time)) * buffer[i][steps] *
                        dt * coefficients[steps];
                }
                for (int j=0; j<steps; ++j){
                    s[i][k] += sqrt(K/(2.0*pi*mu)) * exp(-imu*K*X) *
                        exp( imu*(e[k] + ierg)*(time + (j+1)*dt)) * buffer[i][j] *
                        dt * coefficients[j];
                }
            }
        }
        if (order > 0) {
            buffer[i][steps] = buffer[i][0];
        }
    }
}

void S_Matrix::cross_sections(const def_float & time, std::string & folder) // Derives the cross sections and stores them into a file with appropriate time
{
    char name[50];
    def_float val;
    def_float erg;
    if (ve_channels!=0){
        sprintf_s(name, "LCP/CS/CSVE_t_%.1f.dat", time);
        FILE * file;
        fopen_s(&file,(folder+name).c_str(),"w");
        fprintf(file,"#Cross Sections of vibrational excitations: channels, values\n");
        fprintf(file,"#\t%d, %d\n", ve_channels, size);
        for (int i=0;i<size;++i){
            erg = e[i];
            fprintf(file, "%.12E", erg);
            for (int j=0;j<ve_channels;++j){
                val = pow(abs(s[j][i]),2)*4.0*pow(pi,3)/(2.0*erg);
                fprintf(file, "\t%.12E", val);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
    if (channels>ve_channels) {
        sprintf_s(name, "LCP/CS/CSDA_t_%.1f.dat", time);
        FILE * file;
        fopen_s(&file,(folder+name).c_str(),"w");
        fprintf(file,"#Cross Sections of vibrational excitations: channels, values\n");
        fprintf(file,"#\t%d, %d\n", channels-ve_channels, size);
        for (int i=0;i<size;++i){
            erg = e[i];
            fprintf(file, "%.12E", erg);
            for (int j=ve_channels;j<channels;++j){
                val = pow(abs(s[j][i]),2)*4.0*pow(pi,3)/(2.0*erg);
                fprintf(file, "\t%.12E", val);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

// LCP model methods
void Model_LCP::Initialize(parameters2D& m2dp, parametersMultiGrid& mgp, parametersLCP& LCPp)
{
    folder = m2dp.folder;
// Evolution method
    method = m2dp.evol_par.evolution_c;
// Grids -- Nuclear (main)
    main_grid = LCPp.nuclear_grid;
    grid = femGrid(mgp.gp[main_grid]);
// Electronic grids
    nel_grids = LCPp.nel_grids;
    el_grids = new int[nel_grids];
    for (int i=0; i<nel_grids; i++){
        el_grids[i] = LCPp.electronic_grids[i];
        if (mgp.n - 1 < LCPp.electronic_grids[i]) {
            std::cout << "Error! The index of electronic grid exceeds the bounds of the grids input!" <<std::endl;
            exit(6049940);
        }
    }

    grid_e = new femGrid[nel_grids];
    for (int i=0; i<nel_grids; ++i){
        grid_e[i] = femGrid(mgp.gp[el_grids[i]]);
    }
// Vector setting
    v_zero = gVector(grid);
    vres = gVector(grid);
    eres = gVector(grid);
    gamma = gVector(grid);
// Energy and eigenstates declaration
    erg = dVector(m2dp.evol_par.steps, m2dp.evol_par.e_min, m2dp.evol_par.e_max, false);
    v_zero.function_fill(potentials::V_zero<def_float,def_comp>, m2dp);
    //Hamiltonian_Full<T,Z> Ham(grid,v_zero,m2dp.mu);     // Auxiliary Hamiltonian for determining the eigenvalues of V_0
    //Ham.SetEigen();
    zOperatorF Ham = buildFullHamiltonian(grid, v_zero, m2dp.mu);
    zEigenSystem eSys = Ham.eigen_system();
    channels = LCPp.ve_channels;
    v_states = new gVector[channels];
    Energies = dVector(channels+1);
    for (int k=0; k<channels; ++k){
        v_states[k] = gVector(grid);
        eSys.eigen_vector(v_states[k].body(),k);
        Energies[k] = real(eSys.eigen_value(k));
    }
    InitEnergy = real(eSys.eigen_value(m2dp.init_par.channel));
// Resonant potential setting
    Make_Vres(m2dp);
    Energies[channels] = affinity;

// Building the initial state and projectors for the S-matrix evaluation
    for (int k=0; k<channels; ++k){
        for (int j=0; j<grid.nb(); ++j){
            v_states[k].f(v_states[k].f(j)*sqrt(gamma.f(j)/(2.0*pi)),j);
        }
    }
    psi = v_states[m2dp.init_par.channel];

// Evolution operators declarations
    mu = m2dp.mu;
    dt = m2dp.evol_par.dt;
    time = 0.0;
    method = LCPp.method;
    loop = m2dp.evol_par.loop;
    cn_use = false;
    cheb_use = false;
    if (method == 'c') {
        cheb_order = LCPp.cheb_order;
        //cheb = new Chebyshev1D(grid,vres,mu,dt,cheb_order);
        zOperatorC H(grid);
        H.set_kinetic_term(mu);
        H+=vres;
        cheb = new Chebyshev1D(cn_order, dt, H);
        cheb_use = true;
    } else if (method == 'n'){
        cn_order = LCPp.cn_order;
        //cn = new CrankNicolson1D(grid,vres,mu,dt,cn_order);
        zOperatorC H(grid);
        H.set_kinetic_term(mu);
        H+=vres;
        cn = new CrankNicolson1D(cn_order, dt, H);
        cn_use = true;
    }
    S.Set(LCPp, m2dp.mu, m2dp.init_par.channel);

    opR = zOperatorD(grid, grid.nb(), 0);
    for (int i=0; i<grid.nb(); ++i){
        opR[i] = grid.xr(i);
    }
    norm[0] = abs(psi*psi);
    radius[0] = real(psi*(opR*psi))/norm();
    S.Contribution(psi,v_states,loop);  // Initial contribution
    psi.save_equidistant((folder+"LCP/WF/Psi_t_0.dat").c_str(), 1.6, 5.2, 600);

  // Resonant vibrational states
    n_vibrations_ = 15;
    vibrations_ = new gVector[n_vibrations_];
    zOperatorF Hv(grid);
    Hv.add_kinetic_term(mu);
    gVector hvres(grid);
    for (int i=0; i<grid.nb(); ++i)
        hvres.f(real(vres.f(i)),i);
    Hv += hvres;
    zEigenSystem veSys = Hv.eigen_system();
    for (int i=0; i<n_vibrations_; ++i){
        vibrations_[i] = gVector(grid);
        veSys.eigen_vector(vibrations_[i].body(),i);
    }

    population_file_.open((folder + "vibrational_populations_lcp.txt").c_str());

    init = true;
}

void Model_LCP::Clean()
{
    if (init){
        delete[] v_states;
        delete[] el_grids;
        delete[] grid_e;
        delete[] phi_res;
        delete D;
        if (cn_use) delete cn;
        if (cheb_use) delete cheb;
    }
    if (population_file_.is_open())
        population_file_.close();
    init = false;
}

bool Model_LCP::Read_eigenstates(const char * name, femGrid& g)
{
    bool stat;
    D = new DiscreteStates;
    stat = D->from_file(name,g);
    // TODO CHECK THE PARAMETERS

    if (!stat) { delete D; }
    return stat;
}

bool Model_LCP::Read_Vres()
{
    phi_res = new gVector[grid.nr()];
    bool stat = vres.read_binary((folder+"LCP/bin/vres.bin").c_str(),grid);
    if (stat) stat = eres.read_binary((folder+"LCP/bin/eres.bin").c_str(),grid);
    if (stat) stat = gamma.read_binary((folder+"LCP/bin/gamma.bin").c_str(),grid);
    if (stat) stat = ReadMultiGridVectorBin((folder+"LCP/bin/phi_res.bin").c_str(), phi_res, grid_e[0]);
    if (!stat) {
        delete[] phi_res;
        phi_res = NULL;
    }
    return stat;
}

bool Model_LCP::Save_Vres()
{
    bool stat = vres.save_binary((folder+"LCP/bin/vres.bin").c_str());
    if (stat) stat = eres.save_binary((folder+"LCP/bin/eres.bin").c_str());
    if (stat) stat = gamma.save_binary((folder+"LCP/bin/gamma.bin").c_str());
    if (stat) stat = SaveMultiGridVectorBin((folder+"LCP/bin/phi_res.bin").c_str(), phi_res, grid.nr());
    return stat;
}

void Model_LCP::Make_Vres(parameters2D& m2dp)
{
/*  The procedure  initializing the resonant potential.  The procedure
    tries to read  data from  previous  computations.  If it fails the
    procedure coumputes the poles of the electron Hamiltonian.      */
    std::cout << "Starting the V_res procedure... " << std::endl;
// Initialization of auxiliary parameters
    int nr = grid.nr();
    //int l = m2dp.l;
    //T mu_e = 1.0;
    //T max_eig_e = 1.0;
    //T min_prec = 1.0E-4;
    //T max_prec = 2.0E-3;
    //T prec_e = min_prec;
    def_float dumm;

// Reading or calculating the vibrational states

    std::cout << "Trying to read electronic bound states at R -> infinity from previous run ..." << std::endl;
    // Read - to be completed later
    //   CALL Read_Eigenstates('Data/electron', '', eig_e_DA, grid_e)
    D = new DiscreteStates;
    if (D->read_binary((folder+"LCP/bin/electron_states.bin").c_str(), this->grid_e[0])) {
    // Successfull
        std::cout << "  ... successful " << D->number_of_states() << " states" << std::endl;

    } else {
        delete D;
        std::cout << "  ... not successful. Calculating electronic bound states ... " << std::endl;
    // effective electronic potential Vint(R->inf,r) + J_e (J_e + 1) / 2 mu_e r^2 on the first electronic grid (only real bound states wanted)
        gVector *e_pot = new gVector[nel_grids];
        for (int i=0; i<nel_grids; ++i){
            e_pot[i] = gVector(grid_e[i]);
            e_pot[i].function_fill(potentials::Attached_electron,m2dp);
        }
    // Calculating the electronic states
        D = new DiscreteStates(nel_grids,grid_e,e_pot,1.0);
        if (D->number_of_states() == 0){
            std::cout << "No discrete states were found for the asymptotic state! Aborting computation!" << std::endl;
            exit(340649);
        } else {
            D->save_binary((folder+"LCP/bin/electron_states.bin").c_str());
        }
        delete[] e_pot;
    }

// Electron affinity
    affinity = real(D->get_energy(0));
    phi_a = gVector(grid_e[0]);
    D->get_state(phi_a, 0);

    std::cout << std::endl << "Trying to read V_res from previous run... ";
    if (Read_Vres()) {
        // Success the data were obtained from previous computation
        std::cout << "success!" << std::endl << std::endl ;
    } else {
        // Failure, the data were not found
        std::cout << std::endl << "failed! Calculating V_res... " << std::endl;
        gVector *e_pot = new gVector[nel_grids];     // allocating the electronic grids
        DiscreteStates *W;
        def_float prec = 1.0e-14;
        def_float low_b = affinity;
        def_comp ext;
        int start = grid.nr()-1;
        phi_res = new gVector[grid.nr()];      // Allocating the pointers to the discrete states
        for (int i=start; i>=0; --i){                   // i: position in the nuclear grid
        // setting the effective electronic potential on given grids
            for (int j=0; j<nel_grids; ++j){
                e_pot[j] = gVector(grid_e[j]);
                e_pot[j].function_xy_fill_x(potentials::V_eff_el, grid.xr(i),m2dp);  // fxc
            }
            e_pot[0].save("test.dat");
        // Calculating the discrete states
repeat:
            W = new DiscreteStates(nel_grids,grid_e,e_pot,1.0, prec, low_b - 0.001, low_b + 0.1);
            if (W->number_of_states() == 0){
                std::cout << " at R = "<< grid.xr(i) << " !" << std::endl;
                if (prec < 1e-4){
                    if (i < start - 1) {
                        low_b = real(eres.f(i+1));
                    }
                    prec*=2.0;
                    std::cout << "Adjusting precision to "<< prec << "..." << std::endl;
                    delete W;
                    goto repeat;
                } else {
                    // Evaluation of the interpolated state
                    if (i+2 <= grid.nr()){
                        ext = (vres.f(i + 2) - vres.f(i + 1)) / (grid.xr(i + 2) - grid.xr(i + 1));  // slope
                        ext = ext * (grid.xr(i) - grid.xr(i + 1)) + vres.f(i + 1);                  // extrapolated value
                        vres.f(ext, i);
                        ext -= v_zero.f(i);
                        eres.f(real(ext), i);
                        gamma.f(-2.0 * imag(ext), i);
                        phi_res[i] = gVector(grid_e[0]);
                        phi_res[i] = phi_res[i+1];
                    } else {
                        vres.f(affinity, i);
                        eres.f(affinity, i);
                        gamma.f(0.0, i);
                        // Phi res to be added
                    }
                }
            } else {
                std::cout << " at R = "<< grid.xr(i) << "..." << W->get_energy(0) << std::endl;
                if (W->number_of_states() > 1) {
                    prec/=10.0;
                    std::cout << "Adjusting precision to "<< prec << "..." << std::endl;
                    delete W;
                    goto repeat;
                }
                vres.f(W->get_energy(0) + v_zero.f(i), i);
                eres.f(real(W->get_energy(0)), i);
                dumm = -2.0 * imag(W->get_energy(0));
                gamma.f((dumm>0.0)? dumm: 0.0, i);
                phi_res[i] = gVector(grid_e[0]);
                W->get_state(phi_res[i], 0);
            }
            delete W;
        }
        delete[] e_pot;
        Save_Vres();
    }
    for (int i=grid.nr(); i<grid.nb(); ++i){
        vres.f(vres.f(nr-1),i);
        v_zero.f(v_zero.f(nr-1),i);
    }
    vres.save((folder+"LCP/V_res.dat").c_str());
    v_zero.save((folder+"LCP/V_zero.dat").c_str());
    gamma.save((folder + "LCP/Gamma.dat").c_str());
}

Model_LCP::Model_LCP()
{
    init = false;
    cn_use = false;
    cheb_use = false;
}

Model_LCP::Model_LCP(parameters2D& m2dp, parametersMultiGrid& mgp, parametersLCP& LCPp)
{
    Initialize(m2dp,mgp,LCPp);
}

Model_LCP::~Model_LCP()
{
    Clean();
}

void Model_LCP::MakePhiD_phys(gVector2D& phi)
{
    for (int i=0; i<phi.get_real_ysize(); ++i){
        phi.write_x_section(phi_res[i],i);
    }
}

void Model_LCP::Multistep()
{
    std::cout<<"LCP ";
    for (int i=0; i<loop; ++i){
        if (cheb_use){
            cheb->one_step(psi);
        } else if (cn_use){
            cn->one_step(psi);
        }
        S.Contribution(psi,v_states,i);
        norm << abs(psi*psi);
        radius << real(psi*(opR*psi))/norm();
        if (population_file_.is_open()){
            population_file_ << time + (1.0+i)*dt;
            for (int k=0; k<n_vibrations_;++k){
                def_comp val = psi * vibrations_[k];
                //population_file_ << "\t" << real(val) << "\t" << imag(val) << "\t" << abs(val*val);
                population_file_ << "\t" << abs(val*val);
            }
            population_file_ << std::endl;
        }
        std::cout<<'.';
    }
    S.Close_multistep(time, dt, Energies, InitEnergy, grid.xr(grid.nr()));
    time += loop*dt;

    S.cross_sections(time, folder);

    char* t_name = new char[50];
    sprintf(t_name, "%d", int(time));
    psi.save_equidistant((folder+"LCP/WF/Psi_t_" + t_name + ".dat").c_str(), 1.6, 5.2, 600);
    norm.save_range(0,time, (folder+"LCP/normalisation.dat").c_str());
    radius.save_range(0,time, (folder+"LCP/internuclear_distance.dat").c_str());
    //psi->Save_equidistant((folder+"LCP/PsiEq.dat").c_str(), 50, 0.0, 20.0, 50, 1.8, 2.8);
    std::cout << "Evolved state norm:\t" << norm() << " at time t=" << time << std::endl;
}

gVector2D Model_LCP::MakePhiD() // Make Discrete State (PhiD)
{
    assert(init);
    femGrid2D G(grid_e[0],grid);
    gVector2D PhiD(G);
    for (int i=0; i<grid.nr(); ++i){
        PhiD.write_x_section(phi_res[i],i);
    }
    // Smoothing of the wavefunction psi(r,R) -> psi(r,R) * N(R) * exp(-iD(R)) * f(r)
    // f(r) = 1 - 1 / (1 + exp(-(r-r_d)))
    // r_d = 10.0 a.u.
    int iMax = 0;
    for (int i=grid.nr()-1; i>=0; --i){
        //std::cout << grid.Xr(i) << "... ";
        int nb = PhiD.get_xsize();
        if (i==grid.nr()-1) {
            def_float val = 0.0;
            for (int j=0; j<nb; ++j){
                if (std::abs(PhiD.f(i*nb+j))>val) {
                    iMax = j;
                    val = std::abs(PhiD.f(i*nb+j));

                }
            }
            //std::cout << "Maximum at " << iMax << " with " << val << "... ";
        }
        def_comp val_0 = 1.0;
        if (i<grid.nr()-1) {
            def_comp f = PhiD.f(i*nb + iMax);
            def_float phi = std::arg(f);
            val_0 = exp(-imu*phi);
        }
        //std::cout << val_0 << "... ";
        def_float N = 0.0;
        for (int j=nb-1; j>=0; --j){
            def_comp val = val_0 * PhiD.f(i*nb + j) * (1.0 - 1.0/(1.0 + exp( - (grid_e[0].xr(j) - 10.0))));
            PhiD.f(val , i*nb + j);
            //if (j<grid_e[0].NR())
            N += std::pow(std::abs(PhiD[i*nb + j]/sqrt(grid.w(i))),2) ;
        }
        N = 1.0/std::sqrt(N);
        //std::cout << N << "... " << std::endl;
        for (int j=nb-1; j>=0; --j){
            PhiD[i*nb + j] *= N;
        }
    }
    PhiD.save((folder+"PhiD.dat").c_str());
    PhiD.save_equidistant((folder+"PhiDEq.dat").c_str(), 200, 0.0, 20.0, 400, 0.0, 5.0);

    return PhiD;
}
