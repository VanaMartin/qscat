#include "ModelLCP/ModelLCP.h"
#include "potentials.h"

using std::string;

namespace QSCAT
{

void LCP::ModelLCP::initialize(const pjvalue& parameters)
{
    folder_ = parameters["model"]["folder"].asString() + "/LCP/";

// Evolution method
    method_ = parameters["evolution"]["approximation"].asString();

// Grids -- Nuclear (main)
    const pjvalue& gp = parameters["grids"][ parameters["LCP"]["nuclear_grid"].asString() ];
    grid_ = grid_from_parameters(gp);

// Electronic grids
    const pjarray gpa = parameters["LCP"]["electronic_grids"].asArray();
    nel_grids_ = gpa.size();
    grid_e_ = new femGrid[nel_grids_];
    for (int i=0; i<nel_grids_; ++i){
        grid_e_[i] = grid_from_parameters(parameters["grids"][gpa[i].asString()]);
    }

// Vector setting
    v_zero_ = gVector(grid_);

// Energy and eigenstates declaration
    def_float emin = parameters["cross_sections"]["range"][0].asDouble();
    def_float emax = parameters["cross_sections"]["range"][1].asDouble();
    int esize = int((emax-emin) / parameters["cross_sections"]["dE"].asDouble());
    energies_ = dVector(esize, emin, emax, false);

    fill_grid_vector(v_zero_, parameters["model"]["potential"], MorsePotential );
    v_zero_.save_binary(folder_ + "bin/v_zero.bin");

    zOperatorF Ham(grid_);
    Ham.add_kinetic_term(parameters["model"]["reduced_mass"].asDouble()) += v_zero_;
    zEigenSystem esys = Ham.eigen_system();
    channels_ = parameters["cross_sections"]["testfunctions"]["vibrational_excitation"]["channels"].asInt();
    v_states_ = new gVector[channels_];
    v_energies_ = dVector(channels_+1);
    for (int k=0; k<channels_; ++k){
        v_states_[k] = gVector(grid_);
        esys.eigen_vector(v_states_[k].body(),k);
        v_energies_[k] = real(esys.eigen_value(k));
    }
    init_energy_ = real(esys.eigen_value(parameters["initial_state"]["channel"].asInt()));

// Resonant potential setting
    make_vres(parameters);
    v_energies_[channels_] = affinity_;

// Building the initial state and projectors for the S-matrix evaluation
    for (int k=0; k<channels_; ++k){
        for (int j=0; j<grid_.get_size(); ++j){
            v_states_[k].f(v_states_[k].f(j)*sqrt(gamma_.f(j)/(2.0*pi)),j);
        }
    }
    cout << parameters["initial_state"]["channel"].asInt() << endl;
    psi_ = v_states_[parameters["initial_state"]["channel"].asInt()].copy();

// Evolution operators declarations
    mu_ = parameters["model"]["reduced_mass"].asDouble();
    dt_ = parameters["evolution"]["time_step"].asDouble();
    time_ = 0.0;
    method_ = parameters["evolution"]["approximation"].asString();
    loop_ = parameters["evolution"]["loop_steps"].asInt();
    cn_use_ = false;
    cheb_use_ = false;
    zOperatorC H(grid_);
    H.set_kinetic_term(mu_);
    H += vres_;
    if (method_ == "Chebyshev") {
        cheb_order_ = parameters["evolution"]["order"].asInt();
        //cheb_ = new Chebyshev1D(grid_, vres_, mu_, dt_, cheb_order_);
        cheb_ = new Chebyshev1D(cheb_order_, dt_, H);
        cheb_use_ = true;
    } else if (method_ == "Crank-Nicolson"){
        cn_order_ = parameters["evolution"]["order"].asInt();
        //cn_ = new CrankNicolson1D(grid_, vres_, mu_, dt_, cn_order_);
        cn_ = new CrankNicolson1D(cn_order_, dt_, H);
        cn_use_ = true;
    }
    //S.Set(LCPp, m2dp.mu, m2dp.init_par.channel);

    opR_ = zOperatorD(grid_, grid_.get_size(), 0);
    for (int i=0; i<grid_.get_size(); ++i){
        opR_[i] = grid_.xr(i);
    }
    norm_[0] = abs(psi_*psi_);
    radius_[0] = real(psi_*(opR_*psi_))/norm_();
    S_.contribution(psi_,v_states_,loop_);  // Initial contribution
    psi_.save_equidistant((folder_+"WF/Psi_t_0.dat").c_str(), 1.6, 5.2, 600);
    init_ = true;
}

void LCP::ModelLCP::clean()
{
    if (init_){
        delete[] v_states_;
        delete[] grid_e_;
        delete[] phi_res_;
        delete dstates_;
        if (cn_use_) delete cn_;
        if (cheb_use_) delete cheb_;
    }
    init_ = false;
}

bool LCP::ModelLCP::read_eigenstates(const char* name, femGrid& g)
{
    bool stat;
    dstates_ = new DiscreteStates1D;
    stat = dstates_->from_file(name,g);
    // CHECK THE PARAMETERS -- TO BE ADDED

    if (!stat) { delete dstates_; }
    return stat;
}

bool LCP::ModelLCP::read_vres()
{
    bool stat = vres_.read_binary(folder_+"bin/vres.bin");
    if (stat) stat = (vres_.get_grid() == grid_);
    if (stat) stat = eres_.read_binary(folder_+"bin/eres.bin");
    if (stat) stat = (eres_.get_grid() == grid_);
    if (stat) stat = gamma_.read_binary(folder_+"bin/gamma.bin");
    if (stat) stat = (gamma_.get_grid() == grid_);

    phi_res_ = new gVector[grid_.nr()];
    if (stat) stat = ReadMultiGridVectorBin((folder_+"bin/phi_res.bin").c_str(), phi_res_, grid_e_[0]);
    if (!stat) {
        delete[] phi_res_;
        phi_res_ = NULL;
    }
    return stat;
}

bool LCP::ModelLCP::save_vres()
{
    bool stat = vres_.save_binary((folder_+"bin/vres.bin").c_str());
    if (stat) stat = eres_.save_binary((folder_+"bin/eres.bin").c_str());
    if (stat) stat = gamma_.save_binary((folder_+"bin/gamma.bin").c_str());
    if (stat) stat = SaveMultiGridVectorBin((folder_+"bin/phi_res.bin").c_str(), phi_res_, grid_.nr());
    return stat;
}

void LCP::ModelLCP::make_vres(const pjvalue& parameters)
{
/*  The procedure  initializing the resonant potential.  The procedure
    tries to read  data from  previous  computations.  If it fails the
    procedure coumputes the poles of the electron Hamiltonian.      */
    std::cout << "Starting the V_res procedure... " << std::endl;
// Initialization of auxiliary parameters
    int nr = grid_.nr();
    def_float dumm;

// Reading or calculating the vibrational states

    std::cout << "Trying to read electronic bound states at R -> infinity from previous run ..." << std::endl;
    // Read - to be completed later
    //   CALL Read_Eigenstates('Data/electron', '', eig_e_DA, grid_e)
    dstates_ = new DiscreteStates1D;
    if (dstates_->read_binary((folder_+"bin/electron_states.bin").c_str()) && dstates_->get_grid() == grid_e_[0]) {
    // Successfull
        std::cout << "  ... successful " << dstates_->number_of_states() << " states" << std::endl;

    } else {
        delete dstates_;
        std::cout << "  ... not successful. Calculating electronic bound states ... " << std::endl;
    // effective electronic potential Vint(R->inf,r) + J_e (J_e + 1) / 2 mu_e r^2 on the first electronic grid (only real bound states wanted)
        gVector* e_pot = new gVector[nel_grids_];
        for (int i=0; i<nel_grids_; ++i){
            e_pot[i] = gVector(grid_e_[i]);
            fill_grid_vector(e_pot[i], parameters["model"]["potential"], AsymptoticLambda);
        }
    // Calculating the electronic states
        dstates_ = new DiscreteStates1D(nel_grids_, grid_e_, e_pot, 1.0);
        if (dstates_->number_of_states() == 0){
            std::cout << "No discrete states were found for the asymptotic state! Aborting computation!" << std::endl;
            exit(340649);
        } else {
            dstates_->save_binary((folder_+"bin/electron_states.bin").c_str());
        }
        delete[] e_pot;
    }

// Electron affinity
    affinity_ = real(dstates_->get_energy(0));
    phi_a_ = gVector(grid_e_[0]);
    dstates_->get_state(phi_a_, 0);

    std::cout << std::endl << "Trying to read V_res from previous run... ";
    if (read_vres()) {
        // Success the data were obtained from previous computation
        std::cout << "success!" << std::endl << std::endl;
    } else {
        // Failure, the data were not found
        std::cout << std::endl << "failed! Calculating V_res... " << std::endl;

        vres_ = gVector(grid_);
        eres_ = gVector(grid_);
        gamma_ = gVector(grid_);

        vres_.fill(0.0);
        eres_.fill(0.0);
        gamma_.fill(0.0);

        gVector *e_pot = new gVector[nel_grids_];     // allocating the electronic grids
        DiscreteStates1D *W;
        def_float prec = 1.0e-14;
        def_float low_b = affinity_;
        def_comp ext;
        int start = grid_.nr()-1;
        phi_res_ = new gVector[grid_.nr()];              // Allocating the pointers to the discrete states
        for (int i=start; i>=0; --i){                   // i: position in the nuclear grid
        // setting the effective electronic potential on given grids
            for (int j=0; j<nel_grids_; ++j){
                e_pot[j] = gVector(grid_e_[j]);
                fill_grid_vector_xaxis(e_pot[j],grid_.xr(i),parameters["model"]["potential"], ElectronicLambdaInteraction);
                //e_pot[j].Make_fxc(potentials::V_eff_el,grid.Xr(i),m2dp);
            }
            e_pot[0].save("test.dat");
        // Calculating the discrete states
repeat:
            W = new DiscreteStates1D(nel_grids_, grid_e_, e_pot, 1.0, prec, low_b - 0.001, low_b + 0.1);
            if (W->number_of_states() == 0){
                std::cout << " at R = "<< grid_.xr(i) << " !" << std::endl;
                if (prec < 1e-4){
                    if (i < start - 1) {
                        low_b = real(eres_.f(i+1));
                    }
                    prec*=2.0;
                    std::cout << "Adjusting precision to "<< prec << "..." << std::endl;
                    delete W;
                    goto repeat;
                } else {
                    // Evaluation of the interpolated state
                    if (i+2 <= grid_.nr()){
                        ext = (vres_.f(i + 2) - vres_.f(i + 1)) / (grid_.xr(i + 2) - grid_.xr(i + 1));  // slope
                        ext = ext * (grid_.xr(i) - grid_.xr(i + 1)) + vres_.f(i + 1);                  // extrapolated value
                        vres_.f(ext, i);
                        ext -= v_zero_.f(i);
                        eres_.f(real(ext), i);
                        gamma_.f(-2.0 * imag(ext), i);
                        //phi_res_[i] = gVector(grid_e_[0]);
                        phi_res_[i] = phi_res_[i+1];
                    } else {
                        vres_.f(affinity_, i);
                        eres_.f(affinity_, i);
                        gamma_.f(0.0, i);
                        // Phi res to be added ------ ???? FIXME!
                    }
                }
            } else {
                std::cout << " at R = "<< grid_.xr(i) << "..." << W->get_energy(0) << std::endl;
                if (W->number_of_states() > 1) {
                    prec/=10.0;
                    std::cout << "Adjusting precision to "<< prec << "..." << std::endl;
                    delete W;
                    goto repeat;
                }
                vres_.f(W->get_energy(0) + v_zero_.f(i), i);
                eres_.f(real(W->get_energy(0)), i);
                dumm = -2.0 * imag(W->get_energy(0));
                gamma_.f((dumm>0.0)? dumm: 0.0, i);
                phi_res_[i] = gVector(grid_e_[0]);
                W->get_state(phi_res_[i], 0);
            }
            delete W;
        }
        delete[] e_pot;
        save_vres();
    }
    for (int i=grid_.nr(); i<grid_.nb(); ++i){
        vres_.f(vres_.f(nr-1),i);
        v_zero_.f(v_zero_.f(nr-1),i);
    }
    vres_.save((folder_+"V_res.dat").c_str());
    v_zero_.save((folder_+"V_zero.dat").c_str());
    gamma_.save((folder_ + "Gamma.dat").c_str());
}

LCP::ModelLCP::ModelLCP()
{
    init_ = false;
    cn_use_ = false;
    cheb_use_ = false;
}

LCP::ModelLCP::ModelLCP(const pjvalue& parameters) :
    S_(parameters)
{
    initialize(parameters);
}

LCP::ModelLCP::~ModelLCP()
{
    clean();
}

void LCP::ModelLCP::fill_discrete_state_phys(gVector2D& phi)
{
    for (int i=0; i<phi.get_real_ysize(); ++i){
        phi.write_x_section(phi_res_[i],i);
    }
}

void LCP::ModelLCP::multistep()
{
    std::cout<<"LCP ";
    char t_name[50];
    for (int i=0; i<loop_; ++i){
        if (cheb_use_){
            cheb_->one_step(psi_);
        } else if (cn_use_){
            cn_->one_step(psi_);
        }
        S_.contribution(psi_,v_states_,i);
        norm_ << abs(psi_*psi_);
        radius_ << real(psi_*(opR_*psi_))/norm_();
        std::cout << '.';
        sprintf(t_name, "%d", int(time_ + i*dt_));
        psi_.save_binary(folder_ + "WF/psi_t_" + t_name + ".qbin");
    }
    S_.close_multistep(time_, dt_, v_energies_, init_energy_, grid_.xr(grid_.nr()));
    time_ += loop_*dt_;

    S_.cross_sections(time_, folder_);

    sprintf(t_name, "%d", int(time_));
    psi_.save_equidistant((folder_+"WF/Psi_t_" + t_name + ".dat").c_str(), 1.6, 5.2, 600);
    norm_.save_range(0, time_, (folder_+"normalisation.dat").c_str());
    radius_.save_range(0, time_, (folder_+"internuclear_distance.dat").c_str());
    //psi->Save_equidistant((folder+"LCP/PsiEq.dat").c_str(), 50, 0.0, 20.0, 50, 1.8, 2.8);
    std::cout << "Evolved state norm:\t" << norm_() << " at time t=" << time_ << std::endl;
}

gVector2D LCP::ModelLCP::get_discrete_state(const pjvalue& parameters)
{
    assert(init_);
    femGrid ge = grid_from_parameters(parameters["grids"][parameters["model"]["electronic_grid"].asString()]);

    femGrid2D G(ge, grid_);
    gVector2D PhiD(G);
    for (int i=0; i<grid_.nr(); ++i){
        for (int j=0; j<grid_e_[0].nr(); ++j ) {
            PhiD.f( phi_res_[i].f(j), i, j );
        }
        for (int j=grid_e_[0].nr(); j<PhiD.get_xsize(); ++j)
            PhiD.f(0, i, j);
        //PhiD.write_x_section(phi_res_[i],i);
    }
    // Smoothing of the wavefunction psi(r,R) -> psi(r,R) * N(R) * exp(-iD(R)) * f(r)
    // f(r) = 1 - 1 / (1 + exp(-(r-r_d)))
    // r_d = 10.0 a.u.
    int iMax = 0;
    for (int i=grid_.nr()-1; i>=0; --i){
        //std::cout << grid.Xr(i) << "... ";
        int nb = PhiD.get_xsize();
        if (i==grid_.nr()-1) {
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
        if (i<grid_.nr()-1) {
            def_comp f = PhiD.f(i*nb + iMax);
            def_float phi = std::arg(f);
            val_0 = exp(-imu*phi);
        }
        //std::cout << val_0 << "... ";
        def_float N = 0.0;
        for (int j=nb-1; j>=0; --j){
            def_comp val = val_0 * PhiD.f(i*nb + j) * (1.0 - 1.0/(1.0 + exp( - (ge.xr(j) - 10.0))));
            PhiD.f(val , i*nb + j);
            //if (j<grid_e[0].NR())
            N += std::pow(std::abs(PhiD[i*nb + j]/sqrt(grid_.w(i))),2) ;
        }
        N = 1.0/std::sqrt(N);
        //std::cout << N << "... " << std::endl;
        for (int j=nb-1; j>=0; --j){
            PhiD[i*nb + j] *= N;
        }
    }
    cout << "Saving discrete state" << endl;
    PhiD.save_binary((folder_ + "PhiD.qbin").c_str());
    PhiD.save((folder_+"PhiD.dat").c_str());
    PhiD.save_equidistant((folder_+"PhiDEq.dat").c_str(), 200, 0.0, 20.0, 400, 0.0, 5.0);

    return PhiD;
}

}
