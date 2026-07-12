// Main project file of electron-Molecule-Scattering project
// The project is under construction
// It is expected to contain following modules
// 2D model
// LCP
// NRM

#include <iostream>				// Input/Output library
#include <complex>				// Complex algebra
#include <math.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <omp.h>				// Parallelization library

#ifdef linux
	#include <sec_stream.h>
#endif

#include "common.h"		// Common functions, parameters and classes
#include "bessel.h"		// Solutions to the coulomb potential Schrödinger equation (sph. Bessel functions)
#include "coulomb.h"
#include "blas.h"		// Intel BLAS, SPARSE-BLAS LAPACK and PARDISO wrapper

#include "Arrays.h"
#include "potentials.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"
#include "interface.h"

#include "pjinput.h"
//#include "module_NRM.h"
#include "Model2d.h"

using namespace QSCAT;

enum ComputationMode
{
    MODE_FULL,
    MODE_QDF,
    MODE_ENERGY_THRESHOLD,
    MODE_TIME_INDEPENDENT,
    MODE_FULL_WITH_POPULATIONS
};

/**
 *  Quasi Bound States class
 */

class QBStates
{
    gVector2D **states_;
    gVector **generators_;
    int vibs_;
    int ryds_;
    ofstream *outputs_;
 public:
    QBStates(int vibrations, int rydbergs)
    {
        vibs_ = vibrations;
        ryds_ = rydbergs;
        states_ = NULL;
        outputs_ = NULL;
        generators_ = NULL;
    }
    void init_states(femGrid2D& g)
    {
        assert(!states_);
        states_ = new gVector2D*[vibs_];
        for (int i=0; i<vibs_; ++i) {
            states_[i] = new gVector2D[ryds_];
            for (int j=0; j<ryds_; ++j) {
                states_[i][j] = gVector2D(g);
            }
        }

        assert(!generators_);
        generators_ = new gVector*[ryds_];
        for (int i=0; i<ryds_; ++i) {
            generators_[i] = new gVector[g.get_ysize()];
        }
    }
    void append_generators(zEigenSystem& e, int yi, femGrid& gx)
    {
        for (int i=0; i<ryds_; ++i){
            generators_[i][yi] = gVector(gx);
            e.eigen_vector(generators_[i][yi].body(), i);
        }
    }
    void build_states(zEigenSystem& e, int ri, femGrid& gy)
    {
        if (ri < ryds_) {
            for (int vi=0; vi<vibs_; ++vi) {
                const femGrid2D& g = states_[vi][ri].get_grid();
                gVector aux(g.get_ygrid());
                e.eigen_vector(aux.body(), vi);
                for (int j=0; j<g.get_ysize(); ++j) {
                    for (int i=0; i<g.get_xsize(); ++i)
                         states_[vi][ri][j*g.get_xsize() + i] = generators_[ri][j][i] * aux[j];
                }
            }
        }

    }
    void init_outputs(string folder)
    {
        assert(outputs_ == NULL);
        outputs_ = new ofstream[vibs_];
        for (int i=0; i<vibs_; ++i) {
            char name[50];
            sprintf(name, "QSB/population_v%d.txt", i);
            outputs_[i].open((folder + name).c_str());
        }
    }
    void print_populations(const gVector2D& psi, double time)
    {
        assert(states_);
        assert(outputs_);

        for (int i=0; i<vibs_; ++i) {
            outputs_[i] << time << "\t";
            for (int j=0; j<ryds_; ++j) {
                outputs_[i] << real(states_[i][j] * psi) << "\t";
            }
            outputs_[i] << endl;
        }
    }
    bool save_states(string folder)
    {
        assert(states_);

        bool stat = true;
        char name[50];
        for (int i=0; i<vibs_; ++i) {
            for(int j=0; j<ryds_; ++j) {
                sprintf(name, "/QBS/phi_v%d_r%d.qbin", i, j);
                stat = stat && states_[i][j].save_binary(folder + name);
                if (!stat) {
                    cout << "ERROR ON SAVE: " << folder + name << endl;
                }
            }
        }
        return stat;
    }
    bool read_states(string folder)
    {
        assert(states_ == NULL);

        bool stat = true;
        char name[50];
        states_ = new gVector2D*[vibs_];
        for (int i=0; i<vibs_; ++i) {
            states_[i] = new gVector2D[ryds_];
            for (int j=0; j<ryds_; ++j) {
                sprintf(name, "/QBS/phi_v%d_r%d.qbin", i, j);
                stat = (stat && states_[i][j].read_binary(folder + name));
            }
        }
        if (!stat)
            clean_states();
        return stat;
    }
    void clean_states()
    {
        if (states_) {
            for (int i=0; i<vibs_; ++i) {
                if (states_[i]) {
                    delete[] states_[i];
                }
            }
            delete[] states_;
        }
        states_ = NULL;
    }
    ~QBStates()
    {
        clean_states();
        if (outputs_) {
            for (int i=0; i<vibs_; ++i) {
                outputs_[i].close();
            }
            delete[] outputs_;
        }
        if (generators_){
            for (int i=0; i<vibs_; ++i)
                delete[] generators_[i];
            delete[] generators_;
        }
    }
};

/*

    Basic concept of double two-dimensional model connected with coupling

*/

using std::cout;
using std::endl;
using std::exp;
using std::abs;
using std::sqrt;
using std::string;

def_comp h2p_v_zero(const def_comp& x)
{
    //  V(R) = V0 * [exp(-2*alpha*(R-R0)) - 2*exp(-alpha*(R-R0))] + Q

    def_float V0 = 0.1027;
    //def_float Q =  -0.5;
    def_float Q =  -0.0;
    def_float R0 = 2.0;
    def_float alpha = 0.69;

    return V0 * ( exp( - 2.0 * alpha * (x - R0)) - 2.0 * exp( -alpha * (x - R0)) ) + Q;
}

def_comp h2sigma_eff(const def_comp& x, const def_comp&y)
{
    def_float l = 1.0;

    return h2p_v_zero(y) + QSCAT::p_sigma_potential(x, y) + l*(l+1.0)/ (2.0 * x*x ) - 1.0/x;
}

def_comp h2sigma_el_eff(const def_comp& x)
{
    def_comp l = 1.0;

    return l*(l+1.0)/ (2.0 * std::pow(x,2) ) - 1.0/x + QSCAT::p_sigma_potential(x, 30.0);
}

void compute_vibration_thresholds(bool sp_switch, string input, QBStates *phi=NULL)
{
    cout << "vibrations" << endl;
    int levels = 50;
    double l = 1.0;
    double e_shift;

 // IO
    pjvalue cfg = read_json_file(input);

 // Model
    pjvalue& mdc = cfg["model"];
    def_float mu = mdc["reduced_mass"].asDouble();
    string folder = mdc["folder"].asString();

 // grids:
    femGrid gx = grid_from_parameters(cfg["grids"]["electronic"]);
    femGrid gy = grid_from_parameters(cfg["grids"]["nuclear"]);

    // first determine energy shift
    gVector v_0(gy);
    for (int i=0; i<gy.nb(); ++i)
        v_0.f(h2p_v_zero(gy.xr(i)), i);

    v_0.save_binary(folder + "potentials/vZero.qbin");
    zOperatorF H_0(gy);
    H_0.add_kinetic_term(mu);
    H_0 += v_0;
    zEigenSystem eSys = H_0.eigen_system();
    e_shift = real(eSys.eigen_value(0));

    for (int i=0; i<10; ++i)
        cout << i << "-th state energy " << real(eSys.eigen_value(i)) << endl;
    cout << "V_0 ground state energy: " << e_shift << endl;

    if (false)
    {
        def_comp y = gy.x(gy.nr());
        cout << "Computing electronic spectra for y = " << y << endl;
        gVector pot_e(gx);
        for (int j=0; j<gx.nb(); ++j){
            def_comp x = gx.xr(j);
            def_comp coul = - 1.0 / x + l * (l+1.0) / (2.0 * pow(x, 2));
            //double ppi  = real ( MODEL_2D::p_pi_potential(x,y) );
            def_comp psgm = QSCAT::p_sigma_potential(x,y);
            //pot.f(coul + ppi, j);
            pot_e.f(real(coul + psgm), j);
        }

        zOperatorF H(gx);
        H.add_kinetic_term(1.0);
        H += pot_e;

        zEigenSystem ES = H.eigen_system();
        for (int k=0; k<levels; ++k)
            cout << -0.5/pow(k+2,2) << ": " << ES.eigen_value(k+1) << endl;

    }
    //exit(0);

    //if (true) {

    gVector *pot = new gVector[levels];
    for (int i=0; i<levels; ++i)
        pot[i] = gVector(gy);

    for (int i=0; i<gy.nb(); ++i) {
        def_comp y = gy.x(i);
        cout << "Computing electronic spectra for y = " << y << endl;
        gVector pot_e(gx);
        for (int j=0; j<gx.nb(); ++j){
            def_comp x = gx.x(j);
            def_comp coul = - 1.0 / x + l * (l+1.0) / (2.0 * pow(x, 2));
            //double ppi  = real ( MODEL_2D::p_pi_potential(x,y) );
            def_comp psgm = QSCAT::p_sigma_potential(x,y);
            //pot.f(coul + ppi, j);
            pot_e.f( real(coul + psgm), j);     // HACK want to see
        }

        zOperatorF H(gx);
        H.add_kinetic_term(1.0);
        H += pot_e;

        zEigenSystem ES = H.eigen_system();
        for (int k=0; k<levels; ++k) {
            pot[k].f(real(h2p_v_zero(y)) + ES.eigen_value(k), i);
        }
        if (phi)
            phi->append_generators(ES, i, gx);
    }

    for (int i=0; i<levels; ++i) {
        char name[20];
        sprintf(name, "RydPot_%d.qbin", i);
        pot[i].save_binary(folder + "potentials/" + name);

        if (phi) {
            zOperatorF Hr(gy);
            Hr.add_kinetic_term(mu);
            Hr += pot[i];

            zEigenSystem es = Hr.eigen_system();
            phi->build_states(es, i, gy);
        }
    }

    //}

/*
    QSCAT::parametersMultiGrid gp("input/coupled/grids.txt");
    femGrid gx(gp.gp[0]); // electronic
    femGrid gy(gp.gp[1]); // electronic

    zEigenSystem *ES = new zEigenSystem[levels];
    for (int k=0; k<levels; ++k) {
        zOperatorF H(gy);
        H.add_kinetic_term(mu);
        H += pot[k];
        ES[k] = H.eigen_system();
    }

    FILE* file;
    fopen_s(&file, "output/H2+/vibrational_thresholds.asc", "w");
    fprintf(file,"#Energies of H2+ molecular vibration thresholds\n");
    for (int k=0; k<levels; ++k) {
        for (int i=0; i<gy.nb(); ++i) {
            fprintf(file, "%.12E", real(ES[k].eigen_value(i) - e_shift) );
            if (i<gy.nb()-1) fprintf(file, "\t");
        }
        fprintf(file, "\n");
    }
    fclose(file);

    fopen_s(&file, "output/H2+/neutral_energies.asc", "w");
    fprintf(file,"#Energies of H2+ + e- molecular energy levels\n");
    for (int i=0; i<gy.nb(); ++i) {
        fprintf(file, "%.12E", gy.xr(i) );
        for (int k=0; k<levels; ++k) {
            fprintf(file, "\t");
            fprintf(file, "%.12E", real(pot[k].f(i)) );
        }
        fprintf(file, "\n");
    }
    fclose(file);
*/
}

void compute_2d_noncoupled_model(bool sp_switch, string input, QBStates *phi=NULL)
{
 // IO
    pjvalue cfg = read_json_file(input);

 // Model
    pjvalue& mdc = cfg["model"];
    def_float mu = mdc["reduced_mass"].asDouble();
    string folder = mdc["folder"].asString();

 // grids:
    femGrid gx = grid_from_parameters(cfg["grids"]["electronic"]);
    femGrid gy = grid_from_parameters(cfg["grids"]["nuclear"]);

    femGrid2D g(gx, gy);

    if (phi) {
        if (! phi->read_states(folder)) {
            phi->init_states(g);
            compute_vibration_thresholds(false, input, phi);
            phi->save_states(folder);
        }
        phi->init_outputs(folder);
    }

 // initial state
    pjvalue& icfg = cfg["initial_state"]["wavepacket"];

    def_float gaussX = icfg["position"].asDouble();
    def_float gaussP = icfg["impulse"].asDouble();
    def_float gaussS = icfg["sigma"].asDouble();

 // cross section parametrisation:

    pjvalue& csc = cfg["cross_sections"];
    pjarray& rng = csc["range"].asArray();
    //QSCAT::parameters2D m2dp("input/coupled/model.txt");
    //QSCAT::parametersEvolution& ep(m2dp.evol_par);
    size_t esize = (rng[1].asDouble() - rng[0].asDouble())/csc["dE"].asDouble();
    dVector erg(esize, rng[0].asDouble(), rng[1].asDouble(), false);

 // evolution
    pjvalue& evc = cfg["evolution"];
    def_float dt = evc["time_step"].asDouble();
    size_t innerSpan = evc["loop_steps"].asUInt();
    size_t CNorder = evc["order"].asUInt();

 // Hamiltonian:
    gVector2D pot(g);
    pot.fill(0.0);
    for (int i=0; i<gy.nb(); ++i){
        for (int j=0; j<gx.nb(); ++j){
            pot.f( h2sigma_eff(gx.xr(j), gy.xr(i)), i*gx.nb() + j);
        }
    }

    zOperator2D H(g);
    H.set_kinetic_term(1.0, mu);
    H += pot;

    pot.save_binary(folder + "potential_H2p_sigma.qbin");

    //pot.save(folder + "v_eff.dat");
    //pot.save_equidistant(folder + "v_eff_eq.dat", 100, 1e-10, 10., 100, 1e-10, 10. );

 // Building wavefunction
  // nuclear coordinate:
    gVector Y(gy);
    Y.fill(0.0);
    for (int i=0; i<gy.nb(); ++i) {
        Y.f(h2p_v_zero(gy.x(i)), i);
    }
    zOperatorF Hm(gy);
    Hm.add_kinetic_term(mu);
    Hm += Y;
    zEigenSystem eSys_m = Hm.eigen_system();
    eSys_m.eigen_vector(Y.body(), 0);
    def_float init_erg = real(eSys_m.eigen_value(0));
    cout << "initial state energy: " << init_erg << endl;

    ofstream ergFile;
    ergFile.open((folder + "vibrational_energies.txt").c_str());
    for (int i=0; i<20; ++i) {
        ergFile << real(eSys_m.eigen_value(i)) << endl;
    }
    ergFile.close();
  // electronic coordinate
    gVector X(gx);
    X.fill(0.0);
    for (int i=0; i<gx.nr(); ++i) {
        X.f(zGaussian(gx.xr(i), gaussX, gaussS, gaussP), i);
        //X[i] = psi_out[i];
    }
  // 2D:
    gVector2D psi(g, X, Y);
    psi.save_binary(folder + "WF/psi_t_0.qbin0");

  // fourier coefficients
    zVector ifc(esize);
    for (int i=0; i<esize; ++i) {
        gVector coul(gx);
        coul.fill(0.0);
        def_comp k = sqrt(2.*erg[i]);
        for (int j=0; j<gx.nr(); ++j) {
            coul.f(QSCAT::coulomb::sF_en(gx.x(j), k, -1., 1., 1), j);
        }
        ifc[i] = X * coul;
    }
    ifc.save(erg, (folder + "ifc.dat").c_str());

 // DA testfunction
  // eigensystem
    gVector aux(gx);
    for (int i=0; i<gx.nb(); ++i) {
        aux.f( h2sigma_el_eff(gx.x(i)) , i);
    }
    zOperatorF He(gx);
    He.add_kinetic_term(1.0);
    He += aux;
    zEigenSystem eSys_e = He.eigen_system();

  // testfunction body
    cout << "Building testfunctions " << endl;
    pjvalue& drc = cfg["cross_sections"]["testfunctions"]["dissociative_recombination"];
    QSCAT::TestFunction2d *tannor = new QSCAT::TestFunction2d[drc["channels"].asInt()];
    for (int i=0; i<drc["channels"].asInt(); ++i) {
        cout << "testfunction " << i << endl;
        if (i==0) {
            drc["wavpacket"]["impulse"] = 60.0;  // 52
        } else {
            drc["wavepacket"]["impulse"] = 12.0;
        }
        tannor[i] = QSCAT::TestFunction2d( drc, g, i, eSys_e, init_erg, -1., mu, 1., erg);
        stringstream ss;
        ss << folder << "cf_" << i << ".bin";
        tannor[i].set_output(ss.str().c_str());
        tannor[i] << psi;
    }
  // s matrix
    //MODEL_2D::SMatrix S(ep.steps, innerSpan, 2, tfp.channels, 't');
    zMatrix S(esize, drc["channels"].asInt());

 // Evolution operator
    QSCAT::CrankNicolson2D CN(CNorder, dt, H);

 // evolution:
    def_float t = 0.0;

    bool saving = cfg["evolution"].isMember("saving_step");
    int sstep = (!saving)? 0 : cfg["evolution"]["saving_step"].asUInt();

    ofstream norm;
    norm.open((folder + "/normalization.txt").c_str());
    norm << 0.0 << "\t" << real(psi*psi) << endl;

    for (int loops=0; loops < 1000000; ++loops) {
        for (int inner=0; inner<innerSpan; ++inner) {
            CN.one_step(psi);
            cout << ".";
            cout.flush();
            for (int i=0; i<drc["channels"].asInt(); ++i) {
                tannor[i] << psi;
            }
            if (saving && inner%sstep==0) {
                stringstream ss;
                ss << folder << "WF/psi_t_" << (t + inner*dt) << ".qbin";
                psi.save_binary(ss.str());
            }
            norm << t + inner*dt << "\t" << real(psi*psi) << endl;
            if (phi)
                phi->print_populations(psi, t + inner * dt);
        }
        cout << "loop finished" << endl;
        for (int i=0; i<drc["channels"].asInt(); ++i) {
            tannor[i].contribution(S, i, t,  dt, ifc);
        }
        t += innerSpan * dt;
        S.save(erg, (folder + "SDA.dat").c_str());
        cout << psi*psi << " at t=" << t << endl;
        psi.save_binary(folder + "psi.qbin");
    }
    norm.close();
}

void compute_quantum_defects(bool sp_switch)
{
/*
    QSCAT::parametersMultiGrid gp("input/coupled/grids.txt");
    femGrid g(gp.gp[0]); // electronic

    int samples = 30;
    double R = 10.0;
    int nb = g.nb();
    double mu = 1.0;
    int nMax = 30;
    if (nMax>nb) nMax=nb;

    double l = 1.0;

    gVector pot(g);

    FILE* file;
    fopen_s(&file, "output/H2+/V_el_spectrum.asc", "w");
    fprintf(file,"#Energies of H2+ elcetronic p-wave pi sates energies with respect to R\n");

    FILE* fpot;
    fopen_s(&fpot, "output/H2+/V_0+El.asc", "w");
    fprintf(fpot,"#Energies of H2+ elcetronic p-wave pi sates energies with respect to R\n");

    FILE* defects;
    fopen_s(&defects, "output/H2+/qdf.asc", "w");
    fprintf(defects,"#H2+ quantum defect functions respect to R\n");

    // labels
    fprintf(defects, "#R");
    for (int j=0; j<nMax; ++j) {
        fprintf(defects, "n=%d\t", j+int(l)+1);
    }
    fprintf(defects,"\n");
    // data
    for (int i=0; i<samples; ++i) {
        // Build potential
        double y = R * i / (samples-1);
        cout << "Computing QDF at R = " << y << " ... ";
        cout.flush();
        for (int j=0; j<nb; ++j) {
            double x = g.xr(j);
            double coul = - 1.0 / x + l * (l+1.0) / (2.0 * pow(x, 2));
            if (sp_switch) {
                double ppi  = real ( QSCAT::p_pi_potential(x,y) );
                pot.f(coul + ppi, j);
            } else {
                double psgm = real( QSCAT::p_sigma_potential(x,y));
                pot.f(coul + psgm, j);
            }
        }

        char name[50];
        sprintf(name, "output/H2+/pot_R_%.3f.asc", y);
        pot.save(name);

        zOperatorF H(g);
        H.add_kinetic_term(mu);
        H += pot;

        zEigenSystem ES = H.eigen_system();
        cout << (1.0+l) - 1.0 / sqrt( - 2.0 * real(ES.eigen_value(0))) << endl;

        fprintf(file, "%.12E", i * R / (samples-1) );
        fprintf(fpot, "%.12E", i * R / (samples-1) );
        fprintf(fpot, "\t%.12E", real(h2p_v_zero(y)) );
        fprintf(defects, "%.12E", i * R / (samples-1) );
        for (int j=0; j<nMax; ++j) {
            fprintf(file, "\t%.12E", real(ES.eigen_value(j)) );
            fprintf(fpot, "\t%.12E", real(ES.eigen_value(j) + h2p_v_zero(y)) );
            fprintf(defects, "\t%.12E", (1.0+j+l) - 1.0 / sqrt( - 2.0 * real(ES.eigen_value(j))) );
        }
        fprintf(file, "\n");
        fprintf(fpot, "\n");
        fprintf(defects, "\n");

    }
    fclose(file);
    fclose(fpot);
    fclose(defects);
*/
}

// Equation (H - E) psi = phi
// for given E solves psi
void solveDrivenSchr(zOperator2D H, def_comp E, gVector2D& psi)
{
    H += -E;

    //gVector2D out(psi);
    assert(psi.init());
    H.LU_factorize();
    H.LU_back_substitution(psi);
    assert(psi.init());
}

void time_independent_solution(string input)
{
    pjvalue cfg = read_json_file(input);

 // Model
    pjvalue& mdc = cfg["model"];
    def_float mu = mdc["reduced_mass"].asDouble();
    string folder = mdc["folder"].asString();

 // grids:
    femGrid gx = grid_from_parameters(cfg["grids"]["electronic"]);
    femGrid gy = grid_from_parameters(cfg["grids"]["nuclear"]);

    femGrid2D g(gx, gy);

 // Hamiltonian:
    gVector2D pot(g);
    pot.fill(0.0);
    for (int i=0; i<gy.nb(); ++i){
        for (int j=0; j<gx.nb(); ++j){
            pot.f( h2sigma_eff(abs(gx.x(j)), abs(gy.xr(i))), i*gx.nb() + j);
        }
    }
    zOperator2D H(g);
    H.set_kinetic_term(1.0, mu);
    H += pot;

 // Building wavefunction
    gVector Y(gy);
    Y.fill(0.0);
    for (int i=0; i<gy.nb(); ++i) {
        Y.f(h2p_v_zero(gy.x(i)), i);
    }
    zOperatorF Hm(gy);
    Hm.add_kinetic_term(mu);
    Hm += Y;
    zEigenSystem eSys_m = Hm.eigen_system();
    eSys_m.eigen_vector(Y.body(), 0);
    def_comp init_erg  = eSys_m.eigen_value(0);

    pjvalue& csc = cfg["cross_sections"];
    pjarray& rng = csc["range"].asArray();
    size_t esize = (rng[1].asDouble() - rng[0].asDouble())/csc["dE"].asDouble();
    dVector erg(esize, rng[0].asDouble(), rng[1].asDouble(), false);

 // VDR ... etc.
    gVector pEl(gx);
    //pEl.fill(0.0);
    for (int i=0; i<gx.nb(); ++i) {
        pEl.f(h2sigma_el_eff(gx.x(i)), i);
    }
    zOperatorF Hy(gx);
    Hy.add_kinetic_term(1.0);
    Hy += pEl;
    zEigenSystem eSys_e = Hy.eigen_system();

    //gVector2D VDR(pot);
    zVector VDR(gx.nb() * gy.nb());
    zVector Vint(gx.nb() * gy.nb());
    def_comp val;
    for (blas_int i=0; i<gx.nb(); ++i) {
        for (blas_int j=0; j<gy.nb(); ++j) {
            Vint[gx.nb()*j + i] = QSCAT::p_sigma_potential(gx.x(i), gy.x(j));
            VDR[gx.nb()*j + i] = h2p_v_zero(gy.x(j)) + QSCAT::p_sigma_potential(gx.x(i), gy.x(j)) - QSCAT::p_sigma_potential(gx.x(i), gy.xr(gy.nr()-1));
        }
    }

    int drc = 3;
    gVector chi[drc];
    def_comp edr[drc];
    cout << init_erg;
    for (int i=0; i<drc; ++i) {
        chi[i] = gVector(gx);
        eSys_e.eigen_vector(chi[i].body(), i);
        edr[i] = eSys_e.eigen_value(i);
        cout << " " <<  edr[i];
    }
    cout << endl;

    zMatrix T(erg.get_size(), drc);

    ofstream outfile;
    outfile.open("Tmatrix.txt");
    for (int i=0; i<erg.get_size(); ++i) {
        def_comp e = erg[i];
        outfile << real(e);
        def_comp k = sqrt(2.0 * erg[i]);
        gVector X(gx);
        for (int j=0; j<gx.nb(); ++j) {
            X.f((j<=gx.nr())? QSCAT::coulomb::sF_en(gx.x(j), k, -1., 1., 1.0) : 0.0, j);
        }
        gVector2D psi(g, X, Y);
        gVector2D psi_sc = psi.copy();
        zOperator2D HmE = H.copy();
        HmE += -(e + init_erg);

        psi_sc.body().element_wise_multiplication(Vint);
        HmE.LU_back_substitution(psi_sc);
        psi -= psi_sc;

        psi.body().element_wise_multiplication(VDR);

        cout << erg[i];
        for (int j=0; j<drc; ++j) {
            def_comp E_out = erg[i] + init_erg - edr[j];
            if (real(E_out) < 0) {
                T(i,j) = 0.0;
                cout << " " << T(i,j);
                outfile << " " << real(T(i,j)) << " " << imag(T(i,j));
                continue;
            }
            def_comp p = sqrt( 2.0 * mu * (E_out) );
            gVector FY(gy);
            FY.fill(0.0);
            for (int jj=0; jj<gy.nr(); ++jj) {
                FY.f((jj<gy.nr())? sphBesselJEn(gy.x(jj), p, mu, 0) : 0, jj );
            }
            gVector2D psi_out(g, chi[j], FY);
            T(i,j) = psi_out * psi;
            cout << " " << T(i,j);
            outfile << " " << real(T(i,j)) << " " << imag(T(i,j));
        }
        cout << endl;
        outfile << endl;
    }
    outfile.close();
    T.save_binary("Tmatrix.qbin");

}

int main(int argc, char **argv)
{

    ComputationMode mode = MODE_FULL_WITH_POPULATIONS;
    if (argc<2) {
        cout << "Usage: " << argv[0] << " input.json" << endl;
        return 0;
    }

    //if (false) {
    //    time_independent_solution(argv[1]);
    //    return 0;
    //}

    switch (mode) {
        case MODE_FULL:
            compute_2d_noncoupled_model(false, argv[1]);
            break;

        case MODE_QDF:
            compute_quantum_defects(false);
            break;

        case MODE_ENERGY_THRESHOLD:
            compute_vibration_thresholds(false, argv[1]);
            break;

        case MODE_FULL_WITH_POPULATIONS:
            QBStates phi(10,30);
            compute_2d_noncoupled_model(false, argv[1], &phi);
            break;
    }


    return 0;
}
