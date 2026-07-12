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

#include "arrays.h"
#include "input.h"
#include "potentials.h"
#include "fem_dvr_ecs.h"
#include "fem_dvr_ecs_2D.h"
#include "interface.h"
#include "module_LCP.h"
//#include "module_NRM.h"
#include "model_2D.h"

enum ComputationMode {
    MODE_FULL,
    MODE_QDF,
    MODE_ENERGY_THRESHOLD
};


/*

    Basic concept of double two-dimensional model connected with coupling

*/

using std::cout;
using std::endl;
using std::exp;
using std::abs;
using std::sqrt;

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

    return h2p_v_zero(y) + MODEL_2D::p_sigma_potential(x, y) + l*(l+1.0)/ (2.0 * x*x ) - 1.0/x;
}

def_comp h2sigma_el_eff(const def_comp& x)
{
    def_comp l = 1.0;

    return l*(l+1.0)/ (2.0 * std::pow(x,2) ) - 1.0/x + MODEL_2D::p_sigma_potential(x, 10.0);
}

void compute_2d_noncoupled_model(bool sp_switch)
{
	cout << "++++++++Running version: time_indep with multiple rydbergs, outgoing Bessel++++++++" << endl;
 //constants
    def_float mu = 918.25;	// nuclear mass
	def_float me = 1.0;		// electronic mass
	def_float Q_e = -1.0;	// electronic charge
	int vi = 0;				// vibrational number
	ifstream rydbergs ("input/H2+/rydbergs.txt");
	int ri_min;				// lowest Rydberg state number
	rydbergs >> ri_min;
	int ri_max;				// highest Rydberg state number
	rydbergs >> ri_max;
	def_float Pi=4.0*atan(1.0); // Pi

 // grids:
    parametersMultiGrid gp("input/coupled/grids.txt");
    femGrid gx(gp.gp[0]); // electronic
    femGrid gy(gp.gp[1]); // nuclear
    femGrid gm(gp.gp[2]); // electronic extended real region for Moeller evolution
    femGrid2D g(gx, gy);
	def_float el_inf = gx.xr(gx.nr()-1); // electronic grid real part endpoint
	def_float nu_inf = gy.xr(gy.nr()-1); // nuclear grid real part endpoint

 // Hamiltonian:
    gVector2D pot(g);
    for (int i=0; i<gy.nb(); ++i){
        for (int j=0; j<gx.nb(); ++j){
            pot.f( h2sigma_eff(gx.xr(j), gy.xr(i)), i*gx.nb() + j);	// FIXME
        }
    }
    zOperator2D H(g);
    H.set_kinetic_term(me, mu);
    H += pot;
    pot.save("output/H2+/v_eff.dat");
    pot.save_equidistant("output/H2+/v_eff_eq.dat", 100, 1e-10, 10., 100, 1e-10, 10. );

 // Building wavefunctions
	// nuclear coordinate for initial vibrational state:
	gVector Y_in(gy);
    for (int i=0; i<gy.nb(); ++i) {
        Y_in.f(h2p_v_zero(gy.x(i)), i);
    }
    zOperatorF Hm(gy);
    Hm.add_kinetic_term(mu);
    Hm += Y_in;
    zEigenSystem eSys_m = Hm.eigen_system();
	def_float E_vib = real(eSys_m.eigen_value(vi));
	cout << "initial vibrational state energy: " << E_vib << endl;
    eSys_m.eigen_vector(Y_in.body(), vi);
    for (int i=gy.nr(); i<gy.nb(); ++i) {
        Y_in.f(0.0, i);
    }
	gVector X_in(gx);

	// electronic coordinate for outgoing rydberg state
	gVector X_out(gx);
	for (int i=0; i<gx.nb(); ++i) {
        X_out.f(h2sigma_eff(gx.x(i),nu_inf), i);
	}
    zOperatorF H_el(gx);
    H_el.add_kinetic_term(me);
    H_el += X_out;
    zEigenSystem eSys_el = H_el.eigen_system();
	eSys_el.eigen_vector(X_out.body(), ri_max);		// We'll use this line several more times in the cycle
    def_float E_ryd = real(eSys_el.eigen_value(ri_max));
    cout << "outgoing rydberg state maximum energy: " << E_ryd << endl;
    for (int i=gx.nr(); i<gx.nb(); ++i) {
       X_out.f(0.0, i);
    }
	gVector Y_out(gy);

	// threshold energy
	def_float E_th_max = E_ryd - E_vib;
	def_float E_th = real(eSys_el.eigen_value(ri_min)) - E_vib;	// minimal threshold energy
	if (E_th_max > 0) {
		cout << "the highest threshold energy is: " << E_th_max << endl;
	}
	else {
		cout << "the highest rydberg state energy is less than the vibrational state energy" << endl;
	}

	// interaction potential
	zVector V_int(g.get_size());
	for (int i=0; i<gy.nb(); ++i){
		for (int j=0; j<gx.nb(); ++j){
			V_int[i*gx.nb() + j] = MODEL_2D::p_sigma_potential(gx.x(j), gy.x(i));
		}
	}

	// DR potential: V_DR = V0 + V_int - V_int(r,R->inf)
	zVector V_DR(V_int);
	for (int i=0; i<gy.nb(); ++i){
        for (int j=0; j<gx.nb(); ++j){
            V_DR[i*gx.nb() + j] += h2p_v_zero(gy.x(i)) - MODEL_2D::p_sigma_potential(gx.x(j), nu_inf);
        }
    }
	
	ofstream outfile;
	outfile.open("output/H2+/sigma.txt");
	ofstream Tmatrix;
	Tmatrix.open("output/H2+/tmatrix.txt");

  // Cross section calculation
	def_float sigma = 0.0;
	//def_float sigma_2 = 0.0;
	def_float E_el = 0.0;
	ifstream energ ("input/H2+/energy_grid.txt");
	//setting P_out to lowest rydberg
	eSys_el.eigen_vector(X_out.body(), ri_min);
    E_ryd = real(eSys_el.eigen_value(1));

	while (energ.good()){
		energ >> E_el;
		def_float k_el = sqrt(2.0*me*E_el);

		// electronic coordinate for initial state
		for (int i=0; i<=gx.nr(); ++i) {
			X_in.f(coulomb::sF_en(gx.x(i),k_el,Q_e,me,1.0), i);
		}
		for (int i=gx.nr(); i<gx.nb(); ++i) {
			X_in.f(0.0, i);
		}

		gVector2D psi(g, X_in, Y_in); //the complete initial state

	  // (E-H) P_scattered = V_int P_in
		// (H-E)
		zOperator2D HwE(H);
		def_float E = E_el + E_vib;
		HwE += (-E);

		// V_int * P_in
		gVector2D psi_sc(psi); // this is not yet P_scattered
		psi_sc.body().element_wise_multiplication(V_int);
		// P_scattered (times -1)
		HwE.LU_back_substitution(psi_sc);
		//Psi+ = Psi + P_scattered
		psi -= psi_sc;

		outfile << E_el;
		Tmatrix << E_el;

		for (int ri=ri_min; ri<=ri_max; ++ri){
			if ( E_el > (real(eSys_el.eigen_value(ri)) - E_vib) ){
				gVector2D psi_2(psi);
				def_float E_DR = E_el + E_vib - real(eSys_el.eigen_value(ri));
				def_float K_DR = sqrt(2.0*mu*E_DR);
				eSys_el.eigen_vector(X_out.body(), ri);
				// nuclear coordinate for outgoing state
				for (int i=0; i<=gy.nr(); ++i) {
					Y_out.f(bessel::s_jEn(gy.x(i),K_DR,mu,0), i);
				}
				for (int i=gy.nr(); i<gy.nb(); ++i) {
					if (gy.xr(i)>1000.0) {Y_out.f(0.0, i);}
					else {Y_out.f(bessel::s_jEn(gy.x(i),K_DR,mu,0), i);}
				}
				gVector2D psi_out(g, X_out, Y_out);
				//Psi_out*(V_DR * Psi+)
				psi_2.body().element_wise_multiplication(V_DR);
				def_comp T_DR = psi_out * psi_2;
				sigma = 4.0*Pi*Pi*Pi*abs(T_DR*T_DR)/(k_el*k_el);
				Tmatrix << " " << real(T_DR) << " " << imag(T_DR);
			}
			else{ 
				sigma=0.0;
				Tmatrix << " " << 0 << " " << 0;
			}
			outfile << " " << sigma;
		}

		outfile << endl;
		Tmatrix << endl;
	}


    exit(0);

}

void compute_quantum_defects(bool sp_switch)
{
    parametersMultiGrid gp("input/coupled/grids.txt");
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
                double ppi  = real ( MODEL_2D::p_pi_potential(x,y) );
                pot.f(coul + ppi, j);
            } else {
                double psgm = real(MODEL_2D::p_sigma_potential(x,y));
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
}

void compute_vibration_thresholds(bool sp_switch)
{
    int levels = 40;
    double mu = 918.25;
    double l = 1.0;
    double e_shift;

    parametersMultiGrid gp("input/coupled/grids.txt");
    femGrid gx(gp.gp[0]); // electronic
    femGrid gy(gp.gp[1]); // electronic

    // first determine energy shift
    gVector v_0(gy);
    for (int i=0; i<gy.nb(); ++i)
        v_0.f(h2p_v_zero(gy.xr(i)), i);
    zOperatorF H_0(gy);
    H_0.add_kinetic_term(mu);
    H_0 += v_0;
    zEigenSystem eSys = H_0.eigen_system();
    e_shift = real(eSys.eigen_value(0));

    cout << "V_0 ground state energy: " << e_shift << endl;

    gVector *pot = new gVector[levels];
    for (int i=0; i<levels; ++i)
        pot[i] = gVector(gy);

    for (int i=0; i<gy.nb(); ++i) {
        double y = gy.xr(i);
        cout << "Computing electronic spectra for y = " << y << endl;
        gVector pot_e(gx);
        for (int j=0; j<gx.nb(); ++j){
            double x = gx.xr(j);
            double coul = - 1.0 / x + l * (l+1.0) / (2.0 * pow(x, 2));
            //double ppi  = real ( MODEL_2D::p_pi_potential(x,y) );
            double psgm = real(MODEL_2D::p_sigma_potential(x,y));
            //pot.f(coul + ppi, j);
            pot_e.f(coul + psgm, j);
        }

        zOperatorF H(gx);
        H.add_kinetic_term(1.0);
        H += pot_e;

        zEigenSystem ES = H.eigen_system();
        for (int k=0; k<levels; ++k)
            pot[k].f(h2p_v_zero(y) + ES.eigen_value(k), i);
    }


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
}

int main(int argc, char **argv) {

ComputationMode mode = MODE_FULL;

switch (mode) {
    case MODE_FULL:
        compute_2d_noncoupled_model(false);
        break;

    case MODE_QDF:
        compute_quantum_defects(false);
        break;

    case MODE_ENERGY_THRESHOLD:
        compute_vibration_thresholds(false);
        break;
}
        //model_2D::CoupledModel2D();

        //std::string input_folder;

        //bool isLCP = true;
        //bool is2D  = true;

        //if (argc>1){
        //	input_folder += "input/";
        //	input_folder += argv[1];
        //	input_folder += "/";
        //} else {
        //	input_folder += "input/coupled/";
        //}

    // Parameters input
        //parameters::multi_grid<def_float> gp( (input_folder + "grids.txt").c_str() );



    // 2D declarations

        //LCP::Model_LCP<def_float, def_comp> LCP(m2dp,gp,LCPp);

        //grid_vector_2D<def_float, def_comp> PhiD = LCP.MakePhiD();
        //model_2D::Model M2D(m2dp,gp.gp[0],gp.gp[1], &PhiD);

        //if (M2D.ReadBinary((m2dp.folder + "frame.M2D").c_str())){
        //	std::cout << "The previous run was successfully loaded." << std::endl;
        //}



        //for (int i=0;i<m2dp.evol_par.tcutoff/(m2dp.evol_par.dt*m2dp.evol_par.loop);++i){
        //	if (is2D) M2D.Multistep();
        //	if (isLCP) LCP.Multistep();
        //	//M2D.SaveBinary((m2dp.folder + "frame.M2D").c_str());
        //}

return 0;
}
