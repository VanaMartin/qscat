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
#include "arrays.h"
#include "input.h"
#include "potentials.h"
#include "fem_dvr_ecs.h"
#include "fem_dvr_ecs_2D.h"
#include "interface.h"
#include "module_NRM.h"

/*

	The Non-local Resonance Model (NRM) module
		- at the moment the time independent mode is under construction
		TODO:
			- implement the initialization of fem grid accordingly to the FORTRAN image
			- Compare the interim results for all known steps
			- list the known issues in internal modules 

*/

using namespace std;
using namespace NRM;

// building the two dimensional discrete state wave function 
void NRM::MakePhiD_spec(gVector2D & phi, parameters2D & mp, femGrid2D & g)
{
	femGrid * pge = g.Gx();			// Pointer to the electonic grid
	femGrid * pgn = g.Gy();			// Pointer to the nuclear grid
	gVector potential(*pge);		// Potential
	gVector phi_b(*pge);			// Bound state 

	for (int i=0; i<g.NbY(); ++i){
		potential.Make_fxc(potentials::Lambda_Spec_R,pgn->Xr(i),mp);
		//Hamiltonian_Full<def_float,def_comp> H(*pge, potential, 1.0); 
		//H.SetEigen();
		//H.eigs->GetState(phi_b.a,0);
		zOperatorF H = makeHamiltonianF(*pge, potential, 1.0);
        zEigSystem eSys = H.EigenSystem();
        eSys.GetState(phi_b.a, 0);
        std::cout << "Eigenvalue E = " << eSys.energies[0] << " at R = " << pgn->Xr(i) << std::endl;
		phi.WriteX(phi_b,i);
	}
}

double NRM::MakePhiD_const(gVector2D & phi, parameters2D & mp, femGrid2D & g)
{
	femGrid * pge = g.Gx();			// Pointer to the electonic grid
	//femGrid * pgn = g.Gy();			// Pointer to the nuclear grid
	gVector potential(*pge);		// Potential
	gVector phi_b(*pge);			// Bound state 
	potential.Fill(0.0);

	potential.Make_fx(potentials::Attached_electron,mp);
	//Hamiltonian_Full<def_float,def_comp> H(*pge, potential, 1.0);
	//H.SetEigen();
	//H.eigs->GetState(phi_b.a,0);

    zOperatorF H = makeHamiltonianF(*pge, potential, 1.0);
    zEigSystem eSys = H.EigenSystem();
    eSys.GetState(phi_b.a, 0);

	for (int i=0; i<g.NbY(); ++i){
		phi.WriteX(phi_b,i);
	}

	return real(eSys.energies[0]);
}

// building the two dimensional discrete state wave function 
void MakePhiD_spec(gVector2D & phi, parameters2D & mp, femGrid2D & g)
{
	femGrid * pge = g.Gx();			// Pointer to the electonic grid
	femGrid * pgn = g.Gy();			// Pointer to the nuclear grid
	gVector potential(*pge);		// Potential
	gVector phi_b(*pge);			// Bound state 

	for (int i=0; i<g.NbY(); ++i){
		potential.Make_fxc(potentials::Lambda_Spec_R,pgn->Xr(i),mp);
		//Hamiltonian_Full<def_float,def_comp> H(*pge, potential, 1.0); 
		//H.SetEigen();
		//H.eigs->GetState(phi_b.a,0);
		zOperatorF H = makeHamiltonianF(*pge, potential, 1.0);
        zEigSystem eSys = H.EigenSystem();
        eSys.GetState(phi_b.a, 0);
		std::cout << "Eigenvalue E = " << eSys.energies[0] << " at R = " << pgn->Xr(i) << std::endl;
		phi.WriteX(phi_b,i);
	}
}

void MakePhiD_const(gVector2D & phi, parameters2D & mp, femGrid2D & g)
{
	femGrid * pge = g.Gx();			// Pointer to the electonic grid
	//femGrid * pgn = g.Gy();			// Pointer to the nuclear grid
	gVector potential(*pge);		// Potential
	gVector phi_b(*pge);			// Bound state 
	potential.Fill(0.0);

	potential.Make_fx(potentials::Attached_electron,mp);
	//Hamiltonian_Full<def_float,def_comp> H(*pge, potential, 1.0);
	//H.SetEigen();
	//H.eigs->GetState(phi_b.a,0);
    
    zOperatorF H = makeHamiltonianF(*pge, potential, 1.0);
    zEigSystem eSys = H.EigenSystem();
    eSys.GetState(phi_b.a, 0);
	
    for (int i=0; i<g.NbY(); ++i){
		phi.WriteX(phi_b,i);
	}
}


void sMatrix::Initialize(parametersNRM & p, def_float & MU, int & iChannel, gVector * St, zPolyVector& VDE, femGrid& grid)
{
	if (init) {
		Clean();		// Flush previous data
	}
// Copying constants for further use
	mu = MU;
	init_channel = iChannel;
	order = p.order;
	channels = p.ve_channels + p.da_channels;
	ve_channels = p.ve_channels;
	da_channels = p.da_channels;
	steps = p.steps;
	size = p.e_steps;
// initializing te coefficients once the order has been declared
	MakeCoefficients();
// Creating the dynamically allocated auxiliay variables
	e = new dVector(size, p.e_min, p.e_max, false);
	buffer = new zVector*[channels];
	s = new zVector[channels];
// Allocating and copying the vibrational states
	states = new gVector* [ve_channels];
	for (int i=0; i<ve_channels; ++i){
		states[i] = new gVector[size];
		for (int j=0; j<size; ++j){
			states[i][j].Set(grid);
			(states[i])[j] = VDE.EwSubMul(St[i], j);
		}
	}
	for (int i=0; i<channels; ++i){
		buffer[i] = new zVector[size];
		for (int j=0; j<size; ++j){
			(buffer[i])[j].Set(steps+1);
			(buffer[i])[j].Fill(0.0);
		}
		s[i].Set(size);
		s[i].Fill(0.0);
	}
	init = true;
}
void sMatrix::Clean(){
	if (init) {
		delete e;
		for (int i=0; i<ve_channels; ++i){
			delete[] buffer[i];
			delete[] states[i];
		}
		delete[] buffer;
		delete[] s;
		delete[] states;
		VdE = NULL;
	}
}
void sMatrix::MakeCoefficients()
{
	switch (order){
		case 0:
			coefficients.Set(steps);
			coefficients.Fill(1.0);
			break;
		case 1:
			coefficients.Set(steps+1);
			coefficients.Fill(1.0);
			coefficients[steps-1] = 0.5;
			coefficients[steps] = 0.5;		// The last element has to be added again if the integration continues another time
			break;
		case 2:
			coefficients.Set(steps+1);
			coefficients[0] = 4.0/3.0;
			for (int i=1;i<steps-1;++i){
				if (i%2==1) {
					coefficients[i]= 4.0/3.0;
				} else {
					coefficients[i]= 2.0/3.0;
				}
			}
			coefficients[steps-1] = 1.0/3.0;
			coefficients[steps] = 1.0/3.0;		// The last element has to be added again if the integration continues another time
	}
}
sMatrix::sMatrix()
{
	init=false;
}
sMatrix::sMatrix(parametersNRM & NRMp, def_float & MU, int & iChannel, gVector* St, zPolyVector& VDE, femGrid& grid)
{
	Initialize(NRMp, MU, iChannel, St, VDE, grid);
}
sMatrix::~sMatrix()
{
	Clean();
}
sMatrix & sMatrix::Set(parametersNRM & NRMp, def_float & MU, int & iChannel, gVector* St,  zPolyVector& VDE, femGrid& grid)
{
	Initialize(NRMp, MU, iChannel, St, VDE, grid);
	return *this;
}
void sMatrix::Contribution(zPolyVector * psi, int & i)
{ 
	gVector aux;
	for (int j=0; j<ve_channels; ++j){
		// Building contributions for separate channels
		for (int k=0; k<size; ++k){
			// Computing contributions for separate energies 
			(buffer[j])[k][i] = psi[k].SubMult((states[j])[k],0);
		}
	}
	for (int j=ve_channels; j<channels; ++j){
		for (int k=0; k<size; ++k){
			(buffer[j])[k][i] = psi[k].F(0, psi[0].grid->NR());
		}
	}
}
void sMatrix::CloseMultistep(def_float& time, def_float& dt, dVector& ergs, const def_float& ierg)
{
	for (int i=0; i<ve_channels; ++i){		// Contributions to vibrational excitations
		def_float shift = ierg - ergs[i];
		for (int k=0; k<size; ++k){
			if ((*e)[k] + shift>0) {
				if (order > 0) { // Contribution of the previous loop end
					s[i][k] += (1.0/imu)*exp(imu*((*e)[k] + ierg)*(time))* (buffer[i])[k][steps] *dt*coefficients[steps];	
				}
				for (int j=0; j<steps; ++j){
					s[i][k] += (1.0/imu)*exp(imu*((*e)[k] + ierg)*(time + (j+1)*dt))* (buffer[i])[k][j] *dt*coefficients[j];
				}
			}
			(buffer[i])[k][steps]=(buffer[i])[k][0];
		}
	}
	for (int i=ve_channels; i<channels; ++i){	// Contributions to dissociative attachment
		for (int k=0; k<size; ++k){
			def_float shift = ierg - ergs[i];
			if ((*e)[k] + shift>0) {
				def_float K = sqrt(2.0*mu*((*e)[k]+shift));
				if (order > 0) { // Contribution of the previous loop end
					s[i][k] += sqrt(K/(2*pi*mu)) * exp(-imu*K*R_0) * exp( imu*((*e)[k] + ierg)*(time)) * (buffer[i])[k][steps] * dt * coefficients[steps];
				}
				for (int j=0; j<steps; ++j){
					s[i][k] += sqrt(K/(2*pi*mu)) * exp(-imu*K*R_0) * exp( imu*((*e)[k] + ierg)*(time + (j+1)*dt)) * (buffer[i])[k][j] * dt * coefficients[j];
				}
			}
			(buffer[i])[k][steps] = (buffer[i])[k][0];
		}
		
	}
}
void sMatrix::CrossSections(const def_float & time, std::string & folder)		// Derives the cross sections and stores them into a file with appropriate time 
{
	char name[50];
	def_float val;
	def_float erg;
	//for (int i=0;i<ve_channels;++i){		
	sprintf_s(name, "NRM/CS/CSVE_t=%.1f.dat", time);
	FILE * file;
	fopen_s(&file,(folder+name).c_str(),"w");
	fprintf(file,"#Cross Sections of vibrational excitations: channels, values\n");
	fprintf(file,"#\t%d, %d\n", channels, size);
	for (int i=0;i<size;++i){
		erg = (*e)[i];
		fprintf(file, "%.12E", erg);
		for (int j=0;j<channels;++j){
			val = pow(abs(s[j][i]),2)*4.0*pow(pi,3)/(2.0*erg);
			fprintf(file, "\t%.12E", val);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}		
gVector* sMatrix::GetState(const int& i){ return states[i]; }


// The main model class
void ModelNRM::Initialize(parametersNRM & nrmp, parameters2D & m2dp, parametersMultiGrid & mgp)
{
	
	NRMp = nrmp;
	M2Dp = m2dp;
	MGp  = mgp;

	folder = m2dp.folder;
	grid_e.Set(mgp.gp[0]);
	grid.Set(mgp.gp[1]);
	femGrid2D g2d(grid_e,grid);
	phi_d.Set(&g2d);
	e_affinity = -MakePhiD_const(phi_d,m2dp,g2d);		// Builds the discrete state, and determines the affinity
	phi_d.Save((folder + "NRM/Phi_d.dat").c_str());
// Build the potentials and couplings
	if (	!Vdn.ReadBinary((folder + "NRM/bin/Vdn.bin").c_str(),&grid) || 
			!Vd.ReadBinary((folder + "NRM/bin/Vd.bin").c_str(), &grid)  ||
			!En.ReadBinary((folder + "NRM/bin/En.bin").c_str(), &grid)){
		build_potentials();
	}
	Vd.Save((folder + "NRM/Vd.dat").c_str());
	Vdn.Save((folder + "NRM/Vdn.dat").c_str());
	En.Save((folder + "NRM/En.dat").c_str());
	steps = nrmp.e_steps;
	energies.Set(nrmp.e_steps,m2dp.evol_par.e_min,m2dp.evol_par.e_max, false);
	std::cout << "Energy range: <" << energies[0] << ", " << energies[nrmp.e_steps-1] << ">" << std::endl;
	
	if (!VdE.ReadBinary((folder + "NRM/bin/VdE.bin").c_str(), &grid)){
		build_coupling();
	}
	VdE.Save((folder + "NRM/VdE.dat").c_str());
	VdE.SaveTransposed(energies, (folder + "NRM/VdET.dat").c_str());
// Build Hamiltoninan
	H.Set(Vdn, Vd, grid_e.NB(), m2dp.mu);
// Build the initial vibrational state & outgoing vibrational states
	gVector aux(grid);
	aux.Make_fx(potentials::V_zero,m2dp);
	aux.Save((folder+ "NRM/V_0.dat").c_str());
	
    //Hamiltonian_Full<def_float,def_comp> H_0(grid,aux,m2dp.mu);
	//H_0.SetEigen();
	zOperatorF H0 = makeHamiltonianF(grid, aux, m2dp.mu);
    zEigSystem eSys0 = H0.EigenSystem();

    // Coying Vibrational states energies
	ergs.Set(nrmp.ve_channels + nrmp.da_channels);
	vib_states = new gVector[nrmp.ve_channels];
	ierg = real( eSys0.energies[m2dp.init_par.channel] );
	for (int i=0; i<nrmp.ve_channels; ++i){
		ergs[i] = real( eSys0.energies[i] );
		vib_states[i].Set(grid);
        eSys0.GetState(vib_states[i].a, i);
	}
	// Copying the vibrational states
	
	loop = m2dp.evol_par.loop;
	dt = m2dp.evol_par.dt;
	time = 0.0;
	CN.Set(H,m2dp.mu, m2dp.evol_par.dt, m2dp.evol_par.pade);
	// Build the initial state (For all of the energies separately)
	psi = new zPolyVector[steps];
	for (int i=0;i<steps;++i){
		psi[i].Set(&grid,grid_e.NB());
		for(int j=0;j<grid.NB();++j){
			psi[i].F(aux.F(j)*VdE.F(i,j),0,j);
		}
	}
//	sMatrix & sMatrix::Set(parametersNRM & NRMp, def_float & MU, int & iChannel, gVector* St,  zPolyVector* VDE)
	S.Set(nrmp, m2dp.mu, m2dp.init_par.channel, vib_states, VdE, grid); 

	init = true;
}


// FIXME 
// FORTRAN code comparison failed (Discontinuous values)
void ModelNRM::build_potentials()
{
	// The method creates the Vdn and Vd
	int nb = grid.NB();			// The number of basis elements in the nuclear grid
	int nr = grid.NR();			// Number of real basis elements in the nuclear grid
	int N = grid_e.NB();		// Number of electronic grid basis functions (also discretisation of continuum)	
	
	Vdn.Set(&grid, N);
	Vd.Set(&grid, N);
	En.Set(&grid, N);
	gVector pot(grid_e);		// Vector for storing the potential
	gVector phi_ad(grid_e);		// Auxiliary vector for storing the electronic state
	
	iVector index(N);
	for (int i=nr-1; i>=0; --i){
		std::cout << "Calculating composed potentials Vdn and Vd for R = " << grid.Xr(i) << "..." << std::endl;

		// Build the potential
		pot.Make_fxc(potentials::V_eff_el,grid.Xr(i),M2Dp);

		// Build operators
		zOperatorF Hel(grid_e, grid_e.NB(), 0);
		
		// Add Kinetic Term and potential
		Hel.AddKineticTerm(1.0) += pot;
		
		// Build the Vd_0: <phi_d|H_el|phi_d>
		phi_d.ReadX(phi_ad,i);
		gVector phi_aux(phi_ad);
		Hel *= phi_ad;
		def_comp vzero = potentials::V_zero(grid.Xr(i),M2Dp);
		Vd.F(vzero + phi_aux*phi_ad, 0, i);
		// Builds the Feschbach continuum operator projection
		zOperatorF P(grid_e, grid_e.NB(), 0); 
		
		// Possible FIX
		//for (int q=0; q<grid_e.NB(); ++q) {
		//	for (int r=0; r<grid_e.NB(); ++r) {
		//		P[q*grid_e.NB() + r] = -phi_aux[q]*phi_aux[r];
		//	}
		//	P[q*grid_e.NB() + q] += 1;
		//}
		((P.OuterProduct(phi_aux)) *= -1) += def_comp(1); 

		zOperatorF PHP(P);
		PHP *= Hel; 
		PHP *= P;
		//zOperatorF PHP = P*Hel*P;  // FIXME
		// Build eigensystem of the PHP operator
		eigen_system<def_comp> ES = PHP.EigenSystem();

		if (i==nr-1) {		// Sort the values using the index ordering	
			for (int k=0; k<N; ++k){
				index[k] = k;
			}
		} else {			// Determine by order distance
			iVector oldInd = index;
			for (int k=0; k<N; ++k){
				index[k] = k;
			}

			for (int k=0; k<N; ++k){
				int near = k;			// Considering the k-th level as closest approximation
				def_float min = abs( En.F(k, i+1) - ES.energies[index[k]]);
				for (int l=k+1; l<N; ++l){
					def_float val  = abs(En.F(k,i+1) - ES.energies[index[l]]);
					if (val <= min){
						near = l;
						min  = val;	
					}
				}
				if (near!=k){
					std::swap(index[k],index[near]);
				}
			}

		// Upside down
//			for (int k=N-1; k>=0; --k){
//				int near = k;			// Considering the k-th level as closest approximation
//				def_float min = abs( En.F(k, i+1) - ES.energies[index[k]]);
//				for (int l=k-1; l>=0; --l){
//					def_float val  = abs(En.F(k,i+1) - ES.energies[index[l]]);
//					if (val <= min){
//						near = l;
//						min  = val;	
//					}
//				}
//				if (near!=k){
//					std::swap(index[k],index[near]);
//				}
//			}
			
		} 

		//En.F(ES.energies[index[0]], 0, i);
		for (int j=0; j<N; ++j){
			ES.GetState(phi_aux.a,index[j]);
			Vdn.F(conj(phi_ad*phi_aux), j, i);
			if (j!=0) Vd.F(vzero+ES.energies[index[j]], j, i);
			En.F(ES.energies[index[j]] , j, i); 
		}

		if (i==nr-1){
			for (int j=0; j<N; ++j){
				def_comp En_nr = En.F(j, i);
				def_comp Vdn_nr = Vdn.F(j, i);
				def_comp Vd_nr = Vd.F(j, i);
				for (int k=nb-1; k>i; --k){
					En.F(En_nr, j, k);
					Vdn.F(Vdn_nr, j, k);
					Vd.F(Vd_nr, j, k);
				}
			}
		}
	}
	Vdn.SaveBinary((folder + "NRM/bin/Vdn.bin").c_str());
	Vd.SaveBinary((folder + "NRM/bin/Vd.bin").c_str());
	En.SaveBinary((folder + "NRM/bin/En.bin").c_str());
}

// Calculates coupling for a given energy
gVector ModelNRM::calculate_VdErg(def_float erg)
{
	gVector VdErg(grid);
	//VdE.Set(&grid, steps);
	gVector pot(grid_e);		// Vector for storing the potential
	gVector phi_ad(grid_e);		// Auxiliary vector for storing the electronic state

	gVector phi; 
	phi.SetBessel(grid_e, erg, 1.0, M2Dp.l);

	for (int i=grid.NR()-1; i>=0; --i){
	// Build the potential, crop the complex region
		pot.Make_fxc(potentials::V_eff_el, grid.Xr(i), M2Dp);
		for (int k=grid_e.NR(); k<grid_e.NB(); ++k){
			pot[k] = 0.0;
		}
	// Build the Hamiltonian
		zOperatorF Hel(grid_e, grid_e.NR(), 0);
		Hel.AddKineticTerm(1.0) += pot;
		//def_comp vzero = potentials::V_zero(grid.Xr(i), M2Dp);
	// get the discrete state
		phi_d.ReadX(phi_ad,i);
	// Build projection operator
		zOperatorF P(grid_e, grid_e.NR(), 0);
		(P.OuterProduct(phi_ad) *= -1.0) += def_comp(1.0);
		
		// FIXME
		zOperatorF PHP(P);
		PHP*=Hel;
		PHP*=P;
		//zOperatorF PHP = P*Hel*P;

	// apply the Hamiltonian
		Hel *= phi_ad;
	// Compute the coupling
		PHP -= def_comp(erg);
		gVector aux(phi);
		(PHP *= aux) *= -1.0;
		aux[grid_e.NR()-1] = 0.0;
		PHP.BackSubstitution(aux) += phi;
		VdErg.F(phi_ad*aux,i);
	}
	return VdErg;
}


// Revised in comparison to FORTRAN code 
void ModelNRM::build_coupling()
{
	VdE.Set(&grid, steps);
	gVector pot(grid_e);		// Vector for storing the potential
	gVector phi_ad(grid_e);		// Auxiliary vector for storing the electronic state

	//Hamiltonian_Full<def_comp,def_comp> Hel(grid_e, pot, 1.0);
	gVector * phi = new grid_vector<def_float,def_comp>[steps];
	for (int j=0; j<steps; ++j){
		phi[j].SetBessel(grid_e,energies[j],1.0,M2Dp.l);
	}

	for (int i=grid.NR()-1; i>=0; --i){
		std::cout << "Calculating the coupling VdE at " << grid.Xr(i) << "..." << std::endl;
	// Build the potential, crop the complex region
		pot.Make_fxc(potentials::V_eff_el, grid.Xr(i), M2Dp);
		for (int k=grid_e.NR(); k<grid_e.NB(); ++k){
			pot[k] = 0.0;
		}
	// Build the Hamiltonian
		zOperatorF Hel(grid_e, grid_e.NR(), 0);
		Hel.AddKineticTerm(1.0) += pot;
		//def_comp vzero = potentials::V_zero(grid.Xr(i), M2Dp);
	// get the discrete state
		phi_d.ReadX(phi_ad,i);
	// Build projection operator
		zOperatorF P(grid_e, grid_e.NR(), 0);
		(P.OuterProduct(phi_ad) *= -1.0) += def_comp(1.0);
		
		// FIXME
		zOperatorF PHP(P);
		PHP*=Hel;
		PHP*=P;
		//zOperatorF PHP = P*Hel*P;

	// apply the Hamiltonian
		Hel *= phi_ad;
	// Compute the coupling
		for (int n=0; n<steps; ++n){
			zOperatorF PHP_a(PHP);
			PHP_a -= energies[n];
			gVector aux(phi[n]);
			(PHP_a *= aux) *= -1.0;
			aux[grid_e.NR()-1] = 0.0;
			PHP_a.BackSubstitution(aux) += phi[n];
			VdE.F(phi_ad*aux,n,i);
		}
	}
	delete[] phi;
	VdE.SaveBinary((folder + "NRM/bin/VdE.bin").c_str());
}
ModelNRM::ModelNRM()
{
	init = false;
}
ModelNRM::ModelNRM(parametersNRM & nrmp, parameters2D & m2dp, parametersMultiGrid & mgp)
{
	Initialize(nrmp, m2dp, mgp);
}
void ModelNRM::MultiStep()
{
	zPolyVector *rpsi;
	for (int i=0; i<loop; ++i){
		for (int j=0; j<steps; ++j)	{
			rpsi = &psi[j];
			CN.One_Step(*rpsi);
		}
		S.Contribution(psi,i);
		std::cout << '.'; 
		time += dt;
	}
	S.CloseMultistep(time, dt, ergs, ierg);
	std::cout << std::endl; 
}

int ModelNRM::TimeIndependentSolution()
{

	// Constants declarations

	int nr = grid.NR();
	int nb = grid.NB();
	int nbe = grid_e.NB();
	char name[5];
	

	def_comp val;

	// Determine the initial vibrational state
	gVector aux(grid);
	aux.Make_fx(potentials::V_zero,M2Dp);
	//Hamiltonian_Full<def_float,def_comp> H_0(grid,aux,M2Dp.mu);
	//H_0.SetEigen();
    zOperatorF H0 = makeHamiltonianF(grid, aux, M2Dp.mu);
    zEigSystem eSys0 = H0.EigenSystem();

	for (int ide=0; ide<steps; ++ide){
		def_float erg = real(energies[ide]); //NRMp.e_min + ide * (NRMp.e_max - NRMp.e_min)/(NRMp.e_steps-1);	// Actual kinetic enery of the incoming electron
		def_float E_t = erg + ierg;
		std::cout << "Calculating the dynamics fot E = " << E_t << std::endl;
		
		// Build nonlocal part of the potential
		// F(E,R,R') = Vdn(R) (E - T_N - V_0(R) - E_n(R))^(-1) Vdn(R') 
		zOperatorF FRR(grid, nr, 0); 			// Create a blank operator

		zOperatorF T(grid, nb, 0); 				// Builds blank operator 
		T.AddKineticTerm(M2Dp.mu);  			// Inserts kinetic term 
		for (int i=1; i<nbe; ++i){					// Loops through all of the electronic eigenstates except the lowest (discrete)
			zOperatorF G(grid, nr, 0);				// Blank Green operator
			//G -= T;									// Subtracts the kinetic term
			// Build of the potential
			gVector potential(grid);
			for (int j=0; j<nr; ++j){
				val = E_t; // - potentials::V_zero(grid.Xr(j),M2Dp) - potentials::V_cfg(0.0, M2Dp.mu, grid.Xr(j)) - En.F(i,j);
				potential.F(val ,j);
			}
			for (int j=nr; j<nb; ++j){
				val = E_t; // - En.F(i,nr);
				potential.F(val, j);
			}
			G += potential;
			G.Inverse();				// Adding the potential and compute the matrix inversion

// VDN !!!!
			zOperatorD DVdn(grid, nb, 0);			// Prepare the projection with VdE	
			DVdn = Vdn.GetVector(i);				// Get the projection potential
		
			//G *= DVdn;								// Right multiplitcation
			//((G.Transpose())*=DVdn).Transpose();	// Left multiplication     

			//std::cout << DVdn[0] << std::endl;
			FRR += G;								// Add to the nonlocal operator	
			
			//int start = 0, end = 1;
			//for (int i=start; i<end; ++i){
			//	for (int j=start; j<end; ++j) {
			//		//std::cout << G[i*nb + j] << std::endl;
			//	}
			//}
		}

		// Prepare the dynamics operator: E - T_N - V_d(R) - F(E,R,R')
		zOperatorF Ham(grid, nb, 0);			// Blank operator of Hamiltonian
		gVector potential(grid);				// Builds a blank vector
		Ham.AddKineticTerm(M2Dp.mu) *= -1.0;	// Adds the Kinetic term and multiplies by -1
	
		for (int i=0; i<nb; ++i){
			potential.F(E_t - Vd.F(0, i) - potentials::V_cfg(0.0, M2Dp.mu, grid.Xr(i)), i);
		}
		Ham += potential;
		Ham -= FRR;
	
	// Copying the vibrational states
		gVector ePsi(grid);
		//H_0.EigenState(ePsi, 0);
        eSys0.GetState(ePsi.a, 0);
		for (int i=0; i<nr; ++i){
			ePsi.F(ePsi.F(i)*VdE.F(ide,i), i);
		}
	
//		// Get the right side
//		gVector ePsi = S.GetState(0)[ide];		FIXME
//		// test the zeros
		for (int i=nr; i<nb; ++i){
			ePsi[i] = 0.0;
		}
		// Before back-subst
		sprintf(name, "%f", erg);
		std::cout << "... Norm =" <<  ePsi.Norm();
		ePsi.Save(( folder + "NRM/PsiIn" + name + ".dat").c_str());

		// Solve the equation
		Ham.BackSubstitution(ePsi);

		// After back subst
		std::cout << "... done. Norm of the determined state =" <<  ePsi.Norm() << std::endl;
		ePsi.Save(( folder + "NRM/PsiOut" + name + ".dat").c_str());

	}
	return 0;
}

