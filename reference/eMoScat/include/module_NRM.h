#ifndef __MODULE_NRM__
    #define __MODULE_NRM__

    using namespace ARRAYS;
    using namespace FEM_DVR_ECS;
    using namespace FEM_DVR_ECS_2D;
    using namespace std;

    namespace NRM{


    // building the two dimensional discrete state wave function
        void MakePhiD_spec(gVector2D & phi, parameters2D & mp, femGrid2D & g);
        double MakePhiD_const(gVector2D & phi, parameters2D & mp, femGrid2D & g);

    // S Matrix variable
        class sMatrix {
            bool init;						// Initialization controller
            int order;						// Time integration quadrature order
            int channels;					// Number of channels to be investigated (VE + DA)
            int ve_channels;				// Number of channels to be investigated (VE)
            int da_channels;				// Number of channels to be investigated (DA)
            int steps;						// Number of time qudrature steps
            int size;						// Number of energy discretization points
            int init_channel;				// Incident channel number
            def_float mu;					// Reduced mass
            def_float R_0;
            dVector coefficients;			// Quadrature coeffitients
            dVector * e;					// Energy range distribution
            zVector ** buffer;				// The multistep integration buffer for VE channel
            zVector * s;					// S matrix elements
            gVector ** states;				// The vibrational states
            zPolyVector * VdE;				// Coupling to continuum
        private:
            void Initialize(parametersNRM & p, def_float & MU, int & iChannel, gVector * St, zPolyVector& VDE, femGrid& grid);
            void Clean();
            void MakeCoefficients();
        public:
            sMatrix();
            sMatrix(parametersNRM & NRMp, def_float & MU, int & iChannel, gVector* St, zPolyVector& VDE, femGrid& grid);
            ~sMatrix();
            sMatrix & Set(parametersNRM & NRMp, def_float & MU, int & iChannel, gVector* St,  zPolyVector& VDE, femGrid& grid);
            void Contribution(zPolyVector * psi, int & i);
            void CloseMultistep(def_float& time, def_float& dt, dVector& ergs, const def_float& ierg);
            void CrossSections(const def_float & time, std::string & folder);
            gVector* GetState(const int& i);
        };

    // The main model class
        class ModelNRM{
            bool init;
        // Model configuration parameters
            parametersNRM 		NRMp;
            parameters2D		M2Dp;
            parametersMultiGrid MGp;
        // Dynamics
            femGrid 			grid, grid_e;		// Nuclear and electronic grids
            zPolyVector* 		psi;				// The-poly vectors for all of the incident energy range
            zPolyVector 		Vdn;				// Coupling to continuum for evolution (complex)
            zPolyVector 		Vd;					// diagonal part of the potential
            zPolyVector			En;					// Eigenvalues as potential curves
            zPolyVector 		VdE;				// coupling to continuum for projection (physical)
            gVector2D 			phi_d;				// discrete state (2D)
            gVector*			vib_states;			// Vibrational states of the molecule

            poly_Hamiltonian_RCM<def_float,def_comp> H;			// Hamiltonian for poly-vector
            poly_Crank_Nicolson_RCM<def_float,def_comp> CN;		// Crank-Nicolson operator for evolution (poly-vector)
        // Constants
            def_float e_affinity;						// Electorn affinity
        // Dynamics parameters
            int loop;
            def_float dt;
            def_float time;
            int steps;
            zVector energies;
            dVector ergs;
            def_float ierg;
            sMatrix S;
            std::string folder;					// Output folder
        private:
            void Initialize(parametersNRM & nrmp, parameters2D & m2dp, parametersMultiGrid & mgp);
            void build_potentials();
            void build_coupling();
            gVector calculate_VdErg(def_float erg);
        public:
            ModelNRM();
            ModelNRM(parametersNRM & nrmp, parameters2D & m2dp, parametersMultiGrid & mgp);
            void MultiStep();
            int TimeIndependentSolution();
        };
    }
#endif
