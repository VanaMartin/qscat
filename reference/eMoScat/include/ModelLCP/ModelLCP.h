#ifndef INCLUDE_MODEL_LCP_H_
#define INCLUDE_MODEL_LCP_H_

#include <iostream>             // Input/Output library
#include "picojson/pjson.h"
#include "pjinput.h"

#ifdef linux
    #include <sec_stream.h>
#endif

#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"
#include "interface.h"

#include "ModelLCP/SMatrix.h"

/// Local complex potential model

/// This namespace contains all of the necessary classes definitions for
/// computing local complex approximation of the two dimenstional model
namespace QSCAT
{

namespace LCP
{
/// Main model class
class ModelLCP {
    bool init_;                 //!< Initialization controller
    std::string method_;        //!< Evolution method
// Grids
    //int main_grid_;             //!< Index of the nuclear grid in the input file
    int nel_grids_;             //!< Total of electronic grids
    int *el_grids_;             //!< Indices of electronic grids (0,nel_grids)
    femGrid grid_;              //!< Main grid for the LCP approximation
    femGrid *grid_e_;           //!< Auxiliary grids for the electronic coordinate (for determining the V_res)
// Potentials
    gVector v_zero_;            //!< The molecular potenial: V_0(R)
    gVector vres_;              //!< Resonant potential: V_0(R) + E_res(R) + i/2 * Gamma(R)
    gVector eres_;              //!< Resonant energy: E_res(R)
    gVector gamma_;             //!< Resonance width: Gamma(R)
// Energy problem variables
    dVector energies_;          //!< Energy range on equidistant grid
// Wave functions
    gVector psi_;               //!< The wave function
    int channels_;              //!< Total of surveilled channels
    gVector * v_states_;        //!< Vibrational states
    dVector v_energies_;        //!< Energies of vibrational states
    def_float init_energy_;     //!< Inicident channel energy
// Evolution parameters & operators
    def_float mu_;              //!< The reduced mass
    def_float time_;            //!< Evolution time
    def_float dt_;              //!< Time step
    int loop_;                  //!< Number of steps in one loop
    CrankNicolson1D *cn_;       //!< Crank Nicolson variable
    int cn_order_;              //!< Crank Nicolson approximation order
    bool cn_use_;               //!< Crank Nicolson usage controller
    Chebyshev1D *cheb_;         //!< The Chebyshev variable
    int cheb_order_;            //!< Chebyshev order
    bool cheb_use_;             //!< Chebyshev usage controller
// Discrete state
    gVector * phi_res_;         //!< Resonant wave functions for the fixed nuclei electron problem
    gVector phi_a_;             //!< The aymptotic resonant wave function
    DiscreteStates * dstates_;  //!< Auxiliary pointer to the discrete states variable
    def_float affinity_;        //!< Electron affinity
    SMatrix S_;                 //!< S-Matrix elements variable
    std::string folder_;        //!< Output folder
    dBuffer norm_;              //!< normalization in time buffer
    dBuffer radius_;            //!< mean internuclear distance buffer
    zOperatorD opR_;            //!< radius operator
 // Resonant vibrations
    int n_vibrations_;          //!< Total of Resonant vibrations to be declared
    gVector *vibrations_;       //!< Vibrational wavefunctions in V_res potential
    ofstream population_file_;  //!< Vibrational level populations in time output
 private:
    /// initialization helper
    void initialize(const pjvalue& parameters);
    /// cleanup helper
    void clean();
    /// eigenstates from file reading helper
    bool read_eigenstates(const char * name, femGrid & g);
    /// resonant potential from file reading helper
    bool read_vres();
    /// resonant potential to file savin helper
    bool save_vres();
    /// computation of resonant potential helper
    void make_vres(const pjvalue& parameters);
 public:
    /// default constructor
    ModelLCP();
    /// constructor
    //ModelLCP(parameters2D& m2dp, parametersMultiGrid& mgp, parametersLCP& LCPp);
    ModelLCP(const pjvalue& parameters);
    /// destructor
    ~ModelLCP();
    /// computation of discrete state (two dimensional)
    void fill_discrete_state_phys(gVector2D& phi);
    /// multiple time step
    void multistep();
    /// discrete state getter
    gVector2D get_discrete_state(const pjvalue& parameters);
    //void MakePhiD_phys(gVector2D& phi);
    //void Multistep();
    //gVector2D MakePhiD();
};

} // namespace LCP

} // namespace QSCAT
#endif // INCLUDE_MODEL_LCP_H_
