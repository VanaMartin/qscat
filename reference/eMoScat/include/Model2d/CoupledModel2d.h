#ifndef COUPLED_MODEL_2D_INCLUDE_H_
#define COUPLED_MODEL_2D_INCLUDE_H_

#include "picojson/pjson.h"
#include "pjinput.h"

namespace QSCAT
{
/** \addtogroup Model2d
* @{ */

// Auxiliary
def_comp p_pi_potential(const def_comp& x, const def_comp& y);
def_comp p_sigma_potential(const def_comp& x, const def_comp& y);

/*
// COMPRESSED CASES IMPLEMENTATION
// Test function
class CompressedTestfunction : public TestFunctionInterface2d<gVector2D>
{
    mVector2D body_;           // The full testfunction (for projection use)
    gVector wave_packet_;      // The wave packet in the channel coordinate
 private:
    void initialize( const parametersEvolution& ep,
                     const parametersTestfunction& p, const mGrid2D& grid, int ch,
                     zEigenSystem& Eig);
    void clean();
 public:
    CompressedTestfunction();
    CompressedTestfunction( const parametersEvolution& ep,
                  const parametersTestfunction& p, const mGrid2D& grid, int ch,
                  zEigenSystem& Eig, const def_comp& IE, const def_comp zQ,
                  const def_float& M, const int L, const dVector& En);

    ~CompressedTestfunction();
    void set( const parametersEvolution& ep, const parametersTestfunction& p,
              const mGrid2D& grid, const int & ch, zEigenSystem& Eig,
              const def_comp& IE, const def_comp zQ, const def_float& M,
              const int L, const dVector& En);
};

//The main class for two-dimensional models with coupling
class CoupledModel2D
{
    bool init;                      // Initialization controller
    // Model parametristation
    // parameters
    //parameters2D& p1;               // Reference to the model parametrization
    //parameters2D& p2;
    // grids
    femGrid gx;                     // Electronic coordinate grid
    femGrid gy;                     // Nuclear coordinate grid
    femGrid2D g;                    // Two dimensional grid

    // Dynamics
    //bool Cheb_u, CN_u;              // Chebyshev and Crank-Nicolson usage controller
    def_float time;                 // Evolution time
    def_float dt;                   // Evolutoin time step
    int loop;                       // Number of loop steps
    // potentials
    gVector2D * potential;          // Vector storing the potential
    // state variables
    doubleGVector2D psi;            // Coupled State vector
    gVector2D *psia, *psib;         // Separate state shallow pointer
    //operators
    doubleOperator2D H;             // Main Hamiltonian operator
    doubleCrankNicolson2D CN;       // Crank-Nicolson template overload
    //doubleChebyshev2D CB;         // Chebyshev Evolution approximation

    // Scattering problem
    def_comp init_erg;              // Initial energy
    // testfunctions
    MultiTestfunction *tf;         // Test function in the electronic coordinate

    // Auxiliary variables
    std::string folder;             // Output folder

//            // Projection onto diabatic state statistics
//            gVector2D * PhiD;               // Pointer to the Projection State
//            gVector PsiD;                   // Projected wvefunction onto the discrete state
//            dBuffer norm;                   // Normalistation in a buffer
//            dBuffer radius;                 // Mean internuclear distance in a buffer
//            zOperatorD opR;                 // Operator storing the operation R
private:
    //void Initialize(parameters2D& m2dp, parametersGrid& gpx, parametersGrid& gpy);
    //void InitTestfunctions(parameters2D& p);
    //void InitPotential(parameters2D& p);
    //void Clean();
    //bool save_bin_body(std::ofstream& file);
    //bool read_bin_body(std::ifstream& file);
public:
    CoupledModel2D();   // Simple case (TODO to be removed)
    //CoupledModel2D(parameters2D& m2dp, parametersGrid& gpx, parametersGrid& gpy, gVector2D* phid = NULL);
    //~CoupledModel2D();
    //CoupledModel2D & Set(parameters2D& p, parametersGrid& gpx, parametersGrid& gpy);
    //void DiscreteProjection();
    //void Multistep();
    //bool SaveBinary(const char *name);
    //bool SaveBinary(std::ofstream& file);
    //bool ReadBinary(const char *name);
    //bool ReadBinary(std::ifstream& file);
};
*/

/** @} */
}; // namespace QSCAT

#endif // COUPLED_MODEL_2D_INCLUDE_H

