#ifndef TEST_FUNCTION_INTERFACE_2D_INCLUDE_H_
#define TEST_FUNCTION_INTERFACE_2D_INCLUDE_H_

#include "FemDvrEcs2d.h"

#include "picojson/pjson.h"
#include "pjinput.h"

namespace QSCAT
{
/** \addtogroup Model2d
* @{ */


/// Test function asymptotics type
enum asymptotics_t
{
    BESSEL_TYPE,
    COULOMB_TYPE
};

/// testfunction virtual interface to be used in the rest of the code
/*!
    The virtual interface assembles the common variables and methods for
    all possible test functions
*/
template<typename State>
class TestFunctionInterface2d
{
 protected:
    bool init_;                     //!< Initialization controller
    char axis_;                     //!< Given channel switch (x: electron, y: nuclear)
    int channel_;                   //!< Transversal excitation state number
    int in_out_;                    //!< In/out controller (1: in, (-1): out)
    int size_;                      //!< Number of basis elements of the 2D grid
    int xsize_;                     //!< Number of basis elements of the electron grid
    int ysize_;                     //!< Number of basis elements of the nuclear grid
    def_float mu_x_;                //!< Reduced mass of the X coordinate
    def_float mu_y_;                //!< Reduced mass of the Y coordinate
    int impulse_momentum_;          //!< Impulse momentum quantum number
    def_comp charge_;               //!< Coulombic charge
    def_comp energy_;               //!< Transversal bound state energy
    def_comp initial_energy_;       //!< Initial transversal bound state energy
    def_float energy_shift_;        //!< Actual shift of the incoming particle energy to fit the threshold of the channel (if it is not below zero)
    gVector bound_state_;           //!< The bound state wave function in the second
    zVector fourier_coefficients_;  //!< Associated fourier coefficients on given range
    dVector energies_;              //!< Energy distribution in given range (shifted by subtraction of transversal energies)
    int quad_order_;                //!< Time integration quadrature order
    zBuffer buffer_;                //!< Simple values buffer
    dVector coefficients_;          //!< Integration over time coefficients
    FemDvrEcsGrid2d grid_;          //!< Associated coordinates discretization

    std::ofstream* outfile_;        //!< direct output
    bool opened_;                   //!< output state controller
 public:
    /// Default constructor
    TestFunctionInterface2d();
    /// Destructor
    ~TestFunctionInterface2d();
    /// Copy constructor
    TestFunctionInterface2d(const TestFunctionInterface2d& old);
    /// Swap Operaion
    TestFunctionInterface2d& swap(TestFunctionInterface2d& rhs);
    /// Initialization controller
    bool init() const;
    /// Energy distribution point
    def_float energy(int i) const;
    /// Energy range shifting
    const def_comp& initial_energy() const;
    //const def_comp& fourier_coefficient(int i) const;
    /// Testfunction actual energy with respect to initial shift
    const def_float& energy_shift() const;
    /// Reduced mass in the direction of continuous states
    const def_float& reduced_mass() const;
    /// action on given state and storage result in buffer
    virtual void operator<< (const State& psi) = 0;
    /// Flushing values of time integration to result matrix
    virtual void contribution(zMatrix& S, int idx, def_float t, def_float dt, const zVector& ifc) = 0;
    /// Open the functional output file
    void set_output(const std::string& filename);
    /// Print coefficients
    virtual void print_coefficients(const std::string& path) = 0;
};

/** @} */
}; // namespace QSCAT

#include "Model2d/TestFunctionInterface2d.hpp"

#endif // TEST_FUNCTION_INTERFACE_2D_INCLUDE_H


