#ifndef DIRAC_TEST_FUNCTION_2D_INCLUDE_H_
#define DIRAC_TEST_FUNCTION_2D_INCLUDE_H_

#include "FemDvrEcs2d.h"
#include "Model2d/TestFunctionInterface2d.h"

#include "picojson/pjson.h"
#include "pjinput.h"

namespace QSCAT
{
/** \addtogroup Model2d
* @{ */

/// Dirac delta test function
/*!
    Specifiaction of the virtual interface to modified Tannor & Weeks approach
    with delta-distribution instead of the test-function wave packet
*/
class DiracTestFunction2d : public TestFunctionInterface2d<gVector2D>, public BinaryStorageInterface
{
    int position_;      //!< position of the test function
 private:
    /// Internal cleanup helper
    void clean();

 protected:
    /// Internal save to binary stream helper
    virtual bool save_bin_body(std::ofstream& file) const;
    /// Internal read from bianry stream helper
    virtual bool read_bin_body(std::ifstream& file);

 public:
    /// Default constructor
    DiracTestFunction2d();
    /// Copy constructor
    DiracTestFunction2d(const DiracTestFunction2d& old);
    /// Constructor from pjson
    DiracTestFunction2d( const pjvalue& params, const femGrid2D& g, int ch,
                  zEigenSystem& Eig, def_comp initial_energy, def_comp charge,
                  def_float mass, int impulsemomentum, const dVector& Energy);
    /// Swap operation
    DiracTestFunction2d& swap(DiracTestFunction2d& rhs);
    /// Assignement operation
    DiracTestFunction2d& operator=(DiracTestFunction2d tmp);
    /// Destructor
    ~DiracTestFunction2d();
    /// compute correlation and store in buffer
    virtual void operator<< (const gVector2D& psi);
    /// Flush the buffer to S-matrix
    virtual void contribution(zMatrix& S, int idx, def_float t, def_float dt, const zVector& ifc);
    /// Reduced mass
    const def_float& reduced_mass() const;
    /// Print coefficients
    virtual void print_coefficients(const std::string& path);
};

/** @} */
}; // namespace QSCAT

#endif // DIRAC_TEST_FUNCTION_2D_INCLUDE_H

