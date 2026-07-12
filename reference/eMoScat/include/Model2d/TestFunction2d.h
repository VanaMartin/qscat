#ifndef TEST_FUNCTION_2D_INCLUDE_H_
#define TEST_FUNCTION_2D_INCLUDE_H_

#include "FemDvrEcs2d.h"
#include "Model2d/TestFunctionInterface2d.h"

#include "picojson/pjson.h"
#include "pjinput.h"

namespace QSCAT
{
/** \addtogroup Model2d
* @{ */

/// Test function
/*!
    Specification of the virtual interface to Tannor & Weeks approach
    of testfunction.
*/
class TestFunction2d : public TestFunctionInterface2d<gVector2D>, public BinaryStorageInterface
{
    gVector2D body_;           //!< The full testfunction (for projection use)
    gVector wave_packet_;      //!< The wave packet in the channel coordinate
 private:
    /// Internal cleanup helper
    void clean();

 protected:
    /// Internal save to binary stream helper
    bool save_bin_body(std::ofstream& file) const;
    /// Internal read from bianry stream helper
    bool read_bin_body(std::ifstream& file);

 public:
    /// Default constructor
    TestFunction2d();
    /// Copy constructor
    TestFunction2d(const TestFunction2d& old);
    /// Constructor from pjson
    TestFunction2d( const pjvalue& params, const femGrid2D& g, int ch,
                  zEigenSystem& Eig, def_comp initial_energy, def_comp charge,
                  def_float mass, int impulsemomentum, const dVector& Energy);
    /// Swap operation
    TestFunction2d& swap(TestFunction2d& rhs);
    /// Assignement operation
    TestFunction2d& operator=(TestFunction2d tmp);
    /// Destructor
    ~TestFunction2d();
    /// compute correlation and store in buffer
    virtual void operator<< (const gVector2D& psi);
    /// Flush the buffer to S-matrix
    virtual void contribution(zMatrix& S, int idx, def_float t, def_float dt, const zVector& ifc);
    /// Print coefficients
    virtual void print_coefficients(const std::string& path);
};

/** @} */
}; // namespace QSCAT

#endif // TEST_FUNCTION_2D_INCLUDE_H
