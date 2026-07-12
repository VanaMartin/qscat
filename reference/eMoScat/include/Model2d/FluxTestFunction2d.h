#ifndef FLUX_TEST_FUNCTION_2D_INCLUDE_H_
#define FLUX_TEST_FUNCTION_2D_INCLUDE_H_

#include "FemDvrEcs2d.h"
#include "Model2d/TestFunctionInterface2d.h"

#include "picojson/pjson.h"
#include "pjinput.h"

namespace QSCAT {
/** \addtogroup Model2d
* @{ */

/// Probability flux test function
/*!
    Specification of the virtual interface to the time-energy fourier transform
    of the probability flux projected to outgoing state
*/
class FluxTestFunction2d : public TestFunctionInterface2d<gVector2D>, public BinaryStorageInterface
{
    int element_index_;             //!< Index of the element containing the surface
    int index_start_;               //!< Starting index of the element basis
    int index_end_;                 //!< Ending index of the element basis
    int position_;                  //!< Position of the surface (shifted to the
                                    //!< element border)
    int quadrature_;                //!< Number of basis functions in each element
    mutable gVector projection_;    //!< Auxiliary grid vector for storing the part
                                    //!< of full state (psi) projection necessary for
                                    //!< derivative computation, mutable allows const
                                    //!< qualifier, although values change
    zVector phi_out_;               //!< Outgoing waves for repeated use at boundary
    zVector dphi_out_;              //!< Derivative of the outgoing waves at boundary
    zBuffer derivative_buffer_;     //!< Simple buffer for derivatives
 private:
    /// Internal cleanup hlper
    void clean();

 protected:
    /// Internal save to binary stream helper
    virtual bool save_bin_body(std::ofstream& file) const;
    /// Internal read from binary stream helper
    virtual bool read_bin_body(std::ifstream& file);

 public:
    /// Default constructor
    FluxTestFunction2d();
    /// Destructor
    ~FluxTestFunction2d();
    /// Copy constructor
    FluxTestFunction2d(const FluxTestFunction2d& old);
    /// Constructor from pjson
    FluxTestFunction2d( const pjvalue& params, const femGrid2D& g, int ch,
                  zEigenSystem& Eig, def_comp initial_energy, def_comp charge,
                  def_float mass, int impulsemomentum, const dVector& Energy);
    /// Swap operation
    FluxTestFunction2d& swap(FluxTestFunction2d& rhs);
    /// Assignement operation
    FluxTestFunction2d& operator=(FluxTestFunction2d tmp);
    /// compute flux and store in buffer
    virtual void operator<< (const gVector2D& psi);
    /// Flush the buffer to S-matrix
    virtual void contribution(zMatrix& S, int idx, def_float t, def_float dt, const zVector& ifc);
    /// Print coefficients
    virtual void print_coefficients(const std::string& path);
};

/** @} */
}; // namespace QSCAT

#endif // FLUX_TEST_FUNCTION_2D_INCLUDE_H

