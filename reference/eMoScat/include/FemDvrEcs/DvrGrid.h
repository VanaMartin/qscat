#ifndef INCLUDE_DVR_GRID_H_
#define INCLUDE_DVR_GRID_H_

#include "common.h"
#include "Storage.h"
#include "Arrays/Vector.h"
#include "Arrays/Matrix.h"

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! DVR basis  grid class.
/** The DVR basis consists of Lagrangean interpolation polynomials through
 *  Gaussian quadrature points.  I.e. the integration over segment can be
 *  expressed as \f[ \int_{-1}^{1} f(x) dx \cong \sum_{i=1}^N w_i f_i, \f]
 *  where \f$ w_i \f$ denotes the Gaussian weight factors, \f$ f_i \f$ are the
 *  function values in quadrature points and \f$ N \f$ denotes the total of
 *  basis functions which equals "qudrature". The basis functions themselves
 *  can be expressed as \f[ \phi_j(x) = \prod_{i=1, i\neq j}^{N} \frac{x -
 *  x_i}{x_j - x_i}, \f] where the \f$ x_i \f$ denotes the Gaussian quadrature
 *  points.  The basis functions are orthogonal, i.e. fulfill relations of
 *  orthogonality \f[ \int_{-1}^{1} \phi_i^\dag(x) \phi_j(x) dx = \delta_{ij}.
 *  \f]
 */
class DvrGrid : public Object, public BinaryStorageInterface
{
 protected:
    //! Number of basis elements.
    blas_int quadrature_;
    //! Start point of the grid.
    dfloat x_min_;
    //! End point of the grid.
    dfloat x_max_;
    //! Array of grid points.
    dVector x_;
    //! Array of grid weights.
    dVector weights_;
 protected:
    //! Internal binary saving procedure.
    /** Saves the class structure into binary stream.
     *  @return true on success false on error.
     */
    virtual bool save_bin_body(std::ofstream& file) const;

    //! Internal binary reading procedure.
    /** Reads the class structure from binary stream.
     *  @return true on success false on error.
     */
    virtual bool read_bin_body(std::ifstream& file);

 private:
    //! internal cleanup method
    void clean();

 public:

 // constructors

    //! Default constructor.
    /** Sets all internal variables to default values. */
    DvrGrid();

    //! Basic constructor.
    /** Constructs a dvr basis on segment (-1,1) with given number of basis
     *  functions.
     *  @param  size    Total number of basis functions in discretization grid.
     */
    DvrGrid(blas_int size);

    //! Extended constructor.
    /** Constructs a dvr basis on segment given segment with given number of
     *  basis functions.
     *  @param  start   Starting point of discretized segment.
     *  @param  end     Ending point of discretized segment.
     *  @param  size    Total number of basis functions in discretization grid.
     */
    DvrGrid(dfloat start, dfloat end, blas_int size);

    //! Copy constructor (shallow copy).
    DvrGrid(const DvrGrid& old);

    //! Destructor.
    ~DvrGrid();

 // modifiers

    //! Assignement operator.
    /** @param  tmp     Temporary copy of source grid.
     *  @return Refence onto this updated grid.
     */
    DvrGrid& operator= (DvrGrid tmp);

    //! Swap operation.
    /** Swaps all internal variables between two instantions.
     *  @param  rhs     Grid to be swapped with.
     *  @return Refence onto this updated grid.
     */
    DvrGrid& swap(DvrGrid& rhs);

    //! Deep copy function.
    /** Performs a deep copy of all internal variables.
     *  @return New instantion identical to current one.
     */
    DvrGrid copy() const;

 // accessors

    //! Quadrature point \f$ x_i \f$.
    dfloat x(blas_int i) const;

    //! Quadrature weight \f$ w_i = \f$.
    dfloat w(blas_int i) const;

    //! Qudrature order
    blas_int quadrature() const;
};

//! Lagrange polynomial derivatives.
/** Auxiliary method for generating matrix of derivatives of Lagrangean
 *  polynomials in quadrature points.
 *  @param  nq      Size of the basis.
 *  @param  dLp     Matrix for storing values. Must be of size "nq*nq".
 *  @param  g       Source discretization points and weights.
 */
void DLagPol(int nq, dMatrix& dLp, DvrGrid& g);

/** @} */
}   // namspace QSCAT
# endif //INCLUDE_DVR_GRID_H_
