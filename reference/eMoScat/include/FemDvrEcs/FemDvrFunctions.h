#ifndef INCLUDE_FEM_DVR_FUNCTIONS_H_
#define INCLUDE_FEM_DVR_FUNCTIONS_H_

#include "common.h"
#include "Arrays.h"

//! Some useful functions necessary for building the DVR basis

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! Gamma function
/** @param  x   Point within (0,3).
 *  @return Gamma function in value from range (0,3).
 */
dfloat Gamma0to3(dfloat x);

//! Recurrent relation coefficients, for various orthogonal polynomials.
/** Evaluates coefficients \f$ a_j, b_j \f$ for relation \f[ b_j p_j(x) = (x
 *  -a_j) p_{j-1}(x) - b_{j-1} p_{j-2}(x) \f] for the various classical
 *  (normalized) orthogonal polynomials, and the zero-th moment \f[ \mu = \int
 *  w(x) dx \f] of the given polynomial weight function w(x). Since the
 *  polynomials are orthonormalized, the tridiagonal matrix is guaranteed to be
 *  symmetric.
 *  @note Laguerre and Jacobi polynomials are dependent on Gamma function.
 *  @param  kind    Polynomial type switch: 1-Legendre, 2-Chebyshev 1st kind,
 *                  3-Chebyshec 2nd kind, 3-Hermite, 4-Jacobi, 5-Laguerre.
 *  @param  n       Total number of the polynomials assumed.
 *  @param  alpha   Work value used for Laguerre and Jacobi.
 *  @param  beta    Work value used only for Jacobi.
 *  @param  b       Output vector of values.
 *  @param  a       Output vector of values.
 *  @param  mu      Zeroth moment.
 */
void RecCoef( int kind,
              int n,
              dfloat alpha,
              dfloat beta, dVector& b,
              dVector& a,
              dfloat& mu );

//! Auxiliary functionn for Gaussian Quadrature.
dfloat gbshift(dfloat shift, int m, const dVector& t, const dVector& b);

//! Solving tridiagonal matrix for computing weights.
/** Translation of the algol procedure imtql2, num. math. 12,
 *  377-383(1968) by Martin and Wilkinson, as modified in num. math.  15,
 *  450(1970) by dubrulle. handbook for auto. comp., vol.ii-linear algebra,
 *  241-248(1971).
 *  This function finds the eigenvalues and first components of the
 *  eigenvectors of a symmetric tridiagonal matrix by the implicit QL method,
 *  and is adapted from the eispak routine imtql2.
 *  @param  n   Order of input the matrix.
 *  @param  d   Diagonal elements of the input matrix stored in vector,
 *              on exit overwritten by eigenvalues in ascending order, if an
 *              error exit is made, the eigenvalues are correct but unordered
 *              for indices.
 *  @param  e   Subdiagonal elements of the input matrix stored in vector,
 *              in (0,n-2), (n-1) arbitrary, destroyed by computation.
 *  @param  z   First row of identity matrix, overwritten by first components
 *              of the orthonormal eigenvectors of the symmetric tridiagonal
 *              matrix, if an error exit is made, z contains the eigenvectors
 *              associated with the stored eigenvalues.
 *  @param ierr On exit is set to 0 for normal return, j if the j-th eigenvalue
 *              has not been determined after 30 iterations.
 */
void gbtql2(int n, dVector& d, dVector& e, dVector& z, int& ierr);

//! Gaussian quadrature.
/** Computes abscissas and weights of the n-point quadrature in the
 *  given interval.
 *  @param  kind    Polynomial type switch: 1-Legendre, 2-Chebyshev 1st kind,
 *                  3-Chebyshec 2nd kind, 3-Hermite, 4-Jacobi, 5-Laguerre.
 *  @param  n       Order of the quadrature.
 *  @param  alpha   Used for Laguerre and Jacobi.
 *  @param  beta    Used only for Jacobi.
 *  @param  b       Auxiliary vector;
 *  @param  t       Output vector, on exit contains quadrature points.
 *  @param  w       Output vector, on exit contains weight factors.
 *  @param  x_min   Start of the interval.
 *  @param  x_max   End of the interval.
 *  @param  kpts    Controls the inclusion of endpoints of the interval,
 *                  0-discard, 1-include.
*/
void GaussQuad( int kind,
                int n,
                dfloat alpha,
                dfloat beta,
                dVector& b,
                dVector& t,
                dVector& w,
                dfloat x_min,
                dfloat x_max,
                int kpts );

//! Initialization of Gaussian quadrature grid.
/** Wraps the whole proces with simple input
 *  @param  n   Order of the quadrature.
 *  @param  x   Output vector, on exit contains qudrature points on given
 *              interval.
 *  @param  w   Output vector, on exit contains quadrature weight of given
 *              basis.
 */
void GLo_Quad( int n, dVector& x, dVector& w );

//! Roots of Pade expansion.
/** Computes roots of Pade approximation for Crank-Nicolson variable.
 *  @param  roots   Output vector, on exit contains computed roots.
 *  @param  order   Order of the Pade expansion.
 */
void Pade_Roots(zVector& roots, int order);

//! Extended Romberg algorithm for general equidistant quadrature.
/** Computes weight factors for equidistant quadrature.
 *  @param  order   Order of the quadrature.
 *  @param  size    Size of the quadrature (how many discretiztion points).
 *  @return Vector of weight factors.
 */
dVector equidistant_quadrature(int order, int size);

/** @} */
} // namespace QSCAT
# endif //INCLUDE_FEM_DVR_FUNCTIONS_H_
