#ifndef INCLUDE_KINETIC_ENERGY_H_
#define INCLUDE_KINETIC_ENERGY_H_

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! Kinetic term derived from laplace operator represented by full matrix.
/** The operator
 *
 *  \f[ \hat{T} = - \frac{1}{2\mu} \Delta = -\frac{1}{2\mu} \nabla \cdot \nabla, \f]
 *
 *  where \f$ \mu \f$ stands for mass term and \f$ \Delta \f$ resp. \f$ \nabla
 *  \f$ stands for Laplace resp. Nabla operator.
 *  The derivatives are evaluated from matrix of basis function derivatives on
 *  sample (-1,1) element and rescaled with proper element length (see
 *  FemDvrEcsGrid<T,Z>::dlp() member) The result forms block-diagonal-like
 *  matrix, where blocks are connected with one element
 *
 *  \f[ \hat{T} = \left(
 *      \begin{array}{cccc}
 *          t_{ij}^1 & t_{iq}^1 & 0 & 0 \\
 *          t_{qj}^1 & t_{qq}^1 + t_{00}^2 & t_{j0}^1 & \cdots \\
 *          0 & t_{i0}^2 & t_{ij}^2 & \\
 *          0 & \vdots & & \ddots
 *      \end{array}
 *  \right) \f]
 *
 *  where the element \f$ t_{ij}^k \f$ denotes the kinetic elements for all
 *  values of \f$ i,j \in {1, ... ,q-1} \f$ the specified values for \f$ 0,q
 *  \f$ are denoted in the overlapping row and column.
 *  @param  grid    Coordinate discretization.
 *  @param  mu      Reduced mass.
 *  @return New instantion of Matrix with the result.
 */
zMatrix generateKineticTerm(const FemDvrEcsGrid& grid, dfloat mu);

//! Kinetic term derived from laplace operator represented by row compressed matrix.
/** The operator
 *
 *  \f[ \hat{T} = - \frac{1}{2\mu} \Delta = -\frac{1}{2\mu} \nabla \cdot \nabla, \f]
 *
 *  where \f$ \mu \f$ stands for mass term and \f$ \Delta \f$ resp. \f$ \nabla
 *  \f$ stands for Laplace resp. nabla operator.
 *  The derivatives are evaluated from matrix of basis function derivatives on
 *  sample (-1,1) element and rescaled with proper element length (see
 *  FemDvrEcsGrid<T,Z>::dlp() member) The result forms block-diagonal-like
 *  matrix, where blocks are connected with one element
 *
 *  \f[ \hat{T} = \left(
 *      \begin{array}{cccc}
 *          t_{ij}^1 & t_{iq}^1 &  &  \\
 *          t_{qj}^1 & t_{qq}^1 + t_{00}^2 & t_{j0}^1 & \cdots \\
 *           & t_{i0}^2 & t_{ij}^2 & \\
 *           & \vdots & & \ddots
 *      \end{array}
 *  \right) \f]
 *
 *  where the element \f$ t_{ij}^k \f$ denotes the kinetic elements for all
 *  values of \f$ i,j \in {1, ... ,q-1} \f$ the specified values for \f$ 0,q
 *  \f$ are denoted in the overlapping row and column.
 *  @param  grid    Coordinate discretization.
 *  @param  mu      Reduced mass.
 *  @return New instantion of RowCompressedMatrix with the result.
*/
RowCompressedMatrix<dcomp> generateKineticTermRCM(const FemDvrEcsGrid& grid, dfloat mu);

/** @} */
}   // namespace QSCAT
# endif //INCLUDE_KINETIC_ENERGY_H_
