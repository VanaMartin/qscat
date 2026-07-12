#ifndef INCLUDE_OPERATOR_DIAGONAL_H_
#define INCLUDE_OPERATOR_DIAGONAL_H_

#include "common.h"
#include "Arrays.h"
#include "input.h"

#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */
//! diagonal operator defined on Hilbert space represented by GridVectors
/** The operator is build on top of a specific grid provided to the
 *  constructor.  Operator has zeros everywhere except the diagonal. The values
 *  are stored in a Vector instance and does not necessarily have the same size
 *  as the discrete coordinate representation, i.e. the operator may have
 *  nonzero values only for one line segment on the coordinate. For such case
 *  only values on this segment are stored.
 *  @note The generality of subset operatio might be removed in future since it
 *  does not provide much functionality and complicates the code a lot. If
 *  necessary a derived class shoud be defined.
*/
class OperatorDiagonal : public Object
{
 protected:
    //! Values on the diagonal.
    zVector body_;
    //! Length of the diagonal representation.
    blas_int size_;
    //! Starting index of the representation.
    blas_int start_;
    //! Length of coordinate discretization.
    blas_int max_size_;
    //! Associated grid.
    FemDvrEcsGrid grid_;
 private:
    //! Internal initialization.
    /** Prepares all common settings for the operator representation.
     *  @param  grid    Coordinate discretization.
     *  @param  length  Length of the segment.
     *  @param  start   Index of the starting point.
     */
    void initialize(const FemDvrEcsGrid& grid, blas_int length, blas_int start);

    //! Cleanup helper.
    void clean();

 public:
 // constructors
    //! Default constructor (uninitialized state)
    /** Sets all internal values to uninitialized defaults.*/
    OperatorDiagonal();

    //! Constructor.
    /** Sets all values according to the coordinatie discretization, operator
     *  is defined on whole coordinate.
     *  @param  grid    Coordinate discretization.
     */
    OperatorDiagonal(const FemDvrEcsGrid& grid);

    //! Constructor for given segment of grid.
    /** Sets all values according to the coordinatie discretization, operator
     *  is defined on a continuous segment of given coordinate.
     *  @param  grid    Coordinate discretization, copy will be stored.
     *  @param  length  Length of the segment (total included basis functions).
     *  @param  start   Index of the starting point.
     */
    OperatorDiagonal(const FemDvrEcsGrid& grid, blas_int length, blas_int start);

    //! Copy constructor.
    /** Performs shallow copy operation on all internal members.*/
    OperatorDiagonal(const OperatorDiagonal& old);

    //! Destructor.
    ~OperatorDiagonal();

 // accessors

    //! Coordinate representaion.
    /** @param  i   Index of the element on the diagonal.
     *  @return Constant reference to desired element.
     */
    const dcomp& operator[] (blas_int i) const;

    //! Coordinate representaion.
    /** @param  i   Index of the element on the diagonal.
     *  @return Reference to desired element.
     */
    dcomp& operator[] (blas_int i);

    //! Size operation segment.
    /** @return Total number of the basis functions on which the operator is
     *  defined.
     */
    blas_int get_size() const;

    //! Starting position index.
    /** @returns Shift of the operator starting position (start or zero).*/
    blas_int get_shift() const;

    //! Coordinate discretization grid.
    /** @returns Constant reference on coordinate discretization class.*/
    const FemDvrEcsGrid& get_grid() const;

 // modifiers

    //! Identity operator.
    /** Sets all values on the diagonal equal to one on whole coordinate
     *  (discards the line segment definition).
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& identity();

    //! Swap operation between two instantions.
    /** Swaps all internal variables between two operators.
     *  @param  rhs     Right-hand-side target to swap with.
     *  @returns the reference on this instantion
    */
    OperatorDiagonal& swap(OperatorDiagonal& rhs);

    //! Operator inversion.
    /** Inversion operation on diagonal operator is simply inversion of matrix
     *  elements on the diagonal. If defined on the line segment of the
     *  coordinate, computes inversion on the line segment.
     *  @return Reference on this updated instantion.
    */
    OperatorDiagonal& inverse();

 // operators

    //! Assignement.
    /** Shallow copy and swap operation of the right hand side.
     *  @param  tmp     Temporary shallow copy of right hand side operator.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator=  (OperatorDiagonal tmp);

    //! Assign function values.
    /** Retrieve function values from provided state vector and assign them to
     *  the operator.
     *  @param  rhs     Source values. Must be defined on the same coordinate
     *                  representaion.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator=  (const GridVector& rhs);

    //! Assign simple values.
    /** Sets the values of the diagonal operator from provided vector.
     *  @param  rhs     Source values. Must be of the same size as the operator
     *                  diagonal.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator=  (const zVector& rhs);

    //! Assign constant.
    /** Sets the oparator to multiple of identity operator.
     *  @param  alpha       Scaling factor of identity operator.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator=  (const dcomp& alpha);

    //! Inplace addition of two operators.
    /** Inplace addition of internal representaion of right hand side operator.
     *  @param  rhs     Source operator. Must be defined on the same coordinate
     *                  representation.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator+= (const OperatorDiagonal& rhs);

    //! Inplace addition of values from vector.
    /** Performs inplace addition of right hand side vector elements.
     *  @param  rhs     Source values. Must be of the same size as the operator
     *                  diagonal.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator+= (const zVector& rhs);

    //! Inplace addition of function values from state vector.
    /** Performs inplace addition of function values stored in a state vector
     *  weighted representation.
     *  @param  rhs     Source values. Must be defined on the same coordinate
     *                  representaion.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator+= (const GridVector& rhs);

    //! Inplace addition of a constant.
    /** Adds a constant factpr to all elements of the operator diagonal.
     *  @param  alpha   Addition factor.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator+= (const dcomp& alpha);

    //! Inplace subtraction of values from a vector.
    /** performs inplace subtraction of right hand side vector
     *  @param  rhs     Source values. Must be of the same size as the operator
     *                  diagonal.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator-= (const zVector& rhs);

    //! Inplace subtraction of function values from state vector.
    /** Performs inplace subtraction of function values stored in a state
     *  vector weighted representation.
     *  @param  rhs     Source values. Must be defined on the same coordinate
     *                  representaion.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator-= (const GridVector& rhs);

    //! Inplace subtraction of two operators.
    /** Performs inplace subtraction of internal representaion of right hand
     *  side operator.
     *  @param  rhs     Source operator. Must be defined on the same coordinate
     *                  representation.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator-= (const OperatorDiagonal& rhs);

    //! Inplace subtraction of a constant.
    /** Subtracts constant factor from the operator diagonal elements.
     *  @param  alpha   Subtraction factor.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator-= (const dcomp& alpha);

    //! Inplace multiplication of two diagonal operators.
    /** Performs inplace element wise multiplication of all elements with the
     *  elements of provided operator.
     *  @param  rhs     Source operator. Must be defined on the same coordinate
     *                  representation.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator*= (const OperatorDiagonal& rhs);

    //! Inplace multiplication by a constant.
    /** Multiplies the operator elements by a given constant.
     *  @param  alpha   Multiplication factor.
     *  @return Reference on this updated instantion.
     */
    OperatorDiagonal& operator*= (const dcomp& alpha);

    //! Operator action.
    /** Applicates the operator onto given state vector.
     *  @param  rhs     Source state on which the operator is applied. Must be
     *                  defined on the same coordinate discretization.
     *  @return Oepration result in new state vector.
     */
    GridVector operator* (const GridVector& rhs) const;
};
/** @} */
}   // namespace QSCAT
# endif //INCLUDE_OPERATOR_DIAGONAL_H_
