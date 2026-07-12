#ifndef INCLUDE_OPERATOR_FULL_H_
#define INCLUDE_OPERATOR_FULL_H_

#include "common.h"
#include "Arrays.h"
#include "input.h"

#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"
#include "FemDvrEcs/OperatorDiagonal.h"
#include "FemDvrEcs/KineticEnergy.h"

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! Full representation operator for multi-purpose use
/** Operator defined for given states associated to provided coordinate
 *  discretization via provided FemDvrEcsGrid in its full x-representation.
 *  This representation stores all of the possible matrix elements in dense
 *  matrix (Matrix) internally.
 */
class OperatorFull : public Object
{
 protected:
    //! Values stored in the matrix class.
    zMatrix body_;
    //! Size of the matrix.
    blas_int size_;
    //! Operator starting position
    /** First basis element corresponding to nonzero matrix element.*/
    blas_int start_;
    //! Maximal size of the representation.
    blas_int max_size_;
    //! Associated grid.
    FemDvrEcsGrid grid_;
    //! Transposition controller.
    char transposed_;
    //! LU factorization controller.
    bool factorized_;
 private:
    //! Internal initialization.
    /** Prepares all necessary variables
     *  @param  grid    Associated grid to the operator.
     *  @param  length  Length of the operator segment.
     *  @param  start   Position of the first affected basis function.
     */
    void initialize(const FemDvrEcsGrid& grid, blas_int length, blas_int start);

 public:
  // constructors

    //! Default constructor (uninitialized instantion).
    /** To be used only for temporary uninitialized state.*/
    OperatorFull();

    //! Standard constructor.
    /** Simple constructor, builds the operator on whole coordinate.
     *  @param  grid    Associated coordinate discrete representaion.
     */
    OperatorFull(const FemDvrEcsGrid& grid);

    //! Extended constructor.
    /** Allows to initialize more constrained operator defined only on a part
     *  of the grid.  The resulting operator acts only on basis functions in
     *  given range (typically only real part).
     *  @param  grid    Associated coordinate discrete representaion.
     *  @param  length  Length of the operator segment.
     *  @param  start   Position of the first affected basis function.
     */
    OperatorFull(const FemDvrEcsGrid& grid, blas_int length, blas_int start);

    //! Copy constructor.
    /** Performs a shallow copy of all internal members.
     *  @parma old  Source object.
     */
    OperatorFull(const OperatorFull& old);

    //! Destructor.
    ~OperatorFull();

  // accessors

    //! Flattened array element getter.
    /** Retrives element from the internal matrix representation in flattened
     *  case.
     *  @param  i   Index of the element in flattened array.
     *  @return Constant referece onto the element value.
     */
    const dcomp& operator[] (blas_int i) const;

    //! Flattened array element setter.
    /** Retrives element from the internal matrix representation in flattened
     *  case.
     *  @param  i   Index of the element in flattened array.
     *  @return Referece onto the element value.
     */
    dcomp& operator[] (blas_int i);

    // TODO remove - duplicite
    //! Coordinate discretization.
    /** Retrives the coordinate discretization struture object.
     *  @return Constant referece onto coordinate discretization grid.
     */
    const FemDvrEcsGrid& grid() const;

    //! Operator definition range size.
    /** @return Total number of basis function on which the operator is
     *  defined.
     */
    blas_int get_size() const;

    //! Operator definintion shift.
    /** @return Index of the first basis function in disctretized coordinate
     *  representation on which the operator is defined.
     */
    blas_int get_shift() const;

  // modifiers

    //! Transposition operator.
    /** Switches the internal transposition operator controller.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& transpose();

    //! Set the operator to identity operator form.
    /** Sets all values on diagonal equal to one, all others to zero.  Resets
     *  all controllers to defaults.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& set_identity();

    //! Swap operation between to instantions.
    /** Performs swap operation on all internal variables, thus swaps the
     *  instantions.
     *  @param rhs  Instantion to be swapped with.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& swap(OperatorFull& rhs);

    // TODO check conditions
    //! Set operator to outer product of two states.
    /** Sets the internal representation of the operator as \f[ \phi \otimes
     *  \psi = | \phi \rangle \langle \psi | = \phi \psi^T, \f] where the states
     *  \f$ \phi, \psi \f$ represents the input states respectively.
     *  @param  rhs1    Left product participant state \f$ \phi \f$.
     *  @param  rhs2    Right product participant state \f$ \psi \f$.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& outer_product(const GridVector& rhs1, const GridVector& rhs2);

    //! Inversion operation.
    /** Performs inversion operation on the operator. Internaly inversion of
     *  the full matrix representation, see Matrix::inverse().
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& inverse();

    //! Lower Upper (LU) decomposition.
    /** Performs the LU decomposition of the internal dense matrix
     *  representaion.  Useful for multiple back-substitutions.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& LU_factorize();

    //! Add kinetic term (-Laplace operator).
    /** Subtracts full representaion of Laplace operator to the current state
     *  of the operator, see GenerateKineticTerm method.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& add_kinetic_term(const dfloat& mu);

    //! Coordinate discretization.
    /** Retrives the coordinate discretization struture object.
     *  @return Constant referece onto coordinate discretization grid.
     */
    const FemDvrEcsGrid& get_grid() const { return grid_; }

  // operators

    //! Deep Copy operation.
    /** Performs deep copy of all internal members.
     *  @return A new cpopied instantion.
     */
    OperatorFull copy() const;

    //! Assignement operation.
    /** Performs a shallow copy and swap operation internally.
     *  @param  tmp     Temportary copy of source operator to be assigned.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator=  (OperatorFull tmp);

    //! Assignement operation of diagonal operator.
    /** Sets all values of the right-hand-side to the diagonal, nullifies all
     *  others.
     *  @param  rhs     Source of diagonal values. Must be defined on the same
     *                  coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator=  (const OperatorDiagonal& rhs);

    //! Assignement of state vector to diagonal.
    /** Nullifies the internal representaion and ads the vector function values
     *  to the diagonal.
     *  @param rhs      Source of diagonal values. Must be defined on the same
     *                  coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator=  (const GridVector& rhs);

    //! Assignement of the Vector values to diagonal.
    /** Nullifies the internal representaion and ads the vector
     *  values to the diagonal.
     *  @param rhs  Source of diagonal values. Must have same size as the
     *              coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator=  (const zVector& rhs);

    //! Assigement of constant value to the diagonal.
    /** Assignes all elements of the diagonal to constant value,
     *  nullifies all others.
     *  @param alpha    Value to be assigned to diagonal.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator=  (const dcomp& alpha);


    //! Inplace addition.
    /** Adds all source matrix elements to the current state.
     *  @param rhs      Source of the elements to be added. Must be defined on
     *                  the same coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator+= (const OperatorFull& rhs);

    //! Inplace addition operation of diagonal operator.
    /** Adds values of the right-hand-side to the diagonal to the current state
     *  of the dense operator.
     *  @param rhs  Source of the diagonal values. Must be defined on the
     *              same coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator+= (const OperatorDiagonal& rhs);

    //! Inplace addition of the Vector values to diagonal.
    /** Adds the source values to the diagonal of the current state of dense
     *  operator.
     *  @param  rhs     Source of diagonal values. Must be of the same size as
     *                  the coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator+= (const zVector& rhs);

    //! Inplace addition of the GridVector values to diagonal.
    /** Adds the source function values to the diagonal of the current state of
     *  dense operator.
     *  @param  rhs     Source of diagonal values. Must be defined on the
     *                  coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator+= (const GridVector& rhs);

    //! Inplace addition of constant value to the diagonal.
    /** Adds a constant factor to all elements on the diagonal.
     *  @param  alpha   Scalar factor to be added to diagonal.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator+= (const dcomp& alpha);


    //! Inplace subtraction operation of another dense operator.
    /** Subtracts all source matrix elements from the current state.
     *  @param rhs      Source of the elements to be added. Must be defined on
     *                  the same coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator-= (const OperatorFull& rhs);

    //! Inplace subtraction operation of diagonal operator.
    /** Subtracts values of the right-hand-side from the diagonal of the
     *  current dense operator.
     *  @param rhs  Source of the diagonal values. Must be defined on the
     *              same coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator-= (const OperatorDiagonal& rhs);

    //! Inplace subtraction of values from the diagonal.
    /** Subtracts the source values from the diagonal of the current state of
     *  dense operator.
     *  @param  rhs     Source of diagonal values. Must be of the same size as
     *                  the coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator-= (const zVector& rhs);

    //! Inplace subtraction of the state vector elements from the diagonal.
    /** Subtracts the source function values from the diagonal of the
     *  current state of dense operator.
     *  @param  rhs     Source of diagonal values. Must be defined on the same
     *                  coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator-= (const GridVector& rhs);

    //! Inplace subtraction of constant value from the diagonal.
    /** Subtracts a scalar factor from all elements of the diagonal.
     *  @param  alpha   Scalar factor to be subtracted from the diagonal.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator-= (const dcomp& alpha);


    //! Inplace multiplication by another dense operator.
    /** Performs inplace matrix multiplication of the current matrix elements
     *  with the source dense matrix elements.
     *  @param  rhs     Source values to be multiplied with. Must be defined on
     *                  the same coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator*= (const OperatorFull& rhs);

    //! Inplace multipication by a diagonal operator.
    /** Performs inplace multiplication of each column by a diagonal element of
     *  the source diagonal operator (i.e. matrix multiplication with diagonal
     *  matrix)
     *  @param  rhs     Source of diagonal values. Must be defined on the same
     *                  coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator*= (const OperatorDiagonal& rhs);

    //! Inplace multipication by a constant value.
    /** Multiplies all elements of the operator with constant scaling factor.
     *  @param  rhs     Scaling factor.
     *  @return Reference onto this updated instantion.
     */
    OperatorFull& operator*= (const dcomp& alpha);

    //! Operator application on given state.
    /** Performs matrix vector multiplication internally of source state vector
     *  with the matrix elements of the dense operator.
     *  @param  rhs     Source state vector. Must be defined on the same
     *                  coordinate discretization.
     *  @return New state vector object containing the result.
     */
    GridVector operator* (const GridVector& rhs) const;

  // custom operations

    //! Solution to eigenvalue problem.
    /** Solves the eigenvalue internally via EigenSystem constructor.
     *  @return New eigen system  object containing all eigenvalues and
     *          eigenvectors.
     */
    EigenSystem<dcomp> eigen_system();

    //! Linear problem solution.
    /** Solves the linear problem \f[ A\vec{x} = \vec{b} \f] via LAPACK back
     *  substitution.
     *  @param  rhs     Right hand side state vector \f$\vec{b}\f$. Must be
     *                  defined on the same coordinate discretization. The
     *                  source is overwritten by solution on exit.
     *  @return Reference onto the source vector overwritten by soulution.
     */
    GridVector& back_substitution(GridVector& rhs);

    //! Linear problem solution via LU decomposition.
    /** Solves the linear problem \f[ A\vec{x} = \vec{b} \f] via LU decomposition
     *  (if not already perfomed, the LU decomposition is performed
     *  automatically).
     *  @param  rhs     Right hand side state vector \f$\vec{b}\f$. Must be
     *                  defined on the same coordinate discretization. The
     *                  source is overwritten by solution on exit.
     *  @return Reference onto the source vector overwritten by soulution.
     */
    GridVector& LU_back_substitution(GridVector& rhs);

    //! Linear problem solution via LU decomposition.
    /** Solves the linear problem \f[ A\vec{x} = \vec{b} \f] via LU decomposition
     *  if the LU is already performed otherwise via LAPACK iterative algorithm.
     *  @param  rhs     Right hand side state vector \f$\vec{b}\f$. Must be
     *                  defined on the same coordinate discretization. The
     *                  source is overwritten by solution on exit.
     *  @return Reference onto the source vector overwritten by soulution.
     */
    GridVector& smart_back_substitution(GridVector& rhs);

    //! Extended matrix vector multiplication (wrapper around BLAS gemv).
    /** General fast BLAS method General Matrix Vector multiplication computes
     *  \f[ \vec{y} = \alpha op(A) \vec{x} + \beta \vec{y}, \f] where
     *  \f$\alpha,\beta\f$ are scaling factors \f$A\f$ denotes the operator,
     *  \f$ op() \f$ denotes the possible operation of conjugation (see.
     *  transpose, conjugate), \f$\vec{x}\f$ source vector and \f$\vec{y}\f$ is
     *  source and output vector.
     *  @param  alpha   Scaling factor \f$\alpha\f$.
     *  @param  x       Source vector to applicate the operator on. Must be
     *                  defined on the same coordinate discretization.
     *  @param  beta    Scaling factor \f$\beta\f$.
     *  @param  y       Source vector to be rescaled and to which the operation
     *                  result will be added. Must be defined on the same
     *                  coordiante representation.
     *  @return The refernce on \f$y\f$ state vector with the operation result.
     */
    GridVector& gemv( const dcomp alpha,
                      const GridVector& x,
                      const dcomp beta,
                      GridVector& y );
};

OperatorFull operator* (OperatorFull lhs, OperatorDiagonal& rhs);
OperatorFull operator* (OperatorDiagonal& lhs, OperatorFull rhs);
OperatorFull operator* (OperatorFull lhs, OperatorFull& rhs);
GridVector operator* (GridVector v, OperatorFull& O);

/** @} */
}   // namespace QSCAT
# endif //INCLUDE_OPERATOR_FULL_H_
