#ifndef INCLUDE_OPERATOR_ROW_COMPRESSED_H_
#define INCLUDE_OPERATOR_ROW_COMPRESSED_H_

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! Operator defined on GridVector in row compressed representation.
/** Operator defined for given states associated to provided coordinate
 *  discretization via provided FemDvrEcsGrid in compressed x-representation.
 *  This representation stores all nonzero matrix elements in compressed form
 *  (RowCompressedMatrix) internally.
 */
class OperatorRowCompressed : public Object
{
 protected:
    //! Operator representation in csr format.
    RowCompressedMatrix<dcomp> body_;
    //! Maximal size of the representation.
    blas_int max_size_;
    //! Number of representation rows (outputs).
    blas_int rows_;
    //! Number of representation columns (inputs).
    blas_int columns_;
    //! Representation starting point.
    blas_int start_;
    //! Representation length.
    blas_int size_;
    //! Associated grid.
    FemDvrEcsGrid grid_;
  private:
    //! Internal initialization.
    /** Sets all internal variables, copies the discretization grid,
     *  initializes the internal representation.
     *  @param  grid            Source discretization.
     *  @param  num_nonzeros    Total of nonzero values in the compressed
     *                          representation.
     *  @param  length          Total number of basis function on which the
     *                          operator is defined.
     *  @param  shift           Position of first basis function within the
     *                          representation.
     */
    void initialize( const FemDvrEcsGrid& grid,
                     blas_int num_nonzeros,
                     blas_int length,
                     blas_int shift );

  public:
 // constructors

    //! Default constructor.
    /** Sets all internal variables to uninitialized default.*/
    OperatorRowCompressed();

    //! Simple constructor.
    /** Resulting operator is defined on whole coordinate. The internal
     *  representation is initialized with no non-zero values.
     *  @param  grid    Coordinate discretization.
     */
    OperatorRowCompressed(const FemDvrEcsGrid& grid);

    //! Extended constructor.
    /** Resulting operator in x-representaion is defined on continuous segment
     *  of the coordinate. The internal represetation is initialized with no
     *  non-zero values.
     *  @param  grid    Source coordinate discretization.
     *  @param  length  Lenth of the segment in discrete representation (total
     *                  of basis functions).
     *  @param  shift   Starting postion of the segment in discrete
     *                  representation (position of first basis function).
     */
    OperatorRowCompressed( const FemDvrEcsGrid& grid,
                           blas_int length,
                           blas_int shift );

    //! Full constructor.
    /** Resulting operator in x-representaion is defined on continuous segment
     *  of the coordinate. The internal represetation is initialized with given
     *  total of nonzero values.
     *  @param  grid            Source coordinate discretization.
     *  @param  num_nonzeros    Total of nonzero values in the compressed
     *                          representation.
     *  @param  length          Total number of basis function on which the
     *                          operator is defined.
     *  @param  shift           Position of first basis function within the
     *                          representation.
     */
    OperatorRowCompressed( const FemDvrEcsGrid& grid,
                           blas_int num_nonzeros,
                           blas_int length,
                           blas_int shift );

    //! Copy constructor.
    /** Performs a shallow copy operation on all member variables.*/
    OperatorRowCompressed(const OperatorRowCompressed& old);

    //! Destructor.
    /** Currently empty.*/
    ~OperatorRowCompressed();

 // accessors

    //! Associated representation size.
    /** @return The size of the total number of basis function on which the
     *          operator is defined.
     */
    blas_int get_size() const;

    //! Associated representaion shift.
    /** @return The position of the first basis function on which the operator
     *          is defined, as index in the associated discretized coordiante.
     */
    blas_int get_shift() const;

    //! Associated coordinate discretization.
    /** @return Constant reference onto the associated coordinate
     *  discretization.
     */
    const FemDvrEcsGrid& get_grid() const;

    //! Matrix element.
    /** Retrieves the matrix element of the operator in discretized
     *  x-representaion, i.e.  \f[ A_{ij} =  a_{ij} |x_i\rangle \langle x_j |,
     *  \f] where \f$A\f$ denotes the operator \f$a_{ij}\f$ denotes the matrix
     *  element and \f$|x_i\rangle\f$ is the i-th basis function.
     *  @param  i   Index of the element (row).
     *  @param  j   Index of the element (column).
     *  @return Matrix element value.
     */
    dcomp operator() (int i, int j) const;

 // modifiers

    // TODO rename remove or new constructor?
    //! Set the operator to identity operator form.
    /** Sets the internal sparse matrix to diagonal form with uniform nonzero
     *  value equal to one.
     *  @param  size    Size of the segment, starting with the very first basis
     *                  function.
     *  @param  fem     Associated coordinate discretization grid.
     *  @return Reference on this updated instantion.
     */
    OperatorRowCompressed& set_identity( const blas_int size,
                                        const FemDvrEcsGrid& fem );

    //! Swap operation between to instantions.
    /** Performs swap operation on all internal variables, thus swaps the
     *  instantions.
     *  @param  rhs     Instantion to be swapped with.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& swap(OperatorRowCompressed& rhs);

    //! Set as kinetic term (-Laplace operator).
    /** Overwrites the internal representaion of the operator with kinetic
     *  term, i.e.  \f[ - \frac{1}{2\mu} \frac{\partial^2}{\partial x^2} = -
     *  \frac{1}{2 \mu} \Delta, \f] where \f$\mu\f$ denotes the reduced mass,
     *  see GenerateKineticTermRCM method for more details.
     *  @note The instantion must be previously initialized with proper
     *  coordinate discretization.
     *  @param  mu      Reduced mass.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& set_kinetic_term(const dfloat mu);

    //! Complex conjugate operation.
    /** Performs elementwise multiplication of nonzero values imaginary parts
     *  by factor -1.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& complex_conjugate();

    //! Hermitean conjugation.
    /** Swithces the internal "operation" bool for conjugation operation which
     *  is used in some BLAS methods, thus provide fast HC operation.  Does not
     *  modify the actual values in the sparse matrix representaion.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& conjugate();

 // operators

    //! Deep copy operation.
    /** Simple deep copy operation of all internal members.*/
    OperatorRowCompressed copy() const;


    //! Assignement operation.
    /** Performs shallow copy and swap operation internally.
     *  @param tmp      Temporary copy of source operator.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& operator=  (OperatorRowCompressed tmp);

    //! Addition of function vector to diagonal.
    /** Retrieves actual fucntion values from the given function vector and
     *  adds the to the diagonal elements of the operator.
     *  @param  rhs     Source state vector. Must be defined on the same
     *                  coordinate discretization.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& operator+= (const GridVector& rhs);

    //! Addition of constant to diagonal.
    /** Adds a scalar factor to all diagonal elements.
     *  @param  alpha   Scalar factor to be added to diagonal.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& operator+= (const dcomp& alpha);

    //! Scaling of whole operator by given factor.
    /** Multiplies each nonzero element of the sparse representaion by given
     *  constant value.
     *  @param  alpha   Scaling factor.
     *  @return Reference onto this updated instantion.
     */
    OperatorRowCompressed& operator*= (const dcomp& alpha);

 // custom operations

    //! Extended matrix vector multiplication (wrapper around BLAS gemv).
    /** General fast BLAS method General Matrix Vector multiplication computes
     *  \f[ \vec{y} = \alpha op(A) \vec{x} + \beta \vec{y}, \f] where
     *  \f$\alpha,\beta\f$ are scaling factors \f$A\f$ denotes the operator,
     *  \f$ op() \f$ denotes the possible operation of conjugation (see
     *  conjugate method), \f$\vec{x}\f$ source vector and \f$\vec{y}\f$ is
     *  source and output vector.
     *  @param  alpha   Scaling factor \f$\alpha\f$.
     *  @param  x       Source state vector to applicate the operator on. Must
     *                  be defined on the same coodrinate discretization.
     *  @param  beta    Scaling factor \f$\beta\f$.
     *  @param  y       Source vector to be rescaled and to which the operation
     *                  result will be added. Must be defined on the same
     *                  coordinate discretization.
     *  @return Refernce on the updated state vector \f$y\f$ contianing the
     *          operation result.
     */
    void gemv( const dcomp alpha,
                const GridVector& x,
                const dcomp beta,
                GridVector& y ) const;

    //! Lower Upper (LU) decomposition.
    /** Performs the LU decomposition of the internal sparse matrix
     *  representaion.
     *  @return Reference onto this updated instantion.
     */
    void LU_factorize();

    //! Linear problem solution via LU decomposition.
    /** Solves the linear problem \f[A\vec{x}=\vec{b}\f] via LU decomposition
     *  (if not already perfomed, the LU decomposition is performed
     *  automatically).
     *  @param  rhs     Source right hand side state vector. The source will be
     *                  overwritten by solution on exit. The vector must be
     *                  defined on the same coordinate discretization.
     */
    void LU_back_substitution(GridVector& x);

    //! auxiliary accessor to internal reprezentation
    const RowCompressedMatrix<dcomp>& body() const;
};

/** @} */
}   // namespace QSCAT
# endif //INCLUDE_OPERATOR_ROW_COMPRESSED_H_
