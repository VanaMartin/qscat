#ifndef INCLUDE_QSCAT_MATRIX_H_
#define INCLUDE_QSCAT_MATRIX_H_

#include "blas.h"
#include "common.h"
#include "Storage.h"
#include "ScalarMultiple.h"
#include "Arrays/Vector.h"
#include "Arrays/EigenSystem.h"

namespace QSCAT
{
/** \addtogroup Arrays
* @{ */

// FIXME the Hermitean conjugation controller is ignored in many methods

/// matrix of values with size rows*columns template
/**
 *  Matrix of type <T> stored as a contiguous array of length "rows*columns".
 *  The ordering is  column wise, i.e. element in  m-th row and n-th column is
 *  stored in the  vector under  index of (m*columns_ + n). Therefore the
 *  column  index is contiguous, row index is not. Basic operations are wrapped
 *  with class methods leading to fast call of BLAS operations.
 */
template <typename T>
class Matrix : public Object, public BinaryStorageInterface
{
 protected:
    //! Transposition/conjugation controller.
    char transposed_;
    //! Total of rows.
    blas_int rows_;
    //! Total of columns.
    blas_int columns_;
    //! The elements of the matrix.
    T *array_;
    //! auxiliary pointer to the end of the array.
    T *end_;
    //! Decomposition controller.
    /** If set true, matrix was decomposed via LU decomposition. */
    bool decomposed_;
    //! Permutations of the original matrix.
    /** Auxiliary for pivotization phase of LU decomposition. */
    blas_int *pivots_;
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
    //! Internal initialization procedure.
    /** Allocates internal memory, sets all internal variables.
     *  @param  rows    Total of rows in the matrix, leading dimension.
     *  @param  columns Total of columns in the matrix.
     */
    void initialize(blas_int rows, blas_int columns);

    //! Internal cleaning procedure.
    void clean();
 public:
 // constructors
    //! Standard constructor of matrix "rows$\times$columns".
    /** Initializes the matrix to given shape.
     *  @param  rows    Total of rows in the matrix, leading dimension.
     *  @param  columns Total of columns in the matrix.
     */
    Matrix(blas_int rows, blas_int columns);

    //! Copy constructor (shallow copy).
    /** Performs a shallow copy operation (copies only the pointer and shape
     *  values).
     *  @param  old Instance of the matrix to be copied from.
     */
    Matrix(const Matrix & old);

    //! Default constructor.
    /** Sets all internal variables to defaults. Count as uninitialized state.
     */
    Matrix();

    //! Destructor.
    /** Decrease the refence count, if total remaining refernces equals zero
     * cleans all allocated memory.
     */
    ~Matrix();

 // accessors

    //! Access matrix element.
    /** Retrieves matrix element in giver row and column.
     *  @param  row     From which row.
     *  @param  column  From wich column.
     *  @return Reference to desired element.
     */
    T& operator() (const blas_int &row, const blas_int &column);

    //! Access matrix element.
    /** Retrieves matrix element in giver row and column.
     *  @param  row     From which row.
     *  @param  column  From wich column.
     *  @return Constant reference to desired element.
     */
    const T& operator() (const blas_int &row, const blas_int &column) const;

    //! Access matrix element as in flat array.
    /** Retrieves matrix element in given position index in the flattened
     *  array.
     *  @param  index   Position of desired element.
     *  @return Reference to desired element.
     */
    T& operator[] (const blas_int &index);

    //! Access matrix element as in flat array.
    /** Retrieves matrix element in given position index in the flattened
     *  array.
     *  @param  index   Position of desired element.
     *  @return Constant reference to desired element.
     */
    const T& operator[] (const blas_int &index) const;

    //! Total elements inside the matrix.
    /** The number is equal to the size of the flattened matrix and equals
     *  rows*columns.
     *  @return Total number of elements.
     */
    blas_int get_size() const;

    //! Number of rows.
    blas_int rows() const;

    //! Number of columns.
    blas_int columns() const;

    //! Retrieve a column from the matrix.
    /** Creates a new instance of the Vector class and performs a deep copy of
     *  the matrix elements of the desired column.
     *  @note   This operation is optimized since the matrix is column ordered.
     *  @param  index   Column number, must be less than total number of
     *                  columns.
     *  @return New instantion of Vector class with stored values.
     */
    Vector<T> get_column(blas_int index) const;

    //! Retrieve a row from the matrix.
    /** Creates a new instance of the Vector class and performs a deep copy of
     *  the matrix elements of the desired row.
     *  @note   This operation is NOT optimized since the matrix is column
     *          ordered.
     *  @param  index   Row number, must be less than total number of rows.
     *  @return New instantion of Vector class with stored values.
     */
    Vector<T> get_row(blas_int index) const;

 // modifiers

    //! Fills whole matrix with given value.
    /** @param  constant    Value to be filled with.
     *  @return  Reference on the updated current matrix.
     */
    Matrix& fill(const T &val);

    //! Swaps internal members between instantions.
    /** Performs simple swap operation to all internal memebers.
     *  @param  rhs The object to be swapped with.
     *  @return Reference onto this object swapped with "rhs".
     */
    Matrix& swap(Matrix &rhs);

    //! Deep copy operation.
    /** Performs the actual deep copy operation.
     *  @return A new instance of the same object.
     */
    Matrix copy() const;

    //! Set to identity matrix.
    /** Resets the matrix to a square matrix of given size, fill with zero
     *  everywhere except the matrix diagonal where values are set to one.
     *  @param  size    Size of the newly constructed identity matrix.
     *  @return Reference on this instantion.
     */
    Matrix& set_identity(blas_int size);

    //! LU decomposition.
    /** The LU decomposition consists of three parts. First pivotization
     *  constructs the matrix permutations (via row & column swap) than
     *  computation of lower triangular matrix (L) with ones on diagonal and
     *  last the computation of upper triangular matrix (U). Both L and U are
     *  stored in the original array where the diagonal ones from L are
     *  ommited, the permutation are stored separately in pivots_ member.
     *  @note   This method assumes the matrix to be square.
     *  @return Reference onto this instantion.
     */
    Matrix& LU_factorize();

    //! Hermitean conjugation switch.
    /** Switch the internal controller for Hermitean conjugation. The operation
     *  itself is not performed on the actual data however some methods are
     *  subject to different behaviour if the conjugation operation is switched
     *  to true.
     *  @note   Conjugation is taken in account for methods: gemv, operator*
     *  @note   The operation is very effective (there is no actual changes in
     *  the array).
     *  @return Reference onto this instantion.
     */
    Matrix& conjugate();

    //! Complex conjugations of elements.
    /** This operation is suitable only if one needs the actual complex
     *  conjugated matrix, for hermitean conjugate as input to BLAS routines it
     *  is ineffective to use this method.
     *  @return Reference onto this instantion.
     */
    Matrix& complex_conjugate();

    /// Matrix inversion.
    /** Overwrites the current matrix with computed matrix inversion.
     *  Internally performs the Solve operation for all columns of the
     *  corresponding identity matrix.
     *  @note   This operation is valid only for square matrices.
     *  @return Reference onto this instantion.
     */
    Matrix& inverse();

 // operators

    //! Assignement operation.
    Matrix& operator= (Matrix tmp);

    //! Inplace addition of another matrix.
    /** Performs the inplace addition of righ-hand-side matrix element wise.
     *  @note   The conjugation switch is ignored.
     *  @param  rhs Matrix to be added. Must be of the same shape as this
     *              instantion.
     *  @return Reference onto this instantion.
     */
    Matrix& operator+= (const Matrix &rhs);

    //! Inplace subtraction of another matrix.
    /** Performs the inplace subtraction of righ-hand-side matrix element wise.
     *  @note   The conjugation switch is ignored.
     *  @param  rhs Matrix to be subtracted. Must be of the same shape as this
     *              instantion.
     *  @return Reference onto this instantion.
     */
    Matrix& operator-= (const Matrix &rhs);

    //! Inplace multiplication with another matrix.
    /** Overwrites the matrix with result of matrix matrix product.
     *  @note   This operation respects the conjugation operation on both
     *          product participants.
     *  @param  rhs Matrix to be multiplied with from the right side. Must be
     *              of the compatible shape with respect to the conjugation
     *              operation (see operator* for details).
     *  @return Reference onto this instantion.
     */
    Matrix& operator*= (const Matrix &rhs);

    //! Matrix multiplication.
    /** Computes matrix matrix inner product.
     *  @note   This operation respects the conjugation operation on both
     *          product participants.
     *  @param  rhs Right side prodict participant. The shapes of both matrices
     *              must be compatible with respect to their conjugation
     *              operation: Value of "rhs" resp. "conjugated rhs" dimension
     *              "rows" resp. "columns" must be equal to this intantion
     *              "columns" or "rows" if this instatnion is conjugated.
     *  @return Reference onto this instantion.
     */
    Matrix  operator* (const Matrix &rhs) const;

    //! Multiply by constant.
    /** Performs element wise scaling by constant value.
     *  @param  scalar  Value to scale with.
     *  @return Reference onto this instantion.
     */
    Matrix& operator*= (const T &scalar);

    //! Inserts a vector of values into the matrix.
    /** Copies values from the vector specified either row or column.
     *  @param  index   To which row resp. column the values will be inserted.
     *  @param  source  Vector of values to be inserted.
     *  @param  axis    0 if vector represents a column 1 if vector represents
     *                  a row.
     *  @return Reference onto this instantion.
     */
    Matrix& operator() ( blas_int index,
                         const Vector<T>& source,
                         blas_int axis=0 );

    //! Matrix vector multiplication.
    /** Computes matrix vector multiplication (vector from right side).
     *  @note   This operation respects the cojugation operation.
     *  @param  rhs Source vector values.
     *  @return New Vector instantion contianing result.
     */
    Vector<T> operator* (const Vector<T>& rhs) const;

    //! Scaled matrix vector multiplication.
    /** Computes matrix vector multiplication (vector from right side) scaled
     *  by a multiplicative factor.
     *  @note   This operation respects the cojugation operation.
     *  @param  rhs Source vector and scalar value wrapped in temporary class.
     *  @return New Vector instantion contianing result.
     */
    Vector<T> operator* (ConstScalarMultiple<T, Vector<T> > & rhs) const;

    //! Adds a constant value to all elements on diagonal.
    Matrix& add_to_diagonal(const T& val);

    //! Adds a vector of values of to diagonal (element wise).
    Matrix& add_vector_to_diagonal(const Vector<T>& val);

    //! Back substitution.
    /** Solution to linear problem \f[A\vec{x}=\vec{b}\f].  The LU approach is
     *  considered to be effective if the back substitution will be performed
     *  repeatedly. Use linear_solve otherwise.
     *  @note This method does not take in account the conjugation switch.
     *  @param  rhs Right hand side (\f$\vec{b}\f$) of the linear equation. The
     *              values will be overwritten with the solution.
     *  @return Reference onto this instantion.
     */
    Matrix& LU_back_substitution(Vector<T> & rhs);

    //! Back substitution.
    /** Solution to linear problem \f[A\vec{x}=\vec{b}\f].  Linear solve is
     *  effective for one (or only few) time use. Use LU_back_substitution
     *  otherwise.
     *  @note This method does not take in account the conjugation switch.
     *  @param  rhs Right hand side (\f$\vec{b}\f$) of the linear equation. The
     *              values will be overwritten with the solution.
     *  @return Reference onto upated right hand side vector.
     */
    Vector<T>& linear_solve(Vector<T> & rhs);

    //! Sub-matrix vector multiplication.
    /** Performs a matrix vector multiplication using specified square part of
     *  the matrix.
     *  @note This method does not take in account the conjugation switch.
     *  @param  rhs     Input/output vector for multiplication. Values will be
     *                  overwritten by computed result.
     *  @param  num_rows        Total number of rows in the submatrix.
     *  @param  num_columns     Total number of columns in the submatrix.
     *  @param  row_shift       Index of the submatrix first row.
     *  @param  column_shift    Index ot the submatrix first column.
     *  @return Reference onto upated right hand side vector.
     */
    Vector<T>& sub_matrix_multiplication( Vector<T>& rhs,
                                          blas_int num_rows,
                                          blas_int num_columns,
                                          blas_int row_shift=0,
                                          blas_int column_shift=0 ) const;

    //! Generalized matrix vector multiplication.
    /** Interface to the optimized BLAS level 2 routine (gemv) on the Matrix
     *  and Vector classes. The operation computes \f$ \vec{y} \f$ form the
     *  operation \f[ \vec{y} = \alpha A \vec{x} + \vec{y}, \f] where
     *  \f$\vec{x}, \vec{y}\f$ represents the vectors of source values.
     *  @param  alpha   Scalar value scaling the matrix vector term
     *                  \f$A\vec{x}\f$.
     *  @param  X       Input values for matrix multiplication.
     *  @param  beta    Scalar value scaling the input vector \f$\vec{y}\f$.
     *  @param  Y       Input/output vector of values, overwritten by the
     *                  result.
     *  @return Reference onto updated output vector \f$\vec{y}\f$.
     */
    Vector<T>& gemv( const T alpha,
                     const Vector<T> &X,
                     const T beta,
                     Vector<T> &Y ) const;

    //! Matrix eigen system.
    /** Solution to matrix eigenvalue problem stored in designated class.
     *  @note   Method does not take conjugation switch in account.
     *  @return New instantion of EigenSystem class.
     */
    EigenSystem<T> get_eigen_system() const;

 // storage

    //! Save matrix elements into text file.
    /** The implementation of this method is dependent on the template
     *  specified type.
     *  @param  name    String with file name and relative path.
     */
    void save(const char* name) const;

    //! Save matrix elements into a text file with extra first column.
    /** The implementation of this method is dependent on the template
     *  specified type.
     *  @param  range   First column of the resulting file, size must be equal
     *                  to rows size of the matrix.
     *  @param  name    String with file name and relative path.
     */
    void save(const Vector<def_float>& range, const char* name) const;
};
/** @} */
}   // namespace QSCAT
#endif // INCLUDE_QSCAT_MATRIX_H_
