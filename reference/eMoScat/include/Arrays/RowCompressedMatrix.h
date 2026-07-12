#ifndef INCLUDE_QSCAT_ROW_COMPRESSED_MATRIX_H_
#define INCLUDE_QSCAT_ROW_COMPRESSED_MATRIX_H_

#include "blas.h"
#include "common.h"
#include "Arrays/Vector.h"

namespace QSCAT
{
/** \addtogroup Arrays
* @{ */
// FWD definitions
template<typename T>
class Vector;

// FIXME the Hermitean conjugation controller is ignored in many methods
// FIXME Sparse solver from opensource solution, move handlers to BLAS internals.
//!  Row compressed matrix in (csr format).
/**
 *   Only nonzero values are stored as a contiguous vector (nonzeros_). The
 *   comlumn indices are stored in another array (columns_), finally the amount
 *   of nonzero elements in each row is stored in a separate array
 *   (row_index_). The first element in "row index" array denotes the position
 *   of the first element in the "nonzeros" & "columns" array. The last element
 *   of this array contains the total of nonzero elements, i.e. row_index_[0] =
 *   0; row_index_[num_rows_] = num_nonzeros_;
 */
template <typename T>
class RowCompressedMatrix : public Object, public BinaryStorageInterface
{
 protected:
    //! Factorization controller.
    bool factorized_;
    //! Scaling lock.
    bool locked_;
    //! Conjugation controller.
    char transposed_;
    //! Total number of rows.
    blas_int num_rows_;
    //! Total number of columns.
    blas_int num_columns_;
    //! Total number of nonzero elements.
    blas_int num_nonzeros_;
    //! Column indices of nonzero elemnts.
    blas_int* columns_;
    //! Starting index of elements in given row.
    blas_int* row_index_;
    //! Array of nonzero elements.
    T* nonzeros_;

    //! PARDISO auxiliary variable: Solver internal data address
    _MKL_DSS_HANDLE_t* handle_;
    //! PARDISO auxiliary variable: Parameters for the solver
    blas_int* iparm_;
 private:
    //! Initialization of internal variables.
    /** Allocates necessary space for stroring matrix with "num_nonzeros"
     *  number of nonzero elements
     *  @param  num_rows        Total number of rows in the matrix.
     *  @param  num_columns     Total number of columns in the matrix.
     *  @param  num_nonzeros    Total number of nonzero elements in the matrix.
    */
    void initialize( blas_int num_rows,
                     blas_int num_columns,
                     blas_int num_nonzeros );

    //! Clean internal variables.
    void clean();

 protected:
    //! Internal binary saving procedure.
    /** Saves the class structure into binary stream.
     *  @return true on success false on error.
     */
    virtual bool save_bin_body(std::ofstream &file) const;

    //! Internal binary reading procedure.
    /** Reads the class structure from binary stream.
     *  @return true on success false on error.
     */
    virtual bool read_bin_body(std::ifstream &file);
 public:
 // constructors
    //! Basic constructor.
    /** Allocates necessary space and prepares all variables for insertions
     *  @param  num_rows        Total number of rows in the matrix.
     *  @param  num_columns     Total number of columns in the matrix.
     *  @param  num_nonzeros    Total number of nonzero elements in the matrix.
     */
    RowCompressedMatrix( blas_int num_rows,
                         blas_int num_columns,
                         blas_int num_nonzeros );

    //! Constructor with source pointers.
    /** Allocates necessary space and prepares all variables then inserts
     *  provided data via deep copy.
     *  @param  num_rows            Total number of rows in the matrix.
     *  @param  num_columns         Total number of columns in the matrix.
     *  @param  num_nonzeros        Total number of nonzero elements in the
     *                              matrix.
     *  @param  source_columns      Source array of element column indices.
     *                              Must be at least of size num_nonzeros.
     *  @param  source_row_index    Source array of element row index starting
     *                              position. Must be at least of size
     *                              num_rows+1.
     *  @param  source_nonzeros     Source array of nonzero elements of the
     *                              matrix. Must be at least of size
     *                              num_nonzeros.
     */
    RowCompressedMatrix( blas_int num_rows,
                         blas_int num_columns,
                         blas_int num_nonzeros,
                         const blas_int *source_columns,
                         const blas_int *source_row_index,
                         const T *source_nonzeros );

    //! Copy constructor (shallow copy).
    RowCompressedMatrix(const RowCompressedMatrix& old);

    //! Default constructor.
    /** Sets the internal values to default values. */
    RowCompressedMatrix();

    //! Destructor.
    /** Decrease the refence count, if total remaining refernces equals zero
     * cleans all allocated memory.
     */
    ~RowCompressedMatrix();

    //! Deep copy operation.
    /** Performs actual deep copy.
     *  @return New instance of the matrix.
     */
    RowCompressedMatrix copy() const;

 // accessors

    //! Number of rows.
    blas_int rows() const;

    //! Number of columns.
    blas_int columns() const;

    //! Total number of nonzero elements.
    blas_int num_nonzeros() const;

    //! Access matrix element as in flat array.
    /** Retrieves nonzero element in given position index in the nonzeros
     *  array.
     *  @param  i   Position of desired element.
     *  @return Constant reference to desired element.
     */
    const T& nonzeros(const blas_int& i) const;

    //! Access matrix element as in flat array.
    /** Retrieves matrix element in given position index in the nonzeros
     *  array.
     *  @param  i   Position of desired element.
     *  @return Reference to desired element.
     */
    T& nonzeros(const blas_int& i);

    //! Matrix element.
    /** @param  i   Index of element row.
     *  @param  j   Index of element column.
     *  @return Actual nonzero value if present zero otherwise.
     */
    T get_element(int i, int j) const;

    //! Column index.
    /** @param  i   Index of the nonzero element in the nonzeros array.
     *  @return Constant reference to the column index of the desired element.
     */
    const blas_int& columns(const blas_int& i) const;

    //! Column index.
    /** @param  i   Index of the nonzero element in the nonzeros array.
     *  @return Reference to the column index of the desired element.
     */
    blas_int& columns(const blas_int& i);

    //! Row index start.
    /** @param  i   Index of the desired row.
     *  @return Constant reference to index of the row first nonzero element in
     *          the nonzeros array. If index equals total number of rows retur
     *          total number of nozero elements in sparse representation.
     */
    const blas_int& row_index(const blas_int& i) const;

    //! Row index start.
    /** @param  i   Index of the desired row.
     *  @return Reference to index of the row first nonzero element in the
     *          nonzeros array. If index equals total number of rows retur
     *          total number of nozero elements in sparse representation.
     */
    blas_int& row_index(const blas_int& i);

 // modifiers
    //! Modification lock.
    /** After setting prevents some operations on the matrix. */
    RowCompressedMatrix& lock();

    //! Swap instantces.
    RowCompressedMatrix& swap(RowCompressedMatrix& rhs);

    //! Matrix dimensions expansion.
    /** Expands the number of rows and columns of the matrix into larger a
     *  representation and shift the posiotions of nonzero elements accordingly
     *  if necessary.
     *  @param  new_rows        Total number of rows in new representaion. Must
     *                          be at least of current size.
     *  @param  new_columns     Total number of columns in new representaion.
     *                          Must be at leas of current size.
     *  @param  shift_rows      Starting postion of current matrix first row
     *                          new representation.
     *  @param  shift_columns   Starting postion of the current matrix first
     *                          column in new representaion.
     *  @return Reference to this instantion.
     */
    RowCompressedMatrix<T>& expand( blas_int new_rows,
                                    blas_int new_columns,
                                    blas_int shift_rows,
                                    blas_int shift_columns );

    //! Complex conjugations of elements.
    /** This operation is suitable only if one needs the actual complex
     *  conjugated matrix, for hermitean conjugate as input to BLAS routines it
     *  is ineffective to use this method.
     *  @return Reference onto this instantion.
     */
    RowCompressedMatrix<T>& complex_conjugate();

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
    RowCompressedMatrix<T>& conjugate();

 // operators

    //! Multiply by constant.
    /** Performs element wise scaling by constant value of all nonzero
     *  elements.
     *  @param  scalar  Value to scale with.
     *  @return Reference onto this instantion.
     */
    RowCompressedMatrix& operator*= (const T& alpha);

    //! Assignement operation.
    RowCompressedMatrix& operator= (RowCompressedMatrix tmp);

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
    //Vector<T> operator* (const VectorScalarMultiple<T>& rhs) const;
    Vector<T> operator* (ConstScalarMultiple<T, Vector<T> >& rhs) const;

 // General form of addition operation handling incompatible cases

    //! Generalized matrix vector multiplication.
    /** Interface to the optimized sparse BLAS level 2 routine (gemv) on the
     *  RowCompressedMatrix and Vector classes. The operation computes \f$
     *  \vec{y} \f$ form the operation \f[ \vec{y} = \alpha A \vec{x} +
     *  \vec{y}, \f] where \f$\vec{x}, \vec{y}\f$ represents the vectors of
     *  source values.
     *  @note   This operation respects the cojugation operation.
     *  @param  alpha   Scalar value scaling the matrix vector term
     *                  \f$A\vec{x}\f$.
     *  @param  X       Input values for matrix multiplication.
     *  @param  beta    Scalar value scaling the input vector \f$\vec{y}\f$.
     *  @param  Y       Input/output vector of values, overwritten by the
     *                  result.
     *  @return Reference onto updated output vector \f$\vec{y}\f$.
     */
    Vector<T>& gemv( const T& alpha,
                     const Vector<T>& X,
                     const T& beta,
                     Vector<T>& Y) const;

    //! Inplace scale and add operation.
    /** Compares the internal representation and reshape internal variables to
     *  contain all values.
     *  @note   This operation ignores the cojugation operation.
     *  @param  alpha   Scaling factor.
     *  @param  x       Source matrix in row compressed format to be added to
     *                  current matrix.
     *  @return Reference on updated this instantion.
     */
    RowCompressedMatrix<T>& axpy( const T &alpha,
                                  const RowCompressedMatrix<T>&x);

    //! Inplace scale and add operation.
    /** Adds new nonzero elements or adds to the existing elements. Extends all
     *  internal variables to contain newly added values. Assuming the source
     *  nonzero elements to be compatible via rows and columns size.
     *  @note   This operation respects the cojugation operation.
     *  @param  alpha   Scaling factor.
     *  @param  src_num_nonzeros    Total number of nonzero elements to be
     *                              added.
     *  @param  src_nonzeros        Array of new nonzero elements. Must be of
     *                              size src_num_nonzeros.
     *  @param  src_columns         Array of nonzero elements column indices.
     *                              Must be of size src_num_nonzeros.
     *  @param  src_row_index       Array of nonzero elemnts row starting
     *                              position. Must be of size rows+1.
     *  @return Reference on updated this instantion.
     */
    RowCompressedMatrix<T>& axpy( const T &alpha,
                                  const blas_int src_num_nonzeros,
                                  const T *src_nonzeros,
                                  const blas_int *src_columns,
                                  const blas_int *src_row_index );

 // custom operations

    //! Adds a constant value to all elements on diagonal.
    RowCompressedMatrix& add_to_diagonal(const T& alpha);

    //! Adds a vector of values of to diagonal (element wise).
    RowCompressedMatrix& add_vector_to_diagonal(const Vector<T>& val);

 // storage

    //! LU decomposition.
    /** The LU decomposition consists of three parts. First pivotization
     *  constructs the matrix permutations (via row & column swap) than
     *  computation of lower triangular matrix (L) with ones on diagonal and
     *  last the computation of upper triangular matrix (U). The resulting
     *  values are stored within associated handlers and workspaces.
     *  @note   This method assumes the matrix to be square.
     *  @note   This method does not take in account the conjugation switch.
     *  @note   Only Intel implementation currently support this operation.
     *  @return Reference onto this instantion.
     */
    void LU_factorize();

    //! Back substitution.
    /** Solution to linear problem \f[A\vec{x}=\vec{b}\f].
     *  @note   This method assumes the matrix to be square.
     *  @note   This method does not take in account the conjugation switch.
     *  @note   Only Intel implementation currently support this operation.
     *  @param  rhs Right hand side (\f$\vec{b}\f$) of the linear equation. The
     *              values will be overwritten with the solution.
     */
    void LU_back_substitution(Vector<T>& rhs);

    //! Auxiliary save method
    /**  DEPRECATED - will be removed soon. */
    void save(const char* name) const;
};

//! Tensor sum of two row compressed matrices.
/** Computes \f$C\f$a new instantion as \f$ C = A \oplus B = A \otimes 1 + 1
 *  \otimes B \f$.
 *  @param  A   Source matrix in row compressed format.
 *  @param  B   Source matrix in row compressed format.
 *  @return New instantion of row compressed matrix contiaining the result.
 */
template<typename T>
RowCompressedMatrix<T> TensorSum( const RowCompressedMatrix<T>& A,
                                  const RowCompressedMatrix<T>& B );

/** @} */
}   // namespace QSCAT
#endif // INCLUDE_QSCAT_ROW_COMPRESSED_MATRIX_H_
