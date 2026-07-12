#ifndef INCLUDE_QSCAT_EIGEN_SYSTEM_H_
#define INCLUDE_QSCAT_EIGEN_SYSTEM_H_

#include "blas.h"
#include "common.h"
#include "Storage.h"
#include "Arrays/Vector.h"
#include "Arrays/Matrix.h"

namespace QSCAT
{
/** \addtogroup Arrays
* @{ */
// FWD definitions
template<typename T>
class Vector;
template<typename T>
class Matrix;

/// The eigenvalues and the eigen Vectors of the given matrix
/**
 *  The class contains the solver for eigenvalue problem, i.e.  the \f[
 *  A\vec{x}_i = \lambda_i \vec{x}_i \f] for all \f$i\in(0,N)\f$, where N is
 *  the size of the matrix.
 */
template<typename T>
class EigenSystem : public Object, public BinaryStorageInterface
{
 protected:
    //! Size of the eigensystem
    blas_int size_;
    //! Array of energies stored
    T* eigen_values_;
    //! Array of eigenVectors
    T* eigen_vectors_;

 private:
    //! Internal initialization method.
    /**
     *  Allocates space for eigenvalues and eigenvectors.
     */
    void initialize(blas_int size);

    //! Internal cleaning method.
    /**
     *  Cleans the internal varibles (set as default).
     */
    void clean();

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

 public:
    //! Default constructor.
    EigenSystem();

    //! Constructor.
    /**
     *  Calls internal initialization and solves the eigen value problem
     *  via BLAS procedure.
     *  @param  size    Size of the square matrix to be provided.
     *  @param  source  Square array of values.
     */
    EigenSystem(blas_int size, const T *source);

    //! Copy constructor (shallow copy).
    EigenSystem(const EigenSystem& old);

    //! Destructor.
    /** Decrease the refence count, if total remaining refernces equals zero
     * cleans all allocated memory.
     */
    ~EigenSystem();

    //! Swaps internal members between instantions.
    /** Performs simple swap operation to all internal memebers.
     *  @param  rhs The object to be swapped with.
     *  @return Reference onto this object swapped with "rhs".
     */
    EigenSystem& swap(EigenSystem& rhs);

    //! Assign operator.
    /** Assigns all internals via shallow copy and swap.
     *  @param  tmp Shallow copy of the right-hand-side object.
     *  @return Reference onto this object with copied internals.
     */
    EigenSystem& operator= (EigenSystem tmp);

    //! Deep copy operation.
    /** Performs the actual deep copy operation.
     *  @return A new instance of the same object.
     */
    EigenSystem copy() const;

    //! Size of the source square matrix.
    /** The result is also the size of the eigenvectors.
     *  @return The size.
     */
    blas_int get_size() const;

    // TODO How about a shallow copy?
    //! Store whole eigenvector in the provided vector.
    /** Performs a deep copy of the eigenvector corresponding to index number
     *  into given the vector.
     *  @param  destination     Vector to which the values will be copied.
     *  @param  state           Index of the eigenvector to be copied.
     *                          (must be less than matrix size-1)
     *  @return Reference onto the destination vector.
     */
    Vector<T>& eigen_vector(Vector<T>& destination, blas_int state) const;

    //! Retrieve eigenvector.
    /**
     *  Deep copy of the eigenvector corresponding to given index stored in a
     *  new instantion of the Vector class.
     *  @param  state   The index of the eigenvector
     *  @return New instatnion of the Vector with desired eigenvector.
     */
    Vector<T> eigen_vector(blas_int state) const;

    //! Eigenvalue.
    /** @param  i   Index of the desired eignevalue.
     *  @return The desired value.
     */
    const T& eigen_value(blas_int i) const;

    //! Pointer to the array of eigen values.
    const T* eigen_values_pointer() const { return eigen_values_; }

    //! Pointer to the array eigen vector values.
    /** The array is shaped as a square matrix with eignevectors as columns of
     *  the matrix.
     */
    const T* eigen_vectors_pointer() const { return eigen_vectors_; }
};
/** @} */
}   // namespace QSCAT
#endif // INCLUDE_QSCAT_EIGEN_SYSTEM_H_
