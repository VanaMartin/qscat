#ifndef INCLUDE_QSCAT_VECTOR_H_
#define INCLUDE_QSCAT_VECTOR_H_

#include "blas.h"
#include "common.h"
#include "Storage.h"
#include "Object.h"
#include "ScalarMultiple.h"

namespace QSCAT
{
/** \addtogroup Arrays
* @{  */
//! Vector of given type values.
/** This class wraps some of the most common BLAS level 1 operations.  Some
 * operations are wrapped via opeator overload others are left with general
 * interface for ease of runtime optimization.
 */
template <typename T>
class Vector : public Object, public BinaryStorageInterface
{
 protected:
    //! Size of the vector.
    blas_int size_;
    //! Array of the vector elements.
    T* array_;
    //! Auxiliary pointer.
    /** Points exactly one element beyond array_ end. */
    T* end_;
 private:
    //! Internal initialization procedure.
    /** Allocates internal memory, sets all internal variables.
     *  @param  size    Total number of elements in vector.
     */
    void initialize(blas_int size);

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
 // constructors

    //! Default constructor.
    /** The instantions returns false with init() method. */
    Vector();

    //! Basic constructor.
    /** Creates an instantion of Vector with given size.
     *  @param  size    Total number of elements in vector.
     */
    explicit Vector(blas_int size);

    //! Range constructor.
    /** Creates a vector of equidistant values within given range.
     *  @param  size        Total number of elements in vector.
     *  @param  min         Range start.
     *  @param  max         Range end.
     *  @param  endpoints   If true endpoints are included in range.
     */
    Vector(blas_int size, T min, T max, bool endpoints);

    //! Shallow source constructor.
    /** Creates vector class on top of given pointer. Increases the refernce
     *  counting twice, so the pointer cannot be deallocated.
     */
    Vector(blas_int size, T* source);

    //! Copy constructor (shallow copy).
    /** Performs a shallow copy operation (copies only the pointer and shape
     *  values).
     *  @param  old Instance of the vector to be copied from.
     */
    Vector(const Vector &old);

    //! Destructor.
    /** Decrease the refence count, if total remaining refernces equals zero
     *  cleans all allocated memory.
     */
    ~Vector();

 // accessors (& modifiers)

    //! Size of the vector.
    blas_int get_size() const;

    //  TODO: rename this method to square_norm
    //! Squared norm of the vector.
    /** @return Squared norm of the vector, i.e. result of inner product with
     *          itself.
     */
    T get_norm() const;

    //! Access vector element.
    /** @param  index  Index of desired element.
     *  @return Constant reference to desired element.
     */
    const T& operator[] (blas_int index) const;

    //! Access vector element.
    /** @param  index  Index of desired element.
     *  @return Reference to desired element.
     */
    T& operator[] (blas_int index);

    // TODO: generalize for dimensional interpretation via increment factor
    //! Deep copy of vector part.
    /** Fills the input vector with values from the internal memory.
     *  @param  destination     Destination vector to be overwritten by result.
     *                          Must be of a smaller or equal size.
     *  @param  shift           Offset from the beginning of the vector. Sum of
     *                          offset and destination size must be smaller
     *                          than current vector size.
     *  @return  Reference on the destination Vector.
     */
    Vector& read_sub_vector(Vector& destination, blas_int shift) const;

    // TODO: generalize for dimensional interpretation via increment factor
    //! Deep copy of vector part.
    /** Fills part of the current vector with values from input vector.
     *  @param  source  Vector of source values to be copied from.  Must be of
     *                  a smaller or equal size.
     *  @param  shift   Offset from the beginning of the vector. Sum of offset
     *                  and destination size must be smaller than current
     *                  vector size.
     *  @return  Reference on the updated current vector.
     */
    Vector& write_sub_vector(Vector& source, blas_int shift);

 // modifiers

    //! Fills whole vector with constant value.
    /** @param  constant    Value to be filled with.
     *  @return  Reference on the updated current vector.
     */
    Vector& fill(const T& constant);

    //! Swaps internal members between instantions.
    /** Performs simple swap operation to all internal memebers.
     *  @param  rhs The vector to be swapped with.
     *  @return Reference onto this vector swapped with "rhs".
     */
    Vector& swap(Vector& rhs);

    //! Deep copy operation.
    /** Performs the actual deep copy operation.
     *  @return A new instance of the same object.
     */
    Vector copy() const;

    //! Complex conjugations of elements.
    /** @return Reference onto this instantion. */
    Vector& complex_conjugate();

 // operators

    //! Assignement operation.
    /** @param  tmp     Shallow copy of source vector.
     *  @return Reference on the updated current vector.
     */
    Vector& operator= (Vector tmp);

    //! Inplace addition of scalar.
    /** Performs a fast BLAS addition of one value to whole vector.
     *  @param  scalar  Value to be added.
     *  @return Reference on the updated current vector.
     */
    Vector& operator+= (const T& scalar);

    //! Inplace subtraction of scalar.
    /** Performs a fast BLAS subtraction of one value from whole vector.
     *  @param  scalar  Value to be subtraction.
     *  @return Reference on the updated current vector.
     */
    Vector& operator-= (const T& scalar);

    //! Inplace addition.
    /** Performs fast BLAS addition of a vector to the instantion.
     *  @param  rhs     Vector to be added to the instantion.
     *  @return Reference on the updated current vector.
     */
    Vector& operator+= (const Vector& rhs);

    //! Inplace subtraction.
    /** Performs fast BLAS subtraction of a vector from the instantion.
     *  @param  rhs     Vector to be subtracted from the instantion.
     *  @return Reference on the updated current vector.
     */
    Vector& operator-= (const Vector& rhs);

    //! Multiply by constant.
    /** Performs element wise scaling by constant value.
     *  @param  scalar  Value to scale with.
     *  @return Reference onto this instantion.
     */
    Vector& operator*= (const T& scalar);

    //! Scalar product.
    /** Performs fast scalar product operation.
     *  @param  rhs     Vector to be contracted with.
     *  @return Value of the operation result.
     */
    T operator* (const Vector& rhs) const;

    //! Assign the result of a scalar multiplication.
    /** Assigns a vector scaled by scalar value \f$\alpha \vec{x}f$ in one BLAS
     *  operation.
     *  @param  multiple    Source vector and scalar value wrapped in temporary
     *                      class.
     *  @return  Reference on the updated current vector.
     */
    Vector& operator= (ConstScalarMultiple<T,Vector<T> >& multiple);

    //! Inplace addition of vector scaled with a scalar.
    /** Performs fast BLAS addition of a vector multiplied by scaling factor to
     *  the instantion, i.e. computes \f$\vec{y}\f$ from \f[\vec{y} =
     *  \alpha\vec{x}+\vec{y}\f], given the \f$\alpha \vec{x}\f$ as argument.
     *  @param  multiple    Source vector and scalar value wrapped in temporary
     *                      class.
     *  @return  Reference on the updated current vector.
     */
    Vector& operator+= (ConstScalarMultiple<T,Vector<T> >& multiple);

    //! Inplace subtraction of vector scaled with a scalar.
    /** Performs fast BLAS subtraction of a vector multiplied by scaling factor to
     *  the instantion, i.e. computes \f$\vec{y}\f$ from \f[\vec{y} =
     *  -\alpha\vec{x}+\vec{y}\f] and stores it in this instantion, given the
     *  \f$\alpha \vec{x}\f$ as argument and \f$\vec{y}\f$ from initial state.
     *  @param  multiple    Source vector and scalar value wrapped in temporary
     *                      class.
     *  @return Reference on the updated current vector.
     */
    Vector& operator-= (ConstScalarMultiple<T,Vector<T> >& multiple);

 // custom operations

    //! Generalized inplace addition.
    /** Computes \f$\vec{y}\f$ from \f[\vec{y} = \alpha\vec{x}+\vec{y}\f] and
     *  stores it in this instantion, given the \f$\alpha\f$ and \f$\vec{x}\f$
     *  as arguments and \f$\vec{y}\f$ from initial state.
     *  @param  alpha   Scaling factor \f$\alpha\f$.
     *  @param  rhs     Source vector \f$\vec{x}\f$.
     *  @return Reference on the updated current vector.
     */
    Vector& axpy(const T& alpha, const Vector& rhs);

    //! Sipmlified general inplace addition.
    /** Computes \f$\vec{y}\f$ from \f[\vec{y} = \alpha\vec{x}\f] and
     *  stores it in this instantion, given the \f$\alpha\f$ and \f$\vec{x}\f$
     *  as arguments.
     *  @param  alpha   Scaling factor \f$\alpha\f$.
     *  @param  rhs     Source vector \f$\vec{x}\f$.
     *  @return Reference on the updated current vector.
     */
    Vector& ax(const T& alpha, const Vector& rhs);

    //! Inner product.
    /** Performs inner contraction of two vectors  however the operation is not
     *  positive definite for complex numbers.
     *  @param  rhs     Vector to be contracted with.
     *  @return Value of the operation result.
     */
    T reduction(const Vector& y) const;

    // TODO Investigate inplace approach
    //! Element wise multiplication.
    /** Performs element wise multiplication on two vectors.
     *  @param  rhs     Source vector.
     *  @return Reference on the updated current vector.
     */
    Vector& element_wise_multiplication(const Vector& rhs);

    //  FIXME APPARENTLY BROKEN, NOT PROPERLY DOCUMENTED
    /** Multiply the a subset of current vector elements by elements of source vector.
     *  @param  rhs         Source vector. Its size times striding step plus
     *                      offset must be samller or equal to size of current
     *                      vector.
     *  @param  shift       Offset from current vector start.
     *  @param  increment   Striding step.
     *  @return Reference on the updated current vector.
     */
    Vector& element_wise_sub_multiplication( Vector& rhs,
                                             blas_int shift,
                                             blas_int increment );

    //! Partial assignement.
    /** Performs assignement operation of source vector elements to a subset of
     *  current vector elements.
     *  @param  num_elements    Total of elements to be assigned.
     *  @param  shift           Offset of the first assignment target.
     *  @param  increment       Target striding step.
     *  @param  source          Source vector.
     *  @param  source_shift    Offset of the first element to be assigned.
     *  @param  source_increment    Source striding step.
     *  @return Reference on the updated current vector.
     */
    Vector& partial_assign( blas_int num_elements,
                            blas_int shift,
                            blas_int increment,
                            Vector<T>& source,
                            blas_int source_shift,
                            blas_int source_increment );

    //! Partial scalar product.
    /** Computes scalar product of current vector subset with source vector
     *  subset.
     *  @param  num_elements    Total numbeer of elements in each subset.
     *  @param  pos1            Current vector subset offset.
     *  @param  inc1            Current vector subset striding factor.
     *  @param  P               Source vector.
     *  @param  pos2            Source vector subset offset.
     *  @param  inc2            Source vector subset striding factor.
     *  @return Value of the operation result.
     */
    T partial_dot_product( blas_int num_elements,
                           blas_int pos1,
                           blas_int inc1,
                           const Vector<T>& P,
                           blas_int pos2,
                           blas_int inc2) const;

 // storage methods

    //! Save vector elements into text file.
    /** The implementation of this method is dependent on the template
     *  specified type.
     *  @param  name    String with file name and relative path.
     */
    void save(const char* name) const;

    // TODO unify the methods
    //! Save vector elements into a text file with extra first column.
    /** The implementation of this method is dependent on the template
     *  specified type.
     *  @param  range   First column of the resulting file, size must be equal
     *                  to rows size of the current vector.
     *  @param  name    String with file name and relative path.
     */
    void save(Vector<def_float>& X, const char* name) const;

    // TODO unify the methods
    //! Save vector elements into a text file with extra first column.
    /** The implementation of this method is dependent on the template
     *  specified type.
     *  @param  range   First column of the resulting file, size must be equal
     *                  to rows size of the current vector.
     *  @param  name    String with file name and relative path.
     */
    void save(Vector<dComplex>& X, const char* name) const;

    // Display
    void print() const;
};

//! Scalar multiplication temporary product.
/** Stores the reference on this instantion and the scalar value inside
 *  Vector*scalar multiple as \f$\alpha \vec{x}f$. (see template class
 *  ScalarMultiple)
 *  @param  scalar  Value to scale with.
 *  @return New instantion of temporary wrapping class.
 */
template<typename T>
ScalarMultiple<T, Vector<T> > operator* (const T& scalar, Vector<T>& object);

//! Scalar multiplication temporary product.
/** Stores the reference on this instantion and the scalar value inside
 *  Vector*scalar multiple as \f$\alpha \vec{x}f$. (see template class
 *  ScalarMultiple)
 *  @param  scalar  Value to scale with.
 *  @return New instantion of temporary wrapping class.
 */
template<typename T>
ScalarMultiple<T, Vector<T> > operator* (Vector<T>& object, const T& scalar);

//! Subtraction operation.
/** @return New instantion containing the subtraction result. */
template<typename T>
Vector<T> operator- (Vector<T> lhs, const Vector<T>& rhs)
{
    return lhs -= rhs;
}

//! Addition operation.
/** @return New instantion containing the addition result. */
template<typename T>
Vector<T> operator+ (Vector<T> lhs, const Vector<T>& rhs)
{
    return lhs -= rhs;
}

//! Array of vectors saving to one text file.
/** The vectors are saved as columns inside one file.
 *  @param  N       Length of all vectors.
 *  @param  M       Total of vectors in the array.
 *  @param  X       Array of vectors, must be of size M and each vector of size
 *                  N.
 *  @param  name    String with file name and relative path.
 */
template<typename T>
void SaveMultipleVectors(int N, int M, Vector<T> ** X, const char * name);

/** @}  */
}   // namespace QSCAT
#endif // INCLUDE_QSCAT_VECTOR_H_
