#ifndef INCLUDE_QSCAT_BUFFER_H_
#define INCLUDE_QSCAT_BUFFER_H_

#include "blas.h"
#include "common.h"
#include "Object.h"
#include "Arrays/Vector.h"

namespace QSCAT
{
/** \addtogroup Arrays
* @{ */

/// Simple buffer class.
/**
 *  Automatically increases the size if overflown. Can be interpreted as vector
 *  using ShallowVector class.
 */
template<typename T>
class Buffer : public Object
{
    //! Pointer to the allocated array start.
    T * begin_;
    //! Pointer to the end of the allocated array.
    /** Actually points beyond the last array element. */
    T * end_;
    //! Pointer to the next value to be written.
    T * current_;
    //! Total number of values stored in the array.
    blas_int num_values_;
    //! Size of allocation step.
    blas_int step_;
    //! Actual size of the allocated array.
    /** The size is evaluated in terms of elements. */
    blas_int size_;
 private:
    //! Increase the size of the buffer.
    /**
     *  Private internal method for extending the internal array, the array is
     *  reallocated and copied to new destination.
     */
    void extend();
 public:
    //! Default constructor with size 1024.
    Buffer();
    //! Constructor with explicite initial size.
    explicit Buffer(blas_int size);
    //! Destructor.
    ~Buffer();
    //! Store given value at the end of the row.
    /** If necessary, extend the internal array.
     *  @param value    Value to be stored.
     */
    Buffer& operator<< (const T& value);
    //! Access value at i-th position.
    /** Standard vector-like access.
     *  @param  index   index of the desired element. Must be smaller than the
     *                  total number of stored values.
     *  @return A constant reference to the value.
     */
    const T& operator[] (blas_int index) const;
    //! Access value to i-th position.
    /** Standard vector-like access.
     *  @param  index   index of the desired element. Must be smaller than the
     *                  total number of stored values.
     *  @return A reference to the value.
     */
    T& operator[] (blas_int index);
    //! Retrieve the last stored value.
    /** @return A constant reference onto the last value in the row. */
    const T& operator() (void) const;
    //! Simple saving in text format.
    /** Save values into the desired file in text form.
     *  @param filename     Target file.
     */
    void save(const char* filename) const;
    //! Save values with equidistant x-axis range.
    /** The computed range will be stored in first column of the resulting
     *  file.
     *  @param x            First value of the range (included).
     *  @param y            Last value of the range (included).
     *  @param filename     Target file.
     */
    void save_range(def_float x, def_float y, const char* filename) const;
    //! Total number of elements already stored in the buffer.
    blas_int get_size() const;
    //! Reset the buffer without changing its allocated size.
    Buffer& clear();
    //! Returns a shallow copy of data interpreted as vector.
    const Vector<T> as_vector() const;
};

/** @} */
}   // namespace QSCAT
#endif // INCLUDE_QSCAT_BUFFER_H_
