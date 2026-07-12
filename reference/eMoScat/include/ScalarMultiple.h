#ifndef INCLUDE_QSCAT_VECTOR_SCALAR_MULTIPLE_H_
#define INCLUDE_QSCAT_VECTOR_SCALAR_MULTIPLE_H_

#include "blas.h"
#include "common.h"

namespace QSCAT
{
//! Scalar multiplication temporary wrapper class.
/** Class is designed to postpone multiplication of specified class by wrapping
 *  the class reference together with the scaling factor. Once the assignement
 *  or binary operator is invoked the scaling factor may be used during the
 *  operation.
 */
template<typename T, class C>
class ScalarMultiple
{
    //! Scaling factor.
    T scalar_;
    //! Referece onto the given object.
    C& object_;
 public:
    //! Constructor
    /** @param  scalar  Scaling factor.
     *  @param  object  Reference onto the scaled object.
     */
    ScalarMultiple(T scalar, C& object) : scalar_(scalar), object_(object) {}
    //! Scaling factor.
    /** @return The value of scaling factor. */
    T scalar() const { return scalar_; }
    //! Object constant reference.
    /** @return Constant reference onto the scaled object. */
    const C& object() const { return object_; }
    //! Object reference.
    /** @return Refence onto the scaled object. */
    C& object() { return object_; }
};

//! Scalar multiplication temporary wrapper class.
/** Class is designed to postpone multiplication of specified class by wrapping
 *  the class constant reference together with the scaling factor. Once the
 *  assignement or binary operator is invoked the scaling factor may be used
 *  during the operation.
 */
template<typename T, class C>
class ConstScalarMultiple
{
    //! Scaling factor.
    T scalar_;
    //! Referece onto the given object.
    const C& object_;
 public:
    //! Constructor
    /** @param  scalar  Scaling factor.
     *  @param  object  Reference onto the scaled object.
     */
    ConstScalarMultiple(T scalar, const C& object) :
        scalar_(scalar),
        object_(object)
        {}
    //! Constructor
    /** @param  src  Non-constant refence wrapping class.
     */
    ConstScalarMultiple(ScalarMultiple<T,C>& src) :
        scalar(src.scalar()),
        object_(src.object())
        {}
    //! Scaling factor.
    /** @return The value of scaling factor. */
    T scalar() const { return scalar_; }
    //! Object constant reference.
    /** @return Constant reference onto the scaled object. */
    const C& object() const { return object_; }
};

}   // namespace QSCAT
#endif // INCLUDE_QSCAT_VECTOR_SCALAR_MULTIPLE_H_
