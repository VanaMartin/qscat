#ifndef _INCLUDE_OBJECT_H
#define _INCLUDE_OBJECT_H

namespace QSCAT {
/// Base class defining any object in QSCAT project
/**
 *  Includes the capability of reference counting and initialization
 *  control. The main methods are same for all objects.
 *
 *  If in the future a multithreaded solution should be implemented,
 *  place all common methods here.
 */
class Object {
 private:
    int* ref_;       //!< reference counter
 protected:
    bool* init_;     //!< initialization controller
 protected:
    /// Increase contter and returns updated state of reference counting.
    int incref() const;
    /// Decrease contter and returns updated state of reference counting.
    int decref() const;
    /// internal swap
    void swap(Object& rhs);
    /// Default constructor
    Object();
    /// Copy constructor
    Object(const Object& old);
    /// Destructor
    ~Object();
 public:
    /// Returns actual state of reference counting.
    int ref() const;
    /// Initialization controller
    bool init() const;
};

}   // namespace QSCAT
#endif // _INCLUDE_OBJECT_H
