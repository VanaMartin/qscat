#include <iostream>
#include "Object.h"

namespace QSCAT
{

int Object::incref() const
{
    return ++(*ref_);
}
int Object::decref() const
{
    return --(*ref_);
}
void Object::swap(Object& rhs)
{
    std::swap(init_, rhs.init_);
    std::swap(ref_, rhs.ref_);
}
Object::Object()
{
    init_ = new bool;
    *init_ = false;
    ref_ = new int;
    *ref_ = 1;
}
Object::Object(const Object& old)
{
    ref_ = old.ref_;
    init_ = old.init_;
    incref();
}
Object::~Object()
{
    if (*ref_ == 0) {
        delete ref_;
        delete init_;
    }
}
int Object::ref() const
{
    return *ref_;
}
bool Object::init() const
{
    return *init_;
}

}
