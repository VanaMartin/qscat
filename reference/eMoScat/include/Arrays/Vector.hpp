#include <fstream>
#include <typeinfo>

namespace QSCAT
{

template<typename T>
void Vector<T>::initialize(blas_int size)
{
    assert(size);   // zero sized arrays are deprecated
  //
    size_ = size;
    array_ = new T[size_];
    end_ = array_ + size_;
  // check result
    assert(array_);
    *init_ = true;
}

template<typename T>
bool Vector<T>::save_bin_body(std::ofstream &file) const
{
    assert(init());
    assert(array_);
  //
    if (file.is_open()) {
      // first write the size of the array
        file.write((char*) &size_, sizeof(blas_int));
      // then write the values
        file.write((char*) array_, size_*sizeof(T));
    } else {
        return false;
    }
    return true;
}

template<typename T>
bool Vector<T>::read_bin_body(std::ifstream &file)
{
    if (file.is_open()){
      // first read the size of the array
        blas_int size;
        file.read((char*) &size, sizeof(blas_int));
      // check consistency
        assert(size);

        if (init() && size != size_) {
            delete[] array_;
            *init_ = false;
        }

        if (!init()) {
            initialize(size);
        }
      // check consistency
        assert(size==size_);
        assert(array_);
      // read values
        file.read((char*) array_, size*sizeof(T));
    } else {
        return false;
    }
    return true;
}

// constructors
template<typename T>
Vector<T>::Vector() : Object()
{
    size_ = 0;
    array_ = 0;
    end_ = 0;
    *init_ = false;
}

template<typename T>
Vector<T>::Vector(blas_int size) : Object()
{
    initialize(size);
}

template<typename T>
Vector<T>::Vector(blas_int size, T min, T max, bool endpoints) : Object()
{
    initialize(size);
  // consistency check
    assert(size==size_);

    T* ptr = array_;    // auxiliary pointer
    if (endpoints){
        T step = (max-min) / T(size - 1);
        T value = min;
        while(ptr != end_) {    // includes min, max
            *ptr++ = value;
            value += step;
        }
    } else {
        T step = (max-min) / T(size);
        T value = min;
        while(ptr != end_) {    // excludes min, max
            value += step;
            *ptr++ = value;
        }
    }
}

template<typename T>
Vector<T>::Vector(blas_int size, T* source)
{
    assert(size > 0);
    assert(source);
  //
    size_ = size;
    array_ = source;
    end_ = source + size;
    *init_ = true;
    incref();   // Important, disables deallocation on destruction
}

template<typename T>
Vector<T>::Vector(const Vector<T>& old):
    Object(old),
    size_(old.size_),
    array_(old.array_),
    end_(old.end_)
{}

template<typename T>
Vector<T>::~Vector()
{
    //if (init_){
    if (decref() == 0 && init()){
      // consistency check
        assert(array_);
        delete[] array_;
        array_ = NULL;
    }
    init_ = false;
}


// accessors (& modifiers)

template<typename T>
blas_int Vector<T>::get_size() const          // returns size of the Vector for compatibility checking
{
    return size_;
}

template<typename T>
T Vector<T>::get_norm() const
{
    assert(init());
  //
    T value;
    value = blas::dotproduct(size_, array_, array_);
    return value;
}

template<typename T>
T& Vector<T>::operator[] (blas_int index)
{
    assert(init());
    assert(index >= 0);
    assert(index < size_);
  //
    return *(array_ + index);
}

template<typename T>
const T&  Vector<T>::operator[] (blas_int index) const
{
    assert(init());
    assert(index >= 0);
    assert(index < size_);
  //
    return *(array_ + index);
}

template<typename T>
Vector<T>& Vector<T>::read_sub_vector(Vector& destination, blas_int shift) const
{
    assert(init());              // permitted only on initialized vector
    assert(destination.init());   // only to initialized vector
    assert(size_ >= shift + destination.size_); // check if there is enough allocated space for the operation
  //
    const T* ptr = array_ + shift;
    blas::copy(destination.size_, ptr, destination.array_);
    return destination;
}

template<typename T>
Vector<T>& Vector<T>::write_sub_vector(Vector& source, blas_int shift)
{
    assert(init());              // permitted only on initialized vector
    assert(source.init());  // only to initialized vector
    assert(size_ >= shift + source.size_); // check if there is enough allocated space for the operation
  //
    T* ptr = array_ + shift;
    blas::copy(source.size_, source.array_, ptr);
    return *this;
}

// modifiers
template<typename T>
Vector<T>& Vector<T>::fill(const T& constant)
{
    assert(init());
  //
    array_[0:size_] = constant;
    return *this;
}
template<typename T>
Vector<T> & Vector<T>::swap(Vector& rhs)
{
    Object::swap(rhs);
    std::swap(size_, rhs.size_);
    std::swap(array_, rhs.array_);
    std::swap(end_, rhs.end_);
    return *this;
}
template<typename T>
Vector<T> Vector<T>::copy() const
{
    assert(init());
  //
    Vector<T> out(size_);
    blas::copy(size_, array_, out.array_);
    return out;
}
template<typename T>
Vector<T>& Vector<T>::complex_conjugate()
{
    assert(init());
  //
    blas::conj(size_, array_);
    return *this;
}


// operators
template<typename T>
Vector<T>& Vector<T>::operator= (Vector tmp)
{
    this->swap(tmp);
    return *this;
}

template<typename T>
Vector<T> & Vector<T>::operator+= (const T& scalar)         // adds scalar to all elements int the vector
{
    assert(init());
  //
    blas::axpy(size_, T(1), &scalar, 0, array_, 1);
    return *this;
}

template<typename T>
Vector<T> & Vector<T>::operator-= (const T& scalar)         // subtracts scalar from all elements in the vector
{
    assert(init());
  //
    blas::axpy(size_, T(-1), &scalar, 0, array_, 1);
    return *this;
}

template<typename T>
Vector<T> & Vector<T>::operator+= (const Vector& rhs)       // adds a vector to the present one
{
    assert(init());
    assert(rhs.init());
    assert(size_ == rhs.size_);
  //
    blas::axpy(size_, T(1), rhs.array_, array_);
    return *this;
}

template<typename T>
Vector<T> & Vector<T>::operator-= (const Vector& rhs)       // subtracts a vector from present one
{
    assert(init());
    assert(rhs.init());
    assert(size_ == rhs.size_);
  //
    blas::axpy(size_, T(-1), rhs.array_, array_);
    return *this;
}

template<typename T>
Vector<T> & Vector<T>::operator*= (const T& scalar)         // multiplies all elements in vector by scalar
{
    assert(init());
  //
    blas::scale(size_, array_, scalar);
    return *this;
}

template<typename T>
T Vector<T>::operator* (const Vector & rhs) const
{
    assert(init());
    assert(rhs.init());
    assert(size_ == rhs.size_);
  //
    return blas::dotproduct(size_, array_, rhs.array_);
}

template<typename T>
Vector<T>& Vector<T>::operator= (ConstScalarMultiple<T, Vector<T> >& multiple)
{
    const Vector<T>& rhs = multiple.object();
  //
    assert(rhs.init());
  //
    if (!init() || rhs.size_ != size_) {
        if (init())
            delete[] array_;
        initialize(rhs.size_);
    }

    blas::copy(size_, rhs.array_, array_);
    blas::scale(size_, array_, multiple.scalar());
    return *this;
}

//template<typename T>
//ScalarMultiple<T, Vector<T> > Vector<T>::operator* (const T& scalar)
//{
//    return ScalarMultiple<T, Vector<T> >(scalar, *this);
//}

template<typename T>
Vector<T>& Vector<T>::operator+= (ConstScalarMultiple<T, Vector<T> >& multiple)
{
    const Vector<T>& rhs = multiple.object();
  //
    assert(init());
    assert(rhs.init());
    assert(size_==rhs.size_);
  //
    blas::axpy(size_, multiple.scalar(), rhs.array_, array_);
    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator-= (ConstScalarMultiple<T, Vector<T> >& multiple)
{
    const Vector<T>& rhs = multiple.object();
  //
    assert(init());
    assert(rhs.init());
    assert(size_==rhs.size_);
  //
    blas::axpy(size_, - multiple.scalar(), rhs.array_, array_);
    return *this;
}

// custom operations
template<typename T>
Vector<T> & Vector<T>::axpy(const T& alpha, const Vector& rhs)
{
    assert(init());
    assert(rhs.init());
    assert(size_ == rhs.size_);
  //
    if (alpha!=T(0)) {
        blas::axpy(size_, alpha, rhs.array_, array_);
    } else {
        fill(T(0));
    }
    return *this;
}

template<typename T>
Vector<T> & Vector<T>::ax(const T& alpha, const Vector& rhs)
{
    assert(init());
    assert(rhs.init());
    assert(size_ == rhs.size_);
  //
    blas::copy(size_, rhs.array_, array_);
    blas::scale(size_, array_, alpha);
    return *this;
}

template<typename T>
T Vector<T>::reduction(const Vector<T>& y) const
{
    assert(init());
    assert(y.init());
    assert(size_ == y.size_);
  //
    return blas::reduct(size_, array_, y.array_);
}

template<typename T>
Vector<T> & Vector<T>::element_wise_multiplication(const Vector<T>& rhs)
{
    assert(init());
    assert(rhs.init());
    assert(size_==rhs.size_);
  //
    T* aux = new T[size_];
    blas::ewxy(size_, array_, rhs.array_, aux);
    blas::swap(size_, aux, array_);
    delete[] aux;
    return *this;
}

template<typename T>
Vector<T> & Vector<T>::element_wise_sub_multiplication(Vector<T>& rhs, blas_int shift, blas_int increment)  // Updates the given Vector by multiplying element wise
{
    assert(init());
    assert(rhs.init());
    assert(size_ >= shift + rhs.size_ * increment );
  //
    T *out = new T[rhs.size_];
    blas::subewxy(rhs.size_, increment * rhs.size_, 0, (array_ + shift), rhs.array_, out);
    std::swap(rhs.array_, out);
    delete[] out;
    return rhs;
}

template<typename T>
Vector<T>& Vector<T>::partial_assign(blas_int num_elements, blas_int shift, blas_int increment, Vector<T>& source, blas_int source_shift, blas_int source_increment)
{
    assert(init());
    assert(source.init());
    assert(source.size_ >= source_shift + num_elements * source_increment);
    assert(size_ >= shift + num_elements * increment);
  //
    blas::copy(num_elements, (source.array_ + source_shift), source_increment, (array_ + shift), increment);
    return *this;
}

template<typename T>
T Vector<T>::partial_dot_product(blas_int num_elements, blas_int pos1, blas_int inc1, const Vector<T>& P, blas_int pos2, blas_int inc2) const  // The projector will be conjugated, i.e. <P|X>
{
    return blas::partial_dotproduct(num_elements, &P[pos2], inc2, &array_[pos1], inc1);
}

template<typename T>
ScalarMultiple<T, Vector<T> > operator* (const T& scalar, Vector<T>& object)
{
    return ScalarMultiple<T, Vector<T> >(scalar, object);
}

template<typename T>
ScalarMultiple<T, Vector<T> > operator* (Vector<T>& object, const T& scalar)
{
    return ScalarMultiple<T, Vector<T> >(scalar, object);
}

}
