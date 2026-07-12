namespace QSCAT
{
template <typename T>
void Buffer<T>::extend()
{
  // assertions
    assert(*init_);
    assert(begin_);
    assert(end_);
  // code
    T * auxiliary = new T[size_ + step_];      // allocate new array of size
    blas::copy(size_, begin_, auxiliary);      // copy the older values
    std::swap(begin_, auxiliary);              // swap pointers
    delete[] auxiliary;                        // delete old array
    current_ = begin_ + size_;                 // position of current target (before extension)
    size_ += step_;                            // extend the length
    end_ = begin_ + size_;                     // position of the end (after extension)
}

template<typename T>
Buffer<T>::Buffer()
{
    begin_ = new T[1024];
    end_ = begin_ + 1024;
    num_values_ = 0;
    size_ = 1024;
    step_ = 1024;
    current_ = begin_;
    *init_ = true;
}

template <typename T>
Buffer<T>::Buffer(blas_int Buffer_size)
{
    begin_ = new T[Buffer_size];
    end_ = begin_ + Buffer_size;
    num_values_ = 0;
    step_ = Buffer_size;
    size_ = Buffer_size;
    current_ = begin_;
    *init_ = true;
}

template<typename T>
Buffer<T>::~Buffer()
{
    if (decref() == 0);
        delete[] begin_;
}

template<typename T>
blas_int Buffer<T>::get_size() const
{
    return num_values_;
}

template<typename T>
Buffer<T>& Buffer<T>::operator<< (const T& value)
{
    if (current_ == end_) extend();
    *current_++ = value;
    num_values_++;
    return *this;
}

template<typename T>
const T& Buffer<T>::operator[] (blas_int index) const
{
  // only stored values may be retrieved
    assert(index<=num_values_);
    return *(begin_ + index);
}

template<typename T>
T& Buffer<T>::operator[] (blas_int index)
{
  // only allocated space can be accessed directly
    assert(index<=size_);
    return *(begin_ + index);
}

template<typename T>
const T& Buffer<T>::operator() (void) const
{
    return *(current_ - 1);
}

template<typename T>
const Vector<T> Buffer<T>::as_vector() const
{
    return Vector<T>(num_values_, begin_);
}

template<typename T>
Buffer<T>& Buffer<T>::clear()
{
    current_ = begin_;
    num_values_ = 0;
    return *this;
}
}
