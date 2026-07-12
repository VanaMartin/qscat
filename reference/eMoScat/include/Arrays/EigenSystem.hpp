namespace QSCAT
{
template<typename T>
void EigenSystem<T>::initialize(blas_int size)
{
    assert(size);
  //
    size_ = size;
    eigen_values_ = new T[size];
    eigen_vectors_ = new T[size*size];
    *init_ = true;
}
template<typename T>
void EigenSystem<T>::clean()
{
    if (init()){
        *init_ = false;
        delete[] eigen_values_;
        delete[] eigen_vectors_;
        size_ = 0;
    }
}
template<typename T>
bool EigenSystem<T>::save_bin_body(std::ofstream & file) const
{
    assert(init());
  //
    if(file.is_open()){
        file.write((char*) &size_, sizeof(blas_int));
        file.write((char*) eigen_values_, size_*sizeof(T));
        file.write((char*) eigen_vectors_, size_*size_*sizeof(T));
        return true;
    } else {
        return false;
    }
}
template<typename T>
bool EigenSystem<T>::read_bin_body(std::ifstream &file)
{
    if(file.is_open()){
        blas_int size;
        file.read((char*) &size, sizeof(blas_int));
      //
        assert(size>0);
      //
        if (!init() || size_ != size) {
            if (init()) {
                clean();
            }
            initialize(size);
        }
        file.read((char*) eigen_values_, size_*sizeof(T));
        file.read((char*) eigen_vectors_, size_*size_*sizeof(T));
    } else {
        return false;
    }
    return true;
}
template<typename T>
EigenSystem<T>::EigenSystem() : Object()
{
    size_ = 0;
    eigen_values_ = 0;
    eigen_vectors_ = 0;
    *init_ = false;
}
template<typename T>
EigenSystem<T>::EigenSystem(blas_int size, const T *source) : Object()
{
    assert(size>0);
  //
    initialize(size);
    blas::copy(size*size, source, eigen_vectors_);
    blas::eigen(eigen_vectors_, eigen_values_, size);
}
template<typename T>
EigenSystem<T>::EigenSystem(const EigenSystem & old) : Object(old) //size_(old.size_), init_(old.init_)
{
    size_ = old.size_;
    eigen_values_ = old.eigen_values_;
    eigen_vectors_ = old.eigen_vectors_;

//    if (init()) {
//        assert(size_>0);
//        initialize(size_);
//
//        blas::copy(size_, old.eigen_values_, eigen_values_);
//        blas::copy(size_*size_, old.eigen_vectors_, eigen_vectors_);
//    } else {
//        eigen_values_ = 0;
//        eigen_vectors_ = 0;
//    }
}
template<typename T>
EigenSystem<T>::~EigenSystem()
{
    if(decref()==0)
        clean();
}
template<typename T>
EigenSystem<T>& EigenSystem<T>::swap(EigenSystem<T>& rhs)
{
    Object::swap(rhs);
    std::swap(size_, rhs.size_);
    std::swap(eigen_values_, rhs.eigen_values_);
    std::swap(eigen_vectors_, rhs.eigen_vectors_);
    //std::swap(init_, rhs.init_);
    return *this;
}
template<typename T>
EigenSystem<T>& EigenSystem<T>::operator= (EigenSystem<T> tmp)
{
    swap(tmp);
    return *this;
}
template<typename T>
EigenSystem<T> EigenSystem<T>::copy() const
{
    assert(init());
  //
    EigenSystem<T> out;
    out.initialize(size_);

    blas::copy(size_, eigen_values_, out.eigen_values_);
    blas::copy(size_*size_, eigen_vectors_, out.eigen_vectors_);

    return out;
}

template<typename T>
blas_int EigenSystem<T>::get_size() const
{
    assert(init());
  //
    return size_;
}
template<typename T>
Vector<T>& EigenSystem<T>::eigen_vector(Vector<T>& destination, blas_int state) const
{
    assert(destination.init());
    assert(destination.get_size() == size_);
    assert(state < size_);
  //
    blas::copy(size_, (eigen_vectors_ + state*size_), &destination[0]);
    return destination;
}
template<typename T>
Vector<T> EigenSystem<T>::eigen_vector(blas_int state) const
{
    assert(state < size_);
  //
    Vector<T> out(size_);
    blas::copy(size_, (eigen_vectors_ + state*size_), &out[0]);
    return out;
}
template<typename T>
const T&  EigenSystem<T>::eigen_value(blas_int i) const
{
    assert(init());
    assert(i < size_);
  //
    return eigen_values_[i];
}
} // namespace QSCAT
