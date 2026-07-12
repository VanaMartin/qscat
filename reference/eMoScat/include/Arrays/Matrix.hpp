namespace QSCAT
{
// internals
template<typename T>
bool Matrix<T>::save_bin_body(std::ofstream &file) const
{
    assert(init());
  //
    if (file.is_open()) {
        file.write((char*) &rows_, sizeof(blas_int));
        file.write((char*) &columns_, sizeof(blas_int));
        file.write((char*)  array_, rows_*columns_*sizeof(T));
        return true;
    } else {
        return false;
    }
}
template<typename T>
bool Matrix<T>::read_bin_body(std::ifstream &file)
{
    if (file.is_open()){
        blas_int rows, columns;
        file.read((char*) &rows, sizeof(blas_int));
        file.read((char*) &columns, sizeof(blas_int));
        if (init() || rows_ != rows || columns != columns_) {
            if (init())
                clean();
            initialize(rows, columns);
        }
        file.read((char*) array_, rows_*columns_*sizeof(T));
    } else {
        return false;
    }
    return true;
}
template<typename T>
void Matrix<T>::initialize(blas_int rows, blas_int columns)
{
    assert(rows>0);
    assert(columns>0);
  //
    rows_ = rows;
    columns_ = columns;
    array_ = new T[rows_*columns_];
    end_ = array_ + rows_ * columns_;
    pivots_ = NULL;
    *init_ = true;
    transposed_ = 'N';
    decomposed_ = false;
}
template<typename T>
void Matrix<T>::clean()
{
    if (init()) {
        transposed_ = 'N';
        rows_ = 0;
        columns_ = 0;
        if (decomposed_){
            delete[] pivots_;
        }
        delete[] array_;
        array_ = 0;
        end_ = 0;
    }
    *init_ = false;
    decomposed_ = false;
}

// constructors
template<typename T>
Matrix<T>::Matrix (blas_int rows, blas_int columns) : Object()
{
    assert(rows>0);
    assert(columns>0);
  //
    initialize(rows, columns);
}
template<typename T>
Matrix<T>::Matrix(const Matrix & old):
    Object(old),
    rows_(old.rows_),
    columns_(old.columns_),
    array_(old.array_),
    end_(old.end_),
    pivots_(old.pivots_),
    transposed_(old.transposed_),
    decomposed_(old.decomposed_)
{}
template<typename T>
Matrix<T>::Matrix() : Object()
{
    transposed_ = 'N';
    rows_ = 0;
    columns_ = 0;
    array_ = 0;
    end_ = 0;
    pivots_ = 0;
    *init_ = false;
    decomposed_ = false;
}
template<typename T>
Matrix<T>::~Matrix()
{
    if (decref() == 0)
        clean();
}

// accessors
template <typename T>
blas_int Matrix<T>::get_size() const
{
    return rows_ * columns_;
}
template <typename T>
blas_int Matrix<T>::rows() const
{
    return rows_;
}
template <typename T>
blas_int Matrix<T>::columns() const
{
    return columns_;
}
template <typename T>
T& Matrix<T>::operator[] (const blas_int& index)
{
    assert(index<rows_*columns_);
  //
    return *(array_ + index);
}
template <typename T>
const T& Matrix<T>::operator[] (const blas_int& index) const
{
    assert(index < rows_*columns_);
  //
    return *(array_ + index);
}
template <typename T>
T& Matrix<T>::operator() (const blas_int& row, const blas_int& column)    // NOTE the complex conjugation is not handled
{
    assert(init());
    assert(transposed_!='C');
  //
    if (transposed_ != 'N') {
        assert(transposed_ == 'T' || transposed_ == 'C');
        assert(column < rows_);
        assert(row < columns_);
      //
        return array_[row*rows_ + column];   // Column Major
    }

    assert(row < rows_);
    assert(column < columns_);
  //
    return array_[column*rows_ + row];  // Columns Major
}
template <typename T>
const T& Matrix<T>::operator() (const blas_int& row,const blas_int& column) const     // NOTE the complex conjugation is not handled
{
    assert(init());
    assert(transposed_!='C');
  //
    if (transposed_ != 'N') {
        assert(transposed_ == 'T' || transposed_ == 'C');
        assert(column < rows_);
        assert(row < columns_);
      //
        //return array_[column*columns_ + row];  // Row Major
        return array_[row*rows_ + column];   // Column Major
    }

    assert(row < rows_);
    assert(column < columns_);
  //
    //return array_[row*columns_ + column];  // Row Major
    return array_[column*rows_ + row];  // Columns Major
}
template<typename T>
Vector<T> Matrix<T>::get_row(blas_int index) const
{
    assert(init());
    assert(index < ((transposed_=='N')? rows_ : columns_ ));
  //
    blas_int length = (transposed_=='N')? columns_ : rows_;
    Vector<T> out(length);
    T* begin = array_ + ((transposed_=='N')? index : index*rows_);
    blas_int increment = (transposed_ == 'N')? rows_ : 1;
    blas::copy( length, begin, increment, &out[0], 1);
    return out;
}
template<typename T>
Vector<T> Matrix<T>::get_column(blas_int index) const
{
    assert(init());
    assert(index < ((transposed_=='N')? columns_: rows_ ));
  //
    blas_int length = (transposed_=='N')? rows_ : columns_;
    Vector<T> out(length);
    T* begin = array_ + ((transposed_=='N')? index*rows_ : index);
    blas_int increment = (transposed_ == 'N')? 1 : rows_;
    blas::copy( length, begin, increment, &out[0], 1);
    return out;
}

// modifiers
template<typename T>
Matrix<T>& Matrix<T>::fill(const T& value)
{
    assert(init());
  //
    //blas::copy(rows_*columns_, &value, 0, array_, 1);
    array_[0:rows_*columns_] = value;
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::swap(Matrix<T>& rhs)
{
    Object::swap(rhs);
    std::swap(transposed_, rhs.transposed_);
    std::swap(rows_,rhs.rows_);
    std::swap(columns_,rhs.columns_);
    std::swap(array_, rhs.array_);
    std::swap(end_, rhs.end_);
    std::swap(decomposed_, rhs.decomposed_);
    std::swap(pivots_, rhs.pivots_);
    return *this;
}

//Matrix<T>& Matrix<T>::copy(const Matrix<T>& rhs)
template<typename T>
Matrix<T> Matrix<T>::copy() const
{
    assert(init());
    assert(!decomposed_);    // deprecated case, copying a decomposed matrix (TODO investigate correction)
  //
    Matrix<T> out(rows_, columns_);
    blas::copy(rows_*columns_, array_, out.array_);
    return out;
    //*this = rhs;    // at this level no point of doing otherwise
    //return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::set_identity(blas_int size)
{
    assert(size>0);
  //
    if (init()) clean();
    initialize(size,size);
    fill(T(0));

    T aux = T(1);
    blas::copy(size, &aux, 0, array_, size+1);
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::LU_factorize()
{
    assert(init());
    assert(rows_==columns_);
  //
    if (!decomposed_) {
        pivots_ = new blas_int[rows_];
        blas::lu_factorize(rows_, array_, pivots_);
        decomposed_ = true;
    }
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::conjugate()
{
    assert(init());
  //
    if (transposed_ == 'N')
        transposed_ = 'C';
    else
        transposed_ = 'N';
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::complex_conjugate()
{
    assert(init());
  //
    blas::conj(rows_*columns_, array_);
    return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::inverse()
{
    using std::cout;
    using std::endl;
    assert(init());
    assert(rows_ == columns_);
  //
    Matrix<T> A(rows_, columns_);
    blas::copy(rows_*columns_, array_, A.array_);
    A.LU_factorize();
    Vector<T> x(rows_);
    for (blas_int i=0; i<columns_; ++i) {
        x.fill(0.0);
        x[i] = 1.0;
        A.LU_back_substitution(x);
        blas::copy(rows_, &x[0], 1, &array_[i*rows_], 1);
    }
    return *this;
}


// operators
template <typename T>
Matrix<T>& Matrix<T>::operator= (Matrix<T> tmp)
{
    this->swap(tmp);
    return *this;
}

template <typename T>
Matrix<T> & Matrix<T>::operator+= (const Matrix<T> & rhs)
{
    assert(init());
    assert(rhs.init());
    assert(rows_ == rhs.rows_);
    assert(columns_ == rhs.columns_);
  //
    blas::axpy(rows_*columns_, T(1), rhs.array_, array_);
    return *this;
}

template <typename T>
Matrix<T> & Matrix<T>::operator-= (const Matrix<T> & rhs)
{
    assert(init());
    assert(rhs.init());
    assert(rows_ == rhs.rows_);
    assert(columns_ == rhs.columns_);
  //
    blas::axpy(rows_*columns_, T(-1), rhs.array_, array_);
    return *this;
}
template <typename T>
Matrix<T> Matrix<T>::operator* (const Matrix<T> & rhs) const
{
    assert(init());
    assert(rhs.init());
    assert( ((transposed_=='N')? columns_ : rows_) == ((rhs.transposed_)? rhs.rows_ : rhs.columns_) );
  //
    Matrix<T> tmp( (transposed_=='N')? rows_ : columns_ , (rhs.transposed_=='N')? rhs.columns_ : rhs.rows_ );
    blas::matrix_matrix(transposed_, rhs.transposed_, rows_, rhs.columns_, columns_, T(1), T(0), array_, rhs.array_, tmp.array_);
    return tmp;
}
template <typename T>
Matrix<T>& Matrix<T>::operator*= (const Matrix<T> & rhs)
{
    assert(init());
    assert(rhs.init());
    assert( ((transposed_=='N')? columns_ : rows_) == ((rhs.transposed_)? rhs.rows_ : rhs.columns_) );
  //
    Matrix<T> tmp = (*this) * rhs;
    swap(tmp);
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator*= (const T& scalar)
{
    assert(init());
  //
    blas::scale(rows_*columns_, array_, scalar);
    return *this;
}
template <typename T>   // TODO rename & move to accessors
Matrix<T>& Matrix<T>::operator() (blas_int index, const Vector<T>& source, blas_int axis)
{
    assert(axis==0 || axis==1);
    assert(init());
    assert(source.init());
    assert( (axis==0)? columns_: rows_ > index);
    assert( (axis==0)? rows_ : columns_ == source.get_size() );
  //
    if (axis) {
        blas::copy( columns_, &source[0], 1, (array_ + index), rows_);  // copy row
    } else {
        blas::copy( rows_, &source[0], 1, (array_ + index*rows_), 1);  // copy column
    }
    return *this;
}
template <typename T>
Vector<T> Matrix<T>::operator* (const Vector<T>& rhs) const
{
    assert(init());
    assert(rhs.init());
    assert(rhs.get_size() == (transposed_=='N')? columns_ : rows_ );
  //
    Vector<T> tmp( (transposed_=='N')? rows_ : columns_ );
    blas::matrix_vector(transposed_, rows_, columns_, T(1), T(0), array_, &rhs[0] , &tmp[0]);
    return tmp;
}
template <typename T>
Vector<T> Matrix<T>::operator* (ConstScalarMultiple<T, Vector<T> >& rhs) const
{
  // auxiliary refence holders
    T alpha = rhs.scalar();
    const Vector<T>& source = rhs.object();
  //
    assert(init());
    assert(source.init());
    assert( ((transposed_=='N')? rows_ : columns_) == source.get_size() );
  //
    Vector<T> tmp( (transposed_=='N')? columns_ : rows_);
    blas::matrix_vector(transposed_, rows_, columns_, alpha, T(0), array_, &source[0], &tmp[0]);
    return tmp;
}

// custom operations
template <typename T>
Matrix<T>& Matrix<T>::add_to_diagonal(const T& value)
{
    assert(init());
  //
    blas::axpy(min(rows_,columns_), T(1), &value, 0, array_, rows_+1);
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::add_vector_to_diagonal(const Vector<T>& source)
{
    assert(init());
    assert(source.init());
    assert(source.get_size() == min(rows_, columns_));
  //
    blas::axpy(min(rows_,columns_), T(1), &source[0], 1, array_, rows_+1);
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::LU_back_substitution(Vector<T>& rhs)
{
    assert(init());
    assert(rows_==columns_);
  //
    if (!decomposed_) {
        LU_factorize();
    }
    blas::lu_back_subst(transposed_, rows_, array_, &rhs[0], pivots_);
    return *this;
}
template<typename T>
Vector<T> & Matrix<T>::linear_solve(Vector<T>& rhs)
{
    assert(init());
    assert(!decomposed_);
    assert(rows_ == columns_);
  //
    blas::lapack_solve(rows_, 1, array_, &rhs[0]);
    return rhs;
}

template <typename T>
Vector<T>& Matrix<T>::sub_matrix_multiplication(Vector<T>& rhs, blas_int num_rows , blas_int num_columns, blas_int row_shift, blas_int column_shift) const
{
    assert(init());
    assert(rhs.init());
    assert(transposed_ == 'N');     // NIY
    assert(rhs.get_size() >= rows_);    // Assuming both start in the same point (NIY)
    assert(rhs.get_size() >= columns_);
  // auxiliary variables
    blas_int sub_rows = (num_rows < rows_)? num_rows : rows_;
    blas_int sub_columns = (num_columns < columns_)? num_columns : columns_;
  //
    assert(row_shift + sub_rows <= rows_);
    assert(column_shift + sub_columns <= columns_);
  //
    T* begin = array_ + row_shift + column_shift * rows_;
  //
    Vector<T> tmp(rhs.get_size());
    blas::sub_matrix_vector(transposed_, sub_rows, sub_columns, rows_, T(1),T(0), begin, &rhs[column_shift], &tmp[row_shift]);
    rhs = tmp;
    return rhs;
}
 template <typename T>
Vector<T>& Matrix<T>::gemv(const T alpha, const Vector<T>& source, const T beta, Vector<T>& destination) const
{
    assert(init());
    assert(source.init());
    assert(destination.init());
    assert( ((transposed_ == 'N')? columns_ : rows_) == source.get_size() );
    assert( ((transposed_ == 'N')? rows_ : columns_) == destination.get_size() );
  //
    blas::matrix_vector(transposed_, rows_, columns_, alpha, beta, array_, &source[0], &destination[0]);
    return destination;
}

template <typename T>
EigenSystem<T> Matrix<T>::get_eigen_system() const
{
    assert(init());
    assert(rows_ == columns_);
  //
    EigenSystem<T> out(rows_, array_);
    return out;
}
}   // namespace QSCAT
