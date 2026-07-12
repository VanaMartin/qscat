namespace QSCAT
{

template <typename T>
void RowCompressedMatrix<T>::initialize(blas_int num_rows, blas_int num_columns, blas_int num_nonzeros)
{
    num_rows_ = num_rows;
    num_columns_ = num_columns;
    num_nonzeros_ = num_nonzeros;
    nonzeros_ = new T[num_nonzeros_];
    columns_ = new blas_int[num_nonzeros_];
    row_index_ = new blas_int[num_rows_ + 1];
    transposed_ = 'N';
    factorized_ = false;
    locked_ = false;
    *init_ = true;

    handle_ = NULL;
    iparm_ = NULL;
}

template <typename T>
void RowCompressedMatrix<T>::clean()
{
    if (init()) {
        if (factorized_) {
            blas::Pardiso_clean(handle_, iparm_, num_rows_, nonzeros_, row_index_, columns_);
            delete[] handle_;
            delete[] iparm_;
            factorized_ = false;
        }
        delete[] nonzeros_;
        delete[] columns_;
        delete[] row_index_;
    }
    num_rows_ = 0;
    num_columns_ = 0;
    num_nonzeros_ = 0;
    locked_ = false;
    *init_ = false;
}

template <typename T>
bool RowCompressedMatrix<T>::save_bin_body(std::ofstream &file) const
{
    assert(init());
  //
    if (file.is_open()){
        file.write((char*) &num_rows_, sizeof(blas_int));
        file.write((char*) &num_columns_, sizeof(blas_int));
        file.write((char*) &num_nonzeros_, sizeof(blas_int));
        file.write((char*) columns_, num_nonzeros_*sizeof(blas_int));
        file.write((char*) row_index_, (num_rows_+1)*sizeof(blas_int));
        file.write((char*) nonzeros_, num_nonzeros_*sizeof(T));
        return true;
    } else {
        return false;
    }
}

template <typename T>
bool RowCompressedMatrix<T>::read_bin_body(std::ifstream &file)
{
    if (file.is_open()) {
        blas_int new_rows, new_columns, new_nonzeros;
        file.read((char*) &new_rows, sizeof(blas_int));
        file.read((char*) &new_columns, sizeof(blas_int));
        file.read((char*) &new_nonzeros, sizeof(blas_int));

        if (init()) {
            if (new_nonzeros != num_nonzeros_){
                delete[] columns_;
                delete[] nonzeros_;
                columns_ = new blas_int[new_nonzeros];
                nonzeros_ = new T[new_nonzeros];
                num_nonzeros_ = new_nonzeros;
            }
            if (num_rows_ != new_rows){
                delete[] row_index_;
                row_index_ = new blas_int[new_rows+1];
                num_rows_ = new_rows;
            }
        } else {
            columns_ = new blas_int[new_nonzeros];
            row_index_ = new blas_int[new_rows+1];
            nonzeros_ = new T[new_nonzeros];
            num_rows_ = new_rows;
            num_columns_ = new_columns;
        }

        file.read((char*) columns_, num_nonzeros_*sizeof(blas_int));
        file.read((char*) row_index_, (num_rows_+1)*sizeof(blas_int));
        file.read((char*) nonzeros_, num_nonzeros_*sizeof(T));
        return true;
    } else {
        return false;
    }
}

template <typename T>
RowCompressedMatrix<T>::RowCompressedMatrix(blas_int num_rows, blas_int num_columns, blas_int num_nonzeros) : Object()     // Constructor
{
    assert(num_rows>0);
    assert(num_columns>0);
    assert(num_nonzeros>=0);
  //
    initialize(num_rows, num_columns, num_nonzeros);
}

template <typename T>
RowCompressedMatrix<T>::RowCompressedMatrix(blas_int num_rows, blas_int num_columns, blas_int num_nonzeros, const blas_int *source_columns, const blas_int *source_row_index, const T *source_nonzeros) : Object()
{
    assert(num_rows>0);
    assert(num_columns>0);
    assert(num_nonzeros>=0);
    assert(source_columns);
    assert(source_row_index);
    assert(source_nonzeros);
  //
    initialize(num_rows, num_columns, num_nonzeros);
    blas::copy(num_nonzeros_, source_columns, columns_);
    blas::copy(num_rows_+1, source_row_index, row_index_);
    blas::copy(num_nonzeros_, source_nonzeros, nonzeros_);
}

template <typename T>
RowCompressedMatrix<T>::RowCompressedMatrix(const RowCompressedMatrix & old):
    Object(old),
    num_rows_(old.num_rows_),
    num_columns_(old.num_columns_),
    num_nonzeros_(old.num_nonzeros_),
    nonzeros_(old.nonzeros_),
    columns_(old.columns_),
    row_index_(old.row_index_),
    transposed_(old.transposed_),
    factorized_(old.factorized_),
    locked_(locked_),
    handle_(old.handle_),
    iparm_(old.iparm_)
{}

template <typename T>
RowCompressedMatrix<T>::RowCompressedMatrix()
{
    num_nonzeros_ = 0;
    num_rows_ = 0;
    num_columns_ = 0;
    columns_ = NULL;
    row_index_ = NULL;
    nonzeros_ = NULL;
    handle_ = NULL;
    iparm_ = NULL;
    *init_ = false;
    locked_ = false;
    factorized_ = false;
    transposed_ = 'N';
}

template <typename T>
RowCompressedMatrix<T>::~RowCompressedMatrix()  // Destructor
{
    if(decref() == 0)
        clean();
}

template <typename T>
RowCompressedMatrix<T> RowCompressedMatrix<T>::copy() const
{
    assert(init());
    assert(!locked_);
    assert(!factorized_);
  //
    RowCompressedMatrix out(num_rows_, num_columns_, num_nonzeros_);
    blas::copy(num_nonzeros_, columns_, out.columns_);
    blas::copy(num_rows_+1, row_index_, out.row_index_);
    blas::copy(num_nonzeros_, nonzeros_, out.nonzeros_);
    out.locked_ = false;
    out.factorized_ = false;
    out.handle_ = NULL;
    out.iparm_ = NULL;
    return out;
}

// accessors

template <typename T>
blas_int RowCompressedMatrix<T>::rows() const
{
    assert(init());
  //
    return num_rows_;
}

template <typename T>
blas_int RowCompressedMatrix<T>::columns() const
{
    assert(init());
  //
    return num_columns_;
}
template <typename T>
blas_int RowCompressedMatrix<T>::num_nonzeros() const
{
    return num_nonzeros_;
}
template <typename T>
const T& RowCompressedMatrix<T>::nonzeros(const blas_int& i) const
{
    assert(init());
    assert(i < num_nonzeros_);
  //
    return nonzeros_[i];
}
template <typename T>
T& RowCompressedMatrix<T>::nonzeros(const blas_int& i)
{
    assert(init());
    assert(i < num_nonzeros_);
  //
    return nonzeros_[i];
}
template<typename T>
T RowCompressedMatrix<T>::get_element(int i, int j) const
{
    assert(init());
    assert(i < num_rows_);
    assert(j < num_columns_);
  //
    for (int k=row_index_[i]; k<row_index_[i+1]; ++k) {
        if (j == columns_[k])
            return nonzeros_[k];
    }
    return T(0);
}
template <typename T>
const blas_int& RowCompressedMatrix<T>::columns(const blas_int& i) const
{
    assert(init());
    assert(i < num_nonzeros_);
  //
    return columns_[i];
}
template <typename T>
blas_int& RowCompressedMatrix<T>::columns(const blas_int& i)
{
    assert(init());
    assert(i < num_nonzeros_);
  //
    return columns_[i];
}
template <typename T>
const blas_int& RowCompressedMatrix<T>::row_index(const blas_int& i) const
{
    assert(init());
    assert(i <= num_rows_);
  //
    return row_index_[i];
}
template <typename T>
blas_int& RowCompressedMatrix<T>::row_index(const blas_int& i)
{
    assert(init());
    assert(i <= num_rows_);
  //
    return row_index_[i];
}

// modifiers

template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::lock()
{
    assert(init());
  //
    locked_ = true;
    return *this;
}

template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::swap (RowCompressedMatrix& rhs)
{
    Object::swap(rhs);
    std::swap(factorized_, rhs.factorized_);
    std::swap(transposed_, rhs.transposed_);

    std::swap(num_rows_, rhs.num_rows_);
    std::swap(num_columns_, rhs.num_columns_);
    std::swap(num_nonzeros_, rhs.num_nonzeros_);
    std::swap(columns_, rhs.columns_);
    std::swap(row_index_, rhs.row_index_);
    std::swap(nonzeros_, rhs.nonzeros_);

    std::swap(handle_, rhs.handle_);
    std::swap(iparm_, rhs.iparm_);
    return *this;
}

template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::expand(blas_int new_num_rows, blas_int new_num_columns, blas_int shift_rows, blas_int shift_columns)
{
    assert(init());
    assert(shift_rows>=0);
    assert(shift_columns>=0);
    assert(num_rows_ + shift_rows <= new_num_rows);             // otherwise shrinking is in order
    assert(num_columns_ + shift_columns <= new_num_columns);     // same as above
  //
    if (num_rows_ != new_num_rows) {                            // expanding rows in order
        blas_int *work = new blas_int[new_num_rows+1];
        if (shift_rows) {
            blas_int zero = 0;
            // Zeros for shift
            blas::copy(shift_rows, &zero, 0, work, 1);
            // original values
            blas::copy(num_rows_+1, row_index_, &work[shift_rows]);
            // repeat the last
            if (new_num_rows - num_rows_ - shift_rows)
                blas::copy(new_num_rows - num_rows_ - shift_rows, &work[num_rows_+shift_rows], 0, &work[num_rows_+shift_rows+1], 1);
            std::swap(work, row_index_);
        } else {
            // All original values
            blas::copy(num_rows_+1, row_index_, work);
            // repeat the last
            blas::copy(new_num_rows - num_rows_, &work[num_rows_], 0, &work[num_rows_+1], 1);
            std::swap(work, row_index_);
        }
        delete[] work;
    }
    if (num_columns_!=new_num_columns && shift_columns != 0) {
        // TODO blas method?
        for (blas_int i=0; i<num_nonzeros_; ++i) {
            columns_[i] += shift_columns;
        }
    }
    num_rows_ = new_num_rows;
    num_columns_ = new_num_columns;
    return *this;
}

template<typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::complex_conjugate()
{
    assert(init());
    assert(num_nonzeros_);
  //
    blas::conj(num_nonzeros_, nonzeros_);
    return *this;
}

template<typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::conjugate()
{
    assert(init());
    assert(num_nonzeros_);
  //
    if (transposed_ == 'N')
        transposed_ = 'C';
    else
        transposed_ = 'N';
    return *this;
}

// operators

template <typename T>
RowCompressedMatrix<T> & RowCompressedMatrix<T>::operator*= (const T & alpha)
{
    assert(init());
    assert(num_nonzeros_);
  //
    blas::scale(num_nonzeros_, nonzeros_,alpha);
    return *this;
}

template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::operator= (RowCompressedMatrix tmp)
{
    this->swap(tmp);
    return *this;
}

template <typename T>
Vector<T> RowCompressedMatrix<T>::operator* (const Vector<T>& rhs) const
{
    assert(init());
    assert(rhs.init());
    assert(rhs.get_size() == ((transposed_=='N')? num_columns_:num_rows_));
    assert(!factorized_);   // cannot perform matrix vector operation on factorized matrix
  //
    Vector<T> out((transposed_=='N')? num_rows_:num_columns_);
    gemv(1.0, rhs, 0.0, out);
    return out;
}

template <typename T>
Vector<T> RowCompressedMatrix<T>::operator* (ConstScalarMultiple<T, Vector<T> >& rhs) const
{
    assert(init());
    assert(rhs.object().init());
    assert(rhs.object().get_size() == ((transposed_=='N')? num_columns_:num_rows_));
    assert(!factorized_);   // cannot perform matrix vector operation on factorized matrix
  //
    Vector<T> out((transposed_=='N')? num_rows_:num_columns_);
    gemv(rhs.scalar(), rhs.object(), 0.0, out);
    return out;
}

// general form of operations

template <typename T>
Vector<T>& RowCompressedMatrix<T>::gemv(const T& alpha, const Vector<T>& x, const T& beta, Vector<T>& y) const
{
    assert(init());
    assert(!factorized_);
    assert(x.init());
    assert(y.init());
    assert(x.get_size() == ((transposed_=='N')? num_columns_:num_rows_));
    assert(y.get_size() == ((transposed_=='N')? num_rows_:num_columns_));
  //
    blas::RCmatrix_vector(transposed_, num_rows_, alpha, nonzeros_, row_index_, columns_, &x[0], beta, &y[0]);
    return y;
}

template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::axpy(const T& alpha, const RowCompressedMatrix<T>& x)
{
    assert(init());
    assert(*x.init_);
    assert(num_rows_ == x.num_rows_);
    assert(num_columns_ == x.num_columns_);
  //
    this->axpy(alpha, x.num_nonzeros_, x.nonzeros_, x.columns_, x.row_index_);
    return *this;
}

template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::axpy(const T &alpha, const blas_int src_num_nonzeros, const T *src_nonzeros, const blas_int *src_columns, const blas_int *src_row_index)
{
    assert(init());
    // TODO add additional arguments to assure the rows & columns equivalence
  //
    // First determine number of new nonzero values which are not present in the compressed vector
    blas_int b_size = num_nonzeros_ + src_num_nonzeros;      // buffer size hold all nze and new values
    blas_int idb = 0;                                        // index in buffer

    blas_int *nRI = new blas_int[num_rows_+1];                    // new row index
    T *buffer = new T[b_size];                          // buffer for nonzero values (making enough space for all values to be written)
    blas_int *ibuffer = new blas_int[b_size];                     // buffer for columns

    blas_int idx = 0;                            // index in nze
    blas_int idy = 0;                            // index in vals
    for (blas_int row=0; row<num_rows_; ++row) {
        nRI[row] = idb;
        while ( idx<row_index_[row+1] ) {         // loop through all possible values in given segment of nze
            if ( idy < src_row_index[row+1] ) {     // if there are still some values in vals to be added
                if ( columns_[idx] == src_columns[idy] ) {
                    if (alpha==T(1)) {
                        buffer[idb] = src_nonzeros[idy] + nonzeros_[idx];
                    } else {
                        buffer[idb] = alpha * src_nonzeros[idy] + nonzeros_[idx];
                    }
                    ibuffer[idb] = columns_[idx];
                    idx++;
                    idy++;
                } else {
                    if ( columns_[idx] < src_columns[idy] ) {
                        buffer[idb] = nonzeros_[idx];
                        ibuffer[idb] = columns_[idx];
                        idx++;
                    } else {
                        buffer[idb] = alpha * src_nonzeros[idy];
                        ibuffer[idb] = src_columns[idy];
                        idy++;
                    }
                }
            } else {                        // no more values to be added (simple copy of remaining values is at place)
                buffer[idb] = nonzeros_[idx];
                ibuffer[idb] = columns_[idx];
                idx++;
            }
            idb++;
        }
        while ( idy < src_row_index[row+1] ) {  // if some values are remaining in the vals array for given row
            buffer[idb] = alpha * src_nonzeros[idy];
            ibuffer[idb] = src_columns[idy];
            idy++;
            idb++;
        }
    }
    nRI[num_rows_] = idb;

    if (idb==num_nonzeros_) {           // degenerate case of first kind (all values were added to existing ones)
        blas::copy(num_nonzeros_, buffer, nonzeros_);
    } else if (idb == num_nonzeros_ + src_num_nonzeros) {     // degenerate case of second kind (none of values were overlapping)
        std::swap(buffer, nonzeros_);         // no need to reallocate, simple swap
        std::swap(ibuffer, columns_);
        std::swap(nRI, row_index_);
        num_nonzeros_ = idb;
    } else {                            // non-degenerate case (must reallocate)
        delete[] nonzeros_;
        nonzeros_ = new T[idb];
        delete[] columns_;
        columns_ = new blas_int[idb];
        blas::copy(idb, buffer, nonzeros_);
        blas::copy(idb, ibuffer, columns_);
        std::swap(nRI, row_index_);
        num_nonzeros_ = idb;
    }

    delete[] buffer;
    delete[] ibuffer;
    delete[] nRI;

    return *this;
}

// custom operations

template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::add_to_diagonal(const T & alpha)
{
    assert(init());
  //
    for (blas_int i=0; i<num_rows_; ++i){
        for (blas_int j=row_index_[i]; j<row_index_[i+1]; ++j){
            if (columns_[j] == i){
                nonzeros_[j] += alpha;
            }
        }
    }
    return *this;
}
template <typename T>
RowCompressedMatrix<T>& RowCompressedMatrix<T>::add_vector_to_diagonal(const Vector<T> & values)
{
    assert(init());
    assert(values.init());
  //
    for (blas_int i=0; i<num_rows_; ++i){
        for (blas_int j=row_index_[i]; j<row_index_[i+1]; ++j){
            if (columns_[j] == i){
                nonzeros_[j] += values[i];
            }
        }
    }
    return *this;
}

// INTEL only

template <typename T>
void RowCompressedMatrix<T>::LU_factorize()
{
    assert(init());
    assert(transposed_=='N');
  //
    // TODO move to implementaions, case dependent
    if (!factorized_) {
        handle_ = new _MKL_DSS_HANDLE_t[64];
        iparm_ = new blas_int[64];
        //for (int i=0; i<64; ++i){ // Initialization of solver parameters
        //    iparm_[i] = 0;
        //    handle_[i] = 0;
        //}
        for (blas_int i=0; i<=num_rows_; ++i){     // The PARDISO uses one based naming convention: the values must be shifted before use
            row_index_[i] = row_index_[i]+1;
        }
        for (blas_int i=0; i<num_nonzeros_; ++i){    // The PARDISO uses one based naming convention: the values must be shifted before use
            columns_[i] = columns_[i]+1;
        }
        blas::lu_factorize_RCM(handle_, iparm_, num_rows_, nonzeros_, row_index_, columns_);
        factorized_ = true;
    }
}

template <typename T>
void RowCompressedMatrix<T>::LU_back_substitution(Vector<T>& rhs)
{
    assert(init());
    assert(transposed_=='N');
  //
    if (!factorized_) {
        LU_factorize();
    }
    blas::lu_back_subst_RCM(handle_, iparm_, num_rows_, nonzeros_, row_index_, columns_, &rhs[0]);
}

template<typename T>
RowCompressedMatrix<T> TensorSum(const RowCompressedMatrix<T>& A, const RowCompressedMatrix<T>& B)       // Generates outer tensor sum of two opeators, assuming both are given as tensor product with identity
{
  //  C = sum( c|c><c| )= sum( b|1><1|b><b| + a|a><a|1><1| )
    //  ordering: A - internal loop, B - external loop
    //  i.e. A keeps its shape whilst B is distributed
    // Assertions
    assert(A.init());
    assert(B.init());
    // Code
    // First analyze the matrices and determine the total number of non-zero elements
    blas_int ma, mb, nza, nzb, nnz, nb, row, pos;
    ma = A.rows();                      // Number of basis functions of the X-coordinate discretisation
    mb = B.rows();                      // Number of basis functions of the Y-coordinate discretisation
    nza = A.num_nonzeros();                 // Number of nonzero elements in X-grid kinetic energy
    nzb = B.num_nonzeros();                 // Number of nonzero elements in Y-grid kinetic energy
    nb = ma*mb;                     // Number of basis functions of the two-dimensional discretisation
    nnz = nza*mb + nzb*ma - nb;     // Number of nonzero elements in the new composed matrix (same as the initialization above)

    RowCompressedMatrix<T> C(nb,nb,nnz);        // Initializes the arrays inside the RowCompressedMatrix class

    std::cout << "The number of nonzero elements in 2D kinetic energy: " << nnz << std::endl;

    row = 0;        // Row of the new matrix
    pos = 0;        // Position in the arrays NZE and COLUMNS

    C.row_index(0) = 0; // The starting position of the zeroth row in arrays NZE and COLUMNS

    blas_int posb, posa, offset, col;
    for (blas_int rowb=0; rowb<B.rows(); ++rowb){        // Block column numbering
        for (blas_int rowa=0; rowa<A.rows(); ++rowa){    // Block row numbering
            posb = B.row_index(rowb);                   // Position in the arrays B.nze and B.columns
            for (blas_int colb=0; colb<B.row_index(rowb+1)-B.row_index(rowb); ++colb){  // Loop through all non-zero elements in the rowb-th row of B.nze
                offset = (B.columns(posb))*ma;          // The position of the first column of the current submatrix
                if (rowb == B.columns(posb)){           // Diagonal submatrix case
                    posa = A.row_index(rowa);           // Position in the arrays A.nze and A.columns
                    for (blas_int cola=0; cola<A.row_index(rowa+1)-A.row_index(rowa); ++cola){   // Loop through all non-zero elements in the rowa-th row of A.nze
                        C.nonzeros(pos) = A.nonzeros(posa);         // Non-zero as A term
                        col = offset + A.columns(posa);             // It's column position
                        C.columns(pos) = col;
                        if (col == row) {                           // Pure diagonal case
                            C.nonzeros(pos) += B.nonzeros(posb);    // Ads the B term
                        }
                        posa++;                                     // Next A column
                        pos++;                                      // Next element
                    }
                } else {                                            // Off-diagonal submatrix case
                    C.nonzeros(pos) = B.nonzeros(posb);             // Non-zero as B term
                    C.columns(pos) = offset + rowa;                 // It's column position
                    pos++;                                          // Next element
                }
                posb++;     // Next B column
            }
            row++;
            C.row_index(row) = pos;
        }
    }
    return C;
}

template<typename T>
void RowCompressedMatrix<T>::save(const char* name) const
{
    assert(0);
}
}
