//
//    Wrapper classes for eMoScat ARRAYS objects
//

// * * * * * * * VECTORS * * * * * * * * * * *

#define MAKE_VECTOR_OVERLOADS(name, type) \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Item, name::operator[],  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Eq,   name::operator=,   1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## IA,   name::operator+=,  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## IS,   name::operator-=,  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## MUL,  name::operator*,   1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Save, name::save_binary, 2, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Load, name::read_binary, 2, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Stxt, name::save,        2, 1); \
\
type name ## GetValue (const name& object, blas_int index) { return object[index]; } \
void name ## SetValue (name& object, blas_int index, type value) { object[index] = value; }

MAKE_VECTOR_OVERLOADS(iVector, blas_int)
MAKE_VECTOR_OVERLOADS(dVector, def_float)
MAKE_VECTOR_OVERLOADS(zVector, def_comp)

PyObject* iVectorExport(const iVector& object)
{
    npy_intp dim = object.get_size();
    PyObject* out = PyArray_SimpleNew(1, &dim, NPY_INT64);
    blas_int* dst = (blas_int*) PyArray_GETPTR1((PyArrayObject*) out, 0);

    blas::copy(object.get_size(), &object[0], dst);
    return out;
}
PyObject* dVectorExport(const dVector& object)
{
    npy_intp dim = object.get_size();
    PyObject* out = PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
    double* dst = (double*) PyArray_GETPTR1((PyArrayObject*) out, 0);

    blas::copy(object.get_size(), &object[0], dst);
    return out;
}
PyObject* zVectorExport(const zVector& object)
{
    npy_intp dim = object.get_size();
    PyObject* out = PyArray_SimpleNew(1, &dim, NPY_CDOUBLE);
    dComplex* dst = (dComplex*) PyArray_GETPTR1((PyArrayObject*) out, 0);

    blas::copy(object.get_size(), &object[0], dst);
    return out;
}

// * * * * * * * EIGEN SYSTEMS * * * * * * * *

#define MAKE_EIGENSYSTEM_OVERLOADS(name, type) \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Vector, name::eigen_vector,  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Save,   name::save_binary,   1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Load,   name::read_binary,   1, 1); \
    \
    type name ## GetValue (const name& obj, blas_int i) { return obj.eigen_value(i);  }

MAKE_EIGENSYSTEM_OVERLOADS(dEigenSystem, def_float)
MAKE_EIGENSYSTEM_OVERLOADS(zEigenSystem, def_comp)

PyObject* dEigenSystemExport(const dEigenSystem& object)
{
    PyObject *out = PyTuple_New((Py_ssize_t) 2);
    npy_intp dim[2]; dim[0] = object.get_size(); dim[1] = object.get_size();
    PyObject *arr = PyArray_SimpleNew(2,dim,NPY_DOUBLE);
    PyObject *erg = PyArray_SimpleNew(2,dim,NPY_DOUBLE);

    double *A = (double*) PyArray_GETPTR2((PyArrayObject*) arr, 0, 0);
    double *E = (double*) PyArray_GETPTR1((PyArrayObject*) erg, 0);

    blas::copy(object.get_size(), object.eigen_values_pointer(), E);
    blas::copy(object.get_size()*object.get_size(), object.eigen_vectors_pointer(), A);

    PyTuple_SetItem(out, (Py_ssize_t) 0, erg );
    PyTuple_SetItem(out, (Py_ssize_t) 1, arr );
    return out;
}
PyObject* zEigenSystemExport(const zEigenSystem& object)
{
    PyObject *out = PyTuple_New((Py_ssize_t) 2);
    npy_intp dim[2]; dim[0] = object.get_size(); dim[1] = object.get_size();
    PyObject *arr = PyArray_SimpleNew(2,dim,NPY_CDOUBLE);
    PyObject *erg = PyArray_SimpleNew(1,dim,NPY_CDOUBLE);

    dComplex *A = (dComplex*) PyArray_GETPTR2((PyArrayObject*) arr, 0, 0);
    dComplex *E = (dComplex*) PyArray_GETPTR1((PyArrayObject*) erg, 0);

    blas::copy(object.get_size(), object.eigen_values_pointer(), E);
    blas::copy(object.get_size()*object.get_size(), object.eigen_vectors_pointer(), A);

    PyTuple_SetItem(out, (Py_ssize_t) 0, erg );
    PyTuple_SetItem(out, (Py_ssize_t) 1, arr );
    return out;
}


// * * * * * * * MATRICES  * * * * * * * * * *

#define MAKE_MATRIX_OVERLOADS(name, type) \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## IM, name::operator*=,  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## MUL, name::operator*,  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Save, name::save_binary,  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Load, name::read_binary,  1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Stxt, name::save,        2, 1); \
    \
    type name ## GetItem (const name& obj, const blas_int i) { return obj[i]; } \
    type name ## GetItem2 (const name& obj, const blas_int i, const blas_int j) { return obj(i,j); } \
    void name ## SetItem (name& obj, const blas_int i, const type val) { obj[i] = val; } \
    void name ## SetItem2 (name& obj, const blas_int i, const blas_int j, const type val) { obj(i,j) = val; }

MAKE_MATRIX_OVERLOADS(dMatrix,def_float);
MAKE_MATRIX_OVERLOADS(zMatrix,def_comp);

PyObject* dMatrixExport(const dMatrix& object)
{
    npy_intp dim[2];
    dim[0] = object.rows();
    dim[1] = object.columns();
    PyObject* out = PyArray_SimpleNew(2, dim, NPY_DOUBLE);

    double *A = (double*) PyArray_GETPTR2((PyArrayObject*) out,0,0);
    blas::copy(object.rows()*object.columns(), &object[0], A);
    return out;
}
PyObject* zMatrixExport(const zMatrix& object)
{
    npy_intp dim[2];
    dim[0] = object.rows();
    dim[1] = object.columns();
    PyObject* out = PyArray_SimpleNew(2, dim, NPY_CDOUBLE);

    dComplex *A = (dComplex*) PyArray_GETPTR2((PyArrayObject*) out,0,0);
    blas::copy(object.rows()*object.columns(), &object[0], A);
    return out;
}


// * * * * ROW COMPRESSED MATRICES * * * * * *

#define MAKE_ROW_COMPRESSED_MATRIX_OVERLOADS(name, type) \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Save, name::save_binary,1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Load, name::read_binary,1, 1); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## Cols, name::columns,    1, 0); \
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( name ## MUL, name::operator*,  1, 1); \
    \
    blas_int name ## GetCol (const name& obj, blas_int i) { return obj.columns(i); } \
    void     name ## SetCol (name& obj, blas_int i, blas_int j) { obj.columns(i) = j; } \
    blas_int name ## GetRowId (const name& obj, blas_int i) { return obj.row_index(i); } \
    void     name ## SetRowId (name& obj, blas_int i, blas_int j) { obj.row_index(i) = j; } \
    type     name ## GetNze (const name& obj, blas_int i) { return obj.nonzeros(i); } \
    void     name ## SetNze (name& obj, blas_int i, type j) { obj.nonzeros(i) = j; } \

MAKE_ROW_COMPRESSED_MATRIX_OVERLOADS(dRCMatrix,def_float);
MAKE_ROW_COMPRESSED_MATRIX_OVERLOADS(zRCMatrix,def_comp);

PyObject* zRCMatrixExport(const zRCMatrix& M)
{
    size_t nnz = M.num_nonzeros();

    PyObject* dims = PyTuple_New((Py_ssize_t) 3);
    PyTuple_SetItem( dims, (Py_ssize_t) 0, PyLong_FromLong( M.rows() )  );
    PyTuple_SetItem( dims, (Py_ssize_t) 1, PyLong_FromLong( M.columns() )  );
    PyTuple_SetItem( dims, (Py_ssize_t) 2, PyLong_FromLong( nnz ) );

    npy_intp dim = nnz;
    PyObject* nonzeros = PyArray_SimpleNew(1, &dim, NPY_CDOUBLE);
    dComplex *z = (dComplex*) PyArray_GETPTR1((PyArrayObject*) nonzeros,0);
    blas::copy(nnz, &M.nonzeros(0), z);

    PyObject* columns = PyArray_SimpleNew(1, &dim, NPY_LONGLONG);
    blas_int *bi = (blas_int*) PyArray_GETPTR1((PyArrayObject*) columns,0);
    blas::copy(nnz, &M.columns(0), bi);

    dim = M.rows() + 1;
    PyObject* rowIndex = PyArray_SimpleNew(1, &dim, NPY_LONGLONG);
    bi = (blas_int*) PyArray_GETPTR1((PyArrayObject*) rowIndex,0);
    blas::copy( M.rows() + 1, &M.row_index(0), bi);

    PyObject* out = PyTuple_New((Py_ssize_t) 4);
    PyTuple_SetItem(out, (Py_ssize_t) 0, dims);
    PyTuple_SetItem(out, (Py_ssize_t) 1, nonzeros);
    PyTuple_SetItem(out, (Py_ssize_t) 2, columns);
    PyTuple_SetItem(out, (Py_ssize_t) 3, rowIndex);
    return out;
}
