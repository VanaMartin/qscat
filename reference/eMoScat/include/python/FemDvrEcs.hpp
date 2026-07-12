
 // Grid

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(femGridSave, femGrid::save_binary, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(femGridLoad, femGrid::read_binary, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(femGridDlp, femGrid::dlp, 2, 0);

    #define GRID_METHOD(type,name) \
        type femGrid ## name(const femGrid& obj) { return obj.name(); }
    GRID_METHOD(def_float, x_pos)
    GRID_METHOD(def_float, x_neg)

    #define GRID_METHOD1(type,name) \
        type femGrid ## name(const femGrid& obj, blas_int i) { return obj.name(i); }
    GRID_METHOD1(blas_int, nel)
    GRID_METHOD1(def_float, xr)
    GRID_METHOD1(def_float, ar)
    GRID_METHOD1(def_comp, x)
    GRID_METHOD1(def_comp, w)
    GRID_METHOD1(def_comp, aaz)

    PyObject* femGridExport(const femGrid& g)
    {
        npy_intp dim = g.nb();
        PyObject* out = PyArray_SimpleNew(1,&dim,NPY_COMPLEX128);
        def_comp *o = (def_comp*) PyArray_GETPTR1((PyArrayObject*) out, 0);
        blas::copy(g.nb(), &g.x(0), o );
        return out;
    }
    // TODO from dictionary
    femGrid gridFromString(const std::string& str)
    {
        pjvalue cfg;
        std::string err = picojson::parse(cfg, str);
        return grid_from_parameters(cfg);
    }

 // Grid vector

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVectorBody, gVector::body, 0, 0);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVectorSave, gVector::save_binary, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVectorLoad, gVector::read_binary, 2, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVectorValue, gVector::f, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVectorValue2, gVector::f, 2, 2);

    def_comp gVectorGetItem(const gVector& obj, blas_int i) { return obj[i]; }
    gVector& gVectorSetItem(gVector& obj, blas_int i, def_comp val) { obj[i] = val; return obj; }

    PyObject* gVectorExportSimple(const gVector& obj)
    {
        npy_intp dim = obj.get_size();
        PyObject* out = PyArray_SimpleNew(1,&dim, NPY_CDOUBLE);
        dComplex *o = (dComplex*) PyArray_GETPTR1((PyArrayObject*) out,0);
        for (int i=0; i<obj.get_size(); ++i)
            *o++ = obj.f(i);
        //blas::copy(obj.get_size(), &obj[0], o);
        return out;
    }
    PyObject* gVectorExportRange(const gVector& obj, const double a, const double b, const blas_int s)
    {
        npy_intp dim = s;
        PyObject* F = PyArray_SimpleNew(1, &dim, NPY_CDOUBLE);
        PyObject* X = PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
        dComplex *f= (dComplex*) PyArray_GETPTR1((PyArrayObject*) F, 0);
        double* x = (double*) PyArray_GETPTR1((PyArrayObject*) X, 0);

        for (int i=0; i<s; ++i){
            x[i] = a + (b-a)/double(s-1)*i;
            f[i] = obj.evaluate(x[i]);
        }

        PyObject* out = PyTuple_New((Py_ssize_t) 2);
        PyTuple_SetItem(out, (Py_ssize_t) 0, X);
        PyTuple_SetItem(out, (Py_ssize_t) 1, F);
        return out;
    }

 // Diagonal Operator

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(zOperatorDSet, zOperatorD::operator=, 1, 1);

    def_comp zOperatorDGetItem(const zOperatorD& obj, blas_int i) { return obj[i]; }
    zOperatorD& zOperatorDSetItem(zOperatorD& obj, blas_int i, def_comp val) { obj[i] = val; return obj; }

 // Full Operator

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(zOperatorFSet, zOperatorF::operator=, 1, 1);

    def_comp zOperatorFGetItem(const zOperatorF& obj, blas_int i) { return obj[i]; }
    zOperatorF& zOperatorFSetItem(zOperatorF& obj, blas_int i, def_comp val) { obj[i] = val; return obj; }

 // Row Compressed Operator

    PyObject* rcOperatorExport(const zOperatorC& obj)
    {
        const zRCMatrix& M = obj.body();
        return zRCMatrixExport(M);
    }
