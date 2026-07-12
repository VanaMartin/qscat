
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(femGrid2DSave, femGrid2D::save_binary, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(femGrid2DLoad, femGrid2D::read_binary, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(femGrid2DWeight, femGrid2D::wz, 2, 1);

    #define GRID2D_METHOD1(type,name) \
        type femGrid2D ## name(const femGrid2D& obj, blas_int i) { return obj.name(i); }
    GRID2D_METHOD1(def_float, xr)
    GRID2D_METHOD1(def_float, yr)
    GRID2D_METHOD1(def_comp,  xz)
    GRID2D_METHOD1(def_comp,  yz)
    GRID2D_METHOD1(def_comp,  xwz)
    GRID2D_METHOD1(def_comp,  ywz)

    PyObject* femGrid2DExport(const femGrid2D& obj)
    {
        npy_intp dim = obj.get_xsize();
        PyObject* X = PyArray_SimpleNew(1, &dim, NPY_DOUBLE );
        dim = obj.get_ysize();
        PyObject* Y = PyArray_SimpleNew(1, &dim, NPY_DOUBLE );
        double *x = (double*) PyArray_GETPTR1((PyArrayObject*) X, 0);
        double *y = (double*) PyArray_GETPTR1((PyArrayObject*) X, 0);
        blas::copy(obj.get_xsize(), &(obj.xr(0)), x);
        blas::copy(obj.get_ysize(), &(obj.yr(0)), y);
        PyObject* out = PyTuple_New((Py_ssize_t) 2);
        PyTuple_SetItem(out, (Py_ssize_t) 0, X);
        PyTuple_SetItem(out, (Py_ssize_t) 1, Y);
        return out;
    }

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVector2DSave, gVector2D::save_binary, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVector2DLoad, gVector2D::read_binary, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVector2DValue, gVector::f, 1, 1);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gVector2DValue2, gVector::f, 2, 2);

    def_comp gVector2DGetItem(const gVector2D& obj, blas_int i) { return obj[i]; }
    gVector2D& gVector2DSetItem(gVector2D& obj, blas_int i, def_comp val) { obj[i] = val; return obj; }

    PyObject* gVector2DExportRange(const gVector2D& obj, double xa, double xb, blas_int sx, double ya, double yb, blas_int sy)
    {
      //
        if (!obj.init())
            throw std::runtime_error("wrong input");
      //
        npy_intp dim[2];
        dim[0]=sx; dim[1]=sy;
        PyObject* F = PyArray_SimpleNew(2, dim, NPY_CDOUBLE );
        PyObject* X = PyArray_SimpleNew(1, &dim[0], NPY_DOUBLE );
        PyObject* Y = PyArray_SimpleNew(1, &dim[1], NPY_DOUBLE );

        dComplex *ff = (dComplex*) PyArray_GETPTR1((PyArrayObject*) F, 0);
        double   *x = (double*) PyArray_GETPTR1((PyArrayObject*) X, 0);
        double   *y = (double*) PyArray_GETPTR1((PyArrayObject*) Y, 0);

        for (blas_int j=0; j<sx; ++j){
            x[j] = xa + (xb-xa) / (sx-1) * j;
        }
        for (blas_int i=0; i<sy; ++i){
            y[i] = ya + (yb-ya) / (sy-1) * i;
            for (blas_int j=0; j<sx; ++j) {
                ff[j*sy+i] = obj.evaluate(x[j],y[i]);      // Transform to row ordering
            }
        }
        PyObject* out = PyTuple_New((Py_ssize_t) 3);
        PyTuple_SetItem(out, (Py_ssize_t) 0, X);
        PyTuple_SetItem(out, (Py_ssize_t) 1, Y);
        PyTuple_SetItem(out, (Py_ssize_t) 2, F);
        return out;
    }
    gVector2D& gVector2DSetValues(gVector2D& obj, PyObject* source)
    {
        PyArrayObject* src = PyArray_GETCONTIGUOUS((PyArrayObject*) source);
        npy_intp *dims = PyArray_DIMS(src);
        assert(dims[0]==obj.get_xsize());
        assert(dims[1]==obj.get_ysize());
        dComplex* s = (dComplex*) PyArray_DATA(src);
        for (blas_int i=0; i<dims[1]; ++i){
            for (blas_int j=0; j<dims[0]; ++j){
                obj.f(s[j*dims[1]+i], i*dims[0]+j);       // Row -> column ordering
            }
        }
        return obj;
    }

    // TODO
//    class doubleGridVector2dPy : public FEM_DVR_ECS_2D::DoubleGridVector2D, public wrapper<FEM_DVR_ECS_2D::DoubleGridVector2D>
//    {
//        public:
//            doubleGridVector2dPy()
//                : FEM_DVR_ECS_2D::DoubleGridVector2D() {}
//            //doubleGridVector2dPy(
//    };

    PyObject* rcOperator2dExport(const zOperator2D& obj)
    {
        const zRCMatrix& M = obj.body();
        return zRCMatrixExport(M);
    }

    class EProjector2d : public EquidistantProjector2d
    {
     public:
        EProjector2d(const femGrid2D& g, size_t xs, size_t ys, dfloat xa, dfloat xb, dfloat ya, dfloat yb ) : 
            EquidistantProjector2d(g, xs, ys, xa, xb, ya, yb) {}
        PyObject* Export(const gVector2D& psi) {
            assert(psi.init());
            assert(body_.init());
            assert(values_.init());
         //

            body_.gemv(1.0, psi.body(), 0.0, values_);

            npy_intp dim[2];
            dim[0]=x_samples_; dim[1]=y_samples_;
            PyObject* F = PyArray_SimpleNew(2, dim, NPY_CDOUBLE );
            PyObject* X = PyArray_SimpleNew(1, &dim[0], NPY_DOUBLE );
            PyObject* Y = PyArray_SimpleNew(1, &dim[1], NPY_DOUBLE );

            dComplex *ff = (dComplex*) PyArray_GETPTR1((PyArrayObject*) F, 0);
            double   *x = (double*) PyArray_GETPTR1((PyArrayObject*) X, 0);
            double   *y = (double*) PyArray_GETPTR1((PyArrayObject*) Y, 0);

            blas::copy(x_samples_, &x_coordinate_[0], x);
            blas::copy(y_samples_, &y_coordinate_[0], y);
            //blas::copy(y_samples_*x_samples_, &values_[0], ff);
            for (blas_int i=0; i<y_samples_; ++i){
                y[i] = y_coordinate_[i];
            }
            for (blas_int i=0; i<x_samples_; ++i){
                x[i] = x_coordinate_[i];
            }
            for (blas_int i=0; i<y_samples_; ++i){
                for (blas_int j=0; j<x_samples_; ++j) {
                    ff[j*y_samples_+i] = 0.0;
                    for (int n=body_.row_index(i*x_samples_ + j);
                            n<body_.row_index(i*x_samples_ + j + 1); ++n){
                        ff[j*y_samples_+i] += body_.nonzeros(n) * psi[body_.columns(n)];
                    }
                    if (ff[j*y_samples_+i] != values_[j*y_samples_ + i])      // Transform to row ordering
                        cout << "ERR: " << ff[j*y_samples_+i] << " <> " << values_[i*x_samples_ + j] << endl;
                }
            }

            PyObject* out = PyTuple_New((Py_ssize_t) 3);
            PyTuple_SetItem(out, (Py_ssize_t) 0, X);
            PyTuple_SetItem(out, (Py_ssize_t) 1, Y);
            PyTuple_SetItem(out, (Py_ssize_t) 2, F);
            return out;
        }
    };
