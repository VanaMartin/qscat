    // STD vector class automated wrapper
    class_<std::vector<size_t> >("stdVecUI")
        .def(boost::python::vector_indexing_suite<std::vector<size_t> >())
        .def("size", &std::vector<size_t>::size)
    ;
    class_<std::vector<int> >("stdVecI")
        .def(boost::python::vector_indexing_suite<std::vector<int> >())
        .def("size", &std::vector<size_t>::size)
    ;
    class_<std::vector<unsigned short> >("stdVecSh")
        .def(boost::python::vector_indexing_suite<std::vector<unsigned short> >())
        .def("size", &std::vector<unsigned short>::size)
    ;
    class_<std::vector<float> >("stdVecF")
        .def(boost::python::vector_indexing_suite<std::vector<float> >())
        .def("size", &std::vector<size_t>::size)
    ;
    class_<std::vector<double> >("stdVecD")
        .def(boost::python::vector_indexing_suite<std::vector<double> >())
        .def("size", &std::vector<size_t>::size)
    ;
    class_<std::vector<dComplex> >("stdVecZ")
        .def(boost::python::vector_indexing_suite<std::vector<dComplex> >())
        .def("size", &std::vector<size_t>::size)
    ;

// * * * * * * * VECTORS * * * * * * * * * * *

// this macro expands all necessary methods for all cases of defined vectors (see below)
#define MAKE_VECTOR_CLASS(name, type) \
    class_<name>(#name, init<>()) \
        .def(init<blas_int>()) \
        .def(init<const name&>()) \
        .def("init",        &name::init) \
        .def("size",        &name::get_size) \
        .def("fill",        &name::fill,         return_internal_reference<>()) \
        .def("__getitem__", &name ## GetValue) \
        .def("__setitem__", &name ## SetValue) \
        .def("__iadd__",    static_cast<name& (name::*) (const type&)>(&name::operator+=), return_internal_reference<>(), name ## IA()) \
        .def("__iadd__",    static_cast<name& (name::*) (const name&)>(&name::operator+=), return_internal_reference<>(), name ## IA()) \
        .def("__isub__",    static_cast<name& (name::*) (const type&)>(&name::operator-=), return_internal_reference<>(), name ## IS()) \
        .def("__isub__",    static_cast<name& (name::*) (const name&)>(&name::operator-=), return_internal_reference<>(), name ## IS()) \
        .def("__imul__",    &name::operator*=, return_internal_reference<>()) \
        .def("__mul__",     static_cast<type  (name::*) (const name&) const> (&name::operator*), name ## MUL()) \
        .def("__mul__",     &name::operator*=, return_internal_reference<>()) \
        .def("__rmul__",    &name::operator*=, return_internal_reference<>()) \
        .def("copy",        &name::copy) \
        .def("swap",        &name::swap, return_internal_reference<>()) \
        .def("save",        static_cast<bool (name::*) (const char*) const>(&name::save_binary), name ## Save()) \
        .def("load",        static_cast<bool (name::*) (const char*)>(&name::read_binary), name ## Load()) \
        .def("savetxt",     static_cast<void (name::*) (const char*) const>(&name::save), name ## Stxt()) \
        .def("export",      &name ## Export)
    ;

MAKE_VECTOR_CLASS(iVector, blas_int);
MAKE_VECTOR_CLASS(dVector, def_float);
MAKE_VECTOR_CLASS(zVector, def_comp);


// * * * * * * * EIGEN SYSTEMS * * * * * * * *

#define MAKE_EIGENSYSTEM_CLASS(name,type) \
    class_<name>(#name, init<>()) \
        .def(init<const name&>()) \
        .def("init",        &name::init) \
        .def("size",        &name::get_size) \
        .def("swap",        &name::swap, return_internal_reference<>()) \
        .def("eigenvalue",  &name ## GetValue) \
        .def("save",        static_cast<bool (name::*) (const char*) const> (&name::save_binary), name ## Save()) \
        .def("state",       static_cast<Vector<type> (name::*) (blas_int) const> (&name::eigen_vector), name ## Vector()) \
        .def("export",      &name ## Export) \
    ;

MAKE_EIGENSYSTEM_CLASS(dEigenSystem, def_float);
MAKE_EIGENSYSTEM_CLASS(zEigenSystem, def_comp);


// * * * * * * * MATRICES  * * * * * * * * * *

#define MAKE_MATRIX_CLASS(name, type) \
    class_<name>(#name, init<>()) \
        .def(init<blas_int, blas_int>()) \
        .def(init<const name&>()) \
        .def("init",        &name::init) \
        .def("set_identity", &name::set_identity,    return_internal_reference<>()) \
        .def("fill",        &name::fill,            return_internal_reference<>()) \
        .def("rows",        &name::rows) \
        .def("columns",     &name::columns) \
        .def("__getitem__", &name ## GetItem) \
        .def("__setitem__", &name ## SetItem) \
        .def("__getitem__", &name ## GetItem2) \
        .def("__setitem__", &name ## SetItem2) \
        .def("get_column",  &name::get_column) \
        .def("get_row",     &name::get_row) \
        .def("copy",        &name::copy) \
        .def("swap",        &name::swap, return_internal_reference<>()) \
        .def("__iadd__",    &name::operator+=, return_internal_reference<>()) \
        .def("__isub__",    &name::operator-=, return_internal_reference<>()) \
        .def("__imul__",    static_cast<name& (name::*) (const type&)>(&name::operator*=), return_internal_reference<>(), name ## IM()) \
        .def("__mul__",     static_cast<name (name::*) (const name&) const>(&name::operator*), name ## MUL()) \
        .def("__rmul__",    static_cast<name (name::*) (const name&) const>(&name::operator*), name ## MUL()) \
        .def("__mul__",     static_cast<Vector<type> (name::*) (const Vector<type>&) const>(&name::operator*), name ## MUL()) \
        .def("add_to_diagonal", &name::add_to_diagonal,         return_internal_reference<>()) \
        .def("add_to_diagonal", &name::add_vector_to_diagonal,  return_internal_reference<>()) \
        .def("gemv",            &name::gemv,                    return_internal_reference<>()) \
        .def("LU",              &name::LU_factorize,            return_internal_reference<>()) \
        .def("LU",              &name::LU_back_substitution,    return_internal_reference<>()) \
        .def("solve",           &name::linear_solve,            return_internal_reference<>()) \
        .def("eigen_system",    &name::get_eigen_system) \
        .def("save",        static_cast<bool (name::*) (const char*) const>(&name::save_binary), name ## Save()) \
        .def("load",        static_cast<bool (name::*) (const char*)>(&name::read_binary), name ## Load()) \
        .def("savetxt",     static_cast<void (name::*) (const char*) const>(&name::save), name ## Stxt()) \
        .def("export",      &name ## Export) \
    ;

MAKE_MATRIX_CLASS(dMatrix, def_float)
MAKE_MATRIX_CLASS(zMatrix, def_comp)


// * * * * ROW COMPRESSED MATRICES * * * * * *

#define MAKE_ROW_COMPRESSED_MATRIX_CLASS(name,type) \
    class_<name>(#name, init<>()) \
        .def(init<blas_int, blas_int, blas_int>()) \
        .def(init<const name&>()) \
        .def("init",        &name::init) \
        .def("rows",        &name::rows) \
        .def("columns",     static_cast<blas_int (name::*) () const>(&name::columns), name ## Cols()) \
        .def("nonzeros",    &name::num_nonzeros) \
        .def("row_index",   &name ## GetRowId) \
        .def("row_index",   &name ## SetRowId) \
        .def("columns",     &name ## GetCol) \
        .def("columns",     &name ## SetCol) \
        .def("nonzeros",    &name ## GetNze) \
        .def("nonzeros",    &name ## SetNze) \
        .def("swap",        &name::swap,                        return_internal_reference<>()) \
        .def("add_to_diagonal", &name::add_to_diagonal,         return_internal_reference<>()) \
        .def("add_to_diagonal", &name::add_vector_to_diagonal,  return_internal_reference<>()) \
        .def("__imul__",    &name::operator*=,                  return_internal_reference<>()) \
        .def("__mul__",     static_cast<Vector<type> (name::*) (const Vector<type>&) const>(&name::operator*), name ## MUL()) \
        .def("gemv",        &name::gemv,                    return_internal_reference<>()) \
        .def("LU",          &name::LU_factorize,            return_internal_reference<>()) \
        .def("LU",          &name::LU_back_substitution,    return_internal_reference<>()) \
        .def("save",        static_cast<bool (name::*) (const char*) const>(&name::save_binary), name ## Save()) \
        .def("load",        static_cast<bool (name::*) (const char*)>(&name::read_binary), name ## Load()) \
    ;

MAKE_ROW_COMPRESSED_MATRIX_CLASS(dRCMatrix, def_float)
MAKE_ROW_COMPRESSED_MATRIX_CLASS(zRCMatrix, def_comp)
