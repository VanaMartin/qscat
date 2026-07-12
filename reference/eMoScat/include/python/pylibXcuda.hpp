
    class_<CudaDVector>("CudaDVector", init<int>())
        .def(init<const CudaDVector&>())
        .def("__iadd__", &CudaDVector::iAdd, return_internal_reference<>())
        .def("__isub__", &CudaDVector::iSub, return_internal_reference<>())
        .def("axpy",     &CudaDVector::Axpy, return_internal_reference<>())
        .def("__imul__", &CudaDVector::iMul, return_internal_reference<>())
        .def("__mul__",  &CudaDVector::operator*)
        .def("copy",     &CudaDVector::Copy, return_internal_reference<>())
        .def("swap",     &CudaDVector::Swap, return_internal_reference<>())
        .def(self*=other<double>())
        .def("export",   &CudaDVector::Export)
    ;

    class_<CudaZVector>("CudaZVector", init<int>())
        .def(init<const CudaZVector&>())
        .def("__iadd__", &CudaZVector::iAdd, return_internal_reference<>())
        .def("__isub__", &CudaZVector::iSub, return_internal_reference<>())
        .def("axpy",     &CudaZVector::Axpy, return_internal_reference<>())
        .def("__mul__",  &CudaZVector::operator*)
        .def("copy",     &CudaZVector::Copy, return_internal_reference<>())
        .def("swap",     &CudaZVector::Swap, return_internal_reference<>())
        .def(self*=other<dComplex>())
        .def("export",   &CudaZVector::Export)
    ;

    class_<CudaDMatrix>("CudaDMatrix", init<int, int>())
        .def("axpy",     &CudaDMatrix::axpy, return_internal_reference<>())
        .def("__imul__", &CudaDMatrix::iMul, return_internal_reference<>())
        .def("gemv",     &CudaDMatrix::pyGemv)
    ;

    class_<CudaZMatrix>("CudaZMatrix", init<int, int>())
        .def("axpy",     &CudaZMatrix::axpy, return_internal_reference<>())
        .def("__imul__", &CudaZMatrix::iMul, return_internal_reference<>())
        .def("gemv",     &CudaZMatrix::pyGemv)
    ;

    class_<CudaDRCMatrix>("CudaDRCMatrix", init<int, int, int>())
        .def(init<const CudaDRCMatrix&>())
        .def("__imul__", &CudaDRCMatrix::iMul, return_internal_reference<>())
        .def("gemv",     &CudaDRCMatrix::pyGemv)
    ;

    class_<CudaZRCMatrix>("CudaZRCMatrix", init<int, int, int>())
        .def(init<const CudaZRCMatrix&>())
        .def("__imul__", &CudaZRCMatrix::iMul, return_internal_reference<>())
        .def("gemv",     &CudaZRCMatrix::pyGemv)
    ;

    // Building representations
    def(    "buildGpuRep",
            static_cast<CUDA::CudaVector<double> (*)(const ARRAYS::vector<double>&)>(&HostToGPU<double>),
            DHostToGPU())
    ;
    def(    "buildGpuRep",
            static_cast<CUDA::CudaVector<dcomp> (*)(const ARRAYS::vector<dcomp>&)>(&HostToGPU<dcomp>),
            ZHostToGPU())
    ;

    def(    "buildGpuRep",
            static_cast<CUDA::CudaMatrix<double> (*)(const ARRAYS::matrix<double>&)>(&HostToGPU<double>),
            DHostToGPU())
    ;
    def(    "buildGpuRep",
            static_cast<CUDA::CudaMatrix<dcomp> (*)(const ARRAYS::matrix<dcomp>&)>(&HostToGPU<dcomp>),
            ZHostToGPU())
    ;

    def(    "buildGpuRep",
            static_cast<CUDA::CudaRowCompMatrix<double> (*)(const ARRAYS::RCMatrix<double>&)>(&HostToGPU<double>),
            DHostToGPU())
    ;
    def(    "buildGpuRep",
            static_cast<CUDA::CudaRowCompMatrix<dcomp> (*)(const ARRAYS::RCMatrix<dcomp>&)>(&HostToGPU<dcomp>),
            ZHostToGPU())
    ;

    // Synchronistation Host -> GPU
    def(    "cudaPush",
            static_cast<xCudaStatus (*)(const ARRAYS::vector<double>&, CUDA::CudaVector<double>&)>(&cudaPush),
            cudaPushDPy())
    ;
    def(    "cudaPush",
            static_cast<xCudaStatus (*)(const ARRAYS::vector<dcomp>&, CUDA::CudaVector<dcomp>&)>(&cudaPush),
            cudaPushZPy())
    ;
    def(    "cudaPush",
            static_cast<xCudaStatus (*)(const ARRAYS::matrix<double>&, CUDA::CudaMatrix<double>&)>(&cudaPush),
            cudaPushDPy())
    ;
    def(    "cudaPush",
            static_cast<xCudaStatus (*)(const ARRAYS::matrix<dcomp>&, CUDA::CudaMatrix<dcomp>&)>(&cudaPush),
            cudaPushZPy())
    ;
    def(    "cudaPush",
            static_cast<xCudaStatus (*)(const ARRAYS::RCMatrix<double>&, CUDA::CudaRowCompMatrix<double>&)>(&cudaPush),
            cudaPushDPy())
    ;
    def(    "cudaPush",
            static_cast<xCudaStatus (*)(const ARRAYS::RCMatrix<dcomp>&, CUDA::CudaRowCompMatrix<dcomp>&)>(&cudaPush),
            cudaPushZPy())
    ;
    def(    "cudaPush",
            &cudaPushGPy)
    ;
    // Synchronistation GPU -> Host
    def(    "cudaPull",
            static_cast<xCudaStatus (*)(ARRAYS::vector<double>&, const CUDA::CudaVector<double>&)>(&cudaPull),
            cudaPullDPy())
    ;
    def(    "cudaPull",
            static_cast<xCudaStatus (*)(ARRAYS::vector<dcomp>&, const CUDA::CudaVector<dcomp>&)>(&cudaPull),
            cudaPullZPy())
    ;
    def(    "cudaPull",
            static_cast<xCudaStatus (*)(ARRAYS::matrix<double>&, const CUDA::CudaMatrix<double>&)>(&cudaPull),
            cudaPullDPy())
    ;
    def(    "cudaPull",
            static_cast<xCudaStatus (*)(ARRAYS::matrix<dcomp>&, const CUDA::CudaMatrix<dcomp>&)>(&cudaPull),
            cudaPullZPy())
    ;
    def(    "cudaPull",
            static_cast<xCudaStatus (*)(ARRAYS::RCMatrix<double>&, const CUDA::CudaRowCompMatrix<double>&)>(&cudaPull),
            cudaPullDPy())
    ;
    def(    "cudaPull",
            static_cast<xCudaStatus (*)(ARRAYS::RCMatrix<dcomp>&, const CUDA::CudaRowCompMatrix<dcomp>&)>(&cudaPull),
            cudaPullZPy())
    ;
    def(    "cudaPull",
            //static_cast<xCudaStatus (*)(FEM_DVR_ECS_2D::grid_vector_2D<double, dcomp>&, const CUDA::CudaVector<dcomp>&)>(&cudaPull),
            &cudaPullGPy)
    ;


    // Auxiliary

    enum_<xCudaStatus>("xCudaStatus")
        .value("OK",  STATUS_OK)
        .value("ERROR_SIZE_DIFFERS", ERROR_SIZE_DIFFERS)
        .value("CUDA_ERROR", CUDA_ERROR)
        .value("BLAS_ERROR", BLAS_ERROR)
    ;

    def(    "cudaSync",     &cudaSync);
