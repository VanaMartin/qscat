
    class CudaDVector: public CUDA::CudaVector<double>, public boost::python::wrapper< CUDA::CudaVector<double> >
    {
        public:
            CudaDVector() : CUDA::CudaVector<double>() {}
            CudaDVector(int N) : CUDA::CudaVector<double>(N) {}
            CudaDVector(const CUDA::CudaVector<double>& old) : CUDA::CudaVector<double>(old) {}
            CudaDVector& iAdd(const CUDA::CudaVector<double>& rhs) { *this+=rhs; return *this; }
            CudaDVector& iSub(const CUDA::CudaVector<double>& rhs) { *this-=rhs; return *this; }
            CudaDVector& iMul(const double& alpha) { *this *= alpha; return *this; }
            PyObject* Export()
            {
                numpy_init();
                npy_intp dim = getN();
                double *out=new double[getN()];
                CUDA::check( cudaMemcpy(out, &(*this)[0], getN()*sizeof(double), cudaMemcpyDeviceToHost) );
                return PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE , (void*) out);
            }
    };

    class CudaZVector: public CUDA::CudaVector<dcomp>, public boost::python::wrapper< CUDA::CudaVector<dcomp> >
    {
        public:
            CudaZVector() : CUDA::CudaVector<dcomp>() {}
            CudaZVector(int N) : CUDA::CudaVector<dcomp>(N) {}
            CudaZVector(const CUDA::CudaVector<dcomp>& old) : CUDA::CudaVector<dcomp>(old) {}
            CudaZVector& iAdd(const CUDA::CudaVector<dcomp>& rhs) { *this+=rhs; return *this; }
            CudaZVector& iSub(const CUDA::CudaVector<dcomp>& rhs) { *this-=rhs; return *this; }
            CudaZVector& iMul(const dcomp& alpha) { *this *= alpha; return *this; }
            PyObject* Export()
            {
                numpy_init();
                npy_intp dim = getN();
                dcomp *out=new dcomp[getN()];
                CUDA::check( cudaMemcpy(out, &(*this)[0], getN()*sizeof(dcomp), cudaMemcpyDeviceToHost) );
                return PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE , (void*) out);
            }
    };

    class CudaDMatrix: public CUDA::CudaMatrix<double>, public boost::python::wrapper< CUDA::CudaMatrix<double> >
    {
        public:
            CudaDMatrix() : CUDA::CudaMatrix<double>() {}
            CudaDMatrix(int M, int N) : CUDA::CudaMatrix<double>(M,N) {}
            CudaDMatrix(const CUDA::CudaMatrix<double>& old) : CUDA::CudaMatrix<double>(old) {}
            CudaDMatrix& iMul(const double& alpha) { this->CUDA::CudaArray<double>::operator*=(alpha); return *this; }
            void pyGemv(const double& alpha, const CudaDVector& x, const double& beta, CudaDVector& y) const { this->Gemv(alpha, x, beta, y); }
    };

    class CudaZMatrix: public CUDA::CudaMatrix<dcomp>, public boost::python::wrapper< CUDA::CudaMatrix<dcomp> >
    {
        public:
            CudaZMatrix() : CUDA::CudaMatrix<dcomp>() {}
            CudaZMatrix(int M, int N) : CUDA::CudaMatrix<dcomp>(M,N) {}
            CudaZMatrix(const CUDA::CudaMatrix<dcomp>& old) : CUDA::CudaMatrix<dcomp>(old) {}
            CudaZMatrix& iMul(const dcomp& alpha) { this->CUDA::CudaArray<dcomp>::operator*=(alpha); return *this; }
            void pyGemv(const dcomp& alpha, const CudaZVector& x, const dcomp& beta, CudaZVector& y) const { this->Gemv(alpha, x, beta, y); }
    };

    class CudaDRCMatrix: public CUDA::CudaRowCompMatrix<double>, public boost::python::wrapper< CUDA::CudaRowCompMatrix<double> >
    {
        public:
            CudaDRCMatrix() : CUDA::CudaRowCompMatrix<double>() {}
            CudaDRCMatrix(int M, int N, int NNZ) : CUDA::CudaRowCompMatrix<double>(M,N,NNZ) {}
            CudaDRCMatrix(const CUDA::CudaRowCompMatrix<double>& old) : CUDA::CudaRowCompMatrix<double>(old) {}
            CudaDRCMatrix& iMul(const double& alpha) { this->CUDA::CudaArray<double>::operator*=(alpha); return *this; }
            void pyGemv(const double& alpha, const CudaDVector& x, const double& beta, CudaDVector& y) const { this->Gemv(alpha, x, beta, y); }
    };

    class CudaZRCMatrix: public CUDA::CudaRowCompMatrix<dcomp>, public boost::python::wrapper< CUDA::CudaRowCompMatrix<dcomp> >
    {
        public:
            CudaZRCMatrix() : CUDA::CudaRowCompMatrix<dcomp>() {}
            CudaZRCMatrix(int M, int N, int NNZ) : CUDA::CudaRowCompMatrix<dcomp>(M,N,NNZ) {}
            CudaZRCMatrix(const CUDA::CudaRowCompMatrix<dcomp>& old) : CUDA::CudaRowCompMatrix<dcomp>(old) {}
            CudaZRCMatrix& iMul(const dcomp& alpha) { this->CUDA::CudaArray<dcomp>::operator*=(alpha); return *this; }
            void pyGemv(const dcomp& alpha, const CudaZVector& x, const dcomp& beta, CudaZVector& y) const { this->Gemv(alpha, x, beta, y); }
    };

    void cudaSync() { cudaDeviceSynchronize(); }

    BOOST_PYTHON_FUNCTION_OVERLOADS(DHostToGPU, HostToGPU<double>, 1, 1);
    BOOST_PYTHON_FUNCTION_OVERLOADS(ZHostToGPU, HostToGPU<dcomp>, 1, 1);

    BOOST_PYTHON_FUNCTION_OVERLOADS(cudaPushDPy, cudaPush<double>, 2, 2);
    BOOST_PYTHON_FUNCTION_OVERLOADS(cudaPushZPy, cudaPush<dcomp>, 2, 2);

    BOOST_PYTHON_FUNCTION_OVERLOADS(cudaPullDPy, cudaPull<double>, 2, 2);
    BOOST_PYTHON_FUNCTION_OVERLOADS(cudaPullZPy, cudaPull<dcomp>, 2, 2);

    //BOOST_PYTHON_FUNCTION_OVERLOADS(cudaPullgPy, cudaPull<double,dcomp>, 2, 2);

    xCudaStatus cudaPushGPy(const FEM_DVR_ECS_2D::grid_vector_2D<double,dcomp>& src, CUDA::CudaVector<dcomp>& dst)
    {
        return cudaPush<dcomp>(src.body, dst);
    }
    xCudaStatus cudaPullGPy(FEM_DVR_ECS_2D::grid_vector_2D<double,dcomp>& dst, const CUDA::CudaVector<dcomp>& src)
    {
        return cudaPull<dcomp>(dst.body, src);
    }
