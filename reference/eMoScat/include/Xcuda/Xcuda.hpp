
    enum xCudaStatus
    { 
        STATUS_OK,
        ERROR_SIZE_DIFFERS,
        CUDA_ERROR,
        BLAS_ERROR,
    };

    //    
    //     Building GPU representation
    //    
    template<typename T>
    CUDA::CudaVector<T> HostToGPU(const ARRAYS::vector<T>& src)
    {
        return CUDA::CudaVector<T>(src.GetSize(), &src[0]);
    }

    template<typename T>
    CUDA::CudaMatrix<T> HostToGPU(const ARRAYS::matrix<T>& src)
    {
        return CUDA::CudaMatrix<T>(src.Size(0), src.Size(1), &src[0]);
    }

    template<typename T>
    CUDA::CudaRowCompMatrix<T> HostToGPU(const ARRAYS::RCMatrix<T>& src)
    {
        return CUDA::CudaRowCompMatrix<T>(src.M(), src.N(), src.NNZ(), &src.NZE(0), &src.C(0), &src.RI(0));
    }

    //
    //  Synchronistation : Host To GPU
    //
    //  Note: The arguments on host are first, the CUDA representation
    //
    template<typename T>
    xCudaStatus cudaPush(const ARRAYS::vector<T>& src, CUDA::CudaVector<T>& dst)
    {
        if (dst.getN()!=src.GetSize()) {
            return ERROR_SIZE_DIFFERS;
        }
        dst.set(&src[0]);
        return STATUS_OK;
    }

    template<typename T>
    xCudaStatus cudaPush(const ARRAYS::matrix<T>& src, CUDA::CudaMatrix<T>& dst)
    {
        if (dst.getM()!=src.Size(0) || dst.getN()!=src.Size(1)){
            return ERROR_SIZE_DIFFERS;
        }
        dst.set(&src[0]);
        return STATUS_OK;
    }

    template<typename T>
    xCudaStatus cudaPush(const ARRAYS::RCMatrix<T>& src, CUDA::CudaRowCompMatrix<T>& dst)
    {
        if (dst.getM()!=src.M() || dst.getN()!=src.N() || dst.getNNZ()!=src.NNZ()){
            return ERROR_SIZE_DIFFERS;
        }
        dst.set(&src.NZE(0), &src.C(0), &src.RI(0));
        return STATUS_OK;
    }

    //
    //  Synchronistation : GPU to Host
    //
    //  Note: Same as above, the host arrays are first then the CUDA representation
    //
    template<typename T>
    xCudaStatus cudaPull(ARRAYS::vector<T>& dst, const CUDA::CudaVector<T>& src)
    {
        if (src.getN()!=dst.GetSize()) {
            return ERROR_SIZE_DIFFERS;
        }
        src.get(&dst[0]);
        return STATUS_OK;
    }

    template<typename T>
    xCudaStatus cudaPull(ARRAYS::matrix<T>& dst, const CUDA::CudaMatrix<T>& src)
    {
        if (src.getM()!=dst.Size(0) || src.getN()!=dst.Size(1)){
            return ERROR_SIZE_DIFFERS;
        }
        src.get(&dst[0]);
        return STATUS_OK;
    }

    template<typename T>
    xCudaStatus cudaPull(ARRAYS::RCMatrix<T>& dst, const CUDA::CudaRowCompMatrix<T>& src)
    {
        if (src.getM()!=dst.M() || src.getN()!=dst.N() || src.getNNZ()!=dst.NNZ()){
            return ERROR_SIZE_DIFFERS;
        }
        src.get(&dst.NZE(0), &dst.C(0), &dst.RI(0));
        return STATUS_OK;
    }
   
    // Grid Vector methods
    template<typename T, typename Z>
    xCudaStatus cudaPull(FEM_DVR_ECS_2D::grid_vector_2D<T,Z>& dst, const CUDA::CudaVector<Z>& src) 
    {
        if (src.getN() != dst.NB() ) 
            return ERROR_SIZE_DIFFERS;
        src.get(&dst[0]);
        return STATUS_OK;
    }
