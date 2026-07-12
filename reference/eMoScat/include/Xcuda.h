#ifndef __XCUDA__   // eMoScat internal CUDA extension symbol
    #define __XCUDA__

    enum xCudaStatus;

    //    
    //     Building GPU representation
    //    
    template<typename T>
    CUDA::CudaVector<T> HostToGPU(const ARRAYS::vector<T>& src);

    template<typename T>
    CUDA::CudaMatrix<T> HostToGPU(const ARRAYS::matrix<T>& src);

    template<typename T>
    CUDA::CudaRowCompMatrix<T> HostToGPU(const ARRAYS::RCMatrix<T>& src);

    //
    //  Synchronistation : Host To GPU
    //
    //  Note: The arguments on host are first, the CUDA representation
    //
    template<typename T>
    xCudaStatus cudaPush(const ARRAYS::vector<T>& src, CUDA::CudaVector<T>& dst);

    template<typename T>
    xCudaStatus cudaPush(const ARRAYS::matrix<T>& src, CUDA::CudaMatrix<T>& dst);

    template<typename T>
    xCudaStatus cudaPush(const ARRAYS::RCMatrix<T>& src, CUDA::CudaRowCompMatrix<T>& dst);

    //
    //  Synchronistation : GPU to Host
    //
    //  Note: Same as above, the host arrays are first then the CUDA representation
    //
    template<typename T>
    xCudaStatus cudaPull(ARRAYS::vector<T>& dst, const CUDA::CudaVector<T>& src);

    template<typename T>
    xCudaStatus cudaPull(ARRAYS::matrix<T>& dst, const CUDA::CudaMatrix<T>& src);

    template<typename T>
    xCudaStatus cudaPull(ARRAYS::RCMatrix<T>& dst, const CUDA::CudaRowCompMatrix<T>& src);

    template<typename T, typename Z>
    xCudaStatus cudaPull(FEM_DVR_ECS_2D::grid_vector_2D<T,Z>& dst, const CUDA::CudaVector<Z>& src);

    // Specialized classes for FEM_DVR_ECS

    #include "Xcuda/Xcuda.hpp"

#endif
