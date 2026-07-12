using namespace std;

// This is Highly experimental area

namespace FEM_DVR_ECS_3D {

	//
	//	Three dimensional grid type. The ordering is done by repeating Z * ( Y * X ).
	//
	template<typename T, typename Z>
	class FemDvrEcsGrid3D
    {
        bool init_;
		FemDvrEcsGrid<T,Z> gx_;	// Copy of the x coordinate grid
		FemDvrEcsGrid<T,Z> gy_;	// Copy of the x coordinate grid
		FemDvrEcsGrid<T,Z> gz_;	// Copy of the x coordinate grid
		//T m_x, m_y, m_z;				// reduced masses
		blas_int nx_, ny_, nz_;		// Number of basis elements of the original grids
		blas_int nb_;				// Resulting grid number of basis elements
	public:
        FemDvrEcsGrid3D()
        {
            init_ = false;
            nx_ = 0;
            ny_ = 0;
            nz_ = 0;
            nb_ = 0;
        }
        // FIXME Grid order is opposite to the arrangement of the indices -> should we change it
		FemDvrEcsGrid3D(const FemDvrEcsGrid<T,Z> & grid_x, const FemDvrEcsGrid<T,Z> & grid_y, const FemDvrEcsGrid<T,Z> & grid_z)
            :   gx_(grid_x), gy_(grid_y), gz_(grid_z), nx_(0), ny_(0), nz_(0), nb_(0)
        {
            init_ = true;
			gx_ = grid_x;
			gy_ = grid_y;
			gz_ = grid_z;
			nx_ = grid_x.nb();
			ny_ = grid_y.nb();
			nz_ = grid_z.nb();
			nb_ = nx_*ny_*nz_;
		}
        FemDvrEcsGrid3D(const FemDvrEcsGrid3D& src) :
            init_(src.init_), gx_(src.gx_), gy_(src.gy_), gz_(src.gz_),
            nx_(src.nx_), ny_(src.ny_), nz_(src.nz_), nb_(src.nb_)
        {}
        FemDvrEcsGrid3D& operator= (FemDvrEcsGrid3D tmp)
        {
            return this->swap(tmp);
        }
        FemDvrEcsGrid3D& swap(FemDvrEcsGrid3D<T,Z>& rhs)
        {
            std::swap(init_, rhs.init_);
            gx_.swap(rhs.gx_);
            gy_.swap(rhs.gy_);
            gz_.swap(rhs.gz_);
            std::swap(nx_, rhs.nx_);
            std::swap(ny_, rhs.ny_);
            std::swap(nz_, rhs.nz_);
            std::swap(nb_, rhs.nb_);
            return *this;
        }
		~FemDvrEcsGrid3D() {}
		const blas_int& get_xsize() const { return nx_; }
		const blas_int& get_ysize() const { return ny_; }
		const blas_int& get_zsize() const { return nz_; }
		//const blas_int& nb() const { return nb_; }
		const blas_int& get_size() const { return nb_; }
		const T& xr(const blas_int i) const { return gx_.xr(i); }
		const T& yr(const blas_int i) const { return gy_.xr(i); }
		const T& zr(const blas_int i) const { return gz_.xr(i); }
		const T& x(const blas_int i) const { return gx_.x(i); }
		const T& y(const blas_int i) const { return gy_.x(i); }
		const T& z(const blas_int i) const { return gz_.x(i); }
		//const Z& w(const blas_int i) const { return gx_.w( i%nx_ ) * gy_.w( (i/nx_) % ny_) * gz_.w(i/(nx_*ny_)); }
		Z w(const blas_int i) const
        {
            return gx_.w( i%nx_ ) * gy_.w( (i/nx_) % ny_) * gz_.w(i/(nx_*ny_));
        }
		Z w(const blas_int i, const blas_int j, const blas_int k) const
        {
            return gz_.w(i)*gy_.w(j)*gx_.w(k);
        }
        const Z& wx(const blas_int i) const { return gx_.w(i); }
        const Z& wy(const blas_int i) const { return gy_.w(i); }
        const Z& wz(const blas_int i) const { return gz_.w(i); }
        FemDvrEcsGrid<T,Z>& gx() { return gx_; }
        FemDvrEcsGrid<T,Z>& gy() { return gy_; }
        FemDvrEcsGrid<T,Z>& gz() { return gz_; }
	};

    template<typename T, typename Z>
    class GridVector3D
    {
     protected:
        bool init_;
        FemDvrEcsGrid3D<T,Z> grid_;
        Vector<Z> body_;
     private:
        void initialize(const FemDvrEcsGrid3D<T,Z>& grid)
        {
            grid_ = grid;
            body_ = Vector<Z>(grid.get_size());
            init_ = true;
        }
    public:
        GridVector3D() { init_=false; }
        GridVector3D(const FemDvrEcsGrid3D<T,Z>& grid) { initialize(grid); }
        GridVector3D(const GridVector3D<T,Z>& old) : grid_(old.grid_), body_(old.body_), init_(old.init_) {}
        ~GridVector3D() {}
        GridVector3D& swap(GridVector3D<T,Z>& rhs)
        {
            grid_.swap(rhs.grid_);
            body_.swap(rhs.body_);
            std::swap(init_, rhs.init_);
            return *this;
        }
        GridVector3D& operator= (GridVector3D tmp) { return swap(tmp); }
        GridVector3D& operator+=(const GridVector3D& rhs) { body_ += rhs.body_; return *this; }
        GridVector3D& operator-=(const GridVector3D& rhs) { body_ -= rhs.body_; return *this; }
        GridVector3D& operator*=(const Z& alpha) { body_ *= alpha; return *this; }
        Z operator* (const GridVector3D& rhs) const { return body_ * rhs.body_; }
        GridVector3D& axpy(const Z& alpha, const GridVector3D& x) { body_.axpy(alpha, x.body_); return *this; }
        GridVector3D& ax(const Z& alpha, const GridVector3D& x) { body_.ax(alpha, x.body_); return *this; }
        Z reduction(const GridVector3D& rhs) { return body_.reduction(rhs.body_); }
        GridVector3D& element_wise_multiplication(const GridVector3D& rhs) { body_.element_wise_multiplication(rhs.body_); return *this; }
        bool init() const { return init_;  }
        const Z& operator[] (blas_int i) const
        {
            assert(i < grid_.get_size());
          //
            return body_[i];
        }
        Z& operator[] (blas_int i)
        {
            assert(i < grid_.get_size());
          //
            return body_[i];
        }
        const Z& operator() (blas_int i, blas_int j, blas_int k) const
        {
            assert(i < grid_.get_zsize());
            assert(j < grid_.get_ysize());
            assert(k < grid_.get_xsize());
          //
            return body_[(i*grid_.get_ysize() + j)*grid_.get_xsize() + k];
        }
        Z& operator() (blas_int i, blas_int j, blas_int k)
        {
            assert(i < grid_.get_zsize());
            assert(j < grid_.get_ysize());
            assert(k < grid_.get_xsize());
          //
            return body_[(i*grid_.get_ysize() + j) * grid_.get_xsize() + k];
        }
        Z f(blas_int i) const
        {
            return body_[i] / sqrt(grid_.w(i));
        }
        GridVector3D& f(Z val, blas_int i)
        {
            body_[i] = val*sqrt(grid_.w(i));
            return *this;
        }
        Z f(blas_int i, blas_int j, blas_int k) const
        {
            return (*this)(i,j,k) / sqrt(grid_.w(i,j,k));
            return this->f(u);
        }
        GridVector3D& f(Z val, blas_int i, blas_int j, blas_int k)
        {
            (*this)(i,j,k) = val * sqrt(grid_.w(i,j,k));
            return *this;
        }
        const ARRAYS::Vector<Z>& body() const { return body_; }
        ARRAYS::Vector<Z>& body() { return body_; }
    };

    template<typename T, typename Z>
    class Operator3D
    {
     protected:
        bool init_;
        RowCompressedMatrix<Z> body_;
        FemDvrEcsGrid3D<T,Z> grid_;
     public:
        Operator3D() {}
        Operator3D(const FemDvrEcsGrid3D<T,Z>& grid) : grid_(grid) { body_ = RowCompressedMatrix<Z>(grid_.get_size(), grid_.get_size(), 0); }
        Operator3D& set_kinetic_term(T mu) { set_kinetic_term(mu, mu, mu); return *this; }
        Operator3D& set_kinetic_term(T muX, T muY, T muZ)
        {
            RowCompressedMatrix<Z> Tx = generateKineticTermRCM(grid_.gx(), muX);
            RowCompressedMatrix<Z> Ty = generateKineticTermRCM(grid_.gy(), muY);
            RowCompressedMatrix<Z> Tz = generateKineticTermRCM(grid_.gz(), muZ);

            RowCompressedMatrix<Z> Txy = TensorSum(Tx, Ty);
            body_ = TensorSum(Txy, Tz);

            return *this;
        }
        Z operator() (blas_int i, blas_int j) const { return body_.get_element(i,j); }
        Operator3D& conjugate() { body_.conjugate(); return *this; }
        Operator3D& operator*= (Z alpha) { body_ *= alpha; return *this; }
        Operator3D& operator+= (Z scalar) { body_.add_to_diagonal(scalar); return *this; }
        Operator3D& operator+= (const GridVector3D<T,Z>& rhs) { body_.add_grid_vector_to_diagonal(rhs); return *this; }
        void gemv(Z alpha, const GridVector3D<T,Z>& x, Z beta, GridVector3D<T,Z>& y) { body_.gemv(alpha, x.body(), beta, y.body()); }
        const RowCompressedMatrix<Z>& body() const { return body_; }
    };

}
