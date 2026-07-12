
// POLY GRID VECTORS -- REIMPLEMENTATION REQUIRED

/*

// Special class of  multiple state vectors associated  with one grid. The
// vectors are stored respectively  in one vector class, so it can be used
// in  the  Row  Compressed  Matrix  -  Vector  multiplication   and  back
// substitution.
template<typename T, typename Z>
class PolyGridVector
{
    bool init;                      // Initialization controller
    int base;                       // Number of basis elements in one grid vector
    int n;                          // Number of vector in the multivector
    int length;                     // The total length
public:
    Vector<Z> body;                 // The body of the multivector
    FemDvrEcsGrid * grid;   // Pointer to the appropriate fem-dvr-grid
private:
    void Initialize(FemDvrEcsGrid * g, const int & N);
    void Clean();
    bool save_bin_body(std::ofstream &file);
    bool read_bin_body(std::ifstream &file);
public:
    PolyGridVector();
    PolyGridVector(FemDvrEcsGrid * g, const int & N);
    PolyGridVector(PolyGridVector & old);
    PolyGridVector & Set(FemDvrEcsGrid * g, const int & N);
    void WriteVector(GridVector & V, const int &  i);
    GridVector GetVector(const int i);
    int GetVecLen();
    int GetPolyVecLen();
    int GetVecNum();
    Z Norm();
    Z & operator() (const int N, const int & i);
    const Z & operator() (const int N, const int & i) const;
    GridVector & EwSubMul(GridVector & rhs, const int & n);
    const Z F(const int N, const int & i);
    void F(const Z & val, const int N, const int & i);
    Z operator* (PolyGridVector & C);
    Z SubMult(GridVector& rhs, const int& n);
    bool save_binary(const char * name);
    bool save_binary(std::ofstream &file);
    bool read_binary(const char * name, FemDvrEcsGrid * g);
    bool read_binary(std::ifstream &file, FemDvrEcsGrid * g);
    void Save(const char * filename);
    void SaveTransposed(const Vector<Z>& X, const char* filename);
};

// Forward definitions
template<>
void PolyGridVector<def_float, def_comp>::Save(const char* filename);
template<>
void PolyGridVector<def_float, def_float>::Save(const char* filename);
//template<>
//void PolyGridVector<def_float, def_comp>::SaveTransposed(const Vector<def_float>& X, const char* filename);
template<>
void PolyGridVector<def_float, def_comp>::SaveTransposed(const Vector<def_comp>& X, const char* filename);
template<>
void PolyGridVector<def_float, def_float>::SaveTransposed(const Vector<def_float>& X, const char* filename);

// Special  class of  the  poly  matrix  Hamiltonian  stored  in  the  row
// compressed format.  The order of the matrix  respects the below defined
// order of the potnetial vectors.  The class uses standard RCM operations
// as defined in  the  arrays library.  The structure  of  the poly-vector
// potentials:  N-1 elments for the discretisation of the energy continuum
// range. The zeroth vector part is reserved for the discrete state.
// Vd:     Vd - V_0 + E_0  - ... - V_0 + E_(N-2)
// VdN:    0  -    Vd0     - ... -     Vd(N-2)
template<typename T, typename Z>
class PolyHamiltonian_RCM
{
    bool init;
    int m;
    int n;
    RowCompressedMatrix<Z> body;
private:
    void initialize(PolyGridVector & VdN, PolyGridVector & Vd, const int & N, const T & mu){
        //Kinetic_Energy_RCM<T,Z> KE(*(Vd.grid), mu);
        RowCompressedMatrix<Z> KE = generateKineticTermRCM<T,Z>(*Vd.grid, mu);
    // Auxiliary variables
        assert(N==VdN.GetVecNum());
        int nnz = KE.NNZ();
        int nb =  Vd.grid->NB();
        m = N*nb;
        n = N*nb;

        int index = 0;
        body.Set(m,n,N*nnz + 2*(N-1)*nb);
    // Writing the values to the body
        body.RI(0) = 0;
        for (int i=0; i<N; ++i){            // which row of the block matrix
            for (int k=0; k<nb; ++k){       // which row in the block
                // First makes the copy of the row start position
                body.RI(i*nb + k+1) = body.RI(i*nb + k);
                for (int j=0; j<N; ++j){    // which column of the block matrix
                    if (i==j){
                        // The diagonal block contains the Kinetic term
                        for (int l=KE.RI(k); l<KE.RI(k+1); ++l){
                            if (k == KE.C(l)){
                                // Pure Diagonal contains the potential term
                                body.NZE(index) = Vd.F(j,k) + KE.NZE(l);
                            } else {
                                // Nondiagonal kinetic terms
                                body.NZE(index) = KE.NZE(l);
                            }
                            body.C(index) = j*nb + KE.C(l); // j-th block + column in KE
                            body.RI(i*nb + k+1)++;              // increase of the position of next row first element
                            index++;                            // increase of the global index
                        }
                    } else if (i==0 || j==0){
                        // The non-diagonal term contains only the potential
                        if (i==0){
                            // The first row of blocks
                            body.NZE(index) = VdN.F(j,k);
                        } else {
                            // The first column of blocks
                            body.NZE(index) = conj(VdN.F(j,k));
                        }
                        body.C(index) = j*nb + k;               // j-th block + k-the column
                        body.RI(i*nb + k+1)++;                  // increase of the position of next row first element
                        index++;                                // increase of the global index
                    }
                }
            }
        }
        std::cout << "The NRM Hamiltonian Initiated with total of " << index << " non-zero elements." << std::endl;
        init = true;
    }
    void clean(){
        init = false;

    }
public:
    PolyHamiltonian_RCM(){
        init = false;
    }
    PolyHamiltonian_RCM(PolyGridVector & VdN, PolyGridVector & Vd, const int & N, const T & mu){
        initialize(VdN,Vd,N,mu);
    }
    PolyHamiltonian_RCM(PolyHamiltonian_RCM & old):body(old.body){
        init = old.init;
        m = old.m;
        n = old.n;
    }
    void Set(PolyGridVector & VdN, PolyGridVector & Vd, const int & N, const T & mu){
        initialize(VdN,Vd,N,mu);
    }
    PolyHamiltonian_RCM & operator= (PolyHamiltonian_RCM tmp){
        std::swap(init,tmp.init);
        std::swap(m,tmp.m);
        std::swap(n,tmp.n);
        body.Swap(tmp.body);
        return *this;
    }
    PolyHamiltonian_RCM & operator*= (const Z & alpha){
        body *= alpha;
        return *this;
    }
    PolyGridVector & operator*= (PolyGridVector & rhs){
        body *= rhs.body;
        return rhs;
    }
    PolyHamiltonian_RCM & AddToDiagonal(const Z & alpha){
        body.AddToDiagonal(alpha);
        return *this;
    }
    PolyHamiltonian_RCM & LUFactorize(){
        body.LUFactorize();
        return *this;
    }
    PolyGridVector & BackSubstitution(PolyGridVector & rhs){
        body.LUBackSubstitution(rhs.body);
        return rhs;
    }
};

template<typename T, typename Z>
class PolyCrank_Nicolson_RCM
{
    bool init;
    int order;
    T dt;
    T mu;
    Vector<Z> * roots;
    PolyHamiltonian_RCM<T,Z> * numerators;
    PolyHamiltonian_RCM<T,Z> * denominators;
private:
    void initialize(PolyHamiltonian_RCM<T,Z> & H, const T & MU, const T & DT, const int & O){
        mu = MU;

        //PolyHamiltonian_RCM<T,Z> Ham(H);

        if (O > 19) {
            order = 20;
        } else if (O > 14) {
            order = 15;
        } else if (O > 9) {
            order = 10;
        } else {
            order = O;
        }
        dt = DT;
        roots = new Vector<Z>(order);
        Pade_Roots(*roots,order);
        numerators = new PolyHamiltonian_RCM<T,Z>[order];
        denominators = new PolyHamiltonian_RCM<T,Z>[order];

        PolyHamiltonian_RCM<T,Z> * hlp;
        for (int i=0;i<order;++i) {
        // numerators
            hlp = &(numerators[i]);
            (*hlp) = H;
            (*hlp) *= imu*dt/(*roots)[i];
            (*hlp).AddToDiagonal(1.0);
            hlp = NULL;
        // denominators
            hlp = &(denominators[i]);
            (*hlp) = H;
            (*hlp) *= -imu*dt/conj((*roots)[i]);
            (*hlp).AddToDiagonal(1.0);
            (*hlp).LUFactorize();
            hlp = NULL;
        }
        init = true;
    }
public:
    PolyCrank_Nicolson_RCM(){
        init = false;
    }
    PolyCrank_Nicolson_RCM(PolyHamiltonian_RCM<T,Z> & H, const T & MU, const T & DT, const int & O){
        initialize(H,MU,DT,O);
    }
    PolyCrank_Nicolson_RCM & Set(PolyHamiltonian_RCM<T,Z> & H, const T & MU, const T & DT, const int & O){
        initialize(H,MU,DT,O);
        return *this;
    }
    PolyGridVector & One_Step(PolyGridVector & X){
        PolyHamiltonian_RCM<T,Z> * hlp;
        for (int s=0;s<order;++s){
            hlp = &(numerators[s]);
            (*hlp) *= X;
            hlp = &(denominators[s]);
            (*hlp).BackSubstitution(X);
            hlp = NULL;
        }
        return X;
    }
};

*/


using namespace QSCAT;
/*
Special class of  multiple state vectors associated  with one grid. The
vectors are stored respectively  in one vector class, so it can be used
in  the  Row  Compressed  Matrix  -  Vector  multiplication   and  back
substitution.
*/
template<typename T, typename Z>
void poly_grid_vector<T,Z>::Initialize(fem_dvr_ecs_grid<T,Z> * g, const int & N){
	grid = g;
	base = g->NB();
	n = N;
	length = N*base;
	body.Set(length);
	body.Fill(0.0);
	init = true;
}
template<typename T, typename Z>
void poly_grid_vector<T,Z>::Clean(){
	init = false;
}
template<typename T, typename Z>
bool poly_grid_vector<T,Z>::save_bin_body(std::ofstream &file){
	if (file.is_open()) {
		file.write((char*)&base,sizeof(int));
		file.write((char*)&n,sizeof(int));
		file.write((char*)&length,sizeof(int));
		body.SaveBinary(file);
		return true;
	} else {
		return false;
	}
}
template<typename T, typename Z>
bool poly_grid_vector<T,Z>::read_bin_body(std::ifstream &file){
	bool stat;
	if (file.is_open()){
		file.read((char*) &base, sizeof(int));
		file.read((char*) &n, sizeof(int));
		file.read((char*) &length, sizeof(int));
		stat = body.ReadBinary(file);
		init = stat;
		return stat;
	} else {
		return false;
	}
}
template<typename T, typename Z>
poly_grid_vector<T,Z>::poly_grid_vector(){
	init = false;
}
template<typename T, typename Z>
poly_grid_vector<T,Z>::poly_grid_vector(fem_dvr_ecs_grid<T,Z> * g, const int & N){
	init=false;
	Initialize(g,N);
}
template<typename T, typename Z>
poly_grid_vector<T,Z>::poly_grid_vector(poly_grid_vector<T,Z> & old):body(old.body){
	init = old.init;
	base = old.base;
	n = old.n;
	length = old.length;
	grid = old.grid;
}
template<typename T, typename Z>
poly_grid_vector<T,Z> & poly_grid_vector<T,Z>::Set(fem_dvr_ecs_grid<T,Z> * g, const int & N){
	if (init){
		Clean();
	}
	Initialize(g,N);
	return *this;
}
template<typename T, typename Z>
void poly_grid_vector<T,Z>::WriteVector(grid_vector<T,Z> & V, const int &  i){
	body.writeSubVector(V.a, i*base);
}
template<typename T, typename Z>
grid_vector<T,Z> poly_grid_vector<T,Z>::GetVector(const int i){
	grid_vector<T,Z> out(*grid);
	assert(i<n);
	blas::copy(base, &body[i*base], &out[0]);
	return out;
}
template<typename T, typename Z>
int poly_grid_vector<T,Z>::GetVecLen(){
	return base;
}
template<typename T, typename Z>
int poly_grid_vector<T,Z>::GetPolyVecLen(){
	return length;
}
template<typename T, typename Z>
int poly_grid_vector<T,Z>::GetVecNum(){
	return n;
}
template<typename T, typename Z>
Z poly_grid_vector<T,Z>::Norm(){
	return body*body;
}
template<typename T, typename Z>
Z & poly_grid_vector<T,Z>::operator() (const int N, const int & i) {		// Returns the address of the vector i-th element (weighted)
	assert(N*base + i < length);
	return body(N*base + i);
}
template<typename T, typename Z>
const Z & poly_grid_vector<T,Z>::operator() (const int N, const int & i) const{ // Returns the constant address of the vector i-th element (weighted)
	assert(N*base + i < length);
	return body(N*base + i);
}
template<typename T, typename Z>
grid_vector<T,Z> & poly_grid_vector<T,Z>::EwSubMul(grid_vector<T,Z> & rhs, const int & n){ // Returns the grid vector element wise multiplied by N-th internal vector
	body.EwSubMult(rhs.a, n);
	return rhs;
}
template<typename T, typename Z>
const Z poly_grid_vector<T,Z>::F(const int N, const int & i) {		// Returns the actual value of the function represented by the vector
	assert(N*base + i < length);
	return body[N*base + i]/sqrt((*grid).Wz(i));
}
template<typename T, typename Z>
void poly_grid_vector<T,Z>::F(const Z & val, const int N, const int & i) {	// Writes the actual value of the function to the vector with appropriate weight //grid_vector & F
	assert(N*base + i < length);
	body[N*base + i] = val*sqrt((*grid).Wz(i));
}
template<typename T, typename Z>
Z poly_grid_vector<T,Z>::operator* (poly_grid_vector<T,Z> & C) { // The operator denotes the scalar product of two vectors. The firs argument is taken as hermitean conjugated
	assert( C.getPolyVecLen() == length);
	return body*C.body;
}
template<typename T, typename Z>
Z poly_grid_vector<T,Z>::SubMult(grid_vector<T,Z>& rhs, const int& n){
	return body.SubMult(rhs.a, n);
}
template<typename T, typename Z>
bool poly_grid_vector<T,Z>::SaveBinary(const char * name){
	std::ofstream file;
	file.open(name, std::ios::out | std::ios::binary );
	if (save_bin_body(file)){
		file.close();
		return true;
	} else {
		return false;
	}
}
template<typename T, typename Z>
bool poly_grid_vector<T,Z>::SaveBinary(std::ofstream &file){
	return save_bin_body(file);
}
template<typename T, typename Z>
bool poly_grid_vector<T,Z>::ReadBinary(const char * name, fem_dvr_ecs_grid<T,Z> * g){
	grid = g;
	std::ifstream file;
	file.open(name, std::ios::in | std::ios::binary);
	if (read_bin_body(file)){
		file.close();
		return true;
	} else {
		return false;
	}
}
template<typename T, typename Z>
bool poly_grid_vector<T,Z>::ReadBinary(std::ifstream &file, fem_dvr_ecs_grid<T,Z> * g){
	grid = g;
	return read_bin_body(file);
}

template<typename T, typename Z>
void poly_grid_vector<T,Z>::Save(const char * filename){
std::cout << "The saving implementation for the poly_grid_vector class have not yet been implemented for given typename." << std::endl;
}
