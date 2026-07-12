#ifndef INCLUDE_DOUBLE_GRID_VECTOR_2D_H_
#define INCLUDE_DOUBLE_GRID_VECTOR_2D_H_
namespace QSCAT
{
/// Joined representation of wavefunctions in two (possibly coupled) model spaces
class DoubleGridVector2d
{
    FemDvrEcsGrid2d grid1_;
    FemDvrEcsGrid2d grid2_;
    zVector body_;
    bool init_;
 public:
  // constructors
    DoubleGridVector2d();
    DoubleGridVector2d(FemDvrEcsGrid2d& G);
    DoubleGridVector2d(FemDvrEcsGrid2d& G1, FemDvrEcsGrid2d& G2);
    DoubleGridVector2d(FemDvrEcsGrid2d& G, GridVector2d& psi, blas_int i);
    DoubleGridVector2d(FemDvrEcsGrid2d& G1, FemDvrEcsGrid2d& G2, GridVector2d& psi1, GridVector2d& psi2);
    DoubleGridVector2d(const DoubleGridVector2d& old);

  // access
    bool init() const;
    dcomp& operator[] (blas_int i);
    const dcomp& operator[] (blas_int i) const;
    blas_int get_size1() const;
    blas_int get_size2() const;
    blas_int get_size() const;
    void f(const dcomp& val,const blas_int i);
    dcomp f(const blas_int i) const;
    const zVector& body() const;
    zVector& body();
    const FemDvrEcsGrid2d& get_grid1() const;
    const FemDvrEcsGrid2d& get_grid2() const;

  // modifiers
    DoubleGridVector2d& swap(DoubleGridVector2d &rhs);

  // operators
    dcomp operator* (const DoubleGridVector2d& rhs) const;
    DoubleGridVector2d& operator= (DoubleGridVector2d tmp);
    DoubleGridVector2d& operator+= (const DoubleGridVector2d &rhs);
    DoubleGridVector2d& operator-= (const DoubleGridVector2d &rhs);
    DoubleGridVector2d& operator*= (const dcomp& alpha);

  // custom operators
    DoubleGridVector2d& axpy(const dcomp& alpha, const DoubleGridVector2d &x);
    DoubleGridVector2d& ax(const dcomp& alpha, const DoubleGridVector2d &x);
};
}   //  namespace QSCAT
#endif // INCLUDE_DOUBLE_GRID_VECTOR_2D_H_
