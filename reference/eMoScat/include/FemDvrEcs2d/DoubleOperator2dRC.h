#ifndef INCLUDE_DOUBLE_OPERATOR_2D_H_
#define INCLUDE_DOUBLE_OPERATOR_2D_H_
namespace QSCAT
{
/// Concept of operator on two coupled systems
class DoubleOperator2dRC
{
    FemDvrEcsGrid2d grid1_;
    FemDvrEcsGrid2d grid2_;
    RowCompressedMatrix<dcomp> body_;
    bool init_;
 public:
  // Constructors
    DoubleOperator2dRC();
    DoubleOperator2dRC(FemDvrEcsGrid2d& grid);
    DoubleOperator2dRC(FemDvrEcsGrid2d& grid1, FemDvrEcsGrid2d& grid2);
    DoubleOperator2dRC(const DoubleOperator2dRC& old);

  // accessors
    bool init() const;
    const FemDvrEcsGrid2d& get_grid1() const;
    const FemDvrEcsGrid2d& get_grid2() const;
    blas_int get_size1() const;
    blas_int get_size2() const;
    blas_int get_size() const;

  // modifiers
    DoubleOperator2dRC& add_kinetic_term(const dfloat& mu1, const dfloat& mu2, const dfloat& mu3, const dfloat& mu4);
    DoubleOperator2dRC& add_kinetic_term(const dfloat& mu1, const dfloat& mu2);
    DoubleOperator2dRC& swap(DoubleOperator2dRC& rhs);
    void LU_factorize();

  // operators
    DoubleOperator2dRC& operator*= (const dcomp& alpha);
    DoubleOperator2dRC& operator+= (const dcomp& alpha);
    DoubleOperator2dRC& operator= (DoubleOperator2dRC rhs);

  // custom operations
    void gemv(const dcomp& alpha, const DoubleGridVector2d& x, const dcomp& beta, DoubleGridVector2d& y) const;
    void LU_back_substitution(DoubleGridVector2d&x);

    DoubleOperator2dRC& add_potential(const GridVector2d& p1, const GridVector2d& p2);
    DoubleOperator2dRC& add_coupling(const GridVector2d& c1, const GridVector2d& c2);
};
}   //  namespace QSCAT
#endif // INCLUDE_DOUBLE_OPERATOR_2D_H_
