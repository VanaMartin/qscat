#ifndef INCLUDE_MASK_OPERATOR_2D_H_
#define INCLUDE_MASK_OPERATOR_2D_H_
namespace QSCAT
{
class MaskOperator2d
{
    MaskGrid2d grid_;                   // given grid
    RowCompressedMatrix<dcomp> body_;   // Operator body
    bool init_;
 public:
  // constructors
    MaskOperator2d();
    MaskOperator2d(const MaskGrid2d& grid);
  // accessors
    bool init() const;
    const MaskGrid2d& get_grid() const;
  // modifiers
    MaskOperator2d& set_kinetic_term(dfloat mu_x, dfloat mu_y);
    void LU_factorize();
  // operators
    MaskOperator2d& operator*= (const dcomp& alpha);
    MaskOperator2d& operator+= (const dcomp& alpha);
    MaskOperator2d& operator+= (const MaskVector2d& rhs);

    MaskVector2d operator* (const MaskVector2d& x) const;
    void gemv(const dcomp& alpha, const MaskVector2d& x, const dcomp& beta, MaskVector2d& y) const;
    void LU_back_substitution(MaskVector2d& x);
};
}   //  namespace QSCAT
#endif // INCLUDE_MASK_OPERATOR_2D_H_
