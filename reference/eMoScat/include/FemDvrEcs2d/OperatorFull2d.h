#ifndef INCLUDE_OPERATOR_FULL_2D_H_
#define INCLUDE_OPERATOR_FULL_2D_H_
namespace QSCAT
{
/** \addtogroup FemDvrEcs2d
* @{ */
/** @} */

class Operator2dFull : public Object
{
    FemDvrEcsGrid2d grid_;
    zMatrix body_;
 public:
  // constructors
    Operator2dFull();
    Operator2dFull(FemDvrEcsGrid2d& grid);
    Operator2dFull(const Operator2dFull& old);
  // accessors
    const FemDvrEcsGrid2d& get_grid() const;
  // modifiers
    Operator2dFull& set_kinetic_term(const dfloat& mux, const dfloat& muy);
    void LU_factorize();
  // operators
    Operator2dFull& operator*= (const dcomp& alpha);
    Operator2dFull& operator+= (const dcomp& alpha);
    Operator2dFull& operator+= (const GridVector2d& rhs);
    Operator2dFull& operator+= (const zVector& rhs);
    GridVector2d operator* (const GridVector2d& x) const;
    void gemv(const dcomp& alpha, const GridVector2d& x, const dcomp& beta, GridVector2d& y) const;
    void LU_back_substitution(GridVector2d& x);
};
}   //  namespace QSCAT
#endif // INCLUDE_OPERATOR_FULL_2D_H_
