#ifndef INCLUDE_OPERATOR_ROW_COMPRESSED_2D_H_
#define INCLUDE_OPERATOR_ROW_COMPRESSED_2D_H_
namespace QSCAT
{
/** \addtogroup FemDvrEcs2d
* @{ */

//! General operator class on 2d grid in csr format
class Operator2dRowCompressed : public Object
{
    FemDvrEcsGrid2d grid_;              //!< associated discretization
    RowCompressedMatrix<dcomp> body_;   //!< internal discretized operator representation
 public:
 // constructors

    //! Default constructor
    Operator2dRowCompressed();
    //! Constructor
    Operator2dRowCompressed(FemDvrEcsGrid2d& grid);
    //! Copy constructor
    Operator2dRowCompressed(const Operator2dRowCompressed& old);
    //! Destructor
    ~Operator2dRowCompressed();

    //! Deep copy operation
    Operator2dRowCompressed copy() const;
    //! Swap operation
    Operator2dRowCompressed& swap(Operator2dRowCompressed& rhs);
    //! Assignement operation
    Operator2dRowCompressed& operator= (Operator2dRowCompressed tmp);

 // accessors

    //! Associated 2d discretization
    const FemDvrEcsGrid2d& get_grid() const;
    //! Matrix element in discretized representaion
    dcomp operator() (int i, int j) const;
    //! Whole internal matrix elements sparse representaion
    const RowCompressedMatrix<dcomp>& body() const;

 // modifiers

    //! Set operator to kinetic term
    Operator2dRowCompressed& set_kinetic_term(const dfloat& mux, const dfloat& muy);
    //! Lower Upper factorization via LAPACK
    void LU_factorize();
    //! Store conjugation operation for use in fast BLAS/LAPACK operations
    Operator2dRowCompressed& conjugate();

 // operators

    //! Scaling with scalar constant
    Operator2dRowCompressed& operator*= (const dcomp& alpha);
    //! Inplace addition of constant to operator diagonal
    Operator2dRowCompressed& operator+= (const dcomp& alpha);
    //! Inplace addition of function to operator diagonal
    Operator2dRowCompressed& operator+= (const GridVector2d& rhs);
    //! Inplace addition of values to operator diagonal
    Operator2dRowCompressed& operator+= (const zVector& rhs);
    //! Acting of operator on given state
    GridVector2d operator* (const GridVector2d& x) const;
    //! Scaled operator action on given vector then addition to another scaled function
    void gemv(const dcomp& alpha, const GridVector2d& x, const dcomp& beta, GridVector2d& y) const;
    //! Solution to linear problem
    void LU_back_substitution(GridVector2d& x);
};

/** @} */
}   //  namespace QSCAT
#endif // INCLUDE_OPERATOR_ROW_COMPRESSED_2D_H_
