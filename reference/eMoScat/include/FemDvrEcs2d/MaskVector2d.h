#ifndef INCLUDE_MASK_VECTOR_2D_H_
#define INCLUDE_MASK_VECTOR_2D_H_
namespace QSCAT
{
class MaskVector2d
{
    MaskGrid2d grid_;              // Given Grid containing compression information
    zVector body_;                    // Compressed body
  //
    bool init_;                         // initialization controller
 public:
  // constructors
    MaskVector2d();
    MaskVector2d(const MaskGrid2d& grid);
    MaskVector2d(const MaskGrid2d& grid, const GridVector& x, const GridVector& y);
    MaskVector2d(const MaskVector2d& old);
  // accessors
    bool init() const;
    blas_int get_size() const;
    const MaskGrid2d& get_grid() const;
    const zVector& body() const;
    zVector& body();
    dcomp f(blas_int i) const;
    void f(const dcomp& val,  blas_int i);
    zVector function_values() const;
    //const Z f(blas_int i, blas_int j) const;
    //void f(const Z& val, blas_int i, blas_int j);
  // modifiers
    GridVector2d& get(GridVector2d& dst) const;
    MaskVector2d& set(const GridVector2d& src);
    MaskVector2d& swap(MaskVector2d& rhs);
    MaskVector2d& copy(const MaskVector2d& src);
  // operators
    dcomp operator* (const MaskVector2d& rhs) const;
  // storage
    void save(const char* filename) const;
};
}   //  namespace QSCAT
#endif // INCLUDE_MASK_VECTOR_2D_H_
