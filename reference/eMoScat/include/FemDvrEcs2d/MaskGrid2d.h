#ifndef INCLUDE_MASK_GRID_2D_H_
#define INCLUDE_MASK_GRID_2D_H_
namespace QSCAT
{
class MaskGrid2d
{
    FemDvrEcsGrid2d grid_;
    iMatrix mask_;

    blas_int x_size_;
    blas_int y_size_;
    blas_int full_size_;

    iVector row_index_;             // How many basis functions in each row (as in row compressed matrix)
    iVector columns_;               // which column element

    blas_int x_elements_, y_elements_;       // size controllers
    blas_int x_quadrature_, y_quadrature_;   // element basis size for coordinates
    blas_int total_elements_;                // total of nonzero elements in compressed representation
    blas_int compressed_size_;               // total of nonzero basis functions

    bool init_;
 public:
  // constructors
    MaskGrid2d();
    MaskGrid2d(const FemDvrEcsGrid2d& grid, const iMatrix& mask);
    MaskGrid2d(const MaskGrid2d& old);
  // modifiers
    MaskGrid2d& swap(MaskGrid2d& rhs);
    MaskGrid2d& operator= (MaskGrid2d tmp);
  // accessors
    bool init() const;
    blas_int get_size() const;
    blas_int get_xsize() const;
    blas_int get_ysize() const;
    dcomp wz(blas_int i) const;
    const FemDvrEcsGrid2d& full_grid() const;
    const iVector& get_columns() const;
    const iVector& get_row_index() const;
    const iMatrix& get_mask() const;
  // operators
    bool operator== (const MaskGrid2d& rhs) const;     // FIXME (use hash)
};
}   //  namespace QSCAT
#endif // INCLUDE_MASK_GRID_2D_H_
