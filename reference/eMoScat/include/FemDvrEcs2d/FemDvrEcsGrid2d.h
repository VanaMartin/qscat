#ifndef INCLUDE_FEM_DVR_ECS_GRID_2D_H_
#define INCLUDE_FEM_DVR_ECS_GRID_2D_H_

#include "common.h"
#include "Storage.h"
#include "Arrays/Vector.h"
#include "Arrays/Matrix.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"

namespace QSCAT
{
/** \addtogroup FemDvrEcs2d
* @{ */

//! Two dimensional grid type.
/** The ordering of the grid is given by repeating the X coordinate (as a
 *  column of a matrix) for each point of the Y coordinate. Stored in a matrix,
 *  the elements a(i,j) = a_x(i)*a_y(j)
 */
class FemDvrEcsGrid2d : public Object, public BinaryStorageInterface
{
    //! X coordinate discretisation.
    FemDvrEcsGrid xgrid_;
    //! Y coordinate discretisation.
    FemDvrEcsGrid ygrid_;
    //! X coordinate basis size.
    blas_int xbasis_size_;
    //! Y coordinate basis size.
    blas_int ybasis_size_;
    //! Total of basis functions.
    /** Equals to x_basis_size_ * y_basis_size_ */
    blas_int basis_size_;
 private:
    //! Internal initialization helper.
    void initialize(const FemDvrEcsGrid& gx, const FemDvrEcsGrid& gy);

    //! Internal cleanup helper.
    void clean();

 protected:
    //! Internal save to binary stream helper
    virtual bool save_bin_body(std::ofstream& file) const;
    //! Internal read from binary stream helper
    virtual bool read_bin_body(std::ifstream& file);
 public:

 //constructors

    //! Default constructor
    FemDvrEcsGrid2d();
    //! Constructor
    FemDvrEcsGrid2d(const FemDvrEcsGrid& gx, const FemDvrEcsGrid& gy);
    //! Copy constructor
    FemDvrEcsGrid2d(const FemDvrEcsGrid2d& old);
    //! Destructor
    ~FemDvrEcsGrid2d();

 // accessors

    //! x-coordinate discretization size
    blas_int get_xsize() const;
    //! y-coordinate discretization size
    blas_int get_ysize() const;
    //! 2D-coordinates  discretization size
    blas_int get_size() const;
    //! x-coordinate real value at discretization point
    const dfloat& xr(blas_int i) const;
    //! y-coordinate real value at discretization point
    const dfloat& yr(blas_int i) const;
    //! x-coordinate complex value at discretization point
    const dcomp& xz(blas_int i) const;
    //! y-coordinate complex value at discretization point
    const dcomp& yz(blas_int i) const;
    //! whole x-coordinate discretization grid
    const FemDvrEcsGrid &get_xgrid() const;
    //! whole y-coordinate discretization grid
    const FemDvrEcsGrid &get_ygrid() const;
    //! 2D-discretization weight factor
    dcomp wz(blas_int i, blas_int j) const;
    //! 2D-discretization weight factor (flattened)
    dcomp wz(blas_int i) const;
    //! x-coordinate weight factor at given point
    const dcomp& xwz(blas_int i) const;
    //! y-coordinate weight factor at given point
    const dcomp& ywz(blas_int i) const;
    //! index of x-coordinate element upper bound
    blas_int x_element_end(const dfloat& x) const;
    //! index of y-coordinate element upper bound
    blas_int y_element_end(const dfloat& x) const;
    //! x-coordinate real part discretization size
    blas_int get_real_xsize() const;
    //! y-coordinate real part discretization size
    blas_int get_real_ysize() const;

 // modifiers

    //! Swap operation
    FemDvrEcsGrid2d& swap(FemDvrEcsGrid2d& rhs);

    //! Assignement operator
    FemDvrEcsGrid2d& operator= (FemDvrEcsGrid2d tmp);

    //! Equivalence evaluation
    bool operator== (const FemDvrEcsGrid2d& rhs) const;

    //! Inequivalence evaluation
    bool operator!= (const FemDvrEcsGrid2d& rhs) const;
};

/** @} */
}   //  namespace QSCAT
#endif // INCLUDE_FEM_DVR_ECS_GRID_2D_H_
