#ifndef INCLUDE_GRID_VECTOR_2D_H_
#define INCLUDE_GRID_VECTOR_2D_H_

#include "common.h"
#include "Storage.h"
#include "Arrays/Vector.h"
#include "Arrays/Matrix.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"
#include "FemDvrEcs2d/FemDvrEcsGrid2d.h"

namespace QSCAT
{
/** \addtogroup FemDvrEcs2d
* @{ */

//! The two  dimensional function on the  grid represented as a grid vector type.
/*!
    The ordering  of  the values in  the  internal  arrays  goes as:
    f(x[i],y[j]) = a[nb_x*j + i]
*/
class GridVector2d : public Object, public BinaryStorageInterface
{
 protected:
    FemDvrEcsGrid2d grid_;      //!< 2D coordinate discretization
    zVector body_;              //!< function reprezentation in given basis
 private:
    //! Internal cleanup helper
    void clean();

 protected:
    //! Internal save to binary stream helper
    virtual bool save_bin_body(std::ofstream& file) const;
    //! Internal read from binary stream helper
    virtual bool read_bin_body(std::ifstream& file);

 public:
 // constructors

    //! Default constructor
    GridVector2d();
    //! Simple constructor
    GridVector2d(const FemDvrEcsGrid2d& g);
    //! Constructor from flattened array
    GridVector2d(const FemDvrEcsGrid2d& g, const zVector& unweighted_values);
    //! Constructor as outer product
    GridVector2d(const FemDvrEcsGrid2d& g, const GridVector& x_values, const GridVector& y_values);
    //! Copy constructor
    GridVector2d(const GridVector2d& old);
    //! Destructor
    ~GridVector2d();

 // accessors

    //! x-coordinate real part discretization size
    blas_int get_real_xsize() const;
    //! y-coordinate real part discretization size
    blas_int get_real_ysize() const;
    //! x-coordinate discretization size
    blas_int get_xsize() const;
    //! y-coordinate discretization size
    blas_int get_ysize() const;
    //! 2D-coordinates  discretization size
    blas_int get_size() const;
    //! Whole 2D coordinate system discretization
    const FemDvrEcsGrid2d& get_grid() const;
    //! Representation value from flattened array
    const dcomp& operator[] (blas_int i) const;
    //! Representation value from flattened array
    dcomp& operator[] (blas_int i);
    //! Function value from flattened array
    dcomp f(blas_int i) const;
    //! Function value from flattened array
    void f(const dcomp& val,  blas_int i);
    //! Function value in double index representaion
    dcomp f(blas_int iy, blas_int ix) const;
    //! Function value in double index representaion
    void f(const dcomp& val, blas_int iy, blas_int ix);
    //! All function values arranged in a vector
    zVector function_values() const;
    //! Function inner product associated norm
    dfloat norm() const;
    //! Function value at arbitrary point within discretization range
    dcomp evaluate(const dfloat& x, const dfloat& y) const;
    //! Internal representation
    const zVector& body() const;
    //! Internal representation
    zVector& body();

 // modifiers

    //! Overwrite values along x-axis at some y-axis discretization point
    void write_x_section(GridVector& section, blas_int j);
    //! Retrieve one-dimensional function along x-axis at some y-axis discretization point
    GridVector get_x_section(blas_int j) const;
    //! Swap operation
    GridVector2d& swap(GridVector2d& rhs);
    //! Explicite copy call
    GridVector2d copy() const;
    //! Explicite complex conjugation operation
    GridVector2d& complex_conjugate();
    //! Set whole function as constant value (FIXME : redo as nullify)
    GridVector2d& fill(const dcomp& val);
    //! Fill with given function (TODO: redo via picojson)
    GridVector2d& fill_function(const dcomp (*func)(const dcomp&, const dcomp&,const parameters::model_2D<dfloat>&), const parameters::model_2D<dfloat>& mp);

 // operators

    //! Scaling by a scalar constant
    GridVector2d& operator*= (const dcomp& alpha);
    //! Inner product
    dcomp operator* (const GridVector2d& rhs) const;
    //! Inplace addition
    GridVector2d& operator += (const GridVector2d& rhs);
    //! Inplace subtraction
    GridVector2d& operator -= (const GridVector2d& rhs);
    //! Assignement operator
    GridVector2d& operator = (GridVector2d tmp);
    //! Deffered scaling for next addition or assignement operation
    ScalarMultiple<dcomp, Vector<dcomp> > operator* (const dcomp& alpha);
    //! Assignement of defferred scaled function
    GridVector2d& operator = (ConstScalarMultiple<dcomp, Vector<dcomp> > mlt);
    //! Addition of deffered scaled function
    GridVector2d& operator += (ConstScalarMultiple<dcomp, Vector<dcomp> > mlt);

 // custom functions

    //! Scale and inplace add result
    GridVector2d & axpy(const dcomp& alpha, const GridVector2d& rhs);
    //! Scale and assign result
    GridVector2d & ax(const dcomp& alpha, const GridVector2d& rhs);
    //! Elementwise multiplication sum (inner product without conjugation)
    dcomp reduction(const GridVector2d& y) const;
    //! Function of one dimensional inner products with function
    dcomp line_projection(const GridVector& P, const char& axis, blas_int pos) const;
    //! Partial inner product (along one dimension)
    GridVector contraction(const GridVector2d& P, const char& axis);
    //! Element wise multiplication
    GridVector2d& element_wise_multiplication(const GridVector2d& src);

 // storage

    //! Save into plain txt file
    void save(const char *filename) const;
    //! Save the function evaluated on equidistant grid into plain txt file
    void save_equidistant(const char *name, blas_int x_points, dfloat x_min, dfloat x_max, blas_int y_points, dfloat y_min, dfloat y_max) const;
};

/** @} */
}   //  namespace QSCAT
#endif // INCLUDE_GRID_VECTOR_2D_H_
