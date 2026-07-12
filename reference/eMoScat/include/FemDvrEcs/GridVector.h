#ifndef INCLUDE_GRID_VECTOR_H_
#define INCLUDE_GRID_VECTOR_H_

#include "common.h"
#include "Arrays/Vector.h"

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */
//! State vector class.
/** A classical vector associated to given grid. The values represent
 * discretization of function as a vector of it values multiplied by grid
 * weights: \f[ f_i = f(x_i) \sqrt{w_i}. \f] This representation is equivalent
 * to decomposition to the FEM-DVR-ECS basis, since the basis functions are
 * orthogonal. The \f$ L^2 \f$ integration of the function is then approximated
 * as a sum of vector elements, \f[ \int_{x_{min}}^{x_{max}} |f(x)|^2 dx \cong
 * \sum_{i=0}^n |f(x_i)|^2 w_i = \sum_{i=0}^n (f(x_i) \sqrt{w_i})^* (f(x_i)
 * \sqrt{w_i} ) = \sum_{i=0}^n f_i^* f_i = \vec{f} \cdot \vec{f}, \f] where the
 * \f$ \cdot \f$ denotes vector inner product. The inner products of two
 * GridVectors then represent integrals of the function overlap.
*/
class GridVector : public Object, public BinaryStorageInterface
{
 protected:
    FemDvrEcsGrid grid_;    //!< Associated grid
    zVector vector_;        //!< Vector of values
 protected:
    //! Internal binary saving procedure.
    /** Saves the class structure into binary stream.
     *  @return true on success false on error.
     */
    virtual bool save_bin_body(std::ofstream& file) const;

    //! Internal binary reading procedure.
    /** Reads the class structure from binary stream.
     *  @return true on success false on error.
     */
    virtual bool read_bin_body(std::ifstream& file);
 public:

 // constructors

    //! Default constructor.
    /** Creates uninitialized instance. */
    GridVector();

    //! Constructor, initializes the internal vector.
    /** Determines the size of the vector from provided grid.
     *  @param  grid    Source discretization grid.
    */
    GridVector(const FemDvrEcsGrid& grid);

    //! Basic constructor.
    /** Intializes the internal vector and fills it from given values.  If not
     *  specified differently the values are weighted by square rooot of given
     *  grid weights.
     *  @param  grid    Source discretization grid.
     *  @param  values  Values to be inserted into reprezentation. The values
     *                  vector must be of the same size as the grid number of
     *                  basis functions.
     *  @param  set_weights     If set true, weight factors will be applied.
     */
    GridVector( const FemDvrEcsGrid& grid,
                const zVector& values,
                bool set_weights=true );

    //! Destructor.
    /** Decrease the refence count. */
    ~GridVector();

    //! Copy constructor (shallow copy).
    /** Performs a shallow copy operation on all internal members.
     *  @param  old Instance of the GridVector to be copy from.
     */
    GridVector(const GridVector& old);

 // accessors

    //! Representation size.
    /** @return Size of the vector, i.e. the size of coordinate discretization
     *  basis.
     */
    blas_int get_size() const;

    //! Function representation value access.
    /** Weighted value of the function representation on the discretized
     *  coordinate \f[ \sqrt{w_i} f(x_i) \f], therefore the retrieved reference
     *  expects the weight factor \f$ \sqrt{w_i} \f$ to be included. If the
     *  actual function value is required use the GridVector::f method instead.
     *  @param  i   Index of the discretization point.
     *  @return Refernce onto representation value.
     */
    dcomp& operator[] (blas_int i);

    //! Function representation value access.
    /** Weighted value of the function representation on the discretized
     *  coordinate \f[ \sqrt{w_i} f(x_i) \f], the retrieved value reference
     *  contians the weight factor \f$ \sqrt{w_i} \f$. If the actual function
     *  value is required use the GridVector::f method instead.
     *  @param  i   Index of the discretization point.
     *  @return Constant refernce onto the function value representation.
     */
    const dcomp& operator[] (blas_int i) const;

    //! Actual function value at desired discretization point.
    /** Actual function value \f[ f(x) \f] retrieved from representation value
     *  by division with weight factor \f$ \sqrt{w_i} \f$ before returning.
     *  @param  i   Index of the discretization point.
     *  @return Value of the function at given point \f$ f(x_i) \f$.
     */
    dcomp f(blas_int i) const;

    //! Writes actual function value at desired discretization point.
    /** Actual function value \f[ f(x) \f] is internally multiplied by weight
     *  factor \f$ \sqrt{w_i} \f$ before storing.
     *  @param  val     Actual value of the given function (state) at at given
     *                  discretization point.
     *  @param  i       Index of the discretization point.
     */
    void f(const dcomp& val, blas_int i);

    //! Vector of actual function values.
    /** Internally computes actual function values \f$f(x_i)\f$ of the state by
     *  GridVector::f() method to retrieve all values.
     *  @return Vector with all actual (unweighted) values.
     */
    zVector function_values() const;

    //! Vector of function representation values.
    /** Provides the internal ordered vector of function representation values
     *  \f[ f_i = \sqrt{w_i} f(x_i) \f].
     *  @return Reference onto complex float vector instance with all
     *          representation values.
     */
    zVector& body();

    //! Vector of function representation values.
    /** Provides the internal ordered vector of function representation values
     *  \f[ f_i = \sqrt{w_i} f(x_i) \f].
     *  @return Constant reference onto complex float vector instance with all
     *          representation values.
     */
    const zVector& body() const;

    //! Coordinate discretization grid.
    /** @return Constant reference on the associated FemDvrEcsGrid containing
     *  all information about the coordinate discretization.
     */
    const FemDvrEcsGrid& get_grid() const;

    //! Weight factor.
    /** @param  i   Index of the discretization point.
     *  @return Weight factor associated with the given discretization point.
     */
    const dcomp& w(blas_int i) const;

    //! Coordinate real value.
    /** Real part of discretized coordinate value at given discretization
     *  point.
     *  @param  i   Index of the discretization point.
     *  @return Coordiante value real part.
     */
    const dfloat& xr(blas_int i) const;

    //! Coordinate complex value.
    /** Complex value of discretized coordinate at given discretization point.
     *  @param  i   Index of the discretization point.
     *  @return Complex coordiante value.
     */
    const dcomp& xz(blas_int i) const;

 // modifiers

    // TODO REMOVE
    //! Fill vector with constant value.
    /** DEPRECATED: this method fills all of the function representation values
     *  whith the same constant value, which does not make a lot of sense unless
     *  the value is zero. To be replaced by nullify in the future.
     *  @param  val     Function representation value to be inserted to all
     *                  points.
     *  @return Reference on this instantion
     */
    GridVector& fill(const dcomp& val);

    //TODO redesign method for more general input
    //! Fill with function of model_2D parametrization.
    /** This method evaluates the given function at all discretization points
     *  and stores them multiplied by the weight factor.
     *  @return Reference on this updated instantion.
     */
    //GridVector& function_fill( const dcomp (*func)(const dcomp&, const parameters::model_2D<dfloat> &),
    //                           const parameters::model_2D<dfloat> & mp);

    // TODO redesign method for more general input
    //! Fill with two dimensional function of model_2D parametrization section along x-axis.
    /** This method evaluates the given 2D function at all discretization
     *  points along x-axis (for fixed y) and stores them multiplied by the
     *  weight factor.
     *  @return Reference on this instantion.
     */
    //GridVector& function_xy_fill_x( const dcomp (*func) (const dcomp&, const dcomp&, const parameters::model_2D<dfloat>&),
    //                                const dcomp& y, const parameters::model_2D<dfloat>& mp);

    // TODO redesign method for more general input
    //! Fill with two dimensional function of model_2D parametrization section along y-axis.
    /** This method evaluates the given 2D function at all discretization
     *  points along y-axis (for fixed x) and stores them multiplied by the
     *  weight factor.
     *  @return Reference on this instantion.
     */
    //GridVector& function_xy_fill_y( const dcomp (*func)( const dcomp&, const dcomp&, const parameters::model_2D<dfloat>&),
    //                                const dcomp& x, const parameters::model_2D<dfloat>& mp);
    // TODO DEPRECATED : To be removed from the GridVector namespace
    //! Set as radial part of free spherical problem solution.
    /** Fills the vector with energy normalized spherical Bessel function as a
     * radial part of the 3-dimensional free spherically symmetrical problem
     * \f[ \left( - \frac{1}{2 m r^2} \frac{\partial}{\partial r} \left(r^2
     * \frac{\partial }{\partial r} \right) + \frac{l (l+1)}{2r^2} \right) R(r)
     * = E R(r)  \f].
     *  @param  energy  Complex energy of desired solution.
     *  @param  mass    Effective mass term.
     *  @param  impulse_momentum    Impulse momentum quantum number.
     *  @return Reference on this instantion.
     */
    GridVector& radial_state( const dcomp& energy,
                              const dfloat& mass,
                              int impulse_momentum);

    // TODO DEPRECATED : To be removed from the GridVector namespace
    //! Fills the vector with energy normalized solution of spherical coulomb problem.
    /** Fills the vector with energy normalized spherical Coulomb function as a
     *  radial part of the 3-dimensional spherically symmetrical Coulomb problem
     *  \f[ \left( - \frac{1}{2 m r^2} \frac{\partial}{\partial r} \left(r^2
     *  \frac{\partial }{\partial r} \right) + \frac{l (l+1)}{2r^2} -
     *  \frac{q}{r} \right) R(r) = E R(r)  \f].
     *  @param      energy .. complex float, complex energy of desired solution
     *  @param      mass .. float, mass term
     *  @param      impulse_momentum .. ineger, impulse momentum quantum number
     *  @param      charge .. complex float, reduced charge parameter
     *  @return reference on this instantion
     */
    GridVector& coulomb_state( const dcomp& energy,
                               const dfloat& mass,
                               int momentum,
                               const dcomp& charge );

    //! \f$L_2\f$ norm of stored function.
    /** Computes \f$L_2\f$ norm  given by \f[ \sqrt{ \int_{x_0}^{x_n} |f(x)|^2
     *  dx } \cong \sqrt{ \sum_{i=0}^{n} |f(x_i)|^2 w_i } = \sqrt{
     *  \sum_{i=0}^{n} |f_i|^2 } \f]
     *  @return Computed value of the norm.
     */
    dcomp get_norm() const;

    //! Function derivative at given discretization point.
    /** Evaluates the stored function prime derivative \f[ f'(x_i) =
     *  (\frac{\partial f}{\partial x})(x_i) = \sum_{k=0}^{q} D_{k(i-k_0)}
     *  f(x_{k_0 + k}), \f] where \f$ q\f$ denotes the quadrature size, matrix
     *  \f$ D \f$ denotes basis functions derivative matrix (see
     *  FemDvrEcsGrid::dlp() method for more details), and index \f$ k_0 \f$
     *  denotes the index of the first (connecting) point of element to which
     *  discretization point \f$ x_i \f$ belongs.
     *  @param  pos     Index of the discretization point.
     *  @return Value of the evaluated derivative.
     */
    dcomp derivative(blas_int pos) const;

    //! Funciton value at arbitrary point of the coordinate.
    /** Evaluates the actual function value at given coordinate point (assumed
     *  within coordinate range) as \f[ f(x) = \sum_{k=0}^{q} f(x_{k_0 + k})
     *  \phi_k(x), \f] where \f$ q\f$ denotes the quadrature, index \f$ k_0 \f$
     *  denotes the index of the first (connecting) point of element to which
     *  point \f$ x \f$ belongs and \f$ \phi_k(x) \f$ denotes the function
     *  value of k-th basis function on assumed element.
     *  @param  x   Coordinate value to evalueate function at.
     *  @return Value of evaluated function.
     */
    dcomp evaluate(const dfloat& x) const;

    //! Swap operation between two instantions.
    /** Swaps all internal variables between two instantions.
     *  @param  rsh     Reference on the instantion to swap with.
     *  @return Reference on this updated instantion.
     */
    GridVector& swap(GridVector& rhs);

    //! Deep copy operation.
    /** Explicite call for a deep copy of this instantion. Preforms deep copy
     *  of all internal members.
     *  @return New instantion holding same values and grid original.
     */
    GridVector copy() const;

    //! Complex conjugation.
    /** Performs actual complex conjugation on all elements (elementwise).
     *  @return reference on this updated instantion.
     */
    GridVector& complex_conjugate();

 // operators ==========

    //! Inner product of instantion and right hand side instantion.
    /** Both states must be defined on the same grid. The result is an
     *  approximation of of the inner product of associated Hilbert space. The
     *  approxiamtion is given by the overlap integration \f[ \int_{x_0}^{x_n}
     *  f^*(x) g(x) dx \cong \sum_{i=0}^n w_i f^*(x_i) g(x_i) = \sum_{i=0}^n
     *  f_i g_i = \vec{f} \cdot \vec{g} \f]
     *  @param  C   Right hand side state vector to compute the product with.
     */
    dcomp operator* (const GridVector& C) const;

    //! Inplace scalar multiplication of the state.
    /** Prepares the auxiliary Scalar Multiple class for further use in
     *  computations.
     *  @param  c    Scaling factor.
     *  @return Scalar Multiple object instantion.
     */
    ScalarMultiple<dcomp, Vector<dcomp> > operator* (const dcomp& c);

    //! Assigns the state from right hand size.
    /** Performs shalloe copy & swap internally
     *  @param  tmp  Temporary shallow copy.
     *  @return reference on this updated instantion.
     */
    GridVector& operator=  (GridVector tmp);

    //! Assigns scaled copy stored in right hand side.
    /** Performs fast BLAS axpy operation internally.
     *  @note This method is designed to be called as a result of operator *
     *  between GridVector and scalar.
     *  @param  rhs  Result of state vector * scalar operation to be assigned to
     *              internal representation.
     *  @return Reference onto this updated instantion.
     */
    GridVector& operator=  (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs);

    //! Inplace addition.
    /** Preforms fast BLAS addition operation internally.
     *  @param  rhs      Constant reference state vector to be added to this.
     *  @return Reference onto this updated instantion.
     */
    GridVector& operator+= (const GridVector& rhs);

    //! Inplace addition of scaled vector.
    /** Performs fast BLAS axpy operation internally.
     *  @note This method is designed to be called as a result of operator *
     *  between GridVector and scalar.
     *  @param  rhs     Result of state vector * scalar to be added to internal
     *                  representation.
     *  @return Reference onto this updated instantion.
     */
    GridVector& operator+= (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs);

    //! Inplace subtraction.
    /** Preforms fast BLAS operation internally.
     *  @param  rhs     Constant reference on state vector to be subtracted
     *                  from this.
     *  @return Reference onto this updated instantion.
     */
    GridVector& operator-= (const GridVector& rhs);

    //! Inplace subtraction of scaled vector.
    /** Performs fast BLAS axpy operation internally.
     *  @note This method is designed to be called as a result of operator *
     *  between GridVector and scalar.
     *  @param  rhs     Result of state vector * scalar operation to be
     *                  subtracted from this.
     *  @return Reference onto this updated instantion.
    */
    GridVector& operator-= (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs);

    //! Inplace scaling with a constant.
    /** @return Reference onto this updated instantion. */
    GridVector& operator*= (const dcomp& alpha);

 // custom

    //! Generalized inplace addition.
    /** Wrapper above the BLAS ?axpy function. The operation \f[ \vec{y} =
     *  \alpha \vec{x} + \vec{y} \f] is performed in one call. Extra effective
     *  on CPUs with FMA capability.
     *  @param  alpha   Scaling factor.
     *  @param  rhs     Source state to scaled and added to this.
     *  @return Reference onto this updated instantion.
     */
    GridVector& axpy(const dcomp& alpha, const GridVector& rhs);

    //! Generalized assignement scaled with constant.
    /** @param  alpha   Scaling factor.
     *  @param  rhs     Source state to scaled and assigned to this.
     *  @return Reference onto this updated instantion.
     */
    GridVector& ax(const dcomp& alpha, const GridVector& rhs);

    //! Reduction between two state vectors.
    /** Performs sum of all elements multiplied with each other without complex
     *  conjugation.
     *  @param y    Source state vector to be contracted with.
     *  @return complex float with resulting value
     */
    dcomp reduction(const GridVector& y) const;

    //! Element wise multiplication.
    /** Multiplies each element of the vector with one element of the source
     *  vector. Not quite effective method, use with caution.
     *  @param rhs      Source state vector to be mutliplied with.
     *  @return Reference onto this updated instantion.
     */
    GridVector& element_wise_multiplication(const GridVector& rhs);

 // storage

    //! Save vector into plain text file.
    /** Prints all actual values as real and imaginary part into a file along
     *  the real part of the coordinate. The result is formated for use via
     *  gnuplot or others.
     *  @param name     Destination filename.
     */
    void save(const char* name) const;

    //! Save equidistant sampling on given range to plain text file.
    /** Evaluates function values at all sampling points in given range
     *  (including endpoints) stores the real and imaginary part along with the
     *  associatied coordinate values in plain text file in formatting for use
     *  via gnuplot.
     *  @param  name    Destination file name.
     *  @param  a       Begining of the range.
     *  @param  b       End of the range.
     *  @param  samples Total number of the sampling points including the
     *                  endpoints.
     */
    void save_equidistant(const char* name, const dfloat a, const dfloat b, const blas_int samples) const;
};

//! saves multiple states into one file respectively
bool SaveMultiGridVectorBin(const char* name, GridVector* V, const blas_int& N);
//! saves multiple states into one binary file
bool ReadMultiGridVectorBin(const char* name, GridVector* V, const FemDvrEcsGrid& g);

/** @} */
}   // namespace QSCAT
# endif //INCLUDE_GRID_VECTOR_H_
