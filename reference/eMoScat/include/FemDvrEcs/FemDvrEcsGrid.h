#ifndef INCLUDE_FEM_DVR_ECS_GRID_H_
#define INCLUDE_FEM_DVR_ECS_GRID_H_

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! FEM DVR ECS class containing all information of the grid.
/** The Finite Element Method discretization of coordinate. Each element
 *  contains its separate DVR basis given by DvrGrid class. The elements may
 *  differ by length resulting in different scaling of the Gaussian weight
 *  factors \f$ w_i \f$. The DVR basis representations share their endpoints
 *  hence resulting in Gauss-Lobbato quadrature with weight factor given as \f[
 *  w_{connecting}^{n} = w_N^{n} + w_0^{n+1}, \f] where the \f$ n \f$ denotes
 *  n-th element and \f$ N \f$ denotes the quadrature order. Such formulation
 *  preserves the orthogonality of all basis functions, since outside their
 *  elements they are implicitly assumed to be zero.  The coordinate is than
 *  bend in copmlex plane at some points \f$ x_0^+, x_0^- \f$ where the
 *  coordinate is multiplied by factor \f$ e^{\pm i\theta(x-x_0^{\pm})} \f$,
 *  where \f$ \theta \f$ denotes the angle of the complex bending.
 *
 *  The class contains the DvrGrid class with sample of element basis on range
 *  (-1,1). Proper basis function is then obtained by simple shift and scalig.
 *  The class also contains the derivatives of DVR basis functions on range
 *  (-1,1) defined in all of Gaussian quadrature points stored in a matrix.
*/
class FemDvrEcsGrid : public Object, public BinaryStorageInterface
{
 protected:
    //! Total number of basis functions.
    /** \f[(n_b = t_{nel}*(n_{q} - 1) - 1),\f] where \f$n_q\f$ stands for total
     *  of basis functions on one element and \f$n_{el}\f$ represents total
     *  number of elements.
     */
    blas_int nb_;
    //! Index of the first real point of the grid.
    /** I.e. xr(ix0_neg) = xz(ix0_neg) = x0_neg. */
    blas_int ix0_neg_;
    //! Index of the last real point of the grid.
    /** i.e. xr(ix0_pos) = xz(ix0_pos) = x0_pos. */
    blas_int ix0_pos_;
    //! Total elemnts on the real part of the grid.
    /** For compatibility reasons: nr = ix0_pos (a shortcut if the radial
     *  coordinate (r >= 0) is used only).
     */
    blas_int nr_;
    //! Real (not scaled) left endpoint of the whole grid.
    dfloat x_min_;
    //! Starting point of exterior complex scaling for negative values.
    dfloat x0_neg_;
    //! Starting point of exterior complex scaling for positive values.
    dfloat x0_pos_;
    //! Right hand side complex bending point coordinate value.
    /** For compatibility reasons: R0 = x0_pos (a shortcut if the radial
     *  coordinate (r >= 0) is used only).
     */
    dfloat R0_;
    //! Real (not scaled) right endpoint of the whole grid.
    dfloat x_max_;
    //! Angle of complex scaling in degrees.
    dfloat theta_;
    //! Complex scaling factor exp(i*theta).
    dcomp eit_;
    //! Complex (scaled) left endpoint of the whole grid.
    /** For evaluation of the basis. */
    dcomp z_min_;
    //! Complex (scaled) right endpoint of the whole grid.
    /** For evaluation of the basis. */
    dcomp z_max_;
    //! Numbers of elements for gird parts.
    /** Numbers of elements respectively: Left complex scaled elements, center
     *  real elements and right complex scaled elements.
     */
    iVector nel_;
    //! Vector of complex scaled coordinate values.
    zVector xz_;
    //! Vector of complex weight factors.
    zVector wz_;
    //! Vector of element endpoints real values.
    /** Size equals \f$t_{nel}+1\f$ (tnel+1), useful to find out in which
     *  element a given point x lies.
     */
    dVector ar_;
    //! Weight factor at x0_neg.
    /** Only contribution from right real element. */
    dcomp wx0_neg_;
    //! Weight factor at x0_pos.
    /** Only contribution from left real element. */
    dcomp wx0_pos_;
    //! Total number of elements.
    blas_int tnel_;
    //! Real part of coordinate discretization points.
    /** Vector of \f$n_b\f$ of real parts of coordinate points (useful for
     *  saving functions) for integration on the real interval (x0_neg:x0_pos).
     */
    dVector xr_;
    //! Complex elements lengths
    /** Vector of \f$t_{nel}\f$ values of complex lengths of grid elements. */
    zVector aaz_;
    //! Derivatives of basis functions at grid points.
    dMatrix dLp_;
    //! DvrGrid basis on range (-1,1)
    DvrGrid g1_;

 protected:
    //! Internal binary saving procedure.
    /** Saves the class structure into binary stream.
     *  @return true on success false on error.
     */
    virtual bool save_bin_body(std::ofstream &file) const;

    //! Internal binary reading procedure.
    /** Reads the class structure from binary stream.
     *  @return true on success false on error.
     */
    virtual bool read_bin_body(std::ifstream & file);

 private:
    //! Internal initialization procedure.
    /** Computes values and assignes them to internal variables.
     *  @param  quadrature      Quadrature size.
     *  @param  total_elements  Total number of discretization elements.
     *  @param  elements        Distribution of elements to left complex,
     *                          center real and lef complex parts of the grid.
     *                          Must be at least of size 3, sum of first three
     *                          elements must be equal to total_elements.
     *  @param  aa              Real starting point and lengths of all elements
     *                          respectively. Must be at least of size
     *                          total_elements + 1.
     *  @param  alpha           Angle of the complex bending.
    */
    void initialize( blas_int quadrature,
                     blas_int total_elements,
                     const iVector& elements,
                     const dVector& aa,
                     dfloat alpha );

    //! Internal initialization procedure.
    /** Computes values and assignes them to internal variables.
     *  @param  quadrature      Quadrature size.
     *  @param  total_elements  Total number of discretization elements.
     *  @param  elements        Distribution of elements to left complex,
     *                          center real and lef complex parts of the grid.
     *                          Must be at least of size 3, sum of first three
     *                          elements must be equal to total_elements.
     *  @param  aaz             Complex starting point and lengths of all
     *                          elements respectively. Must be at least of size
     *                          total_elements + 1.
     *  @param  alpha           Angle of the complex bending.
    */
    void initialize( blas_int quadrature,
                     blas_int total_elements,
                     const iVector& elements,
                     const zVector& aaz,
                     dfloat alpha );

 public:

  // Constructors

    //! Default constructor.
    /** Uninitialized instance. */
    FemDvrEcsGrid();

    //! Basic constructor.
    /** Builds the coordinate discrete representation. Invokes internal
     *  initialization method.
     *  @param  quadrature      Quadrature size.
     *  @param  total_elements  Total number of discretization elements.
     *  @param  elements        Distribution of elements to left complex,
     *                          center real and lef complex parts of the grid.
     *                          Must be at least of size 3, sum of first three
     *                          elements must be equal to total_elements.
     *  @param  aa              Real starting point and lengths of all elements
     *                          respectively. Must be at least of size
     *                          total_elements + 1.
     *  @param  alpha           Angle of the complex bending.
     */
    FemDvrEcsGrid( blas_int quadrature,
                   blas_int total_elements,
                   const iVector& elements,
                   const dVector& aa,
                   dfloat alpha );

    //! Basic constructor.
    /** Builds the coordinate discrete representation. Invokes internal
     *  initialization method.
     *  @param  quadrature      Quadrature size.
     *  @param  total_elements  Total number of discretization elements.
     *  @param  elements        Distribution of elements to left complex,
     *                          center real and lef complex parts of the grid.
     *                          Must be at least of size 3, sum of first three
     *                          elements must be equal to total_elements.
     *  @param  aaz             Complex starting point and lengths of all
     *                          elements respectively. Must be at least of size
     *                          total_elements + 1.
     *  @param  alpha           Angle of the complex bending.
     */
    FemDvrEcsGrid( blas_int quadrature,
                   blas_int total_elements,
                   const iVector& elements,
                   const zVector& aaz,
                   dfloat alpha );

    // TODO remove, settle with pjson
    //! Constructor from parametrization structure.
    /** Interprets the parameters given in struct and calls the initialization
     *  procedure.
     *  @param  gp      grid parametrization in appropriate struct.
    */
    FemDvrEcsGrid(parameters::grid<dfloat>& gp);

    //! Copy constructor (shallow copy).
    /** Performs a shallow copy operation of all internal variables.
     *  @param  old Instance of to be copied from.
     */
    FemDvrEcsGrid(const FemDvrEcsGrid& old);

    //! Destructor.
    /** Decrease the refence count. */
    ~FemDvrEcsGrid();

 // operators

    //! Assignement operation.
    /** @param  tmp     Shallow copy of the source.
     *  @return Reference on the updated current instance.
     */
    FemDvrEcsGrid& operator= (FemDvrEcsGrid tmp);

    //! Equivalence determination.
    bool operator== (const FemDvrEcsGrid& rhs) const;

    //! Inequivalence determination.
    bool operator!= (const FemDvrEcsGrid& rhs) const;

 // modifiers

    //! Swaps internal members between instantions.
    /** Performs simple swap operation to all internal memebers.
     *  @param  rhs Instance to be swapped with.
     *  @return Reference onto this instance swapped with "rhs".
     */
    FemDvrEcsGrid& swap(FemDvrEcsGrid& rhs);

    //! Deep copy operation.
    /** Performs the actual deep copy operation.
     *  @return A new instance of the same object.
     */
    FemDvrEcsGrid copy() const;

 // accessors

    //! Real value of coordinate discretization point.
    /** @param  i   Index of desired discretization point. Must be less than
     *              total number of discretization basis functions.
     *  @return Real part of the actual value at discretization point.
     */
    const dfloat& xr(const blas_int& i) const;

    //! Complex value of coordinate discretization point.
    /** @param  i   Index of desired discretization point. Must be less than
     *              total number of discretization basis functions.
     *  @return Actual complex value at discretization point.
     */
    const dcomp& x(const blas_int& i) const;

    //! Weight factor.
    /** @param  i   Index of desired discretization point. Must be less than
     *              total number of discretization basis functions.
     *  @return Complex value of wieght factor corresponding basis function
     *          associated to given discretization point.
     */
    const dcomp& w(const blas_int& i) const;

    //! Elements in discretization regions.
    /** @param  i   Index of discretization region (part), 0: left complex,
     *              1: real center, 2: right complex.
     *  @return Total number of elements in desired region.
     */
    const blas_int& nel(const blas_int& i) const;

    //! Element ending point.
    /** For ease of acces to ending points, complex values of connecting points
     *  on the coordinate. The length of element is then aax(i) - aaz(i-1).
     *  @param  i   Index of desired element. Must be less than total number of
     *              elements.
     *  @return Complex value of element end point.
     */
    const dcomp& aaz(const blas_int& i) const;

    //! Quadrature weight factor.
    /** @param  i   Index of desired weight factor on element sample element
     *              (-1,1).  Must be less than total number of basis functions
     *              on one grid (quadrature order).
     *  @return Weight factor of desired basis function on sample element.
     */
    dfloat wq(const blas_int& i);

    //! Basis functions derivatives on sample element.
    /** Retrieves the derivatives of basis functions on the sample (-1,1) element.
     *  @return Values arranged in square matrix of size given by quadrature.
     */
    const dMatrix& dlp() const;

    //! Derivatives of basis function at one point of sample element.
    /** Retrieves the derivatives of one basis function at one point on the
     *  sample (-1,1) element
     *  @param  i   Index of the basis function to be derived.
     *  @param  j   Index of the discretiztion point where to evaluate
     *              derivative.
     *  @return Value of the derivative at desired point.
    */
    const dfloat& dlp(const blas_int& i, const blas_int& j) const;

    // TODO REMOVE
    //! total number of basis elements (DEPRECATED)
    /**
        returns the total number of discretization points of the whole grid. TODO DEPRECATED: use get_size() instead.
    */
    blas_int nb() const;

    //! Total number of basis functions (discretization points).
    blas_int get_size() const;

    //! Qudrature order.
    blas_int quadrature() const;

    //! Total number of elements.
    blas_int tnel() const;

    //! Number of real part basis functions.
    /** @note: Valid only if left complex region has exactly zero elements.
        @return Total number of basis functions, resp. discretization points on
                real part of the grid.
    */
    blas_int nr() const;

    //! End of coordinate real part.
    /** @return Value of coordinate at positive complex bending (end point if
     *  no positive bending present).
     */
    const dfloat& x_pos() const;

    //! start of coordinate real part
    /** @return Value of coordinate at negative complex bending (start point if
     *  no negative bending present).
     */
    const dfloat& x_neg() const;

    //! Sample DVR basis.
    /** @return Sample DVR basis computed on sample element on (-1,1) range. */
    const DvrGrid& dvr() const;

 // Functions

    //! End of the element.
    /** @param  i   Index of given element.
     *  @return Index of the last basis function defined on given element.
     */
    blas_int get_element_end(const blas_int& i) const;

    //! Start of the element
    /** @param  i   Index of given element.
     *  @return Index of the first basis function defined on given element.
     */
    blas_int get_element_start(const blas_int& i) const;

    //! Index of element.
    /** @param  X   Point on real part of the coordinate.
     *  @return Index of the element containing given point.
     */
    blas_int get_element_index(const dfloat& X) const;

    //! End of element.
    /** @param  X   Point on real part of the coordinate.
     *  @return Index of last basis function of element containing given point.
     */
    blas_int get_element_end_x(const dfloat& X) const;

    //! Start of element.
    /** @param  X   Point on real part of the coordinate.
     *  @return Index of the first basis function of element containing given
     *          point.
     */
    blas_int get_element_start_x(const dfloat& X) const;

    //! Basis function function value.
    /** @param  i       Index of the basis function.
     *  @param  j       Point inside assumed element where the value will be
     *                  evaluated.
     *  @param  start   Index of the first basis function defined on the element.
     *  @param  end     Index of the las basis function defined on the element.
     *  @return Function value of given basis function at given point.
     */
    dcomp basis_function_value( const blas_int& i,
                                const dfloat& x,
                                const blas_int& start,
                                const blas_int& end ) const;

    //! Real part of element ending points.
    /** For ease of acces to ending points.
     *  @param  i   Element index, must be less or equal to total number
     *              elements.
     *  @return Real part of complex values of connecting points on the
     *          coordinate.
    */
    const dfloat& ar(const blas_int& i) const;
};

/** @} */
}   // QSCA# endif //INCLUDE_FEM_DVR_ECS_GRID_H_
# endif //INCLUDE_FEM_DVR_ECS_GRID_H_
