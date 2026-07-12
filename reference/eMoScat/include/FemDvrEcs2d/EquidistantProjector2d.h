#ifndef INCLUDE_EQUIDISTANT_PROJECTOR_2D_H_
#define INCLUDE_EQUIDISTANT_PROJECTOR_2D_H_
namespace QSCAT
{
/** \addtogroup FemDvrEcs2d
* @{ */

//! Projection operator to equidistant grid in csr format
/** General purpose of this operator is to compute the coefficients of the
 *  projection onto the equidistant grid only once and than use the
 *  coefficients repeatedly on different wavefunctions.  In DVR representation
 *  the computation of the function value in arbitrary point is no longer
 *  diagonal and therefore the value is given by value of the normalized basis
 *  function at given point scaled by representation coefficient \f[ f(x) =
 *  \sum_{i=0}^Nf(x_i)b_i(x) = \sum_{i=0}^Nfrac{1}{\sqrt(w_i)}f_ib_i(x), \f] where
 *  \f$b_i(x)\f$ represents the i-th basis function value at point \f$x\f$. By
 *  storing all basis function values in the right order in a sparse matrix we
 *  create a operator which apllied on any state vector returns the actual
 *  values of the desired function at given points. */
class EquidistantProjector2d : public Object
{
 protected:
    //! Associated 2d-coordinate discretization.
    /** Source of all discretization parameters. */
    FemDvrEcsGrid2d grid_;
    //! Operator sparse representation.
    /** Computed values of operation representation in right order. */
    RowCompressedMatrix<dcomp> body_;
    //! Stored input values.
    /** Auxiliary storage space for holding input values.*/
    mutable zVector values_;
    //! X-coordinate range values.
    /** Values of resulting X-coordinate discretization. There will be no
     *  weight factors in the resulting representation. */
    dVector x_coordinate_;
    //! Y-coordinate range values.
    /** Values of resulting Y-coordinate discretization. There will be no
     *  weight factors in the resulting representation. */
    dVector y_coordinate_;
    //! Total X-coordinate samples.
    /** Also represent the size of x_coordinate_ vector.*/
    size_t x_samples_;
    //! Total Y-coordinate samples.
    /** Also represent the size of y_coordinate_ vector.*/
    size_t y_samples_;
    //! X-coordinate starting position.
    /** Starting point of the X-coordinate segment.*/
    dfloat x_start_;
    //! X-coordinate endin position.
    /** Ending point of X-coordinate segment.*/
    dfloat x_end_;
    //! Y-coordinate starting position.
    /** Starting point of the Y-coordinate segment.*/
    dfloat y_start_;
    //! Y-coordinate ending position.
    /** Ending point of Y-coordinate segment.*/
    dfloat y_end_;
    //! Zooming filter.
    /** Auxiliary operator allowing rescaling of some values (2d region).*/
    ZoomFilter zoom_;
 public:
    //! Constructor.
    /** Basic constructor, builds the operator on region given by X,Y ranges
     *  with equidistand distribution of points and including the boundary
     *  points.
     *  @param  grid        Associate coordinate discretization of all future
     *                      inputs to the operator.
     *  @param  x_sampling  Total number of equidistant representation
     *                      X-coordinate discretization.
     *  @param  x_sampling  Total number of equidistant representation
     *                      Y-coordinate discretization.
     *  @param  xa          X-coordinate range starting point.
     *  @param  xb          X-coordinate range ending point.
     *  @param  ya          Y-coordinate range starting point.
     *  @param  yb          Y-coordinate range ending point.
     */
    EquidistantProjector2d( const FemDvrEcsGrid2d& grid, size_t x_sampling,
                            size_t y_sampling, dfloat xa, dfloat xb, dfloat ya,
                            dfloat yb);

    //! Destructor.
    ~EquidistantProjector2d();

    //! Evaluate given state (discard the internal data).
    /** Push input to internal auxiliary storage.
     *  @param  state   Source two-dimensional state vector, must be
     *                  represented on the same XY coordinate discretization
     *                  used to initialize the operator.
     */
    void operator<< (const GridVector2d& state) const;

    //! Export given state to file.
    /** Compute operation result and store it into a plain text file.
     *  @param  state       Source two-dimensional state vector, must be
     *                      represented on the same XY-coordinate.
     *                      discretization used to initialize the operator.
     *  @param  filename    Path to save target. If the file exists, it will be
     *                      overwritten silently.
     */
    void export_state(GridVector2d& state, const char* filename) const;

    //! Export internal data to file.
    /** Store the operation result into a plain text file. The input is taken
     *  from internal auxiliary memory therefore the operator << must be
     *  previously used. Values will be computed at this point and than stored
     *  into given file.
     *  @param  state       Source two-dimensional state vector, must be
     *                      represented on the same XY-coordinate.
     *                      discretization used to initialize the operator.
     *  @param  filename    Path to save target. If the file exists, it will be
     *                      overwritten silently.
     */
    void export_state(const char* filename) const;

    //! Export given state to bitmap via hsv model.
    /** Compute operation result and store it into a plain text file.
     *  @param  state       Source two-dimensional state vector, must be
     *                      represented on the same XY-coordinate.
     *                      discretization used to initialize the operator.
     *  @param  filename    Path to save target. If the file exists, it will be
     *                      overwritten silently.
     *  @param  magnitude   Radius in polar coordinates at which all values are
     *                      mapped to exactly half Lightness.
     *  @param  phase       Phase factor for Hue shift.
     */
    void export_state_hsv( GridVector2d& state, const char *filename,
                           dfloat magnitude, dfloat phase=0.0) const;

    //! Export internal data to bitmap via hsv model.
    /** Compute the HSL model of stored internal state and write it into
     *  specified path as a bitmap.
     *  @param  filename    Path to save target. If the file exists, it will be
     *                      overwritten silently.
     *  @param  magnitude   Radius in polar coordinates at which all values are
     *                      mapped to exactly half Lightness.
     *  @param  phase       Phase factor for Hue shift.
     */
    void export_state_hsv(const char* filename, dfloat magnitude, dfloat phase=0.0) const;

    //! Set internal zoom funciton to given region.
    /** Creates magnitude scaling region where different HSL mapping will take
     *  place.
     *  @param  x_start     Starting point of mapped area in X-coordinate.
     *  @param  x_end       Ending point of mapped area in X-coordinate.
     *  @param  y_start     Starting point of mapped area in Y-coordinate.
     *  @param  y_end       Ending point of mapped area in Y-coordinate.
     *  @param  magnitude   Overriding value of Lightness magnitude (see
     *                      export_state_hsv() for details).
     */
    void set_zoom_filter(dfloat x_start, dfloat x_end, dfloat y_start, dfloat y_end, dfloat magnitude);
};

/** @} */
}   //  namespace QSCAT
#endif // INCLUDE_EQUIDISTANT_PROJECTOR_2D_H_
