#ifndef INCLUDE_ZOOM_FILER_2D_H_
#define INCLUDE_ZOOM_FILER_2D_H_
namespace QSCAT
{
/** \addtogroup FemDvrEcs2d
* @{ */

//! Magnification operation on given range for use in HSV model evaluation.
class ZoomFilter
{
 public:
    blas_int xa, xb;    //!< x-coordinate range
    blas_int ya, yb;    //!< y-coordinate range
    double mag;         //!< maginfication factor
 private:
    bool init_;         //!< initialization controller
 public:
    //! default constructor
    ZoomFilter() : xa(0), xb(0), ya(0), yb(0), mag(1), init_(false) {}
    //! constructor
    ZoomFilter(blas_int Xa, blas_int Xb, blas_int Ya, blas_int Yb, double Mag)
        : xa(Xa), xb(Xb), ya(Ya), yb(Yb), mag(Mag), init_(true)
        {
            assert(xa>=0);
            assert(ya>=0);
            assert(xb>=0);
            assert(yb>=0);
            assert(mag > 0);
        }
    //! initialization check
    bool init() const { return init_; }
    //! Operator assignement
    ZoomFilter& operator=(const ZoomFilter& rhs)
        { xa = rhs.xa; xb = rhs.xb; ya = rhs.ya; yb = rhs.yb; mag = rhs.mag; init_ = rhs.init_; return *this;}
};

/** @} */
}   //  namespace QSCAT
#endif // INCLUDE_ZOOM_FILER_2D_H_
