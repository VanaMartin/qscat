//
//  This file contains the wraping macro for the fortran coulcc library and
//  some necessary functions.
//

#ifndef _QSCAT_BESSEL_INCLUDE_H_
#define _QSCAT_BESSEL_INCLUDE_H_

#include "coulcc.h"

/*! \defgroup Bessel
    Special functions solving either radial or radial part of spherically
    symmetrical problem.
    Naming conventions: PREFIX sph - spherical case
                        SUFFIX En - energy normalized case
*/
namespace QSCAT
{
/** \addgroup Bessel
 * @{
 */

/// Bessel function of first kind
cc_comp BesselJ(const cc_comp & xz, int l);

/// Spherical Bessel function
/*!
    Spherical Bessel function of a complex argument. Calls the external function
    coulcc from fortran file coulcc.lib. The file differs for debug and release
    mode.
    Function: j_l(r)
*/
cc_comp sphBesselJ(const cc_comp & xz, int l);

/// Energy normalized shperical Bessel funciton.
/*!
    The energy normalized spherical Bessel function.
    Function: sqrt(2mk/pi) k y_l(kr)
*/
cc_comp sphBesselJEn(const cc_comp & xz,const cc_comp & k,const double & m,const int & l);

/// Neumann function Y, also known as Bessel of the second kind
cc_comp NeumannY(cc_comp & xz, int l);

/// Spherical Neumann function
cc_comp sphNeumannY(cc_comp & xz, int l);

/// Energy normalized spherical Neumann function.
/*!
    The energy normalized spherical bessel function.
    Function: sqrt(2mk/pi)*k*y_l(kr)
*/
cc_comp sphNeumannYEn(cc_comp & xz, cc_comp & k, double & m, int & l);

/// Spherical Hankel function of the first kind.
/*!
    Spherical  Hankel  function  of the first  kind (h+),  representing the
    outgoing wave.
*/
cc_comp sphHankel1(cc_comp & xz, int l);

/// Energy normalized spherical Hankel function of the first kind.
/*!
    The energy normalized spherical bessel function.
    Function: sqrt(2mk/pi)*k*y_l(kr)
*/
cc_comp sphHankel1En(const cc_comp & xz, const cc_comp & k, const double & m, const int & l);

/// Spherical Hankel function of the second kind.
/*!
    Spherical  Hankel  function  of the second kind (h-),  representing the
    incoming wave.
*/
cc_comp sphHankel2(cc_comp & xz, int l);

/// Energy normalized spherical Hankel function of the second kind.
/*!
    The energy normalized spherical bessel function.
    Function: sqrt(2mk/pi)*k*y_l(kr)
*/
cc_comp sphHankel2En(const cc_comp & xz, const cc_comp & k, const double & m, const int & l);

/** @} */
} // namespace QSCAT
#endif // _QSCAT_BESSEL_INCLUDE_H_
