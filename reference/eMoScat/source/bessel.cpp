//#include <iostream>
#include <fstream>
#include <stdio.h>
#include <complex>
#include <string>
/*
    This file contains the wraping macro for the fortran coulcc library and
    some necessary functions.
*/

#include "bessel.h"

using namespace std;
typedef std::complex<double> cc_comp;

namespace QSCAT
{

cc_comp BesselJ(const cc_comp & xz, int l)
{
    cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
    long int c1, c2, ifail;             // "coulcc" auxiliary variables

    long int i = 1;
    c1 = 3;
    c2 = 1;
    ifail = 0;
    cc_comp z0 = 0.0;
    cc_comp nu = cc_comp(l)-0.5;
    cc_comp xza = xz;
    coulcc(xza, z0, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
    return ZFC*sqrt(2.0*xz/pi_b);
}

cc_comp sphBesselJ(const cc_comp & xz, int l)
{
    cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
    long int c1, c2,  ifail;                // "coulcc" auxiliary variables
    long int i = 1;
    c1 = 3;
    c2 = 1;
    ifail = 0;
    cc_comp z0 = 0.0;
    cc_comp nu = cc_comp(l);
    cc_comp xza = xz;
    coulcc(xza, z0, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
    return ZFC;
}

cc_comp sphBesselJEn(const cc_comp & xz,const cc_comp & k,const double & m,const int & l)
{
    cc_comp xp = k*xz;
    return sqrt(2.0*m*k/pi_b)*sphBesselJ(xp,l)*xz;
}

cc_comp NeumannY(cc_comp & xz, int l)
{
    cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
    long int c1, c2, ifail;             // "coulcc" auxiliary variables
    long int i = 1;
    c1 = 1;
    c2 = 1;
    ifail = 0;
    cc_comp z0 = 0.0;
    cc_comp nu = cc_comp(l)-0.5;
    coulcc(xz, z0, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
    return ZGC*(2.0*xz/pi_b);
}

cc_comp sphNeumannY(cc_comp & xz, int l)
{
    cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
    long int c1, c2, ifail;             // "coulcc" auxiliary variables
    long int i = 1;
    c1 = 1;
    c2 = 1;
    ifail = 0;
    cc_comp z0 = 0.0;
    cc_comp nu = cc_comp(l);
    coulcc(xz, z0, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
    return ZGC;
}

cc_comp sphNeumannYEn(cc_comp & xz, cc_comp & k, double & m, int & l)
{
    cc_comp xp = k*xz;
    return sqrt(2.0*m*k/pi_b)*sphNeumannY(xp,l)*xz;
}

cc_comp sphHankel1(cc_comp & xz, int l)
{
    cc_comp ZFC, ZGC, ZFCP, ZGCP, ZSIG; // "coulcc" auxiliary variables
    long int c1, c2, ifail;             // "coulcc" auxiliary variables
    long int i = 1;
    c1 = 11;
    c2 = 1;
    ifail = 0;
    cc_comp z0 = 0.0;
    cc_comp nu = cc_comp(l);
    coulcc(xz, z0, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
    return ZGC;
}

cc_comp sphHankel1En(const cc_comp & xz, const cc_comp & k, const double & m, const int & l)
{
    cc_comp xp = k*xz;
    return sqrt(2.0*m*k/pi_b)*sphHankel1(xp,l)*xz;
}

cc_comp sphHankel2(cc_comp & xz, int l)
{
    cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
    long int c1, c2, ifail;             // "coulcc" auxiliary variables
    long int i = 1;
    c1 = 21;
    c2 = 1;
    ifail = 0;
    cc_comp z0 = 0.0;
    cc_comp nu = cc_comp(l);
    coulcc(xz, z0, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
    return ZGC;
}

cc_comp sphHankel2En(const cc_comp & xz, const cc_comp & k, const double & m, const int & l)
{
    cc_comp xp = k*xz;
    return sqrt(2.0*m*k/pi_b)*sphHankel2(xp,l)*xz;
}

} // namespace QSCAT
