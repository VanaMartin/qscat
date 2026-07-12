#include <fstream>
#include <stdio.h>
#include <complex>
#include <string>

#include "coulomb.h"

namespace QSCAT
{
    //  Spherical coulomb function of the first type
    cc_comp coulomb::sF(const cc_comp xz, const cc_comp eta, int l)
    {
        cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
        long int c1, c2, ifail;             // "coulcc" auxiliary variables

        long int i = 1;
        c1 = 2;
        c2 = -1;
        ifail = 0;
        cc_comp nu = cc_comp(l);
        cc_comp xz_ = xz, eta_ = eta;       // auxiliary non-const
        coulcc(xz_, eta_, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
        return ZFC;
    }

    //  Spherical coulomb function of the second type
    cc_comp coulomb::sG(const cc_comp xz, const cc_comp eta, int l)
    {
        cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
        long int c1, c2, ifail;             // "coulcc" auxiliary variables

        long int i = 1;
        c1 = 2;
        c2 = -1;
        ifail = 0;
        cc_comp nu = cc_comp(l);
        cc_comp xz_ = xz, eta_ = eta;       // auxiliary non-const
        coulcc(xz_, eta_, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
        return ZGC;
    }

    //  Spherical coulomb function of the second type
    cc_comp coulomb::sH1(const cc_comp xz, const cc_comp eta, int l)
    {
        cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
        long int c1, c2, ifail;             // "coulcc" auxiliary variables

        long int i = 1;
        c1 = 11;
        c2 = 0;
        ifail = 0;
        cc_comp nu = cc_comp(l);
        cc_comp xz_ = xz, eta_ = eta;       // auxiliary non-const
        coulcc(xz_, eta_, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
        return ZFC;
    }

    //  Spherical coulomb function of the second type
    cc_comp coulomb::sH2(const cc_comp xz, const cc_comp eta, int l)
    {
        cc_comp ZFC, ZGC, ZFCP, ZGCP,ZSIG;  // "coulcc" auxiliary variables
        long int c1, c2, ifail;             // "coulcc" auxiliary variables

        long int i = 1;
        c1 = 21;
        c2 = 0;
        ifail = 0;
        cc_comp nu = cc_comp(l);
        cc_comp xz_ = xz, eta_ = eta;       // auxiliary non-const
        coulcc(xz_, eta_, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail);
        return ZGC;
    }

    //  Energy normalized coulomb function of the first type
    cc_comp coulomb::sF_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l)
    {
        cc_comp xk = k*x;
        cc_comp eta = m * z / k;
        return sqrt( 2.0 * m / (pi_b*k) )* coulomb::sF(xk, eta, l) ;
    }

    //  Energy normalized coulomb function of the second type
    cc_comp coulomb::sG_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l)
    {
        cc_comp xk = k*x;
        cc_comp eta = m * z / k;
        return sqrt( 2.0 * m / (k*pi_b) )* coulomb::sG(xk, eta, l) ;
    }

    //  Energy normalized coulomb function of the first type
    cc_comp coulomb::sH1_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l)
    {
        cc_comp xk = k*x;
        cc_comp eta = m * z / k;
        return sqrt(2.0 * m / (pi_b*k) )* coulomb::sH1(xk, eta, l) ;
    }

    //  Energy normalized coulomb function of the second type
    cc_comp coulomb::sH2_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l)
    {
        cc_comp xk = k*x;
        cc_comp eta = m * z / k;
        return sqrt( 2.0 * m / (pi_b*k) )* coulomb::sH2(xk, eta, l) ;
    }
} // namespace QSCAT
