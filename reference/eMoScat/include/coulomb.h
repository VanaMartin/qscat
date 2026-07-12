#ifndef __COULOMB__
    #define __COULOMB__

    #include "coulcc.h"


namespace QSCAT
{
    namespace coulomb
    {
        /// Coulombic functions

        /// The solutions for radial coulombic problem in one dimension.
        /// The naming conventions are similar to bessel namespace.

        ///  Spherical coulomb function of the first type
        cc_comp sF(const cc_comp xz, const cc_comp eta, int l);

        ///  Spherical coulomb function of the second type
        cc_comp sG(const cc_comp xz, const cc_comp eta, int l);

        ///  Spherical coulomb function of the first type
        cc_comp sH1(const cc_comp xz, const cc_comp eta, int l);

        ///  Spherical coulomb function of the second type
        cc_comp sH2(const cc_comp xz, const cc_comp eta, int l);

        ///  Energy normalized coulomb function of the first type
        cc_comp sF_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l);

        ///  Energy normalized coulomb function of the second type
        cc_comp sG_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l);

        ///  Energy normalized coulomb function of the first type
        cc_comp sH1_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l);

        ///  Energy normalized coulomb function of the second type
        cc_comp sH2_en(const cc_comp x, const cc_comp k, const cc_comp z, const double m, const int l);

    }
}

#endif
