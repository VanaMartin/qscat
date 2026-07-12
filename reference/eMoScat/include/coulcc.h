
/*
    Solution of coulomb problem

    laplace f(x) + ( 1 - eta / x - lambda ( lambda + 1 ) / x^2 ) f(x) = 0
    
    extern "C" void coulcc_(
        cc_comp * z,            // X - complex coordinate
        cc_comp * zc,           // eta1 - coulombic potential factor: 0 - for bessel branch, nonzero - coulomb branch
        cc_comp * nu,           // zlmin - sarting value of lambda, rotational quantum number (fractional)
        long int * i,           // number of integer spaced values lambda (mostly =1) 
        cc_comp * ZFC,          // function complex value F (J) 
        cc_comp * ZGC,          // function complex value G (Y)
        cc_comp * ZFCP,         // derivative complex value F' (J')
        cc_comp * ZGCP,         // derivative complex value G' (G')
        cc_comp * ZSIG,         // coulomb phase shift
        long int * c1,          // MODE 
        long int * c2,          // KFN (0,-1 - coulomb), (1 sph Bessel), (2 cyl. Bessel), (3 mod. cyl. Bessel) 
        long int * ifail );     // status function 

*/

#ifndef __COULCC__
    #define __COULCC__
    const double pi_b       = 3.1415926535897932;       // Pi 

    using namespace std;
    typedef std::complex<double> cc_comp;

    // A slightly different definition of the function call for UNIX systems

    #ifdef __cplusplus
        #ifdef linux
            extern "C" void coulcc_(cc_comp * z, cc_comp * zc, cc_comp * nu, long int * i, cc_comp * ZFC, cc_comp * ZGC, cc_comp * ZFCP, cc_comp * ZGCP, cc_comp * ZSIG, long int * c1, long int * c2, long int * ifail );
            #define coulcc(z, zc, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail) coulcc_(&z, &zc, &nu, &i, &ZFC, &ZGC, &ZFCP, &ZGCP, &ZSIG, &c1, &c2, &ifail)
        #else
            extern "C" void COULCC(cc_comp * z, cc_comp * zc, cc_comp * nu, long int * i, cc_comp * ZFC, cc_comp * ZGC, cc_comp * ZFCP, cc_comp * ZGCP, cc_comp * ZSIG, long int * c1, long int * c2, long int * ifail );
            #define coulcc(z, zc, nu, i, ZFC, ZGC, ZFCP, ZGCP, ZSIG, c1, c2, ifail) COULCC(&z, &zc, &nu, &i, &ZFC, &ZGC, &ZFCP, &ZGCP, &ZSIG, &c1, &c2, &ifail)
        #endif  /* linux */
    #endif /* __cplusplus */

#endif
