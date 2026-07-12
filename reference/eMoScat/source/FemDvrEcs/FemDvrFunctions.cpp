#include <cassert>
#include <complex>
#include <math.h>

#include "FemDvrEcs/FemDvrFunctions.h"

using std::abs;
using std::sqrt;

namespace QSCAT
{

dfloat Gamma0to3(dfloat x)  // Function Gamma0to3 returns the value of Euler Gamma function at x in range o (0,3)
{
        dfloat a[18];
        dfloat p, t, val;
        int i;

        a[0] =  1.0;
        a[1] =  0.4227843350984678;
        a[2] =  0.4118403304263672;
        a[3] =  0.0815769192502609;
        a[4] =  0.0742490106800904;
        a[5] = -0.0002669810333484;
        a[6] =  0.0111540360240344;
        a[7] = -0.0028525821446197;
        a[8] =  0.0021036287024598;
        a[9] = -0.0009184843690991;
        a[10] =  0.0004874227944768;
        a[11] = -0.0002347204018919;
        a[12] =  0.0001115339519666;
        a[13] = -0.0000478747983834;
        a[14] =  0.0000175102727179;
        a[15] = -0.0000049203750904;
        a[16] =  0.0000009199156407;
        a[17] = -0.0000000839940496;

        if (x <= 1.0) {
            t = x;
        } else if (x <= 2.0) {
            t = x - 1.0;
        } else {
            t = x - 2.0;
        }

        p = a[17];
        for (i=16;i>=0;i--) {
            p = t * p + a[i];
        }

        if (x <= 1.0) {
            val = p / (x * (x + 1.0));
        } else if (x <= 2.0) {
            val = p / x;
        } else {
            val = p;
        }

        return val;
    }

// Coefficients
void RecCoef(int kind, int n, dfloat alpha, dfloat beta, dVector& b, dVector& a, dfloat& mu)  // Function returns coeficients a(j), b(j) of recurent relation b p(x) = (x-a)p(x) - b p(x)
{
    int nm1;
    dfloat ab, abi, a2b2;

    nm1 = n - 1;
    a.fill(0);
    b.fill(0);
    switch (kind) {
        case(1): // Legendre polynomials p(x) on (-1, +1), w(x) = 1
            {
                mu = dfloat(2);
                a.fill(dfloat(0));
                for (int i=0; i<nm1; i++){
                    b[i] = dfloat(i + 1) / sqrt(dfloat(4) * (i+1)*(i+1)- dfloat(1));
                }
            }
            break;
        case(2): // Chebyshev polynomials of the first kind t(x) on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
            {
                 mu = dfloat(pi);
                 a[0] = dfloat(0);
                 b[0] = dfloat(sqrt2);
                 for (int i=1;i<nm1;i++){
                     a[i] = dfloat(0);
                     b[i] = dfloat(1/2);
                 }
            }
            break;
        case(3): // Chebyshev polynomials of the second kind u(x) on (-1, +1), w(x) =  sqrt(1 - x*x)
            {
                mu = pi / 2.0;
                for (int i=0; i<nm1;i++) {
                    a[i] = 0.0;
                    b[i] = 0.5;
                }
            }
            break;
        case(4): // Hermite polynomials on (-infinity, +infinity), w(x) =  exp(-x**2)
            {
                mu = sqrtpi;
                for (int i=0;i<nm1;i++){
                    a[i] = 0;
                    b[i] = sqrt( (i+1.0)/2.0);
                }
            }
            break;
        case(5): // Jacobi polynomials p(alpha, beta)(x), alpha and beta greater than -1 on (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta
            {
                ab  = alpha + beta;
                abi = 2.0 + ab;
                mu = pow(dfloat(2),dfloat(ab + 1)) * Gamma0to3(alpha + 1.0) * Gamma0to3(beta + 1.0) / Gamma0to3(abi);
                a[0] = (beta - alpha) / abi;
                b[0] = sqrt(4.0 * (1.0 + alpha) * (1.0 + beta) / ((abi + 1.0) * abi * abi));
                a2b2 = beta * beta - alpha * alpha;
                for (int i=1; i<nm1; i++){
                    abi  = 2.0 * (i+1) + ab;
                    a[i] = a2b2 / ((abi - 2.0) * abi);
                    b[i] = sqrt(4.0 * (i+1) * (i+1 + alpha) * (i+1 + beta) * (i+1 + ab) / ((abi * abi - 1) * abi * abi));
                }
                abi  = 2.0* n + ab;
                a[n-1] = a2b2 / ((abi - 2.0) * abi);
            }
            break;
        case(6): // Laguerre polynomials l(alpha)(x), alpha greater than -1 on (0, +infinity), w(x) = exp(-x) * x**alpha
            {
                mu = Gamma0to3(alpha + 1.0);
                for (int i=0; i < nm1; ++i) {
                    a[i] = 2.0 * (i+1.) - 1. + alpha;
                    b[i] = sqrt( (i+1.) * (i + 1. + alpha) );
                }
                a[n-1] = 2. * n - 1. + alpha;
            }
            break;
    }
}

// Auxiliary functionn for Gauss Quad
dfloat gbshift(dfloat shift, int m, const dVector& t, const dVector& b) // Auxiliary function
{
    dfloat al;
    al = t[0] - shift;
    for (int i=1; i<m; i++){
        al = t[i] - shift - pow(b[i-1],2)/al;
    }
    return 1.0/al;
}

// Solving tridiagonal matrix for computing weights
void gbtql2(int n, dVector& d, dVector& e, dVector& z, int& ierr)
{
    int i, j, k, l, m, ii, mml;
    dfloat machep, c, f, b, s, p, g, r;

    /*     ========== machep is a machine dependent parameter specifying
                    the relative precision of floating point arithmetic.
                    machep = 16.0d0**(-13) for long form arithmetic
                    on s360 ========== */
    machep=1.0e-14;

    // d = 0, e != 0

    ierr = 0;
    if (n!=1) {
        e[n-1] = 0.0;
        for (l=1;l<n;l++){
            j = 0;
flag105:    // ========== Look for small sub-diagonal element ==========
            for (m = l;m < n;m++) {
                if (std::abs(e[m-1]) <= machep * (std::abs(d[m-1]) + std::abs(d[m]))) { break; }
            }
            p = d[l-1];
            if (m == l) goto flag240;
            if (j == 30) goto flag1000;
            j = j + 1;

            // ========== form shift ==========
            g = (d[l] - p) / (2.0 * e[l-1]);
            r =  sqrt(g*g+1.0);
            g = d[m-1] - p + e[l-1] / (g +  functions::sign(r,g));
            s = 1.0;
            c = 1.0;
            p = 0.0;
            mml = m - l;

            // ========== for i=m-1 step -1 until l DO -- ==========
            for (ii=1; ii<mml+1; ii++) {
                i = m - ii;
                f = s * e[i-1];
                b = c * e[i-1];
                if ( std::abs(f) < std::abs(g)) { goto flag150;}
                c = g / f;
                r =  sqrt(c*c+1.0);
                e[i] = f * r;
                s = 1.0/r;
                c = c * s;
                goto flag160;
flag150:
                s = f / g;
                r =  sqrt(s*s+1.0);
                e[i] = g * r;
                c = 1.0/r;
                s = s * c;
flag160:
                g = d[i] - p;
                r = (d[i-1] - g)* s + 2.0 * c * b;
                p = s * r;
                d[i] = g + p;
                g = c * r - b;
            // ========== form first component of vector ==========
                f = z[i];
                z[i] = s * z[i-1] + c * f;
                z[i-1] = c * z[i-1] - s * f;
            }
            d[l-1] = d[l-1] - p;
            e[l-1] = g;
            e[m-1] = 0.0;
            goto flag105;
flag240:;
        }
        // ========== order eigenvalues and eigenvectors ==========

        for (ii = 2; ii<=n; ii++){
            i = ii - 1;
            k = i;
            p = d[i-1];
            for (j=ii; j<=n; j++) {
                if (d[j-1] > p) { goto flag260; }
                k = j;
                p = d[j-1];
flag260:;           }
            if (k==i) {goto flag300;}
            d[k-1] = d[i-1];
            d[i-1] = p;
            p = z[i-1];
            z[i-1] = z[k-1];
            z[k-1] = p;
flag300:;       }
        goto flag1001;
        //     ========== set error -- no convergence to an
        //                eigenvalue after 30 iterations ==========
flag1000:       ierr = l;
flag1001:       return;
    }
}

// Gaussian Quadrature
void GaussQuad(int kind, int n, dfloat alpha, dfloat beta, dVector& b, dVector& t, dVector& w, dfloat x_min, dfloat x_max, int kpts)
{
    int nm1 = n - 1;
    int ierr;
    dfloat gam, t1;
    dfloat mu;
    RecCoef(kind,n,alpha,beta,b,t,mu);

    /*
    The matrix of coefficients is assumed to be symmetric.
    The array t contains the diagonal elements,
    the array b the off-diagonal elements.
    Make appropriate changes in the lower right 2 by 2
    submatrix if endpoints are to be part of the grid
    */

    switch (kpts) {
        case(1): {
            t[n] = gbshift(x_min,nm1,t,b) * pow(b[nm1-1],2) + x_min;
        }
        break;
        case(2):
        {
            gam = gbshift(x_min,nm1,t,b);
            t1 = ((x_min - x_max) / (gbshift(x_max,nm1,t,b) - gam));
            b[nm1-1] = sqrt(t1);
            t[n-1] = x_min + gam * t1;
        }
    }

    /*
     Note that the indices of the elements of b run from 0 to n-2
     and thus the value of b(n-1) is arbitrary.
     Now compute the eigenvalues of the symmetric tridiagonal
     matrix, which has been modified as necessary.
     The method used is a QL-type method with origin shifting
    */

    w.fill(0.0);
    w[0] = 1.0;

    gbtql2(n, t, b, w, ierr);

    for (int k=0; k<n; k++) {
        w[k] = mu * w[k] * w[k];
    }
}

// Initialization of Gaussian Quadrature Grid
void GLo_Quad(int n, dVector& x, dVector& w)
{
    int kind = 1;
    int kpts = 2; // endpoints included
    dfloat alpha = dfloat(0), beta = dfloat(0), x_min = dfloat(-1), x_max = dfloat(1);

    dVector scr(n);

    GaussQuad(kind,n,alpha,beta,scr,x,w,x_min,x_max,kpts);
    return;
}

// Roots of pade approximation for Crank-Nicolson variable
void Pade_Roots(zVector& roots, int order)
{
    switch (order) {
        case 1:
            roots[0] = dcomp(2.0,0);
            break;
        case 2:
            roots[0] = dcomp(-3.0, sqrt(3.0));
            roots[1] = conj(roots[0]);
            break;
        case 3:
            roots[0] = dcomp(-4.6443707092521712, 0.0);
            roots[1] = dcomp(-3.6778146453739144, 3.5087619195674433);
            roots[2] = conj(roots[1]);
            break;
        case 4:
            roots[0] = dcomp(-5.7924212056407443, 1.7344682578690075);
            roots[1] = conj(roots[0]);
            roots[2] = dcomp(-4.2075787943592557, 5.3148360837135054);
            roots[3] = conj(roots[2]);
            break;
        case 5:
            roots[0] = dcomp(-7.2934771906592865, 0.0);
            roots[1] = dcomp(-6.7039127983070663, 3.4853228323663954);
            roots[2] = conj(roots[1]);
            roots[3] = dcomp(-4.6493486063632905, 7.1420458406759528);
            roots[4] = conj(roots[3]);
            break;
        case 6:
            roots[0] = dcomp(-8.4967187917267279, 1.7350193464627312);
            roots[1] = conj(roots[0]);
            roots[2] = dcomp(-7.4714167126516293, 5.2525446228942513);
            roots[3] = conj(roots[2]);
            roots[4] = dcomp(-5.0318644956216428, 8.9853459073078851);
            roots[5] = conj(roots[4]);
            break;
        case 7:
            roots[0] = dcomp(-9.9435737170558713, 0.0);
            roots[1] = dcomp(-9.5165810563092579, 3.4785721222610731);
            roots[2] = conj(roots[1]);
            roots[3] = dcomp(-8.1402783272762749, 7.0343480954195063);
            roots[4] = conj(roots[3]);
            roots[5] = dcomp(-5.3713537578865315, 10.841388261433498);
            roots[6] = conj(roots[5]);
            break;
        case 8:
            roots[0] = dcomp(-11.175772086526170, 1.7352288907055729);
            roots[1] = conj(roots[0]);
            roots[2] = dcomp(-10.409681581273764, 5.2323503052850549);
            roots[3] = conj(roots[2]);
            roots[4] = dcomp(-8.7365784344048048, 8.8288850009430782);
            roots[5] = conj(roots[4]);
            roots[6] = dcomp(-5.6779678977952610, 12.707822597209754);
            roots[7] = conj(roots[6]);
            break;
        case 9:
            roots[0] = dcomp(-12.594038363429937, 0.0);
            roots[1] = dcomp(-12.258735808548546, 3.4756967669617250);
            roots[2] = conj(roots[1]);
            roots[3] = dcomp(-11.208843639015563, 6.9963138357721872);
            roots[4] = conj(roots[3]);
            roots[5] = dcomp(-9.2768797743607806, 10.634543350871302);
            roots[6] = conj(roots[5]);
            roots[7] = dcomp(-5.9585215963601425, 14.582927376684364);
            roots[8] = conj(roots[7]);
            break;
        case 10:
            roots[0] = dcomp(-13.844089810854492, 1.7353303909024429);
            roots[1] = conj(roots[0]);
            roots[2] = dcomp(-13.230581930953741, 5.2231358416001798);
            roots[3] = conj(roots[2]);
            roots[4] = dcomp(-11.935056657175572, 8.7698943778838641);
            roots[5] = conj(roots[4]);
            roots[6] = dcomp(-9.7724391337179992, 12.449970964943134);
            roots[7] = conj(roots[6]);
            roots[8] = dcomp(-6.2178324672981964, 16.465398918147175);
            roots[9] = conj(roots[8]);
            break;
        case 15:
            roots[0] = dcomp(-20.546219332644956, 0.0);
            roots[1] = dcomp(-20.341827992880136, 3.4727778389009119);
            roots[2] = conj(roots[1]);
            roots[3] = dcomp(-19.719134456792558, 6.9613424228655327);
            roots[4] = conj(roots[3]);
            roots[5] = dcomp(-18.647198641217939, 10.484517790475234);
            roots[6] = conj(roots[5]);
            roots[7] = dcomp(-17.064918104596681, 14.068787251034091);
            roots[8] = conj(roots[7]);
            roots[9] = dcomp(-14.858793985884308, 17.757965242243031);
            roots[10] = conj(roots[9]);
            roots[11] = dcomp(-11.800303427329295, 21.639998275507146);
            roots[12] = conj(roots[11]);
            roots[13] = dcomp(-7.2947137249766045, 25.959002141520839);
            roots[14] = conj(roots[13]);
            break;
        case 20:
            roots[0] = dcomp(-27.134848566306627, 1.7354725099115966);
            roots[1] = conj(roots[0]);
            roots[2] = dcomp(-26.825194287213204, 5.2108002943589951);
            roots[3] = conj(roots[2]);
            roots[4] = dcomp(-26.197644949154326, 8.6997298235829248);
            roots[5] = conj(roots[4]);
            roots[6] = dcomp(-25.234562633219703, 12.212959740104793);
            roots[7] = conj(roots[6]);
            roots[8] = dcomp(-23.906181604999975, 15.764116868494900);
            roots[9] = conj(roots[8]);
            roots[10] = dcomp(-22.165160667462304, 19.372186483657157);
            roots[11] = conj(roots[10]);
            roots[12] = dcomp(-19.935524957720782, 23.066229457032493);
            roots[13] = conj(roots[12]);
            roots[14] = dcomp(-17.087791453700064, 26.896090546839399);
            roots[15] = conj(roots[14]);
            roots[16] = dcomp(-13.371053756590380, 30.962612375847237);
            roots[17] = conj(roots[16]);
            roots[18] = dcomp(-8.1420371236326346, 35.543738137770912);
            roots[19] = conj(roots[18]);
    }
}

dVector equidistant_quadrature(int order, int size)
{
    assert(order < 10);     // Higher orders are deprecated
  //

    // First build one step vector
    int base = size / std::pow(2,order-1) + 1;

    // First order: Trapezoid rule
    dVector res(base);
    res.fill(1); res[0] = res[res.get_size()-1] = 0.5;
    order--;

    int o = 0;
    while(order--) {
        o++;
        int s = res.get_size();
        dVector work( s * 2 - 1);
        work.fill(0);
        for (int i=0; i<s; ++i) {   // Build double sampled grid
            work[i] += res[i];
            work[s+i-1] += res[i];
        }
        work *= std::pow(4.0, o);   // Romberg
        for (int i=0; i<s; ++i)
            work[2*i] -= 2.0 * res[i];    // previous step has double size of "h"
        work *= 1.0/ ( std::pow(4.0,o) - 1.0 );
        res = work;
    }
    return res;
}

} // namespace QSCAT
