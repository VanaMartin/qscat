// Constants declarations
// The real constants are stored in double

#ifndef INCLUDE_COMMON_H_
#define INCLUDE_COMMON_H_

#include <complex>

#define qfloat double
#define qint long long
#define qcomp std:complex<qfloat>

#ifndef def_float
    #define def_float double        // Invoking constants through def_float()
    #define def_comp std::complex<def_float>
#endif

// Real Constants for easy use

const def_float pi      = 3.1415926535897932;       // Pi
const def_float pi2     = 6.2831853071795865;       // 2*Pi
const def_float sqrtpi  = 1.7724538509055160;       // Sqrt(Pi)
const def_float sqrtpi2 = 2.506628274631000;        // Sqrt(2*Pi)
const def_float sqrt2   = 1.4142135623730950;       // Sqrt(2)
const def_float eu      = 2.7182818284590452;       // e

// Complex number defnition, Intel MKL complex number definition
#if !defined(dComplex)
    #define dComplex  std::complex<double>
#endif


// Globaly defined function for skipping the rest of line & multiple others
void skipline(std::ifstream& file, const int & i);

// Globally defined binary operator functions
//inline comp operator^(const comp& a, const int& e);
//inline double operator^(const double& a, const int& e);

// Common complex constants

const def_comp zone = def_comp(1.0, 0.0);
const def_comp imu = def_comp(0.0, 1.0);
const def_comp zzero = def_comp(0.0, 0.0);


namespace functions {

    bool check_mem_size(const int elements);

    // declarations
    template<typename T, typename Z>
    const T operator^ (const T & v, const Z & k);
    template <class T>
    int sgn(T val);
    template <class T>
    T sign(T x, T y);
    template<typename T, typename Z>
    Z Gaussian(const T & x, const T & x_0, const T & s, const T & p);
    template<typename T>
    T Sine(const T & x, const T & d, const int & j);

    // implemenatations

    template<typename T, typename Z>
    const T operator^ (const T & v, const Z & k)
    {
        return pow(v,k);
    }

    template <class T>
    int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }
    template <class T>
    T sign(T x, T y)
    {
        if (std::abs(y)!=0){
            return x*y/std::abs(y);
        } else {
            return x;
        }
    }
    template<typename Z>
    Z Gaussian(const Z & x, const Z & x_0, const Z & s, const Z & p)  // Complex Gaussian wave packet
    {
        Z core;
        Z val;
        core = exp(-pow((x-x_0),2)/(4*pow(s,2)))/sqrt(sqrt(pi*2)*s);
        val = exp(imu*p*x)*core;
        return val;
    }

    template<typename T>
    T Sine(const T & x, const T & d, const int & j)
    {
        T val = sqrt(2/d)*sin(j*pi*x/d);
        return val;
    }

}

#endif // INCLUDE_COMMON_H_
