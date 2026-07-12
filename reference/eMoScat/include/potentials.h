using namespace std;

#include "pjinput.h"

// FIXME clean the input types

/// Model potentials

/// Contains several potential functions to compose the two-dimensional model
namespace QSCAT
{
namespace potentials
{
// Forward declarations
    template<typename T, typename Z>
    const Z V_zero(const Z & y, const parameters::model_2D<T> & mp);
    template<typename T, typename Z>
    const Z V_int(const Z & x, const Z & y, const parameters::model_2D<T> & mp);
    template<typename T, typename Z>
    const Z V_eff_el(const Z & x, const Z & y, const parameters::model_2D<T> & mp);
    template<typename T, typename Z>
    const Z  Attached_electron(const Z & x,  const parameters::model_2D<T> & mp);
    template<typename T, typename Z>
    const Z V_cfg(const T& J, const T& mu, const Z& X);
    template<typename T, typename Z>
    const Z Lambda_Spec_R(const Z & x, const Z & R, const parameters::model_2D<T> & mp);

    template<typename T, typename Z>
    const Z foo(const T & i, const Z & k)
    {
        std::cout << i*k << " That it was!" << std::endl;
        return i*k;
    }
/// Full 2D potential
    template<typename T, typename Z>
    const Z potential_2D(const Z & x, const Z & y, const parameters::model_2D<T> & mp)
    {
        return V_zero(y,mp) + mp.l*(mp.l + 1.0)/(2.0*x*x) + V_int(x,y,mp);
    }

/// Molecular asymptotic potential (electron far from the molecule)
    template<typename T, typename Z>
    const Z V_zero(const Z & y, const parameters::model_2D<T> & mp)
    {
        return Z(1.0)*mp.D_0 * (exp(-2.0*mp.alpha_0*(y - mp.R_0)) - 2.0*exp(-mp.alpha_0*(y - mp.R_0)));
    }

/// Interaction part of the 2D potential
    template<typename T, typename Z>
    const Z V_int(const Z & x, const Z & y, const parameters::model_2D<T> & mp)
    {
        Z l_0 = (mp.lambda_c - mp.lambda_inf)*(1 + exp(mp.lambda_1*(mp.R_c - mp.R_lambda)));
        return (-1.0)*(mp.lambda_inf + l_0/(1.0 + exp(mp.lambda_1*(y - mp.R_lambda))))*exp(-mp.alpha_c*x*x);
    }

/// Electronic effective potential
    template<typename T, typename Z>
    const Z V_eff_el(const Z & x, const Z & y, const parameters::model_2D<T> & mp)
    {
        Z l = Z(mp.l);
        return l*(l + 1.0)/(2.0*x*x) + V_int(x,y,mp);
    }

    template<typename T, typename Z>
    const Z  Attached_electron(const Z & x,  const parameters::model_2D<T> & mp)
    {
        return -mp.lambda_inf*exp(-mp.alpha_c*pow(x,2)) +  Z(mp.l)*(Z(mp.l) + 1.0)/(2.0*pow(x,2));
    }

    template<typename T, typename Z>
    const Z V_cfg(const T& J, const T& mu, const Z& X)
    {
        if (J == T(0)) return Z(0);
        return  Z(J)*(Z(J)+1.0)/(2.0*pow(X,2));
    }

/// Special potential for obtaining the disctrete state: Parametrically depenedent on R_nuclear
    template<typename T, typename Z>
    const Z Lambda_Spec_R(const Z & x, const Z & R, const parameters::model_2D<T> & mp)
    {
        T R_d = 4.0;    // specified as for lambda spec 2
        T c_d = 1.5;
        T l_s = 25;     // Lambda -inf
        return (1.0 + mp.l)*mp.l/(2.0*x*x) - (mp.lambda_inf + (l_s - mp.lambda_inf)/(1.0 + exp(c_d*(R - R_d))))*exp(-mp.alpha_c*x*x);
    }

    //grid_vector_2D & Make_fxy(const Z (*f)(const Z &, const Z &,const parameters::model_2D<T> &), const parameters::model_2D<T> & mp);

    template<typename T, typename Z>
    void savePotential2D(const char *filename, const Z (*f)(const Z&, const Z&, const parameters::model_2D<T>&),  const parameters::model_2D<T> &mp, int nX, T aX, T bX, int nY, T aY, T bY)
    {
        FILE * file;
        fopen_s(&file,filename,"w");
        fprintf(file,"# Potential saved on the equidistant grid with %d x %d values", nY, nX);
        fprintf(file,"#Coordinate X   \tCoordinate Y   \tReal part of z     \tImaginary part of z\n");
        int xS = (aX==0.0)? 1:0;
        int yS = (aY==0.0)? 1:0;
        for (int i=yS; i<nY+1; ++i){
            T y = (bY - aY)/nY*i + aY;
            for (int j=xS; j<nX+1; ++j){
                T x = (bX - aX)/nX*j + aX;
                Z v = f(x, y, mp);
                fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\n", y, x, real(v), imag(v));
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

/// Molecular asymptotic potential (electron far from the molecule)
    def_comp MorsePotential(const def_comp& y, const Parameters& p);

/// Interaction part of the 2D potential
    def_comp LambdaInteraction(const def_comp& x, const def_comp& y, const Parameters& p);

/// Full 2D potential
    def_comp Neutral2dPotential(const def_comp& x, const def_comp& y, const Parameters& p);

/// Asymptotic interaction
    def_comp AsymptoticLambda(const def_comp& x, const Parameters& p);

/// Effective electronic part (lambda + cent. barrier) of full potential
    def_comp ElectronicLambdaInteraction(const def_comp& x, const def_comp& y, const Parameters& p);

} // namespace QSCAT
