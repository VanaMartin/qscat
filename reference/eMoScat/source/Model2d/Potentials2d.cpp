
#include "common.h"
#include "potentials.h"

namespace QSCAT
{
    def_comp MorsePotential(const def_comp& y, const Parameters& p)
    {
        return p("D_0")*( exp(-2.0*p("alpha_0") * (y - p("R_0"))) - 2.0*exp(-p("alpha_0")*(y - p("R_0"))) );
    }

    def_comp LambdaInteraction(const def_comp& x, const def_comp& y, const Parameters& p)
    {
        def_comp lambda0 = (p("lambda_c") - p("lambda_inf"))*(1 + exp(p("lambda_1")*(p("R_c") - p("R_lambda"))));
        return (-1.0)*(p("lambda_inf") + lambda0/(1.0 + exp(p("lambda_1")*(y - p("R_lambda")))))*exp(-p("alpha_c")*x*x);
    }

    def_comp Neutral2dPotential(const def_comp& x, const def_comp& y, const Parameters& p)
    {
        def_comp l = p("impulsemomentum");
        return MorsePotential(y,p) + l*(l + 1.0)/(2.0*x*x) + LambdaInteraction(x,y,p);
    }

    def_comp AsymptoticLambda(const def_comp& x, const Parameters& p)
    {
        return (-1.0)*p("lambda_inf") * exp(-p("alpha_c")*x*x) + p("impulsemomentum") * (p("impulsemomentum")+1.0) /(2.0 * x * x);
    }

    def_comp ElectronicLambdaInteraction(const def_comp& x, const def_comp& y, const Parameters& p)
    {
        def_comp l = p["impulsemomentum"].asDouble();
        return l*(l + 1.0)/(2.0*x*x) + LambdaInteraction(x,y,p);
    }
}
