#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <complex>
#include <cassert>
#include <math.h>
#include <stdlib.h>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "common.h"
#include "Arrays.h"
#include "input.h"
#include "FemDvrEcs.h"

namespace QSCAT
{
    OperatorFull buildFullHamiltonian(FemDvrEcsGrid& g, GridVector& v, const dfloat& mu)
    {
        OperatorFull H(g);
        return ( H.add_kinetic_term(mu) ) += v;
    }

}   // namespace QSCAT
