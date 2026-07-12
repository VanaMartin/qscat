#ifndef INCLUDE_SMATRIX_LCP_H_
#define INCLUDE_SMATRIX_LCP_H_

#include "picojson/pjson.h"
#include "pjinput.h"

#ifdef linux
    #include <sec_stream.h>
#endif

#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"
#include "interface.h"

namespace QSCAT
{

namespace LCP
{
/// S Matrix variable
class SMatrix {
    bool init_;             //!< Initialization controller
    int order_;             //!< Time integration quadrature order
    int channels_;          //!< Number of channels to be investigated (VE + DA)
    int ve_channels_;       //!< Number of channels to be investigated (VE)
    int da_channels_;       //!< Number of channels to be investigated (DA)
    int steps_;             //!< Number of time qudrature steps
    int size_;              //!< Number of energy discretization points
    int init_channel_;      //!< Incident channel number
    def_float mu_;          //!< Reduced mass
    dVector coefficients_;  //!< Quadrature coeffitients
    dVector energies_;      //!< Energy range distribution
    zVector * buffer_;      //!< Integration qudrature
    zVector * s_;           //!< S matrix elements
 private:
    /// Initialization helper
    void initialize(const pjvalue& parameters);
    /// cleanup helper
    void clean();
    /// building time integration coefficients
    void make_coefficients();
 public:
    /// default constructor
    SMatrix();
    /// constructor
    SMatrix(const pjvalue& parameters);
    //SMatrix(parametersLCP& LCPp, def_float & MU, int & iChannel);
    ~SMatrix();
    /// set contribution to s-matrix elements
    void contribution(gVector & psi, gVector* states, int i);
    /// finish the time multistep (integrate multistep)
    void close_multistep(def_float time, def_float dt, dVector& ergs, const def_float& ierg, const def_float& X);
    /// save cross sections
    void cross_sections(const def_float & time, std::string& folder);
};

} // namespace LCP

} // namespace QSCAT
#endif // INCLUDE_SMATRIX_LCP_H_
