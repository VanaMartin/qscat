#ifndef MULTI_TEST_FUNCTION_2D_INCLUDE_H_
#define MULTI_TEST_FUNCTION_2D_INCLUDE_H_

#include "FemDvrEcs2d.h"
#include "Model2d/TestFunctionInterface2d.h"
#include "Model2d/TestFunction2d.h"
#include "Model2d/DiracTestFunction2d.h"
#include "Model2d/FluxTestFunction2d.h"

#include "picojson/pjson.h"
#include "pjinput.h"

namespace QSCAT
{
/** \addtogroup Model2d
* @{ */

/// Multiple test functions wrapper
class MultiTestFunction2d : public BinaryStorageInterface
{
    bool init_;             //!< Initialization controller
    int channels_;          //!< Number of channels in the assumed coordinate asymptotic region
    int initial_channel_;   //!< The number of the initial channel: (if it is in different coordinate, it is to be set to <0)
    char axis_;             //!< Which coordinate asymptotic region
    char method_;           //!< Which S-matrix derivation methods are to be used
    def_float mu_;          //!< Reduced mass in the given coordinate
    dVector energies_;      //!< Energy distribution in given range (shifted by subtraction of transversal energies)
    TestFunction2d *testfunctions_;               //!< The Tannor&Weeks method test functions in given channel
    DiracTestFunction2d *dirac_testfunctions_;    //!< The T&W with Dirac delta test functions in given channel
    FluxTestFunction2d *flux_testfunctions_;      //!< The probability flux test functions in given channel

    zMatrix tannor_s_matrix_;       //!< Tannor&Weeks approach obtained S-Matrix elements
    zMatrix dirac_s_matrix_;        //!< Modified Tannor&Weeks with delta approach obtained S-Matrix elements
    zMatrix flux_s_matrix_;         //!< Probability flux approach obtained S-Matrix elements

    zVector fourier_coefficients_;              // Initial state fourier coeffitients (the appropriate correlation impulse value may differ for each state, the ifc must be computed for appropriate value)
    dMatrix tannor_cross_sections_; //!< Tannor&Weeks approach obtained cross sections
    dMatrix dirac_cross_sections_;  //!< Modified Tannor&Weeks with delta approach obtained cross sections
    dMatrix flux_cross_sections_;   //!< Probability flux approach obtained cross sections
    std::string folder_;            //!< Output folder

    const pjvalue& p_;
    const pjvalue& tp_;
protected:
    virtual bool save_bin_body(std::ofstream& file) const;
    virtual bool read_bin_body(std::ifstream& file);
public:
    MultiTestFunction2d();
    MultiTestFunction2d(femGrid2D *g, const pjvalue& params, const pjvalue& tfparams, const def_comp& init_erg);
    ~MultiTestFunction2d();
    void step_buffer(gVector2D& Psi, const int& step);
    void close_multistep(const def_float& t, const def_float& dt);
    void cross_sections(const def_float& time);
};

/** @} */
}; // namesapce QSCAT

#endif // MULTI_TEST_FUNCTION_2D_INCLUDE_H

