#ifndef TIME_DEPENDENT_MODEL_2D_INCLUDE_H_
#define TIME_DEPENDENT_MODEL_2D_INCLUDE_H_

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


/// The main class for the two dimensional model of collisions
/*!
    All necessarry operations and preparations are gathered within this class.
*/
class TimeDependentModel2D : public BinaryStorageInterface
{
    bool init_;                         //!< Initialization controller
    parameters2D& params_;              //!< Reference to the model parametrization
    const pjvalue& par_;                //!< picojson parametrization
    femGrid xgrid_;                     //!< Electronic coordinate grid
    femGrid ygrid_;                     //!< Nuclear coordinate grid
    femGrid2D grid_;                    //!< Two dimensional grid
    gVector2D potential_;               //!< Vector storing the potential
    bool xtestfunction_use_;            //!< Test function in electronic coordinate usage controller
    MultiTestFunction2d *xtestfunction_;  //!< Test function in the electronic coordinate
    bool ytestfunction_use_;            //!< Test function in nuclear coordinate usage controller
    MultiTestFunction2d *ytestfunction_;  //!< Test function in the nuclear coordinate
    gVector2D psi_;                     //!< Evolved wave function
    def_comp initial_energy_;           //!< Initial energy
    bool chebyshev_use_;                //!< Chebyshev method usage controller
    Chebyshev2D *chebyshev_;            //!< Pointer to the evolution operator approximated by Chebyshev method
    bool crank_nicolson_use_;           //!< Crank-Nicolson method usage controller
    CrankNicolson2D *crank_nicolson_;   //!< Pointer to the evolution operator approximated by Crank-Nicolson method
    def_float time_;                    //!< Evolution time
    def_float dt_;                      //!< Evolutoin time step
    int loop_size_;                     //!< Number of loop steps
    std::string folder_;                //!< Output folder
    // Projection onto diabatic state statistics
    gVector2D discrete_state_;          //!< Pointer to the Projection State
    gVector psi_discrete_;              //!< Projected wvefunction onto the discrete state
    dBuffer normalization_;             //!< Normalistation in a buffer
    dBuffer radius_;                    //!< Mean internuclear distance in a buffer
    zOperatorD radius_operator_;        //!< Operator storing the operation
    EquidistantProjector2d *ep_;        //!< Equidistant output projection
    EquidistantProjector2d *hsv_;       //!< Equidistant output to HSV-RGB projection
    gVector *vibrations_;               //!< Vibrational wavefunctions in V_res potential
    ofstream population_file_;          //!< Vibrational level populations in time output
 private:
    /// Main internal initialization helper
    //void initialize(parametersGrid& gpx, parametersGrid& gpy);
    /// Testfunction initialization helper
    void initialize_testfunctions();
    /// Potential curves initialization helper
    void initialize_potential();
    /// Cleanup helper
    void clean();

 protected:
    /// Internal save to file helper
    virtual bool save_bin_body(std::ofstream& file) const;
    /// Internal read from file helper
    virtual bool read_bin_body(std::ifstream& file);

 public:
    /// Constructor
    //TimeDependentModel2D(parameters2D& m2dp, parametersGrid& gpx, parametersGrid& gpy, gVector2D* phid = NULL);
    /// Constructor
    TimeDependentModel2D(const pjvalue& parameters, gVector2D* phid = NULL);
    /// Destructor
    ~TimeDependentModel2D();
    /// Reset internals
    TimeDependentModel2D& set(parameters2D& p, parametersGrid& gpx, parametersGrid& gpy);
    /// Perform a discrete projection
    void discrete_projection();
    /// Several time steps of evolution and testfunction correlations
    void multistep();
};

/** @} */
}; // namespace QSCAT
#endif // TIME_DEPENDENT_MODEL_2D_INCLUDE_H

