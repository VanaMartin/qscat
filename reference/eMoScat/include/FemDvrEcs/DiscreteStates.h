#ifndef INCLUDE_DISCRETE_STATES_H_
#define INCLUDE_DISCRETE_STATES_H_

#include "bessel.h"
#include "coulomb.h"
#include "common.h"
#include "Storage.h"
#include "Arrays.h"
#include "input.h"

#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"

#include "FemDvrEcs/OperatorDiagonal.h"
#include "FemDvrEcs/OperatorFull.h"

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! Bound states & resonances in 1D form multiple grids.
/** Class stores all necessary information about the system described by two
 *  different parametrizations, solve the eigen states for both parametrizations
 *  and compare the spectral values.  If the spectral values are indpendent of
 *  the parametrization, the value is considered as discrete state.
 */
class DiscreteStates1D : public Object, public BinaryStorageInterface
{
 protected:
    //! Total number of determined discrete states.
    blas_int num_discrete_states_;
    //! Energy threshold for discrete state determination process.
    dfloat max_energy_;
    //! Machine precision threshold.
    dfloat machine_epsilon_;
    //! Precision threshold.
    /** Convergence boundary condition. */
    dfloat precision_;
    //! Eigenvalues.
    dcomp *spectrum_;
    //! Eigenstates.
    GridVector *states_;
    //! Base coordinate discretization.
    /** In this representation the resulting states will be stored. */
    FemDvrEcsGrid grid_;
 private:
    //! Internal read from file method.
    /** Reads states from external plain text file.
     *  @param  name    Name of source file.
     *  @param  grid    Associated grid. If states differs in discetization,
     *                  method fails.
     *  @return True on success false otherwise.
    */
    bool read(const char *name, const FemDvrEcsGrid& grid);

    //! Cleaning method.
    /** Cleans all internal variables. */
    void clean();

 protected:
    //! Internal binary saving procedure.
    /** Saves the class structure into binary stream.
     *  @return true on success false on error.
     */
    virtual bool save_bin_body(std::ofstream& file) const;

    //! Internal binary reading procedure.
    /** Reads the class structure from binary stream.
     *  @return true on success false on error.
     */
    virtual bool read_bin_body(std::ifstream& file);
 public:
 // constructors

    //! Default constructor.
    /** The result returns false with init() method. */
    DiscreteStates1D();

    //! Basic constructor with simple input.
    /** Initialization with default precision, minimal and maximal eigenvalue
     *  boundaries.
     *  @param  num_grids   Total number of different parametrizations.
     *  @param  grids       Array of different grids. Must be at least of size
     *                      num_grids.
     *  @param  potentials  Array of potential reprezetnations. Must be at
     *                      least of size num_grids.
     *  @param  mass        Effective mass of the model.
     */
    DiscreteStates1D( blas_int num_grids,
                      FemDvrEcsGrid* g,
                      GridVector* potentials,
                      dfloat mass );

    //! constructor with extended input
    /** Initialization with specified precision, minimal and maximal eigenvalue
     *  boundaries.
     *  @param  num_grids   Total number of different parametrizations.
     *  @param  grids       Array of different grids. Must be at least of size
     *                      num_grids.
     *  @param  potentials  Array of potential reprezetnations. Must be at
     *                      least of size num_grids.
     *  @param  mass        Effective mass of the model.
     *  @param  a_prec      Convergence precision threshold.
     *  @param  min_eig     Minimal eigenvalue estimation boundary.
     *  @param  max_eig     Maximal eigenvalue estimation boundary.
     */
    DiscreteStates1D( blas_int num_grids,
                      FemDvrEcsGrid* grids,
                      GridVector* potentials,
                      dfloat mass,
                      dfloat a_prec,
                      dfloat min_eig,
                      dfloat max_eig );

    //! Constructor from plain text file.
    /** Reads the internal variables from the text file, compares with provided
     *  grid, builds the states representations.
     *  @param  name    Name of source file.
     *  @param  grid    Associated grid.
    */
    DiscreteStates1D(char *name, const FemDvrEcsGrid& grid);

    //! Destructor.
    ~DiscreteStates1D();

 // accessors

    // FIXME: DEPRECATED
    //! Retrieve results and store them.
    /** Inserts all results into the provided arrays.
     *  Arguments: E .. complex vector, pointer reference, eigenvalues, size of num_discrete_states
     *             V .. complex vector, two dimensional array,
    */
    void retrieve(zVector*& E, zVector **V);

    //! Obtain discrete state.
    /** Inserts desired discrete state vector into provided GridVector.
     *  @param  destination     Output destination.
     *  @param  i               Index of retrieved state. must be less than
     *                          num_discrete_states.
     */
    void get_state(GridVector& destination, blas_int i);

    //! Obtain 1st discrete state grid reference.
    /** Retrieves a constant reference onto the first state coordinate
     *  discretization.
     */
    const FemDvrEcsGrid& get_grid() const;

    //! Obtian discrete state eigen value.
    /** Retrieves the computed discrete eigen value.
     *  @param  i   Index of the desired eigen value. Must be less than
     *              num_discrete_states.
     *  @return Result value.
     */
    dcomp get_energy(blas_int i);

    //! Total number of determined discrete states.
    blas_int number_of_states();

 // storage

    //! Load states from file.
    /** Loads states from plain text file, compares with provided grid.
     *  @param  name    Name of source file.
     *  @param  grid    Associated grid.
     */
    bool from_file(const char *name, FemDvrEcsGrid& grid);

    //! Store computed states into plain text file.
    /** Saving all computed states into the text file.
     *  @param  name    Name of destination file.
     */
    void save_states(const char *name);
};

/** @} */
} // namespace QSCAT

#endif // INCLUDE_DISCRETE_STATES_H_
