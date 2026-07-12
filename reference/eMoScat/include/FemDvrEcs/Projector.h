#ifndef INCLUDE_PROJECTOR_H_
#define INCLUDE_PROJECTOR_H_

namespace QSCAT
{
/** \addtogroup FemDvrEcs
* @{ */

//! Special class for storing  Feschbach projection operators as symmetric full  matrix class.
/** The  projector is  generalized as  the sum of two operators: \f[ \hat{P} =
 * \alpha |\psi\rangle \langle \psi | + \beta 1. \f] The definition is used  in
 * the multiplication operation.
 */
class Projector : public Object
{
    //! Representaion size.
    blas_int size_;
    //! Internal scaling factor.
    dcomp alpha_;
    //! Internal scaling factor.
    dcomp beta_;
 public:
    //! operator internal representaion via dense matrix
    zMatrix body_;
 private:
    //! Internal initialization.
    /** Initializes the interanal varibles, with default scalar factor values.
     *  @param  psi     Source state gridvector.
     */
    void initialize(GridVector& psi);

    //! Extended initialization.
    /** Initializes the internal variables.
     *  @param  psi     Source state vector.
     *  @param  alpha   Pprojector part scaling factor.
     *  @param  beta    Identity part scaling factor.
     */
    void initialize(GridVector& psi, const dcomp& alpha, const dcomp& beta);
 public:

    //! Default constructor.
    /** Sets internals with purely default values (uninitialized).*/
    Projector();

    //! Constructor.
    /** Sets internals to simplified version: \f$ \alpha = 1, \beta = 0\f$.
     *  @param  psi     Source state for outer product body.
     */
    Projector(GridVector& psi);

    //! Extended constructor.
    /** Sets internals to general projector form.
     *  @param  psi     Source state vector.
     *  @param  alpha   Pprojector part scaling factor.
     *  @param  beta    Identity part scaling factor.
     */
    Projector(GridVector& psi, const dcomp& alpha, const dcomp& beta);

    //! Destructor.
    ~Projector();

    //! Projection operation.
    /** Action of the operator on given state vector.
     *  @param  rhs     Source state vector to operate on. Must be defined on
     *                  the same coordinate discretization as the state vector
     *                  used for projector initialization.
     *  @returns New state vector containing the output.
     */
    GridVector operator* (const GridVector& rhs) const;
};

/** @} */
}   //namespace QSCAT
# endif //INCLUDE_PROJECTOR_H_
