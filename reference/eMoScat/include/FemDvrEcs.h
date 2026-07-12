#ifndef INCLUDE_FEM_DVR_ECS_H_
#define INCLUDE_FEM_DVR_ECS_H_

#include "blas.h"
#include "bessel.h"
#include "common.h"
#include "Arrays.h"
#include "input.h"

// Using own defined namespace is ok.
//using namespace functions;

/// Finite  Elements  Method  (FEM) - Discrete  Variable  Representation  (DVR) - Exterior Complex Scaling (ECS) namespace.

///
/// This library define  basic types and  methods for one dimensional FEM-DVR-ECS
/// representation  of a  function. The DVR grid  specifies representation on one
/// element  of given  length.  Through  one-basis-function  connection  multiple
/// elements  are  used in  FEM approach.  The  ECS allows  extending  the use of
/// Dirichlet  boundary  condition  to a  wide range  of  complex  funcions  with
/// infinite support set.
///
/// The  Grid  Vector  class  combines  the  FEM-DVR-ECS  representation  of  the
/// coordinate  with  contiguous storage  via ARRAYS::Vector<Z>, allowing the use
/// of previously  defined methods via BLAS operation. The function values in the
/// Vector are  multiplied by  the weights  factor overwhelmingly simplifying the
/// integration  process to a  simple vector  dot product.  The evaluation of the
/// function  values at  arbitrary points  within the grid range is iplemented as
/// a member method.
///
/// The operator  classes  defines a representation  of  operations  on the given
/// space of functions.  The Full operator  allows an  arbitrary linear operation
/// on the  function space. The  diagonal and Row Compressed operators carry some
/// basic  restrictions.  The generator of  Quantum  Kinetic  Energy  operator is
/// implemented for the Full and Row Compressed case.
///

/// Some useful functions necessary for DVR basis computation
#include "FemDvrEcs/FemDvrFunctions.h"

/// Dvr basis definitions
#include "FemDvrEcs/DvrGrid.h"

/// Dvr + fem + ecs definitions
#include "FemDvrEcs/FemDvrEcsGrid.h"

/// Representaion of functions on dicretized variable
#include "FemDvrEcs/GridVector.h"

/// Projection operator alpha * |phi><phi| + beta * Id
#include "FemDvrEcs/Projector.h"

/// Simple diagonal operator defined on functions represented by GridVector
#include "FemDvrEcs/OperatorDiagonal.h"

/// Full operator defined on functions represented by GridVector
#include "FemDvrEcs/OperatorFull.h"

/// Operator represented via CSR format defined on functions represented by GridVector
#include "FemDvrEcs/OperatorRowCompressed.h"

/// Kinetic term generators as Full/CSR operators, Hamiltonian build helper
#include "FemDvrEcs/KineticEnergy.h"

/// Discrete states finder class
#include "FemDvrEcs/DiscreteStates.h"

namespace QSCAT
{
OperatorFull buildFullHamiltonian(FemDvrEcsGrid& g, GridVector& v, const dfloat& mu);

}   // namspace QSCAT

#endif // ndef INCLUDE_FEM_DVR_ECS_H_
