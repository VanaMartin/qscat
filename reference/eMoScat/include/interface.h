//
//	Definition of commonly used types for source files.
//

#ifndef __INTERFACE__
    #define __INTERFACE__

    // All headers used below should be defined here
    #include "Arrays.h"
    #include "input.h"
    #include "FemDvrEcs.h"
    #include "FemDvrEcs2d.h"
    #include "crank_nicolson.hpp"
    #include "chebyshev.hpp"
    #include "moeller.hpp"

namespace QSCAT {

    // FEM_DVR_ECS
    typedef FemDvrEcsGrid          femGrid;
    typedef GridVector             gVector;
    typedef OperatorFull           zOperatorF;
    typedef OperatorDiagonal       zOperatorD;
    typedef OperatorRowCompressed  zOperatorC;
    //typedef poly_grid_vector       zPolyVector;
    typedef DiscreteStates1D       DiscreteStates;

    typedef CrankNicolson<def_float, def_comp, zOperatorC, gVector> CrankNicolson1D;
    typedef Chebyshev<def_float, def_comp, zOperatorC, gVector> Chebyshev1D;

    typedef MoellerOperator<def_float, CrankNicolson1D, zOperatorC, gVector> MoellerCN1D;

    // FEM_DVR_ECS_2D
    typedef FemDvrEcsGrid2d             femGrid2D;
    typedef GridVector2d                gVector2D;
    typedef Operator2dRowCompressed     zOperator2D;

    typedef Chebyshev<def_float, def_comp, zOperator2D, gVector2D>       Chebyshev2D;
    typedef CrankNicolson<def_float, def_comp, zOperator2D, gVector2D>   CrankNicolson2D;

    // parameters
    typedef parameters::grid<def_float>          parametersGrid;
    typedef parameters::multi_grid<def_float>    parametersMultiGrid;
    typedef parameters::NRM<def_float>           parametersNRM;
    typedef parameters::LCP<def_float>           parametersLCP;
    typedef parameters::model_2D<def_float>      parameters2D;
    typedef parameters::evolution<def_float>     parametersEvolution;
    typedef parameters::initial_state<def_float> parametersInitState;
    typedef parameters::testfunction<def_float>  parametersTestfunction;

    // TODO: Add this method to all classes
    template<class U>
    void save(const U& X, const std::string& name)
    {
        X.save(name.c_str());
    }
}

    // ARRAYS

    // ARRAYS DOCUMENTATION
    /// \class dBuffer interface.h
    /// \brief floating point type buffer via \link ARRAYS::Buffer \endlink
    /// \copydetails ARRAYS::Buffer
    /// \class zBuffer interface.h
    /// \brief complex floating point type buffer via \link ARRAYS::Buffer \endlink
    /// \copydetails ARRAYS::Buffer
    /// \class iVector interface.h
    /// \brief integer type vector via \link ARRAYS::Vector \endlink
    /// \copydetails ARRAYS::Vector
    /// \class dVector interface.h
    /// \brief floating point type vector via \link ARRAYS::Vector \endlink
    /// \copydetails ARRAYS::Vector
    /// \class zVector interface.h
    /// \brief complex floating point type vector via \link ARRAYS::Vector \endlink
    /// \copydetails ARRAYS::Vector
    /// \class dMatrix interface.h
    /// \brief floating point type matrix via \link ARRAYS::Matrix \endlink
    /// \copydetails ARRAYS::Matrix
    /// \class zMatrix interface.h
    /// \brief complex floating point type matrix via \link ARRAYS::Matrix \endlink
    /// \copydetails ARRAYS::Matrix

    /// \class zRCMatrix interface.h
    /// \brief complex floating point type row compressed matrix via \link ARRAYS::RowCompressedMatrix \endlink
    /// \copydetails ARRAYS::RowCompressedMatrix

    /// \class dEigenSystem interface.h
    /// \brief floating point eigensystem via \link ARRAYS::EigenSystem \endlink
    /// \copydetails ARRAYS::EigenSystem
    /// \class zEigenSystem interface.h
    /// \brief complex floating point eigensystem via \link ARRAYS::EigenSystem \endlink
    /// \copydetails ARRAYS::EigenSystem



    // FEM_DVR_ECS DOCUMENTATION
    /// \class femGrid interface.h
    /// \brief FEM-DVR-ECS grid on of floating point variable via \link FEM_DVR_ECS::FemDvrEcsGrid \endlink
    /// \class gVector interface.h
    /// \brief general complex function state vector on given coordinate grid via \link FEM_DVR_ECS::GridVector \endlink
    /// \class zOperatorF interface.h
    /// \brief general complex operator full representation of operator on given coordinate grid via\link FEM_DVR_ECS::OperatorFull \endlink
    /// \class zOperatorD interface.h
    /// \brief diagonal complex operator on given coordinate grid via \link FEM_DVR_ECS::OperatorDiagonal \endlink
    /// \class zOperatorC interface.h
    /// \brief general complex operator in row compressed format on given coordinate grid via \link FEM_DVR_ECS::OperatorRowCompressed \endlink
    /// \class DiscreteStates interface.h
    /// \brief Schrodinger equation discrete spectrum state evaluation for complex wave functions via \link FEM_DVR_ECS::DiscreteStates1D \endlink

    /// \class CrankNicolson1D interface.h
    /// \brief one-dimensional case of evolution operator approximation for complex wavefunctions via \link CrankNicolson \endlink
    /// \class Chebyshev1D interface.h
    /// \brief one-dimensional case of evolution operator approximation for complex wavefunctions via \link Chebyshev \endlink
    /// \class MoellerCN1D interface.h
    /// \brief one-dimensional case of Moeller operator approximation for comlex wavefunctions \link MoellerOperator \endlink

    /// \class femGrid2D interface.h
    /// \brief two-dimensional FEM-DVR-ECS grid on of floating point variable via \link FEM_DVR_ECS_2D::FemDvrEcsGrid2D \endlink
    /// \class gVector2D interface.h
    /// \brief general two-dimensional complex function state vector on given coordinate grid via  \link FEM_DVR_ECS_2D::GridVector2D \endlink
    /// \class zOperator2D interface.h
    /// \brief general complex operator on two-dimensional states in row compressed format on given coordinate grid via  \link FEM_DVR_ECS_2D::Operator2DRowCompressed \endlink

    /// \class Chebyshev2D interface.h
    /// \brief two-dimensional case of evolution operator approximation for complex wavefunctions via  \link Chebyshev \endlink
    /// \class CrankNicolson2D interface.h
    /// \brief two-dimensional case of evolution operator approximation for complex wavefunctions via  \link CrankNicolson \endlink

    // - MASK2D - row compressed version of FEM_DVR_ECS_2D
//    typedef FEM_DVR_ECS_2D::MaskGrid2D                 mGrid2D;
//    typedef FEM_DVR_ECS_2D::MaskVector2D               mVector2D;
//    typedef FEM_DVR_ECS_2D::MaskOperator2D             mOperator2D;

    /// \class mGrid2D interface.h
    /// \brief specialized case of two-dimensional coordinate discretization grid with element usage mask \link FEM_DVR_ECS_2D::MaskGrid2D \endlink
    /// \class mVector2D interface.h
    /// \brief specialized case of two-dimensional state vector associated to masked grid \link FEM_DVR_ECS_2D::MaskVector2D \endlink
    /// \class mOperator2D interface.h
    /// \brief specialized case of operator on two-dimensional state vectors on masked grid in row compressed format  \link FEM_DVR_ECS_2D::MaskOperator2D \endlink

//    typedef CrankNicolson<def_float, def_comp, mOperator2D, mVector2D>      mCrankNicolson2D;

    /// \class mCrankNicolson2D interface.h
    /// \brief two-dimensional case of evolution operator approximation for masked complex wavefunctions via \link CrankNicolson \endlink

    // coupled 2D
//    typedef FEM_DVR_ECS_2D::DoubleGridVector2D  doubleGVector2D;
//    typedef FEM_DVR_ECS_2D::DoubleOperator2DRC  doubleOperator2D;
    // shallow connection
//    typedef FEM_DVR_ECS_2D::ShallowGridVector2D sgVector2D;
//    typedef CrankNicolson<def_float, def_comp, doubleOperator2D, doubleGVector2D> doubleCrankNicolson2D;


    // \class  interface.h
    // \brief  \link  \endlink
    // \class  interface.h
    // \brief  \link  \endlink
    // \class  interface.h
    // \brief  \link  \endlink
    // \class  interface.h
    // \brief  \link  \endlink
    // \class  interface.h
    // \brief  \link  \endlink
    // \class  interface.h
    // \brief  \link  \endlink
    // \class  interface.h
    // \brief  \link  \endlink
    // \class  interface.h
    // \brief  \link  \endlink

    // functions
    #define zGaussian(x,xo,s,p) functions::Gaussian<def_comp>(x,xo,s,p)
    #define dSine(x,l,j) functions::Sine<def_float>(x,l,j)

    /// \fn zGaussian interface.h
    /// \brief  \link functions::Gaussian \endlink
    /// \fn dSine interface.h
    /// \brief  \link functions::Sine \endlink

    // potentials

    // Helper methods
    #define makeHamiltonianF(g, v, m) FEM_DVR_ECS::buildFullHamiltonian<def_float, def_comp>(g, v, m)

    /// \fn makeHamiltonianF interface.h
    /// \brief  \link FEM_DVR_ECS::buildFullHamiltonian \endlink


#endif
