#ifndef INCLUDE_FEM_DVR_ECS_2D_H_
#define INCLUDE_FEM_DVR_ECS_2D_H_
#include <string.h>
#include <vector>


/// Two dimensional state and operator representations associated to give discretisation.

/// Based on the one dimensional grids, the two dimentsional package is a straightforward
/// generalization of the one-dimensional representation.

#include "FemDvrEcs2d/FemDvrEcsGrid2d.h"
#include "FemDvrEcs2d/GridVector2d.h"
#include "FemDvrEcs2d/ShallowGridVector2d.h"
#include "FemDvrEcs2d/OperatorRowCompressed2d.h"
#include "FemDvrEcs2d/OperatorFull2d.h"
#include "FemDvrEcs2d/ZoomFilter.h"
#include "FemDvrEcs2d/EquidistantProjector2d.h"

//#include "FemDvrEcs2d/DoubleGridVector2d.h"
//#include "FemDvrEcs2d/DoubleOperator2dRC.h"

// EXPERIMENTAL SECTION - MASK VECTORS AND OPERATORS
//#include "FemDvrEcs2d/MaskGrid2d.h"
//#include "FemDvrEcs2d/MaskVector2d.h"
//#include "FemDvrEcs2d/MaskOperator2d.h"

namespace QSCAT
{
// EXPERIMENT Preconditioner
template<typename T, typename Z>
class PreconditionerOperator2d
{
    FemDvrEcsGrid2d grid_;
    Matrix<Z> x_kinetic_term_;
    Matrix<Z> x_kinetic_term_inversion_;
 public:
    PreconditionerOperator2d(const FemDvrEcsGrid2d& grid, T mu_x) : grid_(grid)
    {
        x_kinetic_term_ = generateKineticTerm<T,Z>(grid_.get_xgrid(), mu_x);
        x_kinetic_term_inversion_ = x_kinetic_term_;
        x_kinetic_term_inversion_.inverse();
    }
    void gemv(T alpha, GridVector2d& x, T beta, GridVector2d& y) const
    {
        blas_int shift;
        for (blas_int i=0; i<grid_.get_ysize(); ++i) {
            shift = i * grid_.get_xsize();
            ShallowVector<Z> p1(grid_.get_xsize(), &x[shift]);
            ShallowVector<Z> p2(grid_.get_xsize(), &y[shift]);
            x_kinetic_term_.gemv(alpha, p1, beta, p2);
        }
    }
    void igemv(T alpha, GridVector2d& x, T beta, GridVector2d& y) const
    {
        blas_int shift;
        for (blas_int i=0; i<grid_.get_ysize(); ++i){
            shift = i * grid_.get_xsize();
            ShallowVector<Z> p1(grid_.get_xsize(), &x[shift]);
            ShallowVector<Z> p2(grid_.get_xsize(), &y[shift]);
            x_kinetic_term_inversion_.gemv(alpha, p1, beta, p2);
        }
    }
    //void conjugate()
    //{
    //    x_kinetic_term_.conjugate();
    //    x_kinetic_term_inversion_.conjugate();
    //}
};

} // QSCAT

#endif // INCLUDE_FEM_DVR_ECS_2D_H_
