#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <complex>
#include <string>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "bessel.h"
#include "common.h"
#include "blas.h"
#include "Arrays.h"
#include "input.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d/FemDvrEcsGrid2d.h"
#include "FemDvrEcs2d/GridVector2d.h"
#include "FemDvrEcs2d/ShallowGridVector2d.h"
#include "FemDvrEcs2d/OperatorFull2d.h"

namespace QSCAT
{
  // constuctors
    Operator2dFull::Operator2dFull()
    {
        *init_ = false;
    }
    Operator2dFull::Operator2dFull(FemDvrEcsGrid2d& grid)
    {
        assert(grid.init());
      //
        grid_ = grid;
        *init_ = true;
        incref();
    }

    Operator2dFull::Operator2dFull(const Operator2dFull& old) :
        Object(old),
        grid_(old.grid_),
        body_(old.body_)
    {}

  // asccessors
    const FemDvrEcsGrid2d& Operator2dFull::get_grid() const
    {
        assert(init());
      //
        return grid_;
    }

  // modifiers
    Operator2dFull& Operator2dFull::set_kinetic_term(const dfloat& mux, const dfloat& muy)
    {
        assert(init());
      //
        zMatrix ke_x = generateKineticTerm(grid_.get_xgrid(), mux);
        zMatrix ke_y = generateKineticTerm(grid_.get_ygrid(), muy);
        // FIXME FIXME FIXME
        assert(0);
        // body_ = TensorSum(ke_x, ke_y);
        return *this;
    }
    void Operator2dFull::LU_factorize()
    {
        assert(init());
        assert(body_.init());
      //
        body_.LU_factorize();
    }

  // operators
    Operator2dFull& Operator2dFull::operator*= (const dcomp& alpha)
    {
        assert(init());
      //
        body_*=alpha;
        return *this;
    }
    Operator2dFull& Operator2dFull::operator+= (const dcomp& alpha)
    {
        assert(init());
        assert(body_.init());
      //
        body_.add_to_diagonal(alpha);
        return *this;
    }
    Operator2dFull& Operator2dFull::operator+= (const GridVector2d& rhs)
    {
        assert(init());
        assert(body_.init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        assert(0); // NIY FIXME
        //body_.add_grid_vector_to_diagonal(rhs);
        return *this;
    }
    Operator2dFull& Operator2dFull::operator+= (const zVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_.get_size() == rhs.get_size());
      //
        body_.add_vector_to_diagonal(rhs);
        return *this;
    }
    GridVector2d Operator2dFull::operator* (const GridVector2d& x) const
    {
        assert(init());
        assert(x.init());
        assert(grid_ == x.get_grid());
      //
        GridVector2d y(x);
        gemv(1.0, x, 0.0, y);
        return y;
    }

  // custom operations
    void Operator2dFull::gemv(const dcomp& alpha, const GridVector2d& x, const dcomp& beta, GridVector2d& y) const
    {
        assert(init());
        assert(x.init());
        assert(y.init());
        assert(grid_ == x.get_grid());
        assert(grid_ == y.get_grid());
      //
        body_.gemv(alpha, x.body(), beta, y.body());
    }
    void Operator2dFull::LU_back_substitution(GridVector2d& x)
    {
        assert(init());
        assert(x.init());
        assert(grid_ == x.get_grid());
      //
        body_.LU_back_substitution(x.body());
    }
}   // namespace QSCAT
