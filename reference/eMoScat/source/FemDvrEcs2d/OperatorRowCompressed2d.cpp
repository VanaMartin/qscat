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
#include "FemDvrEcs2d/OperatorRowCompressed2d.h"

namespace QSCAT
{
  // constuctors
    Operator2dRowCompressed::Operator2dRowCompressed() : Object()
    {
        *init_ = false;
    }
    Operator2dRowCompressed::Operator2dRowCompressed(FemDvrEcsGrid2d& grid) : Object(), grid_(grid)
    {
        assert(grid.init());
      //
        *init_ = true;
    }

    Operator2dRowCompressed::Operator2dRowCompressed(const Operator2dRowCompressed& old) :
        Object(old),
        grid_(old.grid_),
        body_(old.body_)
    {}
    Operator2dRowCompressed::~Operator2dRowCompressed()
    {
        decref();
    }

    Operator2dRowCompressed Operator2dRowCompressed::copy() const
    {
        Operator2dRowCompressed out;
        out.grid_ = grid_;
        out.body_ = body_.copy();
        *out.init_ = *init_;
        return out;
    }
    Operator2dRowCompressed& Operator2dRowCompressed::swap(Operator2dRowCompressed& rhs)
    {
        Object::swap(rhs);
        grid_.swap(rhs.grid_);
        body_.swap(rhs.body_);
        return *this;
    }
    Operator2dRowCompressed& Operator2dRowCompressed::operator= (Operator2dRowCompressed tmp)
    {
        swap(tmp);
        return *this;
    }


  // asccessors
    const FemDvrEcsGrid2d& Operator2dRowCompressed::get_grid() const
    {
        assert(init());
      //
        return grid_;
    }
    dcomp Operator2dRowCompressed::operator() (int i, int j) const
    {
        assert(init());
        assert(i < grid_.get_size());
        assert(j < grid_.get_size());
      //
        return body_.get_element(i,j);
    }
    const RowCompressedMatrix<dcomp>& Operator2dRowCompressed::body() const
    {
        return body_;
    }
  // modifiers
    Operator2dRowCompressed& Operator2dRowCompressed::set_kinetic_term(const dfloat& mux, const dfloat& muy)
    {
        assert(init());
      //
        RowCompressedMatrix<dcomp> ke_x = generateKineticTermRCM(grid_.get_xgrid(), mux);
        RowCompressedMatrix<dcomp> ke_y = generateKineticTermRCM(grid_.get_ygrid(), muy);
        body_ = TensorSum(ke_x, ke_y);
        return *this;
    }
    void Operator2dRowCompressed::LU_factorize()
    {
        assert(init());
        assert(body_.init());
      //
        body_.LU_factorize();
    }
    Operator2dRowCompressed& Operator2dRowCompressed::conjugate()
    {
        assert(init());
      //
        body_.conjugate();
        return *this;
    }

  // operators
    Operator2dRowCompressed& Operator2dRowCompressed::operator*= (const dcomp& alpha)
    {
        assert(init());
      //
        body_*=alpha;
        return *this;
    }
    Operator2dRowCompressed& Operator2dRowCompressed::operator+= (const dcomp& alpha)
    {
        assert(init());
        assert(body_.init());
      //
        body_.add_to_diagonal(alpha);
        return *this;
    }
    Operator2dRowCompressed& Operator2dRowCompressed::operator+= (const GridVector2d& rhs)
    {
        assert(init());
        assert(body_.init());
        assert(rhs.init());
        assert(grid_ == rhs.get_grid());
      //
        body_.add_vector_to_diagonal(rhs.function_values());
        return *this;
    }
    Operator2dRowCompressed& Operator2dRowCompressed::operator+= (const zVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_.get_size() == rhs.get_size());
      //
        body_.add_vector_to_diagonal(rhs);
        return *this;
    }
    GridVector2d Operator2dRowCompressed::operator* (const GridVector2d& x) const
    {
        assert(init());
        assert(x.init());
        assert(grid_ == x.get_grid());
      //
        GridVector2d y = x.copy();
        gemv(1.0, x, 0.0, y);
        return y;
    }

  // custom operations
    void Operator2dRowCompressed::gemv(const dcomp& alpha, const GridVector2d& x, const dcomp& beta, GridVector2d& y) const
    {
        assert(init());
        assert(x.init());
        assert(y.init());
        assert(grid_ == x.get_grid());
        assert(grid_ == y.get_grid());
      //
        body_.gemv(alpha, x.body(), beta, y.body());
    }
    void Operator2dRowCompressed::LU_back_substitution(GridVector2d& x)
    {
        assert(init());
        assert(x.init());
        assert(grid_ == x.get_grid());
      //
        body_.LU_back_substitution(x.body());
    }
}   // namespace QSCAT
