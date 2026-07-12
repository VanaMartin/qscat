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
#include "FemDvrEcs2d.h"

namespace QSCAT
{
  // Constructors
    DoubleOperator2dRC::DoubleOperator2dRC() :
        grid1_(),
        grid2_(),
        body_()
    {
        init_ = false;
    }
    DoubleOperator2dRC::DoubleOperator2dRC(FemDvrEcsGrid2d& grid) :
        grid1_(grid),
        grid2_(grid),
        body_()
    {
        init_ = true;
    }
    DoubleOperator2dRC::DoubleOperator2dRC(FemDvrEcsGrid2d& grid1, FemDvrEcsGrid2d& grid2) :
        grid1_(grid1),
        grid2_(grid2),
        body_()
    {
        init_ = true;
    }
    DoubleOperator2dRC::DoubleOperator2dRC(const DoubleOperator2dRC& old) :
        grid1_(old.grid1_),
        grid2_(old.grid2_),
        body_(old.body_),
        init_(old.init_)
    {}

  // accessors
    bool DoubleOperator2dRC::init() const
    {
        return init_;
    }
    const FemDvrEcsGrid2d& DoubleOperator2dRC::get_grid1() const
    {
        assert(init_);
      //
        return grid1_;
    }
    const FemDvrEcsGrid2d& DoubleOperator2dRC::get_grid2() const
    {
        assert(init_);
      //
        return grid2_;
    }
    blas_int DoubleOperator2dRC::get_size1() const
    {
        assert(init_);
      //
        return grid1_.get_size();
    }
    blas_int DoubleOperator2dRC::get_size2() const
    {
        assert(init_);
      //
        return grid2_.get_size();
    }
    blas_int DoubleOperator2dRC::get_size() const
    {
        assert(init_);
      //
        return get_size1() + get_size2();
    }

  // modifiers
    DoubleOperator2dRC& DoubleOperator2dRC::add_kinetic_term(const dfloat& mu1, const dfloat& mu2, const dfloat& mu3, const dfloat& mu4)
    {
        assert(init_);
      //
        blas_int M = get_size();

        RowCompressedMatrix<dcomp> ke1 = generateKineticTermRCM(grid1_.get_xgrid(), mu1);
        RowCompressedMatrix<dcomp> ke2 = generateKineticTermRCM(grid1_.get_ygrid(), mu2);
        RowCompressedMatrix<dcomp> KE1 = ARRAYS::TensorSum(ke1, ke2);

        ke1 = generateKineticTermRCM(grid2_.get_xgrid(), mu3);
        ke2 = generateKineticTermRCM(grid2_.get_ygrid(), mu4);
        RowCompressedMatrix<dcomp> KE2 = ARRAYS::TensorSum(ke1, ke2);

        KE1.expand(M, M, 0, 0);
        KE2.expand(M, M, get_size1(), get_size1());

        KE1.axpy(dcomp(1.0), KE2);

        if (body_.init()) {
            body_.axpy(dcomp(1.0), KE1);
        } else {
            body_.swap(KE1);
        }
        return *this;
    }
    DoubleOperator2dRC& DoubleOperator2dRC::add_kinetic_term(const dfloat& mu1, const dfloat& mu2)
    {
        assert(init_);
      //
        return add_kinetic_term(mu1, mu2, mu1, mu2);
    }
    DoubleOperator2dRC& DoubleOperator2dRC::swap(DoubleOperator2dRC& rhs)
    {
        grid1_.swap(rhs.grid1_);
        grid2_.swap(rhs.grid2_);
        body_.swap(rhs.body_);
        std::swap(init_, rhs.init_);
        return *this;
    }
    void DoubleOperator2dRC::LU_factorize()
    {
        assert(init_);
        assert(body_.init());
      //
        body_.LU_factorize();
    }

  // operators
    DoubleOperator2dRC& DoubleOperator2dRC::operator*= (const dcomp& alpha)
    {
        assert(init_);
      //
        body_*=alpha;
        return *this;
    }
    DoubleOperator2dRC& DoubleOperator2dRC::operator+= (const dcomp& alpha)
    {
        assert(init_);
      //
        body_.add_to_diagonal(alpha);
        return *this;
    }
    DoubleOperator2dRC& DoubleOperator2dRC::operator= (DoubleOperator2dRC rhs)
    {
        return this->swap(rhs);
    }
  // custom operations
    void DoubleOperator2dRC::gemv(const dcomp& alpha, const DoubleGridVector2d& x, const dcomp& beta, DoubleGridVector2d& y) const
    {
        assert(init_);
        assert(body_.init());
        assert(x.init());
        assert(y.init());
        assert(grid1_ == x.get_grid1());
        assert(grid2_ == x.get_grid2());
        assert(grid1_ == y.get_grid1());
        assert(grid2_ == y.get_grid2());
      //
        body_.gemv(alpha, x.body(), beta, y.body());
    }
    void DoubleOperator2dRC::LU_back_substitution(DoubleGridVector2d&x)
    {
        assert(init_);
        assert(body_.init());
        assert(grid1_ == x.get_grid1());
        assert(grid2_ == x.get_grid2());
      //
        body_.LU_back_substitution(x.body());
    }
    DoubleOperator2dRC& DoubleOperator2dRC::add_potential(const GridVector2d& p1, const GridVector2d& p2)
    {
        // TODO : separate cases with AddVectorToDiagonal
        assert(init_);
        assert(grid1_ == p1.get_grid());
        assert(grid2_ == p2.get_grid());
        assert(body_.init());
        // Build diagonal representation of the vectors
        zVector vals(get_size());       // vector of non zero values, size of the diagonal representation
        iVector cols(get_size());             // vector of columns indice, same size
        iVector inds(get_size()+1);           // vector of row indices, same size + 1
        for (blas_int j=0; j<p1.get_size(); ++j) {       // First grid diagonal
            vals[j] = p1.f(j);                      // Function value
            cols[j] = j;                            // column
            inds[j] = j;                            // row index value
        }
        for (blas_int j=get_size1(); j<get_size(); ++j) {   // Second grid diagonal
            vals[j] = p2.f(j-get_size1());              // Function value
            cols[j] = j;                            // column
            inds[j] = j;                            // row index
        }
        inds[get_size()] = get_size();
        body_.axpy(1.0, blas_int(get_size()), &vals[0], &cols[0], &inds[0]);
        return *this;
    }
    DoubleOperator2dRC& DoubleOperator2dRC::add_coupling(const GridVector2d &c1, const GridVector2d &c2)
    {
        // TODO: add asymmetrical cases
        assert(init_);
        assert(body_.init());
        assert(get_size1() == get_size2());   // Symmetrical case only
        assert(get_size1() == c1.get_size());
        assert(get_size2() == c2.get_size());
      //
        blas_int nb = get_size1();
        zVector vals(2*nb);
        iVector cols(2*nb);
        iVector inds(2*nb + 1);
        for (blas_int j=0; j<nb; ++j) {
            vals[j] = c1.f(j);    // Function value
            cols[j] = j + nb;
            inds[j] = j;
        }
        for (blas_int j=0; j<nb; ++j) {
            vals[j+nb] = c2.f(j);    // Function value
            cols[j+nb] = j;
            inds[j+nb] = j+nb;
        }
        body_.axpy(1.0, get_size(), &vals[0], &cols[0], &inds[0]);
        return *this;
    }
}   // namespace QSCAT
