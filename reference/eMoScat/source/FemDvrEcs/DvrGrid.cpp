#include <cassert>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "Arrays.h"
#include "input.h"

#include "FemDvrEcs/FemDvrFunctions.h"
#include "FemDvrEcs/DvrGrid.h"

namespace QSCAT
{

    bool DvrGrid::save_bin_body(std::ofstream & file) const
    {
        assert(init());
      //
        if (file.is_open()){
            file.write((char*) &quadrature_, sizeof(blas_int));
            file.write((char*) &x_min_, sizeof(dfloat));
            file.write((char*) &x_max_, sizeof(dfloat));
            if(!x_.save_binary(file)) goto save_break;
            if(!weights_.save_binary(file)) goto save_break;
            return true;
        }
    save_break:
        return false;
    }
    bool DvrGrid::read_bin_body(std::ifstream & file)
    {
        if(file.is_open()){
            file.read((char*) &quadrature_, sizeof(blas_int));
            file.read((char*) &x_min_, sizeof(dfloat));
            file.read((char*) &x_max_, sizeof(dfloat));
            if(!x_.read_binary(file)) goto read_break;
            if(!weights_.read_binary(file)) goto read_break;
            *init_ = true;
            return true;
        }
    read_break:
        clean();
        return false;
    }
    void DvrGrid::clean()
    {
        x_min_ = 0;
        x_max_ = 0;
        quadrature_ = 0;
        x_ = dVector();
        weights_ = dVector();
        *init_ = false;
    }
    DvrGrid::DvrGrid() : Object()            // Constructor of uninitialized DvrGrid
    {
        *init_ = false;
        x_min_ = 0;
        x_max_ = 0;
        quadrature_ = 0;
    }
    DvrGrid::DvrGrid(blas_int quadrature) : Object()           // Constructor of DVR basis grid and weights on range -1,1
    {
        *init_ = true;
        quadrature_ = quadrature;
        x_min_ = dfloat(-1);
        x_max_ = dfloat(1);
        x_ = dVector(quadrature);
        weights_ = dVector(quadrature);
        GLo_Quad(quadrature, x_, weights_);
    }
    DvrGrid::DvrGrid(dfloat start, dfloat end, blas_int quadrature) : Object()     // Constructor: creates DVR basis grid and weights on given range.
    {
        assert(quadrature>1);
        assert(end - start > 0);
      //
        *init_ = true;
        quadrature_ = quadrature;
        x_min_ = start;
        x_max_ = end;
        x_ = dVector(quadrature);
        weights_ = dVector(quadrature);
        GLo_Quad(quadrature_, x_, weights_);
        for (blas_int i=0; i<quadrature; i++) {
            x_[i] = 0.5 * (x_min_ + x_max_ + (x_max_ - x_min_) * x_[i]);
            weights_[i] = 0.5 * (x_max_ - x_min_) * weights_[i];
        }
    }
    DvrGrid::DvrGrid(const DvrGrid& old):
        Object(old),
        quadrature_(old.quadrature_),
        x_min_(old.x_min_),
        x_max_(old.x_max_),
        x_(old.x_),
        weights_(old.weights_)
    {}
    DvrGrid & DvrGrid::operator= (DvrGrid tmp)
    {
        this->swap(tmp);
        return *this;
    }
    DvrGrid& DvrGrid::swap (DvrGrid& rhs)
    {
        Object::swap(rhs);
        std::swap(x_min_, rhs.x_min_);
        std::swap(x_max_, rhs.x_max_);
        std::swap(quadrature_, rhs.quadrature_);
        x_.swap(rhs.x_);
        weights_.swap(rhs.weights_);
        return *this;
    }
    DvrGrid DvrGrid::copy () const
    {
        assert(init());
      //
        DvrGrid out = DvrGrid();
        *out.init_ = *init_;
        out.quadrature_ = quadrature_;
        out.x_min_ = x_min_;
        out.x_max_ = x_max_;
        out.x_ = x_.copy();
        out.weights_ = weights_.copy();
        return out;
    }
    DvrGrid::~DvrGrid()            // Deconstructor
    {
        if(decref()==0)
            clean();
    }
    dfloat DvrGrid::x(blas_int i) const
    {
        assert(i<=quadrature_);
      //
        return x_[i];
    }
    dfloat DvrGrid::w(blas_int i) const
    {
        assert(i<=quadrature_);
      //
        return weights_[i];
    }
    blas_int DvrGrid::quadrature() const
    {
        return quadrature_;
    }

    void DLagPol(int nq, dMatrix& dLp, DvrGrid& g)
    {
        dfloat hlp;
        for (int i=0; i<nq; i++) {
        // Diagonal terms
            dLp(i,i) = 0.0;
            for (int k=0; k<nq; k++) {
                if (i != k) {
                    hlp = dLp(i,i) + 1.0/(g.x(i) - g.x(k));
                    dLp(i,i) = hlp;
                }
            }
        // Non-diagonal terms: dLP(j,i) = - 1 / (dLP(i,j) * (x1(i) - x1(j))^2)
            for (int j=i+1; j<nq; j++) {
                hlp = 1.0;
                for (int k=0;k<nq;k++) {
                    if (k!=i && k!=j) {
                        hlp *= (g.x(j) - g.x(k))/(g.x(i) - g.x(k));
                    }
                }
                dLp(i,j) = hlp/(g.x(i)- g.x(j));
                dLp(j,i) = 1.0/(hlp*(g.x(j)- g.x(i)));
            }
        }
    }
} // namespace QSCAT
