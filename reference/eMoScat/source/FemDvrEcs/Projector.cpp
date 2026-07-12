#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <complex>
#include <cassert>
#include <math.h>
#include <stdlib.h>

#ifdef linux
    #include <sec_stream.h>
#endif

#include "bessel.h"
#include "coulomb.h"
#include "common.h"
#include "Arrays.h"
#include "input.h"

#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"
#include "FemDvrEcs/Projector.h"

namespace QSCAT
{
    void Projector::initialize(GridVector& psi)
    {
        alpha_ = dcomp(1);
        beta_ = dcomp(0);
        size_ = psi.get_size();
        body_ = zMatrix(size_,size_);
        for (blas_int i=0; i<size_; ++i){
            for (blas_int j=0; j<size_; ++j){
                body_(i,j) = conj(psi[i])*psi[j];
            }
        }
        *init_ = true;
    }
    void Projector::initialize(GridVector& psi, const dcomp& alpha, const dcomp& beta)
    {
        alpha_ = alpha;
        beta_ = beta;
        size_ = psi.get_size();
        body_ = zMatrix(size_, size_);
        for (blas_int i=0; i<size_; ++i){
            for (blas_int j=0; j<size_; ++j){
                body_(i,j) = alpha_ * conj(psi[i]) * psi[j];
            }
            body_(i,i) += beta_;
        }
        *init_ = true;
    }
    Projector::Projector() : Object()
    {
        *init_ = false;
    }
    Projector::Projector(GridVector& psi) : Object()
    {
        initialize(psi);
    }
    Projector::Projector(GridVector& psi, const dcomp& alpha, const dcomp& beta) : Object()
    {
        initialize(psi,alpha,beta);
    }
    Projector::~Projector()
    {
        decref();
    }
    GridVector Projector::operator* (const GridVector& rhs) const
    {
        assert(init());
        assert(rhs.init());
        assert(rhs.get_size() == size_);
      //
        GridVector out = rhs.copy();
        out.body() =  body_ * rhs.body();
        return out;
    }
}   // namespace QSCAT
