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

#include "common.h"
#include "Arrays.h"
#include "input.h"
#include "bessel.h"
#include "coulomb.h"

#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"
#include "FemDvrEcs/GridVector.h"

namespace QSCAT
{

  // internals

    bool GridVector::save_bin_body(std::ofstream & file) const
    {
        assert(init());
      //
        bool stat = file.is_open();
        if (stat) stat = grid_.save_binary(file);
        if (stat) stat = vector_.save_binary(file);
        return stat;
    }

    bool GridVector::read_bin_body(std::ifstream & file)
    {
        bool stat = file.is_open();
        FemDvrEcsGrid g;
        if (stat) stat = g.read_binary(file);
        if (!init() || !(g == grid_)) grid_ = g;
        if (stat) stat = vector_.read_binary(file);
        if (stat) *init_ = true;
        return stat;
    }

  // constructors

    GridVector::GridVector() : Object()
    {
        *init_ = false;
    }
    GridVector::GridVector(const FemDvrEcsGrid& grid): vector_(grid.nb()), grid_(grid), Object()
    {
        assert(grid.init());
      //
        *init_ = true;
    }
    GridVector::GridVector(const FemDvrEcsGrid& grid, const zVector& values, bool set_weights) : grid_(grid), Object()
    {
        assert(grid.init());
        assert(values.init());
        assert(grid.nb() == values.get_size());
      //
        *init_ = true;
        if (set_weights) {
            vector_ = zVector(grid.nb());
            for (blas_int i=0; i<grid.nb(); ++i){
                vector_[i] = values[i] * sqrt(w(i));                // TODO: add intel procedure
            }
        } else {
            vector_ = values.copy();
        }
    }
    GridVector::~GridVector()
    {
        decref();
    }
    GridVector::GridVector(const GridVector& old):
        Object(old),
        grid_(old.grid_),
        vector_(old.vector_)
    {}

    // accessors

    blas_int GridVector::get_size() const
    {
        return (init())? grid_.nb() : 0;
    }
    dcomp& GridVector::operator[] (blas_int i)
    {
        assert(init());
        assert(i < get_size());
      //
        return vector_[i];
    }
    const dcomp& GridVector::operator[] (blas_int i) const
    {
        assert(init());
        assert(i < get_size());
      //
        return vector_[i];
    }
    dcomp GridVector::f(blas_int i) const
    {
        assert(init());
        assert(i < get_size());
      //
        return vector_[i]/sqrt(w(i));
    }
    void GridVector::f(const dcomp& value, blas_int i)
    {
        assert(init());
        assert(i < get_size());
      //
        vector_[i] = value * sqrt(w(i));
    }
    zVector GridVector::function_values() const
    {
        assert(init());
      //
        zVector out(get_size());
        for (blas_int i=0; i<get_size(); ++i)
            out[i] = f(i);
        return out;
    }
    zVector& GridVector::body()
    {
        assert(init());
      //
        return vector_;
    }
    const zVector& GridVector::body() const
    {
        assert(init());
      //
        return vector_;
    }
    const FemDvrEcsGrid& GridVector::get_grid() const
    {
        return grid_;
    }
    const dcomp& GridVector::w(blas_int i) const
    {
        assert(init());
        assert(grid_.init());
        assert(i < grid_.nb());
      //
        return grid_.w(i);
    }
    const dfloat& GridVector::xr(blas_int i) const
    {
        assert(init());
        assert(grid_.init());
        assert(i < get_size());
      //
        return grid_.xr(i);
    }
    const dcomp& GridVector::xz(blas_int i) const
    {
        assert(init());
        assert(grid_.init());
        assert(i < get_size());
      //
        return grid_.x(i);
    }

 // modifiers

    GridVector& GridVector::fill(const dcomp& value)
    {
        assert(init());
      //
        vector_.fill(value);
        return *this;
    }

// FIXME REMOVE
    /*
    GridVector& GridVector::function_fill(const dcomp (*func)(const dcomp&, const parameters::model_2D<dfloat>&), const parameters::model_2D<dfloat>& mp) // Fills the vector with function values
    {
        assert(init());
        assert(func!=NULL);
      //
        for (blas_int i=0; i<get_size(); ++i){
            vector_[i] = func(xz(i), mp) * sqrt(w(i));
        }
        return *this;
    }
    GridVector& GridVector::function_xy_fill_x(const dcomp (*func)(const dcomp&, const dcomp&, const parameters::model_2D<dfloat>&), const dcomp& y, const parameters::model_2D<dfloat>& mp) // Fills the vector with function values of two dimensional function with one variable fixed (i.e. the second)
    {
        assert(init());
        assert(func!=NULL);
      //
        for (blas_int i=0; i<get_size(); ++i){
            vector_[i] = func(xz(i), y, mp) * sqrt(w(i));
        }
        return *this;
    }
    GridVector& GridVector::function_xy_fill_y(const dcomp (*func)(const dcomp&, const dcomp&,const parameters::model_2D<dfloat> &), const dcomp& x, const parameters::model_2D<dfloat> & mp) // Fills the vector with function values of two dimensional function with one variable fixed (i.e. the first)
    {
        assert(init());
        assert(func!=NULL);
      //
        for (blas_int i=0; i<get_size(); ++i){
            vector_[i] = func(x, xz(i), mp) * sqrt(w(i));
        }
        return *this;
    }
    */

    GridVector& GridVector::radial_state(const dcomp& energy, const dfloat& mass, int impulse_momentum)
    {
        assert(init());
      //
        for (blas_int i=0; i<get_size(); ++i){
            vector_[i] = sphBesselJEn(xz(i), sqrt(2.*mass*energy), mass, impulse_momentum) * sqrt(w(i));
        }
        return *this;
    }
    GridVector& GridVector::coulomb_state(const dcomp& energy, const dfloat& mass, int impulse_momentum, const dcomp& charge)
    {
        assert(init());
      //
        for (blas_int i=0; i<get_size(); ++i){
            vector_[i] = coulomb::sF_en(xz(i), sqrt(2.*mass*energy), charge, mass, impulse_momentum) * sqrt(w(i));
        }
        return *this;
    }
    dcomp GridVector::get_norm() const
    {
        assert(init());
      //
        return sqrt(vector_*vector_);
    }
    dcomp GridVector::derivative(blas_int pos) const
    {
        assert(init());
        assert(pos < get_size());
      //
        blas_int i, element, start, end;
        blas_int quadrature = grid_.quadrature();
        dcomp length, out;
    // First determining the element index and the position within the element
        if ( (pos + 1) % (quadrature-1) != 0){
            i = (pos + 1) % (quadrature-1);
            element = (pos + 1 - i)/(quadrature-1);
        } else {
            i = quadrature-1;
            element = (pos+1)/(quadrature-1);
        }

        length = grid_.ar(element) - grid_.ar(element-1); // !
        start = 0;
        end =  quadrature;
        if (element==1) { start++; }
        if (element==grid_.tnel()-1) { end--; }

        out = 0.0;

        for (blas_int k=start; k<end; ++k){
            out += grid_.dlp(k,i) * vector_[(element-1)*(quadrature-1) + k - 1] / sqrt(w((element-1)*(quadrature-1) + k - 1));
        }
        return out/length*2.0;
    }
 // FIXME evaluation on endpoints
    dcomp GridVector::evaluate(const dfloat& X) const
    {
        assert(init());
        assert(X > xr(0));
        assert(X < xr(get_size()-1));
      //
        blas_int el_x = grid_.get_element_index(X);
        blas_int el_xs = grid_.get_element_start(el_x);
        blas_int el_xe = grid_.get_element_end(el_x);
        zVector hlp_x(el_xe - el_xs + 1);
        for (blas_int i=el_xs; i<=el_xe; ++i){
            hlp_x[i - el_xs] = grid_.basis_function_value(i,X,el_xs,el_xe);
        }
        dcomp fx = 0.0;
        for (blas_int j=el_xs; j<=el_xe; ++j){
            fx += vector_[j] * hlp_x[j-el_xs];
        }
        return fx;
    }
    GridVector& GridVector::swap(GridVector& rhs)
    {
        Object::swap(rhs);
        vector_.swap(rhs.vector_);
        grid_.swap(rhs.grid_);
        return *this;
    }
    GridVector GridVector::copy() const
    {
        GridVector out;
        *out.init_ = *init_;
        out.grid_ = grid_.copy();
        out.vector_ = vector_.copy();
        return out;
    }
    GridVector&  GridVector::complex_conjugate()
    {
        assert(init());
      //
        vector_.complex_conjugate();
        return *this;
    }

// operators

    dcomp GridVector::operator* (const GridVector& rhs) const
    {
        assert(init());
        assert(rhs.init());
        assert(get_size() == rhs.get_size());
      //
        return vector_ * rhs.vector_;
    }
    ScalarMultiple<dcomp, Vector<dcomp> > GridVector::operator* (const dcomp& scalar)
    {
        assert(init());
      //
        return vector_ * scalar;
    }
    GridVector& GridVector::operator= (GridVector tmp)
    {
        swap(tmp);
        return *this;
    }
    GridVector& GridVector::operator= (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs)
    {
        vector_ = rhs;
        return *this;
    }
    GridVector& GridVector::operator+= (const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(get_size() == rhs.get_size());
      //
        vector_ += rhs.vector_;
        return *this;
    }
    GridVector& GridVector::operator+= (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs)
    {
        assert(init());
        assert(rhs.object().init());
        assert(get_size() == rhs.object().get_size());
      //
        vector_ += rhs;
        return *this;
    }
    GridVector& GridVector::operator-= (const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(get_size() == rhs.get_size());
      //
        vector_ -= rhs.vector_;
        return *this;
    }
    GridVector& GridVector::operator-= (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs)
    {
        assert(init());
        assert(rhs.object().init());
        assert(get_size() == rhs.object().get_size());
      //
        vector_ -= rhs;
        return *this;
    }
    GridVector& GridVector::operator*= (const dcomp& alpha)
    {
        assert(init());
      //
        vector_ *= alpha;
        return *this;
    }

 // custom

    GridVector& GridVector::axpy(const dcomp& alpha, const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(get_size() == rhs.get_size());
      //
        vector_.axpy(alpha,rhs.vector_);
        return *this;
    }
    GridVector & GridVector::ax(const dcomp& alpha, const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(get_size() == rhs.get_size());
      //
        vector_.ax(alpha,rhs.vector_);
        return *this;
    }
    dcomp GridVector::reduction(const GridVector& y) const
    {
        assert(init());
        assert(y.init());
        assert(get_size() == y.get_size());
      //
        return vector_.reduction(y.vector_);
    }
    GridVector & GridVector::element_wise_multiplication(const GridVector& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(get_size() == rhs.get_size());
      //
        vector_.element_wise_multiplication(rhs.vector_);
        return *this;
    }

 // Storage

    bool SaveMultiGridVectorBin(const char* name, GridVector* V, const blas_int& N)
    {
        bool stat, ostat;
        std::ofstream file;
        file.open(name, std::ios::out | std::ios::binary);
        ostat = file.is_open();
        stat = ostat;
        if (stat) file.write((char*) &N, sizeof(blas_int));
        if (stat) stat = (V[0].get_grid()).save_binary(file);
        for (blas_int i=0; i<N; ++i){
            if (stat) stat = V[i].body().save_binary(file);
        }
        if (ostat) file.close();
        return stat;
    }
    bool ReadMultiGridVectorBin(const char* name, GridVector* V, const FemDvrEcsGrid& g)
    {
        bool stat, ostat;
        blas_int N;
        std::ifstream file;
        file.open(name, std::ios::in | std::ios::binary);
        ostat = file.is_open();
        stat = ostat;
        if (stat) file.read((char*) &N, sizeof(blas_int));
        FemDvrEcsGrid G;
        if (stat) stat = G.read_binary(file);
        if (stat) stat = (G==g);
        for (blas_int i=0; i<N; ++i){
            V[i] = GridVector(g);
            if (stat) stat = V[i].body().read_binary(file);
        }
        if (ostat) file.close();
        return stat;
    }

    void GridVector::save(const char*name) const
    {
        assert(init());
      //
        blas_int size = get_size();
        FILE * file;
        def_comp v;
        fopen_s(&file,name,"w");
        fprintf(file,"#Function on FEM-DVR-ECS grid of %lld basis functions.\n", size);
        fprintf(file,"#Coordinate X   \tReal part of z     \tImaginary part of z    \tAbsolute value z squared\n");
        for (blas_int i=0; i<size;++i){
            v = f(i);
            fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\n", xr(i), real(v), imag(v), std::pow(std::abs(v),2)); // 
        }
        fclose(file);
    }
    void GridVector::save_equidistant(const char*name, const def_float a, const def_float b, const blas_int samples) const
    {
        assert(init());
      //
        FILE * file;
        Vector<def_comp> V(samples);
        V.fill(0.0);
        Vector<def_float> X(samples, a, b, true);

        blas_int start = (a == xr(0))? 1:0;
        blas_int end = ( b == xr(grid_.nb()-1))? samples-1:samples;
        for (blas_int i=start; i<end; ++i){
            V[i] = evaluate( X[i] );
        }

        fopen_s(&file,name,"w");
        fprintf(file,"#Function on FEM-DVR-ECS grid of %lld basis functions.\n", grid_.nb());
        fprintf(file,"#Coordinate X   \tReal part of z     \tImaginary part of z    \tAbsolute value z squared\n");
        for (blas_int i=0;i<samples;++i){
            fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\n", X[i], real(V[i]), imag(V[i]), std::pow(std::abs(V[i]),2)); // 
        }
        fclose(file);
    }
} // namespace QSCAT
