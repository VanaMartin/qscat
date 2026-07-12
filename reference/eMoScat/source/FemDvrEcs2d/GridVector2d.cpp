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

namespace QSCAT
{
  //  The two  dimensional function on the  grid represented as a grid vector
  //  type.  The ordering  of  the values in  the  internal  arrays  goes as:
  //  f(x[i],y[j]) = a[nb_x*j + i]

  // internals
    void GridVector2d::clean()
    {
        *init_ = false;
    }
    bool GridVector2d::save_bin_body(std::ofstream& file) const
    {
        assert(init());
        bool stat = file.is_open();
        if (stat) stat = grid_.save_binary(file);
        if (stat) stat = body_.save_binary(file);
        return stat;
    }
    bool GridVector2d::read_bin_body(std::ifstream& file)
    {
        bool stat = file.is_open();
        if (stat) stat = grid_.read_binary(file);
        if (stat) stat = body_.read_binary(file);
        *init_ = true;
        return stat;
    }

  // constructors
    GridVector2d::GridVector2d() : Object()
    {
        clean();
        incref();
        *init_ = false;
    }
    GridVector2d::GridVector2d(const FemDvrEcsGrid2d& grid) : Object(), grid_(grid)
    {
        assert(grid.init());
      //
        *init_ = true;
        body_ = zVector(grid_.get_size());
        incref();
    }
    GridVector2d::GridVector2d(const FemDvrEcsGrid2d& grid, const zVector& unweighted_values) : Object(), grid_(grid)
    {
        assert(grid.init());
        assert(unweighted_values.init());
        assert(unweighted_values.get_size() == grid.get_size());
      //
        body_ = zVector(grid_.get_size());

        blas_int index;
        for (blas_int j=0;j<grid_.get_ysize();++j){
            for (blas_int i=0;i<grid_.get_xsize();++i){
                index = j*grid_.get_xsize() + i;
                body_[index] = unweighted_values[index] * sqrt(grid_.wz(i,j));
            }
        }
        *init_ = true;
        incref();
    }
    GridVector2d::GridVector2d(const FemDvrEcsGrid2d& grid, const GridVector& x_values, const GridVector& y_values) : Object(),  grid_(grid)
    {
        assert(grid.init());
        assert(grid.get_xgrid() == x_values.get_grid());
        assert(grid.get_ygrid() == y_values.get_grid());
      //
        body_ = zVector(grid_.get_size());
        for (blas_int j=0; j<grid_.get_ysize(); ++j){
            for (blas_int i=0; i<grid_.get_xsize(); ++i){
                body_[j*grid_.get_xsize() + i] = x_values[i] * y_values[j];
            }
        }
        *init_ = true;
        incref();
    }
    GridVector2d::GridVector2d(const GridVector2d& old) :
        Object(old),
        grid_(old.grid_),
        body_(old.body_)
    {}
    GridVector2d::~GridVector2d()
    {
        decref();
    }

  // accessors
    blas_int GridVector2d::get_real_xsize() const
    {
        assert(init());
      //
        return grid_.get_real_xsize();
    }
    blas_int GridVector2d::get_real_ysize() const
    {
        assert(init());
      //
        return grid_.get_real_ysize();
    }
    blas_int GridVector2d::get_xsize() const
    {
        assert(init());
      //
        return grid_.get_xsize();
    }
    blas_int GridVector2d::get_ysize() const
    {
        assert(init());
      //
        return grid_.get_ysize();
    }
    blas_int GridVector2d::get_size() const
    {
        assert(init());
      //
        return grid_.get_size();
    }
    const FemDvrEcsGrid2d& GridVector2d::get_grid() const
    {
        assert(init());
      //
        return grid_;
    }
    const dcomp& GridVector2d::operator[] (blas_int i) const
    {
        assert(init());
        assert(i < get_size());
      //
        return body_[i];
    }
    dcomp& GridVector2d::operator[] (blas_int i)
    {
        assert(init());
        assert(i < get_size());
      //
        return body_[i];
    }
    dcomp GridVector2d::f(blas_int i) const
    {
        assert(init());
        assert(i < get_size());
      //
        return body_[i]/sqrt(grid_.wz(i));
    }
    void GridVector2d::f(const dcomp& val,  blas_int i)
    {
        assert(init());
        assert(i < get_size());
      //
        body_[i] = val * std::sqrt(grid_.wz(i));
    }
    /// \brief Two dimensional accessor
    /// \detail The order of indeces respects tho column wise ordering of storage
    /// NOTE: the order of indeces is opposite to the order of grid creation (FIXME?)
    dcomp GridVector2d::f(blas_int iy, blas_int ix) const
    {
        assert(init());
        assert(ix < get_xsize());
        assert(iy < get_ysize());
      //
        return f(iy * get_xsize() + ix);
    }
    /// \brief Two dimensional accessor
    /// \detail The order of indeces respects tho column wise ordering of storage
    /// NOTE: the order of indeces is opposite to the order of grid creation (FIXME?)
    void GridVector2d::f(const dcomp& val, blas_int iy, blas_int ix)
    {
        assert(init());
        assert(ix < get_xsize());
        assert(iy < get_ysize());
      //
        f(val, iy * get_xsize() + ix);
    }
    zVector GridVector2d::function_values() const
    {
        assert(init());
      //
        zVector out(get_size());
        for (blas_int i=0; i<get_size(); ++i)
            out[i] = f(i);
        return out;
    }
    dfloat GridVector2d::norm() const
    {
        assert(init());
      //
        return sqrt( abs(body_*body_) );
    }
    dcomp GridVector2d::evaluate(const dfloat& x, const dfloat& y) const
    {
        const FemDvrEcsGrid *p;   // Auxiliary pointer;
        p = &grid_.get_xgrid();
        blas_int el_x = p->get_element_index(x);
        blas_int el_xs = p->get_element_start(el_x);
        blas_int el_xe = p->get_element_end(el_x);
        p = &grid_.get_ygrid();
        blas_int el_y = p->get_element_index(y);
        blas_int el_ys = p->get_element_start(el_y);
        blas_int el_ye = p->get_element_end(el_y);

        zVector hlp_x(el_xe - el_xs + 1);
        zVector hlp_y(el_ye - el_ys + 1);
        p = &grid_.get_xgrid();
        for (blas_int i=el_xs; i<=el_xe; ++i){
            hlp_x[i - el_xs] = p->basis_function_value(i,x,el_xs,el_xe);
        }
        p = &grid_.get_ygrid();
        for (blas_int i=el_ys; i<=el_ye; ++i){
            hlp_y[i - el_ys] = p->basis_function_value(i,y,el_ys,el_ye);
        }

        dcomp fx = 0.0;
        for (blas_int i=el_ys; i<=el_ye; ++i){
            for (blas_int j=el_xs; j<=el_xe; ++j){
                if (i < get_ysize() && j < get_xsize())
                    fx += body_[i*get_xsize() + j] * hlp_y[i-el_ys] * hlp_x[j-el_xs];
            }
        }
        return fx;
    }
    const zVector& GridVector2d::body() const
    {
        assert(init());
      //
        return body_;
    }
    zVector& GridVector2d::body()
    {
        assert(this->init());
      //
        return body_;
    }

  // modifiers
    void GridVector2d::write_x_section(GridVector& section, blas_int j) // Writes the column of values: j*grid_x.NB() to (j+1)*grid_x.NB() - 1
    {
        assert(init());
        assert(grid_.get_xgrid() == section.get_grid());
      //
        for(blas_int i=0; i<get_xsize(); ++i){
            body_[j*get_xsize() + i] = section[i] * sqrt(grid_.ywz(j));
        }
    }
    GridVector GridVector2d::get_x_section(blas_int j) const     // Retrieves a column of values: j*grid_x.NB() to (j+1)*grid_x.NB() - 1
    {
        assert(init());
      //
        GridVector out(grid_.get_xgrid());
        for(blas_int i=0; i<get_xsize(); ++i){
            out[i] = body_[j*get_xsize() + i] / sqrt(grid_.ywz(j));
        }
        return out;
    }
    GridVector2d& GridVector2d::swap(GridVector2d& rhs)
    {
        Object::swap(rhs);
        body_.swap(rhs.body_);
        grid_.swap(rhs.grid_);
        return *this;
    }
    GridVector2d GridVector2d::copy() const
    {
      //
        GridVector2d out;
        out.grid_ = grid_;
        out.body_ = body_.copy();
        *out.init_ = *init_;
        return out;
    }
    GridVector2d& GridVector2d::complex_conjugate()
    {
        assert(init());
      //
        body_.complex_conjugate();
        return *this;
    }
    GridVector2d& GridVector2d::fill (const dcomp& val)
    {
        assert(init());
      //
        body_.fill(val);
        return *this;
    }
    GridVector2d& GridVector2d::fill_function(const dcomp (*func)(const dcomp&, const dcomp&,const parameters::model_2D<dfloat>&), const parameters::model_2D<dfloat>& mp) // Fills the vector with function values of two dimensional function with one variable fixed (i.e. the first)
    {
        assert(init());
      //
        for (blas_int j=0; j<get_ysize(); ++j){
            for (blas_int i=0; i<get_xsize(); ++i){
                body_[j*get_xsize() + i] = func(grid_.xz(i), grid_.yz(j), mp) * sqrt(grid_.wz(j*get_xsize() + i));
            }
        }
        return *this;
    }

  // operators
    GridVector2d& GridVector2d::operator*= (const dcomp & alpha)
    {
        assert(init());
      //
        body_ *= alpha;
        return *this;
    }
    dcomp GridVector2d::operator* (const GridVector2d& rhs) const
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        return body_ * rhs.body_;
    }
    GridVector2d& GridVector2d::operator+= (const GridVector2d& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        body_ += rhs.body_;
        return *this;
    }
    GridVector2d& GridVector2d::operator-= (const GridVector2d& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        body_ -= rhs.body_;
        return *this;
    }
    GridVector2d& GridVector2d::operator= (GridVector2d tmp)
    {
        swap(tmp);
        return *this;
    }
    ScalarMultiple<dcomp, Vector<dcomp> > GridVector2d::operator* (const dcomp& alpha)
    {
        assert(init());
      //
        return body_ * alpha;
    }
    GridVector2d& GridVector2d::operator= (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs)
    {
        assert(init());                                  // cannot initialize withou grid
        assert(get_size() == rhs.object().get_size());  // incompatible size
      //
        body_ = rhs;
        return *this;
    }
    GridVector2d& GridVector2d::operator+= (ConstScalarMultiple<dcomp, Vector<dcomp> > rhs)
    {
        assert(init());
        assert(get_size() == rhs.object().get_size());
      //
        body_ += rhs;
        return *this;
    }
    GridVector2d& GridVector2d::element_wise_multiplication(const GridVector2d& src)  // TODO: isn't there more effective way?
    {
        assert(init());
        assert(src.init());
        assert(grid_ == src.grid_);
      //
        for (blas_int i=0; i<grid_.get_size(); ++i)
            body_[i] *= src.f(i);
        return *this;
    }

  // custom functions
    GridVector2d& GridVector2d::axpy(const dcomp& alpha, const GridVector2d& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        body_.axpy(alpha, rhs.body_);
        return *this;
    }
    GridVector2d& GridVector2d::ax(const dcomp& alpha, const GridVector2d& rhs)
    {
        assert(init());
        assert(rhs.init());
        assert(grid_ == rhs.grid_);
      //
        body_.ax(alpha, rhs.body_);
        return *this;
    }
    dcomp GridVector2d::reduction(const GridVector2d& y) const
    {
        assert(init());
        assert(*y.init_);
        assert(grid_ == y.grid_);
      //
        return body_.reduction(y.body_);
    }
    dcomp GridVector2d::line_projection(const GridVector& P, const char &axis, blas_int pos) const   // Projection on deltafunction <P|X>, where P = delta(x).chi_j(y)
    {
        assert(init());
        assert( axis == 'x' || axis == 'y' );
        assert( P.get_size() == ((axis=='x')? get_xsize() : get_ysize()) );
      //
        if (axis == 'x'){   // NOTE: The axis specifies the direction of the projection line - IT IS PERPENDICULAR TO THE CHANNEL AXIS!
            return body_.partial_dot_product(P.get_size(), pos * get_xsize(), 1, P.body(), 0 ,1) / sqrt(grid_.ywz(pos));
        } else if (axis == 'y'){
            return body_.partial_dot_product(P.get_size(), pos, get_xsize(), P.body(), 0, 1) / sqrt(grid_.xwz(pos));
        }
        return 0.0;
    }
    GridVector GridVector2d::contraction(const GridVector2d& P, const char& axis)    // axis specifies direction of contraction
    {
        assert(*P.init_);
        assert(init());
        assert(grid_ == P.grid_);
        assert(axis == 'x' || axis == 'y');
      //
        GridVector out;
        if (axis=='x'){
            out = GridVector(grid_.get_ygrid());
            for (blas_int j=0; j<get_ysize(); ++j){
                blas_int pos = j*get_xsize();
                out[j] = body_.partial_dot_product(get_xsize(), pos, 1, P.body(), pos, 1) / sqrt(grid_.ywz(j));
            }
        } else {
            out = GridVector(grid_.get_xgrid());;
            for (blas_int j=0; j<get_xsize(); ++j){
                blas_int inc = get_xsize();
                out[j] = body_.partial_dot_product(get_ysize(), j, inc, P.body(), j, inc) / sqrt(grid_.xwz(j));
            }
        }
        return out;
    }

    // Shallow copy of grid vector
    ShallowGridVector2d::ShallowGridVector2d(FemDvrEcsGrid2d& grid, dcomp *src)
        : GridVector2d()
    {
        this->grid_ = grid;
        Vector<dcomp> X(grid.get_size(), src);
        this->body_.swap(X);
    }
    ShallowGridVector2d::~ShallowGridVector2d()
    {
        //ShallowVector<dcomp> X; // = static_cast<ARRAYS::shallow_vector<dcomp>& >(this->body);
        //this->body_.swap(X);
        //X.clear();
    }


    void GridVector2d::save(const char*filename) const
    {
        assert(init());
      //
        FILE * file;
        def_comp v;
        fopen_s(&file,filename,"w");
        fprintf(file,"#Function on FEM-DVR-ECS grid of %lld basis functions.\n", get_size());
        fprintf(file,"#Coordinate X   \tCoordinate Y   \tReal part of z     \tImaginary part of z\n");
        for (blas_int i=0; i<get_ysize(); ++i){
            for (blas_int j=0; j<get_xsize(); ++j){
                v = f(i*get_xsize() + j);
                fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\n", grid_.yr(i), grid_.xr(j), real(v), imag(v));
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
    void GridVector2d::save_equidistant(const char *name, blas_int x_points, dfloat x_min, dfloat x_max, blas_int y_points, dfloat y_min, dfloat y_max) const
    {
        assert(init());
      //
        FILE * file;
        def_comp v;
        fopen_s(&file,name,"w");

    // Auxiliary variables declarations
        dVector Yr(y_points+1,y_min,y_max, true);
        dVector Xr(x_points+1,x_min,x_max, true);


        blas_int x_start, x_end, y_start, y_end;
        (x_min <= grid_.xr(0))? x_start = 1: x_start = 0;
        (x_max >= grid_.xr(get_xsize()-1))? x_end = x_points-1: x_end = x_points;
        (y_min <= grid_.yr(0))? y_start = 1: y_start = 0;
        (y_max >= grid_.yr(get_ysize()-1))? y_end = y_points-1: y_end = y_points;
        //fprintf(file,"#Function on 2d range (%f,%f)x(%f,%f) enumerated on %dx%d points.\n", x_min, x_max, y_min, y_max, x_end - x_start, y_end - y_start);
        //fprintf(file,"#Coordinate X   \tCoordinate Y   \tReal part of z     \tImaginary part of z\n");
        for (blas_int i=0;i<=y_points;++i){
            if (i<y_start || i > y_end){
                for (blas_int j=0;j<=x_points;++j){
                    fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\n", Yr[i], Xr[j], 0.0, 0.0, 0.0, 0.0);
                }
            } else {
                for (blas_int j=0;j<=x_points;++j){
                    if (j<x_start || j > x_end){
                        fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\n", Yr[i], Xr[j], 0.0, 0.0, 0.0, 0.0);
                    } else {
                        v = evaluate(Xr[j],Yr[i]);
                        fprintf(file, "%.12E\t%.12E\t%.12E\t%.12E\t%.12E\t%.12E\n",Yr[i],Xr[j], real(v), imag(v), pow(abs(v),2), arg(v));
                    }
                }

            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
} // namespace QSCAT
