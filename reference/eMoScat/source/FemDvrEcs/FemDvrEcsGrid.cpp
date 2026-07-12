#include <cassert>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "Arrays.h"
#include "input.h"

#include "FemDvrEcs/DvrGrid.h"
#include "FemDvrEcs/FemDvrEcsGrid.h"

namespace QSCAT
{

    void FemDvrEcsGrid::initialize(blas_int quadrature, blas_int total_elements, const iVector& elements, const dVector& aa, dfloat alpha)
    {
    // auxiliary variables
        blas_int nall;
        blas_int i;
    // initialization controller
        *init_ = true;
    // build dvr grid on -1,1
        g1_ = DvrGrid(quadrature);
    // build matrix of derivatives
        dLp_ = dMatrix(quadrature,quadrature);
        DLagPol(quadrature, dLp_, g1_);
        tnel_ = total_elements;
        nel_ = iVector(3);
        nel_[0] = elements[0];
        nel_[1] = elements[1];
        nel_[2] = elements[2];
        theta_ = alpha;

        nall = tnel_ * (quadrature - 1) + 1;     // number of all grid points including x_min and x_max
        nb_ = nall - 2;                  // number of basis functions = nall - 2 (functions are assumed to be zero at the ends)

        ix0_neg_ = nel_[0] * (quadrature - 1);
        ix0_pos_ = (nel_[0] + nel_[1]) * (quadrature - 1);

        if (nel_[0] == 0) {
            ix0_neg_ = 0;        // to avoid problems if these
        }
        if (nel_[2] == 0) {
            ix0_pos_ = nb_;       // indices are used
        }
        nr_ = ix0_pos_;
        eit_ = exp(imu * theta_ * pi / def_float(180.0));

    // real endpoints of elements
        ar_ = dVector(tnel_+1);
        ar_[0] = aa[0];
        for (i=1; i<=tnel_; ++i) {
            ar_[i] = ar_[i-1] + aa[i];
        }
        x_min_ = ar_[0];
        x0_neg_ = ar_[nel_[0]];
        x0_pos_ = ar_[nel_[0] + nel_[1]];
        R0_ = x0_pos_;
        x_max_ = ar_[tnel_];

    // scaled lengths and endpoints of elements
        aaz_ = zVector(tnel_+1);
        for (i=0; i<=tnel_; ++i){
            aaz_[i] = dcomp(aa[i]);
        }
        zVector az(tnel_+1);
        for (i=0;i<=tnel_;++i){
            az[i] = ar_[i];
        }
        az[0] = dcomp(x0_neg_) + dcomp((x_min_ - x0_neg_) * eit_);
        aaz_[0] = aaz_[0]*eit_;
        for (i=1; i<=nel_[0]; ++i){
            aaz_[i] = aaz_[i]*eit_;
            az[i] = az[i-1] + aaz_[i];
        }
        for (i=nel_[0] + nel_[1] + 1; i<=tnel_; ++i){
            aaz_[i] = aaz_[i]*eit_;
            az[i] = az[i-1] + aaz_[i];
        }
        z_min_ = az[0];
        z_max_ = az[tnel_];

    // full grid including endpoints stored in auxiliary arrays starting from index 0
        zVector all_wz(nall);
        zVector all_xz(nall);
        dVector all_xr(nall);
        blas_int ii;
        dcomp hz;
        all_wz.fill(dComplex(0.0));
        for (i=0; i<tnel_; ++i){
            ii = i * (quadrature - 1);          // index of the first point of ith-element
            hz = def_float(0.5) * aaz_[i+1];
            for (blas_int j=0; j<quadrature; ++j){
                all_xz[ii + j] = hz * g1_.x(j) + hz + az[i];        // shift of the grid to [az(i-1),az(i)]
                all_wz[ii + j] = all_wz[ii + j] + hz * g1_.w(j);    // scaling of weight to [az(i-1),az(i)]
                                                                        // at points where one element ends and another begins,
                                                                        // weights are symbolically wz = w_i(nq) + w_i+1(1)
                all_xr[ii + j] = def_float(0.5) * (aa[i+1] * g1_.x(j) + aa[i+1]) + ar_[i];
                if (i == nel_[0]) { wx0_neg_ = hz * g1_.w(1); } // 1 or 0?
                if (i == nel_[0] + nel_[1]) { wx0_pos_ = hz * g1_.w(quadrature-1);}
            }
        }
    // actual grid w/o endpoints (arrays starting from index 1)
        xz_ = zVector(nb_);
        wz_ = zVector(nb_);
        xr_ = dVector(nb_);
        for (i=0;i<nb_;++i) {
            xz_[i] = all_xz[i+1];
            wz_[i] = all_wz[i+1];
            xr_[i] = all_xr[i+1];
        }

        *init_ = true;
    }

    void FemDvrEcsGrid::initialize(blas_int quadrature, blas_int total_elements, const iVector& elements, const zVector& aaz, dfloat alpha)
    {
    // auxiliary variables
        blas_int nall;
        blas_int i;
    // initialization controller
        *init_ = true;
    // build dvr grid on -1,1
        g1_ = DvrGrid(quadrature);
    // build matrix of derivatives
        dLp_ = dMatrix(quadrature,quadrature);
        DLagPol(quadrature, dLp_, g1_);
        tnel_ = total_elements;
        nel_ = iVector(3);
        nel_[0] = elements[0];
        nel_[1] = elements[1];
        nel_[2] = elements[2];
        theta_ = alpha;

        nall = tnel_ * (quadrature - 1) + 1;     // number of all grid points including x_min and x_max
        nb_ = nall - 2;                  // number of basis functions = nall - 2 (functions are assumed to be zero at the ends)

        ix0_neg_ = nel_[0] * (quadrature - 1);
        ix0_pos_ = (nel_[0] + nel_[1]) * (quadrature - 1);

        if (nel_[0] == 0) {
            ix0_neg_ = 0;        // to avoid problems if these
        }
        if (nel_[2] == 0) {
            ix0_pos_ = nb_;       // indices are used
        }
        nr_ = ix0_pos_;
        eit_ = exp(imu * theta_ * pi / def_float(180.0));

    // scaled lengths and endpoints of elements
        aaz_ = zVector(tnel_+1);
        for (i=0; i<=tnel_; ++i){
            aaz_[i] = aaz[i];
        }

        zVector az(tnel_+1);
        dVector aa(tnel_+1);

        az[0] = aaz_[0];
        for (i=1;i<=nel_[0];++i){
            az[i] = az[i-1] + aaz_[i];
            aa[i] = abs(aaz_[i]);
        }

      // kill the cumulated imaginary part
        assert(imag(az[nel_[0]]) < 1e-15);
        az[nel_[0]] = real(az[nel_[0]]);

        for (i=nel_[0]+1;i<=nel_[1];++i){
            assert(imag(aaz_[i]) == 0);
            az[i] = az[i-1] + aaz_[i];
            aa[i] = abs(aaz_[i]);
        }

        for (i=nel_[0] + nel_[1] + 1; i<=tnel_; ++i) {
            az[i] = az[i-1] + aaz_[i];
            aa[i] = abs(aaz_[i]);
        }

        //az[0] = dcomp(x0_neg_) + dcomp((x_min_ - x0_neg_) * eit_);
        //aaz_[0] = aaz_[0]*eit_;
        //for (i=1; i<=nel_[0]; ++i){
        //    aaz_[i] = aaz_[i]*eit_;
        //    az[i] = az[i-1] + aaz_[i];
        //}
        //for (i=nel_[0] + nel_[1] + 1; i<=tnel_; ++i){
        //    aaz_[i] = aaz_[i]*eit_;
        //    az[i] = az[i-1] + aaz_[i];
        //}
        z_min_ = az[0];
        z_max_ = az[tnel_];


    // real endpoints of elements
        ar_ = dVector(tnel_+1);
        //ar_[0] = aa[0];
        //for (i=1; i<=tnel_; ++i) {
        //    ar_[i] = ar_[i-1] + abs(aaz[i]);
        //}

        // values should have real part only
        for (i=nel_[0]+1; i<=nel_[0] + nel_[1]; ++i)
            ar_[i] = real(az[i]);

        // left complex
        for (i=nel_[0]-1; i>=0; --i)
            ar_[i] = ar_[i+1] - abs(aaz_[i+1]);
        // right complex
        for (i=nel_[0]+nel_[1]+1; i<=tnel_; ++i)
            ar_[i] = ar_[i-1] + abs(aaz_[i]);

        aa[0] = ar_[0];
        x_min_ = ar_[0];
        x0_neg_ = ar_[nel_[0]];
        x0_pos_ = ar_[nel_[0] + nel_[1]];
        R0_ = x0_pos_;
        x_max_ = ar_[tnel_];

    // full grid including endpoints stored in auxiliary arrays starting from index 0
        zVector all_wz(nall);
        zVector all_xz(nall);
        dVector all_xr(nall);
        blas_int ii;
        dcomp hz;
        all_wz.fill(dComplex(0.0));
        for (i=0; i<tnel_; ++i){
            ii = i * (quadrature - 1);          // index of the first point of ith-element
            hz = def_float(0.5) * aaz_[i+1];
            for (blas_int j=0; j<quadrature; ++j){
                all_xz[ii + j] = hz * g1_.x(j) + hz + az[i];        // shift of the grid to [az(i-1),az(i)]
                all_wz[ii + j] = all_wz[ii + j] + hz * g1_.w(j);    // scaling of weight to [az(i-1),az(i)]
                                                                        // at points where one element ends and another begins,
                                                                        // weights are symbolically wz = w_i(nq) + w_i+1(1)
                all_xr[ii + j] = def_float(0.5) * (aa[i+1] * g1_.x(j) + aa[i+1]) + ar_[i];
                if (i == nel_[0]) { wx0_neg_ = hz * g1_.w(1); } // 1 or 0?
                if (i == nel_[0] + nel_[1]) { wx0_pos_ = hz * g1_.w(quadrature-1);}
            }
        }
    // actual grid w/o endpoints (arrays starting from index 1)
        xz_ = zVector(nb_);
        wz_ = zVector(nb_);
        xr_ = dVector(nb_);
        for (i=0;i<nb_;++i) {
            xz_[i] = all_xz[i+1];
            wz_[i] = all_wz[i+1];
            xr_[i] = all_xr[i+1];
        }

        *init_ = true;
    }

    bool FemDvrEcsGrid::save_bin_body(std::ofstream &file) const
    {
        if(file.is_open()){
            file.write((char*) &nb_, sizeof(blas_int));
            file.write((char*) &ix0_neg_, sizeof(blas_int));
            file.write((char*) &ix0_pos_, sizeof(blas_int));
            file.write((char*) &nr_, sizeof(blas_int));
            file.write((char*) &x_min_, sizeof(dfloat));
            file.write((char*) &x0_neg_, sizeof(dfloat));
            file.write((char*) &x0_pos_, sizeof(dfloat));
            file.write((char*) &R0_, sizeof(dfloat));
            file.write((char*) &x_max_, sizeof(dfloat));
            file.write((char*) &theta_, sizeof(dfloat));
            file.write((char*) &eit_, sizeof(dcomp));
            file.write((char*) &z_min_, sizeof(dcomp));
            file.write((char*) &z_max_, sizeof(dcomp));

            if(!nel_.save_binary(file)) goto save_break;
            if(!xz_.save_binary(file)) goto save_break;
            if(!xr_.save_binary(file)) goto save_break;
            if(!wz_.save_binary(file)) goto save_break;
            if(!ar_.save_binary(file)) goto save_break;

            file.write((char*) &wx0_neg_, sizeof(dcomp));
            file.write((char*) &wx0_pos_, sizeof(dcomp));
            file.write((char*) &tnel_, sizeof(blas_int));

            if(!aaz_.save_binary(file)) goto save_break;
            if(!g1_.save_binary(file)) goto save_break;
            if(!dLp_.save_binary(file)) goto save_break;
            return true;
        }
    save_break:
        assert(0);      // should not happen
        return false;
    }
    bool FemDvrEcsGrid::read_bin_body(std::ifstream & file)
    {
        if (file.is_open()){
            file.read((char*) &nb_, sizeof(blas_int));
            file.read((char*) &ix0_neg_, sizeof(blas_int));
            file.read((char*) &ix0_pos_, sizeof(blas_int));
            file.read((char*) &nr_, sizeof(blas_int));
            file.read((char*) &x_min_, sizeof(dfloat));
            file.read((char*) &x0_neg_, sizeof(dfloat));
            file.read((char*) &x0_pos_, sizeof(dfloat));
            file.read((char*) &R0_, sizeof(dfloat));
            file.read((char*) &x_max_, sizeof(dfloat));
            file.read((char*) &theta_, sizeof(dfloat));
            file.read((char*) &eit_, sizeof(dcomp));
            file.read((char*) &z_min_, sizeof(dcomp));
            file.read((char*) &z_max_, sizeof(dcomp));

            if(!nel_.read_binary(file)) goto read_break;
            if(!xz_.read_binary(file)) goto read_break;
            if(!xr_.read_binary(file)) goto read_break;
            if(!wz_.read_binary(file)) goto read_break;
            if(!ar_.read_binary(file)) goto read_break;

            file.read((char*) &wx0_neg_, sizeof(dcomp));
            file.read((char*) &wx0_pos_, sizeof(dcomp));
            file.read((char*) &tnel_, sizeof(blas_int));

            if(!aaz_.read_binary(file)) goto read_break;
            if(!g1_.read_binary(file)) goto read_break;
            if(!dLp_.read_binary(file)) goto read_break;

            *init_ = true;
            return true;
        }
    read_break:
        return false;
    }

    // Constructors

    FemDvrEcsGrid::FemDvrEcsGrid() : Object()
    {
        nb_ = 0;
        nr_ = 0;
        ix0_neg_ = 0;
        ix0_pos_ = 0;
        x_min_ = dfloat(0);
        x0_neg_ = dfloat(0);
        x0_pos_ = dfloat(0);
        R0_ = dfloat(0);
        x_max_ = dfloat(0);
        theta_ = dfloat(0);
        eit_ = dcomp(0);
        z_min_ = dcomp(0);
        z_max_ = dcomp(0);
        tnel_ = 0;
        wx0_neg_ = dcomp(0);
        wx0_pos_ = dcomp(0);
    }

    FemDvrEcsGrid::FemDvrEcsGrid(blas_int quadrature, blas_int total_elements, const iVector& elements, const dVector& aa, dfloat alpha) : Object()
    {
        assert(total_elements > 1);
        assert(elements.get_size() == 3);
        assert(quadrature > 1);
      //
        initialize(quadrature, total_elements, elements, aa, alpha);
    }
    FemDvrEcsGrid::FemDvrEcsGrid(blas_int quadrature, blas_int total_elements, const iVector& elements, const zVector& aaz, dfloat alpha) : Object()
    {
        assert(total_elements > 1);
        assert(elements.get_size() == 3);
        assert(quadrature > 1);
      //
        initialize(quadrature, total_elements, elements, aaz, alpha);
    }
    FemDvrEcsGrid::FemDvrEcsGrid(parameters::grid<dfloat> &gp) : Object()
    {
        assert(gp.tnel > 1);
        assert(gp.nq > 1);
      //
        Vector<blas_int> aux(3);
        aux[0] = gp.nel[0];
        aux[1] = gp.nel[1];
        aux[2] = gp.nel[2];
        initialize(gp.nq, gp.tnel, aux, *gp.aa, gp.theta);
    }
    FemDvrEcsGrid::FemDvrEcsGrid(const FemDvrEcsGrid &old):
        Object(old),
        nb_(old.nb_),
        ix0_neg_(old.ix0_neg_),
        ix0_pos_(old.ix0_pos_),
        nr_(old.nr_),
        x_min_(old.x_min_),
        x0_neg_(old.x0_neg_),
        x0_pos_(old.x0_pos_),
        R0_(old.R0_),
        x_max_(old.x_max_),
        theta_(old.theta_),
        eit_(old.eit_),
        z_min_(old.z_min_),
        z_max_(old.z_max_),
        wx0_neg_(old.wx0_neg_),
        wx0_pos_(old.wx0_pos_),
        tnel_(old.tnel_),
        xz_(old.xz_),
        xr_(old.xr_),
        wz_(old.wz_),
        nel_(old.nel_),
        ar_(old.ar_),
        aaz_(old.aaz_),
        g1_(old.g1_),
        dLp_(old.dLp_)
    {}
    FemDvrEcsGrid::~FemDvrEcsGrid()
    {
        decref();
    }

 // operators

    FemDvrEcsGrid& FemDvrEcsGrid::operator= (FemDvrEcsGrid tmp)
    {
        this->swap(tmp);
        return *this;
    }
    bool FemDvrEcsGrid::operator== (const FemDvrEcsGrid& rhs) const
    {
        bool error = true;
        if(nb_!=rhs.nb_) error = false;
        if(ix0_neg_!=rhs.ix0_neg_) error = false;
        if(ix0_pos_!=rhs.ix0_pos_) error = false;
        if(nr_!=rhs.nr_) error = false;
        if(x_min_!=rhs.x_min_) error = false;
        if(x0_neg_!=rhs.x0_neg_) error = false;
        if(x0_pos_!=rhs.x0_pos_) error = false;
        if(R0_!=rhs.R0_) error = false;
        if(x_max_!=rhs.x_max_) error = false;
        if(theta_!=rhs.theta_) error = false;
        if(eit_!=rhs.eit_) error = false;
        if(z_min_!=rhs.z_min_) error = false;
        if(z_max_!=rhs.z_max_) error = false;
        if(g1_.quadrature() != rhs.g1_.quadrature()) error = false;
        return error;
    }
    bool FemDvrEcsGrid::operator!= (const  FemDvrEcsGrid& rhs) const
    {
        return !(*this==rhs);
    }

 // modifiers

    FemDvrEcsGrid& FemDvrEcsGrid::swap(FemDvrEcsGrid& rhs)
    {
        Object::swap(rhs);
        std::swap(nb_,       rhs.nb_);
        std::swap(ix0_neg_,  rhs.ix0_neg_);
        std::swap(ix0_pos_,  rhs.ix0_pos_);
        std::swap(nr_,       rhs.nr_);
        std::swap(x_min_,    rhs.x_min_);
        std::swap(x0_neg_,   rhs.x0_neg_);
        std::swap(x0_pos_,   rhs.x0_pos_);
        std::swap(R0_,       rhs.R0_);
        std::swap(x_max_,    rhs.x_max_);
        std::swap(theta_,    rhs.theta_);
        std::swap(eit_,      rhs.eit_);
        std::swap(z_min_,    rhs.z_min_);
        std::swap(z_max_,    rhs.z_max_);
        xz_.swap(rhs.xz_);
        xr_.swap(rhs.xr_);
        wz_.swap(rhs.wz_);
        nel_.swap(rhs.nel_);
        ar_.swap(rhs.ar_);
        std::swap(wx0_neg_,  rhs.wx0_neg_);
        std::swap(wx0_pos_,  rhs.wx0_pos_);
        std::swap(tnel_,     rhs.tnel_);
        aaz_.swap(rhs.aaz_);
        g1_.swap(rhs.g1_);
        dLp_.swap(rhs.dLp_);
        return *this;
    }
    FemDvrEcsGrid FemDvrEcsGrid::copy() const
    {
        FemDvrEcsGrid out;
        *out.init_ = init();
        out.nb_ = nb_;
        out.ix0_neg_ = ix0_neg_;
        out.ix0_pos_ = ix0_pos_;
        out.nr_ = nr_;
        out.x_min_ = x_min_;
        out.x0_neg_ = x0_neg_;
        out.x0_pos_ = x0_pos_;
        out.R0_ = R0_;
        out.x_max_ = x_max_;
        out.theta_ = theta_;
        out.eit_ = eit_;
        out.z_min_ = z_min_;
        out.z_max_ = z_max_;
        out.xz_ = xz_.copy();
        out.xr_ = xr_.copy();
        out.wz_ = wz_.copy();
        out.nel_ = nel_.copy();
        out.ar_ = ar_.copy();
        out.wx0_neg_ = wx0_neg_;
        out.wx0_pos_ = wx0_pos_;
        out.tnel_ = tnel_;
        out.aaz_ = aaz_.copy();
        out.g1_ = g1_.copy();
        out.dLp_ = dLp_.copy();
        out.incref();
        return *this;
    }

 // acessors

    const dfloat& FemDvrEcsGrid::xr(const blas_int& i) const
    {
        assert(i<nb_);
      //
        return xr_[i];
    }
    const dcomp & FemDvrEcsGrid::x(const blas_int& i) const
    {
        assert(i<nb_);
      //
        return xz_[i];
    }
    const dcomp & FemDvrEcsGrid::w(const blas_int& i) const
    {
        assert(i<nb_);
      //
        return wz_[i];
    }
    const blas_int & FemDvrEcsGrid::nel(const blas_int& i) const
    {
        assert(i<3);
      //
        return nel_[i];
    }
    const dcomp& FemDvrEcsGrid::aaz(const blas_int& i) const
    {
        assert(i<tnel_+1);
      //
        return aaz_[i];
    }
    dfloat FemDvrEcsGrid::wq(const blas_int& i)
    {
        assert(i < quadrature());
        return g1_.w(i);
    }
    const dMatrix& FemDvrEcsGrid::dlp() const
    {
        return dLp_;
    }
    const dfloat& FemDvrEcsGrid::dlp(const blas_int& i, const blas_int& j) const
    {
        return dLp_(i,j);
    }
    blas_int FemDvrEcsGrid::nb() const
    {
        assert(init());
      //
        return nb_;
    }
    blas_int FemDvrEcsGrid::get_size() const
    {
        assert(init());
      //
        return nb_;
    }
    blas_int FemDvrEcsGrid::quadrature() const
    {
        return g1_.quadrature();
    }
    blas_int FemDvrEcsGrid::tnel() const
    {
        return tnel_;
    }
    blas_int FemDvrEcsGrid::nr() const
    {
        return nr_;
    }
    const dfloat& FemDvrEcsGrid::x_pos() const
    {
        return x0_pos_;
    }
    const dfloat& FemDvrEcsGrid::x_neg() const
    {
        return x0_neg_;
    }
    const DvrGrid& FemDvrEcsGrid::dvr() const
    {
        assert(init());
      //
        return g1_;
    }

    // Functions

    blas_int FemDvrEcsGrid::get_element_end(const blas_int& i) const
    {
        if (i >= tnel_){
            return nb_-1;
        } else {
            return i * (g1_.quadrature() - 1) - 1;
        }
    }
    blas_int FemDvrEcsGrid::get_element_start(const blas_int& i) const
    {
        if (i!=1){
            return (i-1)*(g1_.quadrature() - 1) - 1;
        } else {
            return 0;
        }
    }
    blas_int FemDvrEcsGrid::get_element_index(const dfloat& X) const           // If X in (ar_(i-1), ar_i) returns i
    {
        for (blas_int i=0; i<tnel_; ++i){
            if (X < ar_[i]){
                return i;
            }
        }
        return tnel_;
    }
    blas_int FemDvrEcsGrid::get_element_end_x(const dfloat& X) const    // Determines the index of the higher endpoint of element where X belongs
    {
        blas_int i = get_element_index(X);
        return get_element_end(i);
    }
    blas_int FemDvrEcsGrid::get_element_start_x(const dfloat& X) const  // Determines the index of the higher endpoint of element where X belongs
    {
        blas_int i = get_element_index(X);
        return get_element_start(i);
    }
    dcomp FemDvrEcsGrid::basis_function_value(const blas_int& i, const dfloat& x, const blas_int& start, const blas_int& end) const    // Value of the basis function at point x in its element
    {
    // Applying the complex scaling if necessary
        dcomp z;
        if (x < x0_neg_) {
            z = x0_neg_ + (x - x0_neg_) * eit_;
        } else if(x <= x0_pos_) {
            z = dcomp(x);
        } else {
            z = x0_pos_ + (x - x0_pos_)*eit_;
        }

        dcomp lag = 1.0;
        for (blas_int l=start;l<=end;++l){
            if (l!=i){
                lag *= (z-xz_[l])/(xz_[i] - xz_[l]);
            }
        }
        if (start == 0) {
            lag *= (z - z_min_)/(xz_[i] - z_min_);
        }
        if (end == nb_ - 1){
            lag *= (z - z_max_)/(xz_[i] - z_max_);
        }
        return lag/sqrt(wz_[i]);
    }
    const dfloat& FemDvrEcsGrid::ar(const blas_int& i) const
    {
        assert(i<tnel_+1);
        return ar_[i];
    }
} // namespace QSCAT
