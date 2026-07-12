#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>

#include "picojson/pjson.h"
#include "common.h"
#include "Arrays.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"
#include "pjinput.h"

namespace QSCAT
{
using std::exp;
using std::cout;
using std::endl;

Parameters::Parameters(const pjvalue& ref) : ref_(ref) {}
const pjvalue& Parameters::operator[] (const std::string& key) const { return ref_[key]; }
double Parameters::operator() (const std::string& key) const { return ref_[key].asDouble(); }

pjvalue read_json_file(const std::string& filename)
{

    ifstream source;
    source.open(filename.c_str());
    if (! source.is_open())
      throw invalid_argument(string("Cannot open source json file ") + filename);

    string content((istreambuf_iterator<char>(source)), istreambuf_iterator<char>());
    // see http://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring for the first argument's parentheses

    source.close();

    pjvalue res;
    string err;
    picojson::parse(res, content.begin(), content.end(), &err);
    if (! err.empty())
      throw invalid_argument(string("Cannot parse json file ") + filename + " reason: " + err);

    return res;
}

int parse_points(dBuffer& aa, const pjvalue& s, def_float& pos)
{
    const pjarray& a = s["lengths"].asArray();
    const pjarray& p = s["points"].asArray();
    int elements = 0;
    for (int i=0; i<p.size(); ++i) {
        while(pos < p[i].asDouble()) {
            if (p[i].asDouble() - pos >= a[i].asDouble() ) {
                aa << a[i].asDouble();
                pos += a[i].asDouble();
                elements++;
            } else {
                if (p[i].asDouble() - pos > 0.1 * a[i].asDouble() ) {
                    aa << p[i].asDouble() - pos;
                    pos = p[i].asDouble();
                    elements++;
                } else { // increase by at max 10%
                    aa[aa.get_size() - 1 ] += p[i].asDouble() - pos;
                    pos = p[i].asDouble();
                }
            }
        }
    }
    return elements;
}

int parse_uniform_increase(dBuffer& aa, const pjvalue& s, def_float& pos)
{

    int elements = s["elements"].asInt();
    int count = 0;

    const pjvalue& inc = s["increment"];

    def_float base;
    if ( inc["base"].asString() == "last" ) {
        base = aa[ aa.get_size() - 1 ];
    } else {
        base = inc["base"].asDouble();
    }

    // Skip first num of element increase (same size as base)
    if (inc.isMember("skip")) {
        for (int i=0; i<inc["skip"].asInt(); i++) {
            aa << base;
        }
        count += inc["skip"].asInt();
    }
    elements -= count;

    def_float max = std::numeric_limits<def_float>::max();
    if (inc.isMember("max"))
        max = inc["max"].asDouble();
    def_float sum = 0;
    // The rest is scaled by given function
    if ( inc["type"].asString() == "exp" ) {
        def_float alpha = inc["alpha"].asDouble();
        for (int i=0; i<elements; ++i) {
            if ((sum += base*exp(alpha*i)) > max )
                break;
            aa << base * exp(alpha * i);
            count++;
        }
    } else if ( inc["type"].asString() == "mult" ) {
        def_float alpha = inc["alpha"].asDouble();
        def_float mult = 1.0;
        for (int i=0; i<elements; ++i) {
            if ((sum += base*mult) > max )
                break;
            aa << base * mult;
            mult *= alpha;
            count++;
        }
    }
    return count;
}

FemDvrEcsGrid grid_from_parameters(const pjvalue& src)
{
    int nq = src["dvr_order"].asInt();
    int total=0;
    Vector<blas_int> elements(3);
    def_float pos=0.0, angle[2];

    dBuffer aa;

    // add total negative complex values
    if (src.isMember("complex_negative")) {
        // NIY
    } else {
        elements[0] = 0;
        angle[0] = 0;
    }

 // real elements
    if (src.isMember("real")) {
        const pjvalue& s = src["real"];
        if (total==0) {
            if (s.isMember("start"))
                pos = s["start"].asDouble();
            aa << pos;
        }
        if (s["type"].asString() == "points") {
            elements[1] = parse_points(aa, s, pos);
        } else if (s["type"].asString() == "uniform_increment") {
            elements[1] = parse_uniform_increase(aa, s, pos);
        }
        total += elements[1];
    } else {
        elements[1] = 0;
    }

 // complex elements
    if (src.isMember("complex_positive")) {
        const pjvalue& s = src["complex_positive"];

        angle[1] = s["angle"].asDouble();
        if (s["type"].asString() == "points") {
            elements[2] = parse_points(aa, s, pos);
        } else if (s["type"].asString() == "uniform_increment") {
            elements[2] = parse_uniform_increase(aa, s, pos);
        }
        total += elements[2];
    } else {
        elements[2] = 0;
        angle[1] = 0;
    }

    if (    (src.isMember("complex_negative") && src["complex_negative"].isMember("special")) ||
            (src.isMember("complex_positive") && src["complex_positive"].isMember("special")) ) {

        zBuffer aaz;
        for (int i=0; i<=total; ++i){
            if (i > elements[0] + elements[1]+1) {
                def_float ang = pi * angle[1] / 180 * ( (i - elements[0] - elements[1] ) * 1.0 / (elements[2]-1)  );
                def_comp eit = exp( imu * ang );
                aaz << aa[i] * eit;
            } else {
                aaz << aa[i];
            }
        }

        return FemDvrEcsGrid(nq, total, elements, aaz.as_vector(), angle[1]);
    }

    return FemDvrEcsGrid(nq, total, elements, aa.as_vector(), angle[1]);
}

GridVector& fill_grid_vector(GridVector& dst, const pjvalue& p, def_comp (*func)(const def_comp&, const Parameters&))
{
    assert(dst.init());
    assert(func!=NULL);
  //
    for (blas_int i=0; i<dst.get_size(); ++i){
        dst.f(func(dst.xz(i), p),i);
    }
    return dst;
}

GridVector& fill_grid_vector_xaxis(GridVector& dst, const def_comp& y, const pjvalue& p, def_comp (*func)(const def_comp&, const def_comp&, const Parameters&))
{
    assert(dst.init());
    assert(func!=NULL);
  //
    for (blas_int i=0; i<dst.get_size(); ++i){
        dst.f(func(dst.xz(i), y, p),i);
    }
    return dst;
}

GridVector& fill_grid_vector_yaxis(GridVector& dst, const def_comp& x, const pjvalue& p, def_comp (*func)( const def_comp&, const def_comp&, const Parameters&))
{
    assert(dst.init());
    assert(func!=NULL);
  //
    for (blas_int i=0; i<dst.get_size(); ++i){
        dst.f(func(x, dst.xz(i), p),i);
    }
    return dst;
}

GridVector2d& fill_grid_vector_2d(GridVector2d& dst, const pjvalue& p, def_comp (*func)(const def_comp&, const def_comp&,const Parameters&))
{
    assert(dst.init());
  //
    const FemDvrEcsGrid2d& g = dst.get_grid();

    for (blas_int j=0; j<dst.get_ysize(); ++j){
        for (blas_int i=0; i<dst.get_xsize(); ++i){
            dst.f( func(g.xz(i), g.yz(j), p), j*dst.get_xsize() + i);
        }
    }
    return dst;
}

def_comp Gaussian(const def_comp& x, const Parameters& p)
{
    def_comp core = exp(-pow((x-p("position")),2)/(4*pow(p("sigma"),2)))/sqrt(sqrt(pi*2)*p("sigma"));
    def_comp val = exp(imu*p("impulse")*x)*core;
    return val;
}

} // QSCAT
