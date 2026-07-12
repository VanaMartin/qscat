#ifndef INCLUDE_PJINPUT_H_
#define INCLUDE_PJINPUT_H_

#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>

#include "picojson/pjson.h"
#include "common.h"
#include "Arrays.h"
#include "FemDvrEcs.h"
#include "FemDvrEcs2d.h"

#include "bessel.h"
#include "coulomb.h"

namespace QSCAT
{
/// Shortcut for the ease of use
typedef picojson::value pjvalue;

/// Shortcut for the ease of use
typedef picojson::array pjarray;

class Parameters
{
    const pjvalue& ref_;
 public:
    Parameters(const pjvalue& ref);
    def_float operator()(const std::string& key) const;
    const pjvalue& operator[](const std::string& key) const;
};

/// Simple reading one file helper
/*!
    This function provides the reading of one json file.
    There are no assumptions on the structure of the file, except json formatting.
*/
pjvalue read_json_file(const std::string& filename);

/// Interpreting the points (trace) grid segment definition.
/*!
    This function provides the interpretation of the points structure
    of the grid segment definition. The points definition assumes the
    elements on the segment to be divided into parts of equidistant
    elements. This definition type expects two arrays: "points" and
    "lengths". The point array marks the postion of borders of equidistant
    parts, the first is provided separately. The lengths then determines the
    length of each element in the parts respectively. NOTE: It is up to the
    provider of the source file to assure that \f[ (p_{i} - p_{i-1} ) = N l_i \f],
    where \f$p_i\f$ denotes the i-th point, \f$l_i\f$ denotes the length of each
    element in give part and \f$N\f$ is the total number of elements in the part.
    If the length is not consistent, the actual postion of the borders may be
    shifted by the over size of the sum (LAZY APPROACH).
    inputs: aa .. default float Buffer reference, destination storing the length of each
    element respectively.
        s .. picojson value, source of paramentrization.
        pos .. default float, starting point of segment discretization
    returns total number of elements.
*/
int parse_points(dBuffer& aa, const pjvalue& s, def_float& pos);

/// Interpreting the increasing grid segment definition.
/*!
    This function provides the interpretation of the points structure
    of the grid segment definition. The points definition assumes the
    elements on the segment to be geometrically increasing its size via
    two operation modes "exp" resp. "mult", resulting in size increment
    as: \f[ l_i = b e^{\alpha i}, \f] resp. \f[ l_i = b m^i, \f]
    where \f$b\f$ dentoes the base length of the first element which is
    providided as parameter "base" either as a number or as string 'last'
    setting the value to be retrieved from previous segment. The maximal
    length of the segment may be provided as "max", the maximum of elements
    must be provided as "elements".
    inputs: aa .. default float Buffer reference, destination storing the length of each
    element respectively.
        s .. picojson value, source of paramentrization.
        pos .. default float, starting point of segment discretization
    returns actual total number of elements.
*/
int parse_uniform_increase(dBuffer& aa, const pjvalue& s, def_float& pos);

/// Interpreting the json node as parametrization of FEM-DVR-ECS grid
/*!
    This function builds all necessary artefacts for building the actual
    FemDvrEcsGrid class. Currently only Gauss-Lobatto quadrature is
    supported.
    input: src .. picojson value, source parametrization;
    returns the FemDvrEcsGrid instantion built accoriding to parametrization.
    TODO: generalize for more complex cases
*/
FemDvrEcsGrid grid_from_parameters(const pjvalue& src);

/// fill with function parametrized by picojson
/*!
    this method evaluates the given function at all discretization points and stores
    them multiplied by the weight factor.
    returns reference on this instantion
*/
GridVector& fill_grid_vector(GridVector& dst, const pjvalue& p, def_comp (*func)(const def_comp&, const Parameters&));

/// fill with two dimensional function along x-axis parametrized by picojson.
/*!
    this method evaluates the given 2D function at all discretization points along x-axis (for fixed y) and stores
    them multiplied by the weight factor.
    returns reference on this instantion
*/
GridVector& fill_grid_vector_xaxis(GridVector& dst, const def_comp& y, const pjvalue& p, def_comp (*func) (const def_comp&, const def_comp&, const Parameters&));

/// fill with two dimensional function along y-axis parametrized by picojson.
/*!
    this method evaluates the given 2D function at all discretization points along y-axis (for fixed x) and stores
    them multiplied by the weight factor.
    returns reference on this instantion
*/
GridVector& fill_grid_vector_yaxis(GridVector& dst, const def_comp& x, const pjvalue& p, def_comp (*func)( const def_comp&, const def_comp&, const Parameters&));

/// Fill with given function
GridVector2d& fill_grid_vector_2d(GridVector2d& dst, const pjvalue& p, def_comp (*func)(const def_comp&, const def_comp&,const Parameters&));

/// Gaussian wavepacket from picojson
def_comp Gaussian(const def_comp& x, const Parameters& p);

} // QSCAT

#endif // INCLUDE_PJINPUT_H_
