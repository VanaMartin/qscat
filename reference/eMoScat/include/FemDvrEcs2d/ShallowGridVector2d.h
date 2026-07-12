#ifndef INCLUDE_SHALLOW_GRID_VECTOR_2D_H_
#define INCLUDE_SHALLOW_GRID_VECTOR_2D_H_
namespace QSCAT
{
/** \addtogroup FemDvrEcs2d
* @{ */

//! Simple shallow wrapper for grid vector
class ShallowGridVector2d : public GridVector2d
{
 public:
    //! Constructor
    ShallowGridVector2d(FemDvrEcsGrid2d& grid, dcomp *src);
    //! Destructor
    ~ShallowGridVector2d();
};

/** @} */
}   //  namespace QSCAT
#endif // INCLUDE_SHALLOW_GRID_VECTOR_2D_H_
