/** @brief Class to serve as abstract generic mesh class
 * @author Aditya Kashi
 * @date Feb 24, 2016
 */

#ifndef __AMESHBASE_H
#define __AMESHBASE_H 1

#ifndef __ACONSTANTS_H
#include <aconstants.h>
#endif

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

namespace amc {

/// Possible mesh types
enum MeshType
{
	TRIANGULAR,
	QUADRANGULAR,
	TETRAHEDRAL,
	HEXAHEDRAL,
	PRISMATIC,
	HYBRID,
	HYBRID_PRISM_TET
};

/// Abstract generic mesh class
class Mesh
{
protected:
	MeshType meshtype;							///< A number describing the [kind of mesh](@ref MeshType). This is purely for information.
	int ndim;									///< Dimensionality of mesh - 1D, 2D or 3D
	amc_int npoin;								///< Number of points in the mesh
	amc_int nelem;								///< Number of cells/elements in the mesh
	amc_int nbface;								///< Number of boundary faces
	amat::Matrix<amc_real> coords;				///< Coordinates of the points
	std::vector<std::vector<amc_int>> inpoel;	///< Element-node interconnectivity matrix
	std::vector<std::vector<amc_int>> bface;	///< boundary face data - points and labels

public:
	amc_real gcoords(amc_int ipoin, int idim);
	amc_int ginpoel(amc_int ielem, int inode);
	amc_int gbface(amc_int iface, int inofa);
};

class Mesh3d : public Mesh
{
protected:
	amc_int nbedge;								///< Number of boundary edges
};

} // end namespace

#endif
