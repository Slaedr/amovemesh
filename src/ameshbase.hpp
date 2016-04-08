/** @brief Class to serve as abstract generic mesh class
 * @author Aditya Kashi
 * @date Feb 24, 2016
 */

#ifndef __AMESHBASE_H
#define __AMESHBASE_H 1

#ifndef __ACONSTANTS_H
#include <aconstants.h>
#endif

#ifndef __AMATRIXT_H
#include <amatrixt.hpp>
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
template<int ndim> class Mesh
{
protected:
	MeshType meshtype;							///< A number describing the [kind of mesh](@ref MeshType). This is purely for information.
	amc_int npoin;								///< Number of points in the mesh
	amc_int nelem;								///< Number of cells/elements in the mesh
	amc_int nbface;								///< Number of boundary faces
	amat::Matrix<amc_real> coords;				///< Coordinates of the points

public:
	amc_real gcoords(amc_int ipoin, int idim);
	amc_int ginpoel(amc_int ielem, int inode);
	amc_int gbface(amc_int iface, int inofa);
};

class Mesh3d : public Mesh<3>
{
protected:
	amc_int nbedge;								///< Number of boundary edges
	constexpr ndim;
public:
	Mesh3d() : ndim(3);
	{ }
};

} // end namespace

#endif
