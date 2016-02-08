/** \brief Hybrid DG-elasticity method for mesh movement
 */

#ifndef __AMESH2DH_H
#include <amesh2dh.hpp>
#endif

#ifndef __ADGM_H
#include <adgm.hpp>
#endif

#ifndef __ALINELAST_P2_STIFFENED_H
#include <alinelast_p2_stiffened.hpp>
#endif

namespace amc {

class DGhybrid
{
	/// input mesh
	UMesh2dh& m;
	/// "Background mesh", which is a coarse mesh generated from [input mesh](@ref m)
	UMesh2d& bm;
	
	/// Delaunay graph mapping context
	DGmove dgm;
	/// Stiffened linear elasticity context
	LinElast_P2 linm;

public:
	DGhybrid(UMesh2dh& mesh, UMesh2d& backmesh) : m(mesh), bm(backmesh);
	void get_backmesh_points();
};

}		// end namespace
