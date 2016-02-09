/** \brief Hybrid DG-elasticity method for mesh movement
 *  \author Aditya Kashi
 *  \date Feb 8, 2016
 */

#ifndef __AMESH2DH_H
#include <amesh2dh.hpp>
#endif

#ifndef __ADGM_H
#include <adgm.hpp>
#endif

#ifndef __ALINELAST_P1_H
#include <alinelast_p1.hpp>
#endif

namespace amc {

class DGhybrid
{
	/// input mesh (linear)
	UMesh2dh& m;
	/// quadratic straight version of the input mesh
	UMesh2dh& mq;
	/// "Background mesh", which is a coarse mesh generated from [input mesh](@ref m)
	UMesh2d bm;
	
	/// Delaunay graph mapping context
	DGmove dgm;
	/// Stiffened linear elasticity context
	LinElastP1 linm;
	
	/// Background mesh points
	amat::Matrix<double> backpoints;
	/// Number of background mesh points
	int nbackp;
	/// number of layers at which to take backmesh points
	const int nlayers;
	/// indices of points in the required layer
	std::vector<int> layerpoints;

public:
	DGhybrid(UMesh2dh& mesh, UMesh2dh& qmesh, const int num_layers);
	void compute_backmesh_points();
	void move();
};

DGhybrid::DGhybrid(UMesh2dh& mesh, UMesh2dh& qmesh, const int num_layers) : m(mesh), mq(qmesh), nlayers(num_layers)
{
}

/// Generate a list of points to use for the background mesh by advancing though [layers](@ref nlayers)
/** Currently cannot handle intersecting layer fronts from different unconnected boundaries.
 */
void DGhybrid::compute_backmesh_points()
{
	std::vector<int> laypo(m.gnpoin(),0);
	std::vector<int> layel(m.gnelem(),0);
	int ib, iel, j, ip, idim, nbpoinq = 0;
	int morder = 2; // for quadratic mesh

	// mark the boundary points, 0th layer
	for(ib = 0; ib < m.gnface(); ib++)
	{
		for(j = 0; j < m.gnnofa(); j++)
			laypo[m.gbface(ib,j)] = 1;
	}
	
	std::cout << "DGhybrid: compute_backmesh_points(): Number of boundary points = " << nbpoin << std::endl;

	// get number of boundary points in the original mesh
	for(ip = 0; ip < m.gnpoin(); ip++)
		nbpoinq += laypo[ip];
	// finally, add number of high-order points
	nbpoinq += m.gnface() * morder;

	// for each layer, mark elements containing marked points, and then mark all points of these elements
	for(int ilayer = 0; ilayer < nlayers; ilayer++)
	{
		// mark elements surrounding marked points, and further mark points of these elements
		for(int ip = 0; ip < m.gnpoin(); ip++)
			if(laypo[ip] == 1)
				for(int iel = m.gesup_p(ip); iel < m.gesup_p(ip+1); iel++)
				{
					ele = m.gesup(iel);
					layel[ele] = 1;
					
					// indices of points in the final layer go into layerpoints
					if(ilayer == nlayers-1)
						for(j = 0; j < m.gnnode(ele); j++)
							if(laypo[m.ginpoel(ele,j)] != 1)
								layerpoints.push_back(m.ginpoel(iel,j));
					
					// mark points of this element
					for(j = 0; j < m.gnnode(ele); j++)
						laypo[m.ginpoel(ele,j)] = 1;
				}
	}

	std::cout << "DGhybrid: compute_backmesh_points(): Found points in layer " << nlayers << std::endl;
	
	nbackp = nbpoin + layerpoints.size();
	backpoints.setup(nbackp,m.gndim());

	// add boundary points of the high-order mesh
	for(ib = 0; ib < mq.gnface(); ib++)
	{
		for(j = 0; j < mq.gnnofa(); j++)
			for(idim = 0; idim < mq.gndim(); idim++)
				backpoints(mq.gbface(ib,j), idim) = mq.gcoords(mq.gbface(ib,j),idim);
	}
	// add layer points of the linear mesh
	for(ip = 0; ip < layerpoints.size(); ip++)
		for(idim = 0; idim < m.gndim(); idim++)
			backpoints(nbpoin+ip,idim) = m.gcoords(layerpoints[ip],idim);

	// now generate the background mesh
}

}		// end namespace
