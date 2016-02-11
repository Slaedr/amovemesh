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
	/// Prescribed movement of each node in the quadratic mesh
	const amat::Matrix<double>& b_motion_q;
	/// "Background mesh", which is a coarse mesh generated from [input mesh](@ref m)
	UMesh2d bm;
	
	/// Delaunay graph mapping context
	DGmove dgm;
	/// Stiffened linear elasticity context
	LinElastP1 linm;
	
	/// Background mesh points
	amat::Matrix<double> backpoints;
	/// Background mesh motion
	amat::Matrix<double> motion_b;
	/// Interior points of the input quadratic mesh
	amat::Matrix<double> inpoints_q;
	/// Number of background mesh points
	int nbackp;
	/// Number of boundary nodes in quadratic mesh
	int nbpoin_q;
	/// number of layers at which to take backmesh points
	const int nlayers;
	/// indices of points in the required layer
	std::vector<int> layerpoints;

	double lambda;				///< Lame elasticity constant 1
	double mu;					///< Lame elasticity constant 2

public:
	DGhybrid(UMesh2dh& mesh, UMesh2dh& qmesh, const amat::Matrix<double>& boundary_motion_quadratic, const int num_layers, const double young, const double nu);
	void compute_backmesh_points();
	void generate_backmesh_and_compute_displacements();
};

DGhybrid::DGhybrid(UMesh2dh& mesh, UMesh2dh& qmesh, const amat::Matrix<double>& b_motion_quadratic, const int num_layers, const double young, const double nu) 
	: m(mesh), mq(qmesh), b_motion_q(b_motion_quadratic), nlayers(num_layers)
{
	lambda = nu*young/((1.0+nu)*(1.0-2.0*nu));
	mu = young/(2.0*(1.0+nu));
}

/// Generate a list of points to use for the background mesh by advancing though [layers](@ref nlayers)
/** \todo Currently cannot handle intersecting layer fronts from different unconnected boundaries.
 */
void DGhybrid::compute_backmesh_points()
{
	std::vector<int> laypo(m.gnpoin(),0);
	std::vector<int> layel(m.gnelem(),0);
	int ib, iel, j, ip, idim; 
	nbpoin_q = 0;

	int morder = 2; // for quadratic mesh ************

	// mark the boundary points, 0th layer
	for(ib = 0; ib < m.gnface(); ib++)
	{
		for(j = 0; j < m.gnnofa(); j++)
			laypo[m.gbface(ib,j)] = 1;
	}
	
	std::cout << "DGhybrid: compute_backmesh_points(): Number of boundary points = " << nbpoin << std::endl;

	// get number of boundary points in the original mesh
	for(ip = 0; ip < m.gnpoin(); ip++)
		nbpoin_q += laypo[ip];
	// finally, add number of high-order points
	nbpoin_q += m.gnface() * morder;

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
	
	nbackp = nbpoin_q + layerpoints.size();
	backpoints.setup(nbackp,m.gndim());
	b_motion.setup(nbackp,m.gndim());

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
}

/// Generates the background mesh and computes displacements of its nodes using linear elasticity
void DGhybrid::generate_backmesh_and_compute_displacements()
{
	// make a list of interior points of the quadratic mesh
	
	std::vector<int> onboun(mq.gnpoin(),0);
	int ipoin, ib, ilp, idim, ninpoin_q = 0, k = 0;
	
	for(ib = 0; ib < mq.gnface(); ib++)
		for(ilp = 0; ilp < mq.gnnofa(); ilp++)
			onboun[mq.gbface(ib,ilp)] = 1;

	for(ipoin = 0; ipoin < mq.gnpoin(); ipoin++)
		ninpoin_q += onboun[ipoin];

	ninpoin_q = mq.gnpoin() - ninpoin_q;

	inpoints_q.setup(ninpoin_q, mq.gndim());

	for(ipoin = 0; ipoin < mq.gnpoin(); ipoin++)
		if(!onboun[ipoin])
		{
			for(idim = 0; idim < mq.gndim(); idim++)
				inpoints_q(k,idim) = mq.gcoords(ipoin,idim);
			k++;
		}

	// setup DGM and get back mesh
	dgm.setup(mq.gndim(), &inpoints_q, &backpoints, &motion_b);
	dgm.generateDG();
	bm = dgm.getDelaunayGraph();

	motion_b.setup(nbackp,m.gndim());
	motion_b.zeros();
	std::vector<int> cflags(nbackp,0);

	// get displacements of the boundary points of the quadratic mesh
	for(ib = 0; ib < mq.gnface(); ib++)
	{
		for(j = 0; j < mq.gnnofa(); j++)
			for(idim = 0; idim < mq.gndim(); idim++)
			{
				motion_b(mq.gbface(ib,j), idim) = b_motion_q.get(mq.gbface(ib,j),idim);
				cflags[mq.gbface(ib,j),idim] = 1;
			}
	}

	linm.setup(bm, lambda, mu);
	linm.assembleStiffnessMatrix();
	linm.assembleLoadVector();
	linm.dirichletBC_points(cflags, motion_b);
}

}		// end namespace
