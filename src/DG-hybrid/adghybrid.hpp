/** \brief Hybrid DG-elasticity method for curved mesh generation
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

#define __ADGHYBRID_H 1

namespace amc {

/// A hybrid DG-elasticity method for curved mesh generation via mesh movement involving a "background" mesh
/** A background mesh is first constructed using points from the original mesh by a layer-advancing scheme;
 * the user must specify the number of layers to move out from the boundary, and points of this layer are chosen, along with the boundary points, to form the background mesh.
 * This background mesh is moved by linear elasticity, and the rest of the points are moved by using this background mesh as a Delaunay graph for DGM.
 */
class DGhybrid
{
	/// input mesh (linear)
	const UMesh2dh* m;
	/// quadratic straight version of the input mesh
	UMesh2dh* mq;
	/// Prescribed movement of each node in the quadratic mesh
	const amat::Matrix<double>* b_motion_q;
	/// "Background mesh", which is a coarse mesh generated from [input mesh](@ref m)
	UMesh2d bm;
	
	/// Delaunay graph mapping context
	DGmove dgm;
	/// linear elasticity context
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
	int nlayers;
	/// indices of points in the required layer
	std::vector<int> layerpoints;
	/// For each point in the [quadratic mesh](@ref mq), this contains 1 or 0, depending on whether it is a boundary point.
	std::vector<int> bounflag_q;

	double lambda;				///< Lame` elasticity constant 1
	double mu;					///< Lame` elasticity constant 2
	double tol;
	int maxiter;
	std::string solver;

public:
	DGhybrid(UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* boundary_motion_quadratic, const int num_layers, 
			const double young, const double nu, const double tol, const int maxiter, const std::string solver);
	void setup(UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* boundary_motion_quadratic, const int num_layers, 
			const double young, const double nu, const double tol, const int maxiter, const std::string solver);
	void compute_backmesh_points();
	void generate_backmesh_and_compute_displacements();
	void movemesh();
};

DGhybrid::DGhybrid(UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* b_motion_quadratic, const int num_layers, const double young, const double nu, 
		const double toler, const int max_iter, const std::string _solver) 
	: m(mesh), mq(qmesh), b_motion_q(b_motion_quadratic), nlayers(num_layers), tol(toler), maxiter(max_iter), solver(_solver)
{
	lambda = nu*young/((1.0+nu)*(1.0-2.0*nu));
	mu = young/(2.0*(1.0+nu));
	
	// flag the boundary points of the quadratic mesh
	bounflag_q.resize(mq->gnpoin(),0);
	int i,j;
	for(i = 0; i < mq->gnface(); i++)
		for(j = 0; j < mq->gnnofa(); j++)
			bounflag_q[mq->gbface(i,j)] = 1;
}
	
void DGhybrid::setup(UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* boundary_motion_quadratic, const int num_layers, 
		const double young, const double nu, const double toler, const int max_iter, const std::string _solver)
{
	m = mesh;
	mq = qmesh;
	b_motion_q = boundary_motion_quadratic;
	nlayers = num_layers;
	tol = toler;
	maxiter = max_iter;
	solver = _solver;
	lambda = nu*young/((1.0+nu)*(1.0-2.0*nu));
	mu = young/(2.0*(1.0+nu));

	// flag the boundary points of the quadratic mesh
	bounflag_q.resize(mq->gnpoin(),0);
	int i,j;
	for(i = 0; i < mq->gnface(); i++)
		for(j = 0; j < mq->gnnofa(); j++)
			bounflag_q[mq->gbface(i,j)] = 1;
}

/// Generate a list of points to use for the background mesh by advancing though [layers](@ref nlayers)
/** \todo Currently cannot handle intersecting layer fronts from different unconnected boundaries.
 */
void DGhybrid::compute_backmesh_points()
{
	std::vector<int> laypo(m->gnpoin(),0);
	std::vector<int> layel(m->gnelem(),0);
	int ib, iel, j, ip, idim, ele;
	nbpoin_q = 0;

	int morder = 2; // for quadratic mesh ************

	// mark the boundary points, 0th layer
	for(ib = 0; ib < m->gnface(); ib++)
	{
		for(j = 0; j < m->gnnofa(); j++)
			laypo[m->gbface(ib,j)] = 1;
	}

	//debug:
	/*for(j = 0; j < m->gnpoin(); j++)
		std::cout << laypo[j] << " ";
	std::cout << std::endl;*/

	// get number of boundary points in the original linear mesh
	for(ip = 0; ip < m->gnpoin(); ip++)
		nbpoin_q += laypo[ip];
	// add number of high-order points to get number of boundary points in the quadratic mesh
	nbpoin_q += m->gnface() * (morder-1);
	
	std::cout << "DGhybrid: compute_backmesh_points(): Number of boundary points = " << nbpoin_q << std::endl;

	// for each layer, mark elements containing marked points, and then mark all points of these elements
	for(int ilayer = 0; ilayer < nlayers; ilayer++)
	{
		// mark elements surrounding marked points, and further mark points of these elements
		for(int ip = 0; ip < m->gnpoin(); ip++)
			if(laypo[ip] == 1)
				for(int iel = m->gesup_p(ip); iel < m->gesup_p(ip+1); iel++)
				{
					ele = m->gesup(iel);
					layel[ele] = 1;
					
					// indices of points in the final layer go into layerpoints
					if(ilayer == nlayers-1)
						for(j = 0; j < m->gnnode(ele); j++)
							if(laypo[m->ginpoel(ele,j)] != 1)
								layerpoints.push_back(m->ginpoel(ele,j));
					
					// mark points of this element
					for(j = 0; j < m->gnnode(ele); j++)
						laypo[m->ginpoel(ele,j)] = 1;
				}
	}

	std::cout << "DGhybrid: compute_backmesh_points(): Found " << layerpoints.size() << " points in layer " << nlayers << std::endl;
	
	nbackp = nbpoin_q + layerpoints.size();
	backpoints.setup(nbackp,m->gndim());

	// add boundary points of the high-order mesh
	int k = 0;
	for(ip = 0; ip < mq->gnpoin(); ip++)
		if(bounflag_q.at(ip) == 1)
		{
			for(idim = 0; idim < mq->gndim(); idim++)
				backpoints(k,idim) = mq->gcoords(ip,idim);
			k++;
		}
	if(k != nbpoin_q) std::cout << "DGhybrid: compute_backmesh_points(): Error in getting the points!" << std::endl;

	// add layer points of the linear mesh
	for(ip = 0; ip < layerpoints.size(); ip++)
		for(idim = 0; idim < m->gndim(); idim++)
			backpoints(nbpoin_q+ip,idim) = m->gcoords(layerpoints[ip],idim);

	std::cout << "DGhybrid: compute_backmesh_points(): Backmesh has " << nbackp << " points." << std::endl;
}

/// Generates the background mesh and computes displacements of its nodes using linear elasticity
void DGhybrid::generate_backmesh_and_compute_displacements()
{
	// make a list of interior points of the quadratic mesh
	
	std::vector<int> onboun(mq->gnpoin(),0);
	int ipoin, ib, ilp, idim, ninpoin_q = 0, k = 0, j;
	
	for(ib = 0; ib < mq->gnface(); ib++)
		for(ilp = 0; ilp < mq->gnnofa(); ilp++)
			onboun[mq->gbface(ib,ilp)] = 1;

	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		ninpoin_q += onboun[ipoin];

	ninpoin_q = mq->gnpoin() - ninpoin_q;

	inpoints_q.setup(ninpoin_q, mq->gndim());

	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		if(!onboun[ipoin])
		{
			for(idim = 0; idim < mq->gndim(); idim++)
				inpoints_q(k,idim) = mq->gcoords(ipoin,idim);
			k++;
		}

	// setup DGM and get back mesh
	dgm.setup(mq->gndim(), &inpoints_q, &backpoints, &motion_b);
	dgm.generateDG();
	bm = dgm.getDelaunayGraph();
	std::cout << "DGhybrid: Back mesh has " << bm.gnpoin() << " points, " << bm.gnelem() << " elements." << std::endl;
	bm.writeGmsh2("testdg.msh");

	motion_b.setup(nbackp,m->gndim());
	motion_b.zeros();

	// cflags contains 1 if the corresponding backmesh point has a Dirichlet BC
	std::vector<int> cflags(nbackp,0);

	// get displacements of the boundary points of the quadratic mesh
	for(ib = 0; ib < mq->gnface(); ib++)
	{
		for(j = 0; j < mq->gnnofa(); j++)
			for(idim = 0; idim < mq->gndim(); idim++)
			{
				motion_b(mq->gbface(ib,j), idim) = b_motion_q->get(mq->gbface(ib,j),idim);
				cflags[mq->gbface(ib,j),idim] = 1;
			}
	}

	// setup and solve the elasticity equations to get displacement of the background mesh
	std::cout << "Starting linelast" << std::endl;	
	linm.setup(&bm, lambda, mu);
	linm.assembleStiffnessMatrix();
	linm.assembleLoadVector();
	linm.dirichletBC_points(cflags, motion_b);

	amat::SpMatrix A = linm.stiffnessMatrix();
	amat::Matrix<double> b = linm.loadVector();
	amat::Matrix<double> xb(2*nbackp,1);
	xb.zeros();

	/// \todo add a switch to change solver
	xb = sparseCG_d(&A, b, xb, tol, maxiter);

	for(int i = 0; i < nbackp; i++)
		for(idim = 0; idim < mq->gndim(); idim++)
			motion_b(nbpoin_q+i, idim) = xb.get(i+idim*nbpoin_q);
}

void DGhybrid::movemesh()
{
	int ibp, idim;
	
	// now solve Delaunay -- may need to set up once again - CHECK
	dgm.movemesh();

	inpoints_q = dgm.getInteriorPoints();
	for(ibp = 0; ibp < nbackp; ibp++)
		for(idim = 0; idim < mq->gndim(); idim++)
			backpoints(ibp,idim) += motion_b.get(ibp,idim);

	bm.setcoords(&backpoints);
	bm.writeGmsh2("testdg_moved.msh");
}

}		// end namespace
