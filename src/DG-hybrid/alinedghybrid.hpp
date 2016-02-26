/** \brief Hybrid DG-elasticity method for curved mesh generation
 *  \author Aditya Kashi
 *  \date Feb 8, 2016
 */

#ifndef __AMESH2DHYBRID_H
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
	/** Ordering of the quadratic mesh's boundary points in the back-mesh is decided by bounflag_q.
	 */
	std::vector<int> bounflag_q;

	double lambda;				///< Lame` elasticity constant 1
	double mu;					///< Lame` elasticity constant 2
	double tol;
	int maxiter;
	std::string solver;

public:
	DGhybrid() { }
	DGhybrid(const UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* boundary_motion_quadratic, const int num_layers, 
			const double young, const double nu, const double tol, const int maxiter, const std::string solver);
	void setup(const UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* boundary_motion_quadratic, const int num_layers, 
			const double young, const double nu, const double tol, const int maxiter, const std::string solver);
	void compute_backmesh_points();
	void generate_backmesh_and_compute_displacements();
	void movemesh();
};

DGhybrid::DGhybrid(const UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* b_motion_quadratic, const int num_layers, const double young, const double nu, 
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
	
void DGhybrid::setup(const UMesh2dh* mesh, UMesh2dh* qmesh, const amat::Matrix<double>* boundary_motion_quadratic, const int num_layers, 
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
	std::vector<int> prevlaypo(m->gnpoin(),0);
	std::vector<int> curlaypo(m->gnpoin(),0);
	std::vector<int> layel(m->gnelem(),0);
	int ib, iel, j, ip, idim, ele;
	nbpoin_q = 0;

	// mark the boundary points of linear mesh, 0th layer
	for(ib = 0; ib < m->gnface(); ib++)
	{
		for(j = 0; j < m->gnnofa(); j++)
			curlaypo[m->gbface(ib,j)] = 1;
	}

	// get number of boundary points in the quadratic mesh
	for(ip = 0; ip < mq->gnpoin(); ip++)
		nbpoin_q += bounflag_q[ip];
	
	std::cout << "DGhybrid: compute_backmesh_points(): Number of boundary points = " << nbpoin_q << ", number of layers = " << nlayers << std::endl;

	// for each layer, mark elements containing marked points, and then mark all points of these elements
	for(int ilayer = 0; ilayer < nlayers; ilayer++)
	{
		for(ip = 0; ip < m->gnpoin(); ip++)
		{
			prevlaypo[ip] = curlaypo[ip];
			curlaypo[ip] = 0;
		}

		// mark elements surrounding marked points, and further mark points of these elements
		for(ip = 0; ip < m->gnpoin(); ip++)
		{
			if(prevlaypo[ip] == 1)
			{
				for(int iel = m->gesup_p(ip); iel < m->gesup_p(ip+1); iel++)
				{
					ele = m->gesup(iel);
					if(layel[ele] == 0)
					{
						layel[ele] = 1;
						
						for(j = 0; j < m->gnnode(ele); j++)
							if(prevlaypo[m->ginpoel(ele,j)] != 1)
								curlaypo[m->ginpoel(ele,j)] = 1;
					}
				}
			}
		}
	}

	// add points marked in curlaypo to layerpoints
	for(ip = 0; ip < m->gnpoin(); ip++)
		if(curlaypo[ip] == 1)
			layerpoints.push_back(ip);

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
	
	int ip, ipoin, ib, ilp, idim, ninpoin_q = 0, k = 0, j;

	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		ninpoin_q += bounflag_q[ipoin];

	ninpoin_q = mq->gnpoin() - ninpoin_q;

	inpoints_q.setup(ninpoin_q, mq->gndim());

	k = 0;
	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		if(!bounflag_q[ipoin])
		{
			for(idim = 0; idim < mq->gndim(); idim++)
				inpoints_q(k,idim) = mq->gcoords(ipoin,idim);
			k++;
		}
	std::cout << "DGhybrid: generate_backmesh_and_compute_displacements(): No. of interior points to move = " << inpoints_q.rows() << std::endl;

	// setup DGM and get backmesh
	
	dgm.setup(mq->gndim(), &inpoints_q, &backpoints, &motion_b);
	dgm.generateDG();
	bm = dgm.getDelaunayGraph();
	std::cout << "DGhybrid: Back mesh has " << bm.gnpoin() << " points, " << bm.gnelem() << " elements." << std::endl;
	bm.writeGmsh2("testdg.msh");

	// prepare input for linear elasticity problem
	
	motion_b.setup(nbackp,m->gndim());
	motion_b.zeros();

	// cflags contains 1 if the corresponding backmesh point has a Dirichlet BC
	std::vector<int> cflags(nbackp,0);
	
	// get displacements of the boundary points of the quadratic mesh
	k = 0;
	for(ip = 0; ip < mq->gnpoin(); ip++)
	{
		if(bounflag_q[ip] == 1)
		{
			for(idim = 0; idim < mq->gndim(); idim++)
			{
				motion_b(k, idim) = b_motion_q->get(ip,idim);
				cflags[k] = 1;
			}
			k++;
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
	amat::Matrix<double> x(2*nbackp,1);
	xb.zeros();

	// TODO: add a switch to change solver
	x = sparseCG_d(&A, b, xb, tol, maxiter);

	for(int i = 0; i < nbackp; i++)
		for(idim = 0; idim < mq->gndim(); idim++)
			motion_b(i, idim) = x.get(i+idim*nbackp);
}

void DGhybrid::movemesh()
{
	int ip, ibp, idim;
	
	// now solve Delaunay -- may need to set up once again - CHECK
	dgm.movemesh();

	// re-assemble coords for the quadratic mesh
	amat::Matrix<double> mqcoords(mq->gnpoin(), mq->gndim());
	amat::Matrix<double> intpointsq = dgm.getInteriorPoints();

	int ipoin, k = 0;
	// get interior points
	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		if(!bounflag_q[ipoin])
		{
			for(idim = 0; idim < mq->gndim(); idim++)
				mqcoords(ipoin,idim) = inpoints_q.get(k,idim);
			k++;
		}
	
	// update the backmesh coords
	for(ibp = 0; ibp < nbackp; ibp++)
		for(idim = 0; idim < mq->gndim(); idim++)
			backpoints(ibp,idim) += motion_b.get(ibp,idim);

	// get boundary points from backpoints
	k = 0;
	for(ip = 0; ip < mq->gnpoin(); ip++)
		if(bounflag_q.at(ip) == 1)
		{
			for(idim = 0; idim < mq->gndim(); idim++)
				mqcoords(ip,idim) = backpoints(k,idim);
			k++;
		}
	mq->setcoords(&mqcoords);
	
	bm.setcoords(&backpoints);
	bm.writeGmsh2("testdg_moved.msh");
}

}		// end namespace
