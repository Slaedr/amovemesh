/** @brief Class to convert a 3D linear mesh into quadratic mesh.
 * @author Aditya Kashi
 * @date March 19, 2016
 * 
 * Notes:
 * We could set the support radius automatically based on some function of the curvature of the boundary.
 */

#ifndef __AGEOMETRYH_H
	#include <ageometry3d.hpp>
#endif

#ifndef __ARBF_H
	#include <arbf.hpp>
#endif

#define __ACURVEDMESHGEN3D_H 1

namespace amc {

typedef RBFmove Meshmove;

/** Class to generate curved mesh from a linear mesh using cubic spline reconstruction and one of the mesh movement techniques. */

class CurvedMeshGen
{
	UMesh* m;						///< Data about the original linear mesh. We need this to compute spline reconstruction of the boundary.
	UMesh* mq;						///< Data of the corresponding (straight-faced) quadratic mesh
	Meshmove* mmv;					///< Pointer to parent class for the mesh-movement classes, such RBF, DGM or linear elasticity.
	BoundaryReconstruction* br;		///< Object to reconstruct the boundary using cubic splines.
	int degree;						///< Degree of the generated mesh (2 or 3; currently only 2 is supported)
	
	double tol;						///< Tolerance for linear solver used for computing spline coefficients.
	int maxiter;					///< Maximum number of iterations for linear solver used to compute spline coefficients.
	int rbfchoice;					///< Parameters for mesh movement - the type of RBF to be used, if applicable
	amc_real supportradius;			///< Parameters for mesh movement - the support radius to be used, if applicable
	int nummovesteps;				///< Number of steps in which to accomplish the total mesh movement.
	string rbfsolver;				///< string describing the method to use for solving the RBF equations

	amc_int nbounpoin;						///< Number if boundary points.
	amc_int ninpoin;						///< Number of interior points.
	amat::Matrix<amc_real> disps;			///< Displacement of midpoint of each face
	amat::Matrix<amc_real> boundisps;		///< Displacement at each boundary point of the quadratic mesh, computed using [disps](@ref disps).
	amat::Matrix<amc_real> bounpoints;
	amat::Matrix<amc_real> inpoints;
	amat::Matrix<amc_int> bflagg;			///< This flag is true if the corresponding mesh node lies on a boundary.
	amat::Matrix<amc_int> toRec;			///< This flag is true if a boundary face is to be reconstructed.

public:
	void setup(UMesh* mesh, UMesh* meshq, Meshmove* mmove, int num_parts, std::vector<std::vector<int>> boundarymarkers, double angle_threshold, double toler, int maxitera, int rbf_choice, amc_real support_radius, int rbf_steps, string rbf_solver);

	~CurvedMeshGen();

	void compute_boundary_displacements();

	void generate_curved_mesh();
};

void CurvedMeshGen::setup(UMesh* mesh, UMesh* meshq, Meshmove* mmove, int num_parts, std::vector<std::vector<int>> boundarymarkers, double angle_threshold, double toler, int maxitera, int rbf_choice, amc_real support_radius, int rbf_steps, string rbf_solver)
{
	degree = 2;
	
	m = mesh;
	mq = meshq;
	mmv = mmove;
	br = new VertexCenteredBoundaryReconstruction(m, degree);
	tol = toler;
	maxiter = maxitera;
	rbfchoice = rbf_choice; supportradius = support_radius;
	nummovesteps = rbf_steps;
	rbfsolver = rbf_solver;
	disps.setup(m->gnface(),m->gndim());
	disps.zeros();
	bflagg.setup(mq->gnpoin(),1);

	// demarcate which faces are to be reconstructed
	toRec.setup(m->gnface(),1);
	toRec.zeros();
	for(amc_int iface = 0; iface < m->gnface(); iface++)
		for(int i = 0; i < boundarymarkers.size(); i++)
			for(int j = 0; j < boundarymarkers[i].size(); j++)
				if(m->gbface(iface,m->gnnofa()) == boundarymarkers[i][j])
					toRec(iface) = 1;
}

CurvedMeshGen::~CurvedMeshGen()
{
	delete br;
}

/// Computes displacement of midpoint of each face.
void CurvedMeshGen::compute_boundary_displacements()
{
	br.preprocess();
	br.solve();

	// get coords of midpoints of each boundary edge
	amat::Matrix<amc_real> edgemidpoints(m->gnface(),m->gndim());
	for(amc_int iface = 0; iface < m->gnface(); iface++)
		for(int idim = 0; idim < m->gndim(); idim++)
		{
			amc_real sum = 0;
			for(int inode = 0; inode < m->gnnofa(); inode++)
				sum += m->gcoords(m->gbface(iface,inode),idim);
			facemidpoints(iface,idim) = sum/m->gnnofa();
		}
	
	double uh = 0.5;
	for(amc_int iface = 0; iface < m->gnface(); iface++)
	{
		// first check if iface was reconstructed!
		if(toRec(iface))
			for(int idim = 0; idim < m->gndim(); idim++)
				disps(iface,idim) = br.getcoords(iface,idim,uh) - facemidpoints.get(iface,idim);
	}

	/// We do not need the linear mesh once we have the displacements of the faces' midpoints.
}

/** Uses the previously computed displacements of the face midpoints to curve the mesh.
*/
void CurvedMeshGen::generate_curved_mesh()
{
	/** 
	Note that this function works with the straight quadratic mesh.
	We assume that the face numberings of the linear mesh and the quadratic mesh are the same.
	*/
	
	/// Get a vector of displacements for each node of the quadratic mesh
	amat::Matrix<amc_real> allpoint_disps(mq->gnpoin(),mq->gndim());
	allpoint_disps.zeros();
	for(amc_int iface = 0; iface < mq->gnface(); iface++)
	{
		for(int idim = 0; idim < mq->gndim(); idim++)
			allpoint_disps(mq->gbface(iface,mq->gnnofa()-1), idim) = disps(iface,idim);
	}
	
	// first get bflag
	bflagg.zeros();
	for(amc_int iface = 0; iface < mq->gnface(); iface++)
	{
		for(int inode = 0; inode < mq->gnnofa(); inode++)
			bflagg(mq->gbface(iface,inode)) = 1;
	}

	nbounpoin = 0;
	for(amc_int i = 0; i < mq->gnpoin(); i++)
		nbounpoin += bflagg(i);

	ninpoin = mq->gnpoin()-nbounpoin;
	std::cout << "CurvedMeshGen: generate_curved_mesh(): Number of boundary points in quadratic mesh = " << nbounpoin << std::endl;
	std::cout << "CurvedMeshGen: generate_curved_mesh(): Number of interior points in quadratic mesh = " << ninpoin << std::endl;
	bounpoints.setup(nbounpoin,mq->gndim());
	boundisps.setup(nbounpoin,mq->gndim());
	inpoints.setup(ninpoin,mq->gndim());
	
	///We divide mesh nodes into boundary points and interior points. We also populate boundisp so that it holds the displacement of each boundary point.
	amc_int k = 0, l = 0;
	for(amc_int ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		if(bflagg(ipoin))
		{
			for(int idim = 0; idim < mq->gndim(); idim++){
				bounpoints(k,idim) = mq->gcoords(ipoin,idim);
				boundisps(k,idim) = allpoint_disps(ipoin,idim);
			}
			k++;
		}
		else
		{
			for(int idim = 0; idim < mq->gndim(); idim++)	
				inpoints(l,idim) = mq->gcoords(ipoin,idim);
			l++;
		}
	
	/// We now have all we need to call the mesh-movement functions and generate the curved mesh.
	//Call RBF functions here

	mmv->setup(&inpoints, &bounpoints, &boundisps, rbfchoice, supportradius, nummovesteps, tol, maxiter, rbfsolver);
	mmv->move();

	bounpoints = mmv->getBoundaryPoints();
	inpoints = mmv->getInteriorPoints();

	/// Finally, we reassemble the coord array for the curved mesh using the mesh-mover's computed values.
	//Get coord array of curved mesh
	amat::Matrix<amc_real> newcoords(mq->gnpoin(),mq->gndim());
	k = 0; l = 0;
	for(amc_int ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
	{
		if(bflagg(ipoin)) {
			for(int idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = bounpoints(k,idim);
			k++;
		}
		else {
			for(int idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = inpoints(l,idim);
			l++;
		}
	}

	// set it in mesh mq
	mq->setcoords(&newcoords);
}

// ------------ end --------------------
}
