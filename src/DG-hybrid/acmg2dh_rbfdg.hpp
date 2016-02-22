/** @brief Class to convert hybrid linear mesh into quadratic mesh using hybrid RBF-DG method (not to be confused with DG-RBF)
 * @author Aditya Kashi
 * @date November 2, 2015
 * 
 * Notes:
 * \todo We could set the support radius automatically based on some function of (1) the curvature of the boundary and/or (2) shape of elements.
 */

#ifndef __AGEOMETRYH_H
	#include <ageometryh.hpp>
#endif

#ifndef __ARBFDGHYBRID_H
	#include "arbfdghybrid.hpp"
#endif

#define __ACMG2DH_RBFDG_H 1

using namespace std;
using namespace amat;
using namespace amc;

/** \brief Class to generate curved mesh from a linear mesh using cubic spline reconstruction and RBF-DGM hybrid mesh movement technique. */
class Curvedmeshgen2d
{
	UMesh2dh* m;				///< Data about the original linear mesh. We need this to compute spline reconstruction of the boundary.
	UMesh2dh* mq;					///< Data of the corresponding (straight-faced) quadratic mesh
	DGhybrid* mmv;					///< the mesh-movement class DG hybrid
	BoundaryReconstruction2d br;	///< Object to reconstruct the boundary using cubic splines.
	
	double spltol;						///< Tolerance for linear solver used for computing spline coefficients.
	double splmaxiter;					///< Maximum number of iterations for linear solver used to compute spline coefficients.
	int rbfchoice;					///< Parameters for mesh movement - the type of RBF to be used, if applicable
	double supportradius;			///< Parameters for mesh movement - the support radius to be used, if applicable
	int nummovesteps;				///< Number of steps in which to accomplish the total mesh movement.
	string rbfsolver;				///< string describing the method to use for solving the RBF equations
	double mmtol;					///< Tolerance to use while solving for mesh movement
	int mmmaxiter;					///< Maximum iterations for mesh movement solver
	string mmsolver;				///< String describing the solver to use for mesh movement - "LU" or "CG"
	double suprad;					///< Support radius for RBF mesh movement
	int nlayers;					///< Number of layers to advance from the boundary in order to select points for the background mesh in the DG-hybrid method
	bool allocmove;					///< Stores whether the [mesh movement context](@ref mmv) has been allocated

	int nbounpoin;					///< Number if boundary points.
	int ninpoin;					///< Number of interior points.
	Matrix<double> disps;			///< Displacement of midpoint of each face
	Matrix<double> boundisps;		///< Displacement at each boundary point of the quadratic mesh, computed using [disps](@ref disps).
	Matrix<double> bounpoints;
	Matrix<double> inpoints;
	Matrix<int> bflagg;				///< This flag is true if the corresponding mesh node lies on a boundary.
	Matrix<int> toRec;				///< This flag is true if a boundary face is to be reconstructed.

public:
	Curvedmeshgen2d()
	{
		allocmove = false;
	}

	~Curvedmeshgen2d()
	{
		if(allocmove)
			delete mmv;
	}

	void setup(UMesh2dh* const mesh, UMesh2dh* const meshq, const int num_parts, const vector<vector<int>> boundarymarkers, const double angle_threshold, 
			const double spl_tol, const int spl_maxiter, const double le_toler, const int le_maxiter, const string le_solver, const double support_radius, const int num_layers);

	void compute_boundary_displacements();

	void generate_curved_mesh();
};

void Curvedmeshgen2d::setup(UMesh2dh* const mesh, UMesh2dh* const meshq, const int num_parts, const vector<vector<int>> boundarymarkers, const double angle_threshold, 
		const double spl_tol, const int spl_maxiter, const double le_toler, const int le_maxiter, const string le_solver, const double supp_rad, const int num_layers)
{
	m = mesh;
	mq = meshq;
	br.setup(m, num_parts, boundarymarkers, angle_threshold);
	spltol = spl_tol;
	splmaxiter = spl_maxiter;
	mmtol = le_toler;
	mmmaxiter = le_maxiter;
	mmsolver = le_solver;
	suprad = supp_rad;
	nlayers = num_layers;
	
	disps.setup(m->gnface(),m->gndim());
	disps.zeros();
	bflagg.setup(mq->gnpoin(),1);

	mmv = new DGhybrid();
	allocmove = true;
	
	// demarcate which faces are to be reconstructed
	toRec.setup(m->gnface(),1);
	toRec.zeros();
	for(int iface = 0; iface < m->gnface(); iface++)
		for(int i = 0; i < boundarymarkers.size(); i++)
			for(int j = 0; j < boundarymarkers[i].size(); j++)
				if(m->gbface(iface,m->gnnofa()) == boundarymarkers[i][j])
					toRec(iface) = 1;
}

/// Computes displacement of midpoint of each face.
void Curvedmeshgen2d::compute_boundary_displacements()
{
	br.preprocess();
	br.detect_corners();
	br.split_parts();
	br.compute_splines(spltol,splmaxiter);

	// get coords of midpoints of each boundary face
	Matrix<double> facemidpoints(m->gnface(),m->gndim());
	for(int iface = 0; iface < m->gnface(); iface++)
		for(int idim = 0; idim < m->gndim(); idim++)
		{
			double sum = 0;
			for(int inode = 0; inode < m->gnnofa(); inode++)
				sum += m->gcoords(m->gbface(iface,inode),idim);
			facemidpoints(iface,idim) = sum/m->gnnofa();
		}
	
	double uh = 0.5;																	// for quadratic mesh

	for(int iface = 0; iface < m->gnface(); iface++)
	{
		// first check if iface was reconstructed!
		if(toRec(iface))
			for(int idim = 0; idim < m->gndim(); idim++)
				disps(iface,idim) = br.getcoords(iface,idim,uh) - facemidpoints.get(iface,idim);
	}

	/// We do not need the linear mesh once we have the displacements of the faces' midpoints.
}

/// Uses the previously computed displacements of the face midpoints to curve the mesh
void Curvedmeshgen2d::generate_curved_mesh()
{
	/** Note that this function works with the straight quadratic mesh.
	 * We assume that the face numberings of the linear mesh and the quadratic mesh are the same.
	 */
	
	/// Get a vector of displacements for each node of the quadratic mesh
	Matrix<double> allpoint_disps(mq->gnpoin(),mq->gndim());
	allpoint_disps.zeros();
	for(int iface = 0; iface < mq->gnface(); iface++)
	{
		for(int idim = 0; idim < mq->gndim(); idim++)
			allpoint_disps(mq->gbface(iface,mq->gnnofa()-1), idim) = disps(iface,idim);
	}
	
	// first get bflag
	bflagg.zeros();
	for(int iface = 0; iface < mq->gnface(); iface++)
	{
		for(int inode = 0; inode < mq->gnnofa(); inode++)
			bflagg(mq->gbface(iface,inode)) = 1;
	}

	nbounpoin = 0;
	for(int i = 0; i < mq->gnpoin(); i++)
		nbounpoin += bflagg(i);

	ninpoin = mq->gnpoin()-nbounpoin;
	cout << "Curvedmeshgen2d: generate_curved_mesh(): Number of boundary points in quadratic mesh = " << nbounpoin << endl;
	cout << "Curvedmeshgen2d: generate_curved_mesh(): Number of interior points in quadratic mesh = " << ninpoin << endl;
	bounpoints.setup(nbounpoin,mq->gndim());
	boundisps.setup(nbounpoin,mq->gndim());
	inpoints.setup(ninpoin,mq->gndim());
	
	///We divide mesh nodes into boundary points and interior points. We also populate boundisp so that it holds the displacement of each boundary point.
	int k = 0, l = 0;
	for(int ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
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

	mmv->setup(m, mq, &allpoint_disps, nlayers, suprad, mmtol, mmmaxiter, mmsolver);
	mmv->compute_backmesh_points();
	mmv->generate_backmesh_and_compute_displacements();
	mmv->movemesh();

	/*bounpoints = mmv->getBoundaryPoints();
	inpoints = mmv->getInteriorPoints();

	/// Finally, we reassemble the coord array for the curved mesh using the mesh-mover's computed values.
	//Get coord array of curved mesh
	Matrix<double> newcoords(mq->gnpoin(),mq->gndim());
	k = 0; l = 0;
	for(int ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
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
	mq->setcoords(&newcoords);*/
}

// ------------ end --------------------
