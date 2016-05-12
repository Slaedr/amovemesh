/** @brief Class to convert hybrid linear mesh into quadratic mesh.
 * @author Aditya Kashi
 * @date November 2, 2015
 * 
 * Notes:
 * \todo We could set the support radius automatically based on some function of the curvature of the boundary.
 */

#ifndef __AGEOMETRYH_H
	#include <ageometryh.hpp>
#endif

#ifndef __ABOUNDARYINFLUENCEDISTANCE_H
	#include <aboundaryinfluencedistance.hpp>
#endif

#ifndef __ARBF_SR_H
	#include <arbf_sr.hpp>
#endif

#define __ACURVEDMESHGEN2DH_SR_H 1

using namespace std;
using namespace amat;
using namespace amc;

typedef RBFmove Meshmove;

/** Class to generate curved mesh from a linear mesh using cubic spline reconstruction and one of the mesh movement techniques. */

class Curvedmeshgen2d
{
	UMesh2dh* m;						///< Data about the original linear mesh. We need this to compute spline reconstruction of the boundary.
	UMesh2dh* mq;					///< Data of the corresponding (straight-faced) quadratic mesh
	Meshmove* mmv;					///< Pointer to parent class for the mesh-movement classes, such RBF, DGM or linear elasticity.
	BoundaryReconstruction2d br;	///< Object to reconstruct the boundary using cubic splines.
	double spltol;					///< Tolerance for spline solver
	int splmaxiter;					///< Max iterations for spline solver
	
	double tol;							///< Tolerance for linear solver used for computing spline coefficients.
	int maxiter;						///< Maximum number of iterations for linear solver used to compute spline coefficients.
	int rbfchoice;						///< Parameters for mesh movement - the type of RBF to be used, if applicable
	Matrix<amc_real> supportradius;	///< Support radius to be used for each point in the quadratic mesh (only values for boundary nodes are used)
	int nummovesteps;					///< Number of steps in which to accomplish the total mesh movement.
	string rbfsolver;					///< string describing the method to use for solving the RBF equations

	int nbounpoin;					///< Number if boundary points.
	int ninpoin;					///< Number of interior points.
	Matrix<double> disps;			///< Displacement of midpoint of each face
	Matrix<double> boundisps;		///< Displacement at each boundary point of the quadratic mesh, computed using [disps](@ref disps).
	Matrix<double> bounpoints;
	Matrix<double> inpoints;
	Matrix<int> bflagg;				///< This flag is true if the corresponding mesh node lies on a boundary.
	Matrix<int> toRec;				///< This flag is true if a boundary face is to be reconstructed.

public:
	void setup(UMesh2dh* mesh, UMesh2dh* meshq, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, 
			double _spltol, int _splmaxiter, double toler, int maxitera, int rbf_choice, int rbf_steps, string rbf_solver);

	void compute_boundary_displacements();

	void generate_curved_mesh();
};

void Curvedmeshgen2d::setup(UMesh2dh* mesh, UMesh2dh* meshq, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, double _spltol, int _splmaxiter, double toler, int maxitera, 
		int rbf_choice, int rbf_steps, string rbf_solver)
{
	m = mesh;
	mq = meshq;
	br.setup(m, num_parts, boundarymarkers, angle_threshold);
	spltol = _spltol;
	splmaxiter = _splmaxiter;
	tol = toler;
	maxiter = maxitera;
	rbfchoice = rbf_choice;
	supportradius.setup(mq->gnpoin(),1);
	nummovesteps = rbf_steps;
	rbfsolver = rbf_solver;
	disps.setup(m->gnface(),m->gndim());
	disps.zeros();
	bflagg.setup(mq->gnpoin(),1);
	
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
	
	double uh = 0.5;
	for(int iface = 0; iface < m->gnface(); iface++)
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
void Curvedmeshgen2d::generate_curved_mesh()
{
	/** 
	Note that this function works with the straight quadratic mesh.
	We assume that the face numberings of the linear mesh and the quadratic mesh are the same.
	*/
	
	// Get scaling parameter for mesh
	/*amc_real scale, mindist = 1.0, dist;
	int idim;
	for(int iface = 0; iface < mq->gnface(); iface++)
	{
		dist = 0;
		for(idim = 0; idim < 2; idim++)
			dist += (mq->gcoords(mq->gbface(iface,1),idim)-mq->gcoords(mq->gbface(iface,0),idim)) * (mq->gcoords(mq->gbface(iface,1),idim)-mq->gcoords(mq->gbface(iface,0),idim));
		dist = sqrt(dist);
		if(dist < mindist) mindist = dist;
	}
	std::cout << "Curvedmeshgen2d: generate_curved_mesh(): Minimum distance between 2 boundary points is " << mindist << std::endl;
	scale = 1.0/pow(mindist,0.5);*/
	
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

	// compute support radii
	boundaryInfluenceDist2D(mq, &allpoint_disps, &supportradius);

	nbounpoin = 0;
	for(int i = 0; i < mq->gnpoin(); i++)
		nbounpoin += bflagg(i);

	ninpoin = mq->gnpoin()-nbounpoin;
	cout << "Curvedmeshgen2d: generate_curved_mesh(): Number of boundary points in quadratic mesh = " << nbounpoin << endl;
	cout << "Curvedmeshgen2d: generate_curved_mesh(): Number of interior points in quadratic mesh = " << ninpoin << endl;
	bounpoints.setup(nbounpoin,mq->gndim());
	boundisps.setup(nbounpoin,mq->gndim());
	inpoints.setup(ninpoin,mq->gndim());
	Matrix<amc_real> srad(nbounpoin,1);
	
	///We divide mesh nodes into boundary points and interior points. We also populate boundisp so that it holds the displacement of each boundary point.
	int k = 0, l = 0, idim;
	amc_int ipoin;
	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		if(bflagg(ipoin))
		{
			for(idim = 0; idim < mq->gndim(); idim++){
				bounpoints(k,idim) = mq->gcoords(ipoin,idim);
				boundisps(k,idim) = allpoint_disps(ipoin,idim);
				srad(k) = supportradius.get(ipoin);
			}
			k++;
		}
		else
		{
			for(idim = 0; idim < mq->gndim(); idim++)	
				inpoints(l,idim) = mq->gcoords(ipoin,idim);
			l++;
		}
	
	// get avg support radius (for debug)
	amc_real ssum = 0, smin = 1.0;
	for(ipoin = 0; ipoin < nbounpoin; ipoin++)
		ssum += srad.get(ipoin);
	ssum /= nbounpoin;
	cout << "CurvedMeshGen2d: generate_curved_mesh(): Average support radius = " << ssum << endl;
	
	for(ipoin = 0; ipoin < nbounpoin; ipoin++)
		if(smin > srad.get(ipoin))
			smin = srad.get(ipoin);
	cout << "CurvedMeshGen2d: generate_curved_mesh(): Min support radius = " << smin << endl;
	
	for(ipoin = 0; ipoin < nbounpoin; ipoin++)
		srad(ipoin) = ssum;
	
	/*// before calling RBF, scale everything
	for(int ipoin = 0; ipoin < nbounpoin; ipoin++)
		for(idim = 0; idim < 2; idim++)
		{
			bounpoints(ipoin,idim) *= scale;
			boundisps(ipoin,idim) *= scale;
		}
	for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		for(idim = 0; idim < 2; idim++)
			inpoints(ipoin,idim) *= scale;*/
	
	/// We now have all we need to call the mesh-movement functions and generate the curved mesh.
	//Call RBF functions here

	mmv = new RBFmove(&inpoints, &bounpoints, &boundisps, rbfchoice, &srad, nummovesteps, tol, maxiter, rbfsolver);
	mmv->move();

	/*bounpoints = mmv->getBoundaryPoints();
	inpoints = mmv->getInteriorPoints();*/
	
	/*// scale coords back down
	for(int ipoin = 0; ipoin < nbounpoin; ipoin++)
		for(idim = 0; idim < 2; idim++)
		{
			bounpoints(ipoin,idim) /= scale;
			boundisps(ipoin,idim) /= scale;
		}
	for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		for(idim = 0; idim < 2; idim++)
			inpoints(ipoin,idim) /= scale;*/

	/// Finally, we reassemble the coord array for the curved mesh using the mesh-mover's computed values.
	Matrix<double> newcoords(mq->gnpoin(),mq->gndim());
	k = 0; l = 0;
	for(int ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
	{
		if(bflagg(ipoin)) {
			for(int idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = bounpoints.get(k,idim);
			k++;
		}
		else {
			for(int idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = inpoints.get(l,idim);
			l++;
		}
	}

	// set it in mesh mq
	mq->setcoords(&newcoords);

	delete mmv;
}

// ------------ end --------------------
