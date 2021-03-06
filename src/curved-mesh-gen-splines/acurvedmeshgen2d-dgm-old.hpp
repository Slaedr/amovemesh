/** @brief Class to convert linear mesh into quadratic curved mesh using Delaunay graph mapping, by triangulating only high-order boundary elements.
 * @author Aditya Kashi
 * @date January 12, 2015
*/

#ifndef __AGEOMETRY_H
	#include <ageometry.hpp>
#endif

#ifndef __ADGM_H
	#include <adgm.hpp>
#endif

#define __ACURVEDMESHGEN2D_H 1

using namespace std;
using namespace amat;
using namespace acfd;

typedef DGmove Meshmove;

/** Class to generate curved mesh from a linear mesh using cubic spline reconstruction and one of the mesh movement techniques. */

class Curvedmeshgen2d
{
	UMesh2d* m;						///< Data about the original linear mesh. We need this to compute spline reconstruction of the boundary.
	UMesh2d* mq;					///< Data of the corresponding (straight-faced) quadratic mesh
	Meshmove* mmv;					///< Pointer to parent class for the mesh-movement classes, such RBF, DGM or linear elasticity.
	BoundaryReconstruction2d br;	///< Object to reconstruct the boundary using cubic splines.
	
	double tol;						///< Tolerance for linear solver used for computing spline coefficients.
	double maxiter;					///< Maximum number of iterations for linear solver used to compute spline coefficients.

	int nbounpoin;					///< Number if boundary points.
	int ninpoin;					///< Number of interior points.
	Matrix<double> disps;			///< Displacement of midpoint of each face
	Matrix<double> boundisps;		///< Displacement at each boundary point of the quadratic mesh, computed using [disps](@ref disps).
	Matrix<double> bounpoints;
	Matrix<double> inpoints;
	Matrix<int> bflagg;				///< This flag is true if the corresponding mesh node lies on a boundary.
	Matrix<int> hflagg;				///< This flag is true if the corresponding node is a high-order node.
	Matrix<int> toRec;				///< This flag is true if a boundary face is to be reconstructed.

public:
	void setup(UMesh2d* mesh, UMesh2d* meshq, Meshmove* mmove, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, double _tol, int _maxiter);

	void compute_boundary_displacements();

	void generate_curved_mesh();
};

/** tol and maxiter are for the spline computation.
 */
void Curvedmeshgen2d::setup(UMesh2d* mesh, UMesh2d* meshq, Meshmove* mmove, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, double _tol, int _maxiter)
{
	m = mesh;
	mq = meshq;
	mmv = mmove;
	br.setup(m, num_parts, boundarymarkers, angle_threshold);
	tol = _tol;
	maxiter = _maxiter;
	
	disps.setup(m->gnface(),m->gndim());
	disps.zeros();
	bflagg.setup(mq->gnpoin(),1);
	hflagg.setup(mq->gnpoin(),1);
	
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
	br.compute_splines(tol,maxiter);

	// get coords of midpoints of each face
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

/// Uses the previously computed displacements of the face midpoints to curve the mesh.
void Curvedmeshgen2d::generate_curved_mesh()
{
	/** 
	This function works with the straight quadratic mesh.
	\note We assume that the face numberings of the linear mesh and the quadratic mesh are the same.
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

	// differentiate high-order nodes from others. In 2D, number of vertices always equals the number of faces.
	hflagg.zeros();
	for(int i = 0; i < mq->gnelem(); i++)
	{
		for(int jnode = mq->gnfael(); jnode < mq->gnnode(); jnode++)
			hflagg(mq->ginpoel(i,jnode)) = 1;
	}

	nbounpoin = 0;
	for(int i = 0; i < mq->gnpoin(); i++)
		nbounpoin += bflagg(i);

	int ndgpoin;		// number of points to be triangulated
	ndgpoin = 0;
	for(int i = 0; i < mq->gnpoin(); i++)
		if(hflagg(i) && bflagg(i))
			ndgpoin++;

	ninpoin = mq->gnpoin()-nbounpoin;
	cout << "Curvedmeshgen2d: generate_curved_mesh(): Number of boundary points in quadratic mesh = " << nbounpoin << endl;
	cout << "Curvedmeshgen2d: generate_curved_mesh(): Number of boundary high order points in quadratic mesh = " << ndgpoin << endl;
	cout << "Curvedmeshgen2d: generate_curved_mesh(): Number of interior points in quadratic mesh = " << ninpoin << endl;
	bounpoints.setup(ndgpoin,mq->gndim());
	boundisps.setup(ndgpoin,mq->gndim());
	inpoints.setup(ninpoin,mq->gndim());
	
	/// We divide mesh nodes into high-order boundary points and interior points. We also populate boundisp so that it holds the displacement of each high-order boundary point.
	int k = 0, l = 0;
	for(int ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		if(bflagg(ipoin) && hflagg(ipoin))
		{
			for(int idim = 0; idim < mq->gndim(); idim++){
				bounpoints(k,idim) = mq->gcoords(ipoin,idim);
				boundisps(k,idim) = allpoint_disps(ipoin,idim);
			}
			k++;
		}
		else if(!bflagg(ipoin))		// regular (non-high-order) boundary nodes are not moved
		{
			for(int idim = 0; idim < mq->gndim(); idim++)	
				inpoints(l,idim) = mq->gcoords(ipoin,idim);
			l++;
		}
	
	/// We now have all we need to call the mesh-movement functions and generate the curved mesh.

	mmv->setup(2, &inpoints, &bounpoints, &boundisps);
	mmv->move();
	mmv->dg.writeGmsh2("dg.msh");

	bounpoints = mmv->getBoundaryPoints();
	inpoints = mmv->getInteriorPoints();

	/// Finally, we reassemble the coord array for the curved mesh using the mesh-mover's computed values.
	//Get coord array of curved mesh
	Matrix<double> newcoords(mq->gnpoin(),mq->gndim());
	newcoords = *(mq->getcoords());
	/*for(int i = 0; i < mq.gnpoin(); i++)
		for(int j = 0; j < mq->gndim(); j++)
			newcoords(i,j) = mq->gcoords(i,j);*/

	k = 0; l = 0;
	for(int ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
	{
		if(bflagg(ipoin) && hflagg(ipoin)) {
			for(int idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = bounpoints(k,idim);
			k++;
		}
		else if(!bflagg(ipoin))
		{
			for(int idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = inpoints(l,idim);
			l++;
		}
	}

	// set it in mesh mq
	mq->setcoords(&newcoords);
}

// ------------ end --------------------
