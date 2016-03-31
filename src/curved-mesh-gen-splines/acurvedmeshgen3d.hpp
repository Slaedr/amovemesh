/** @brief Class to convert a 3D linear mesh into quadratic mesh.
 * @author Aditya Kashi
 * @date March 19, 2016
 * 
 * Notes:
 * We could set the support radius automatically based on some function of the curvature of the boundary.
 */

#ifndef __AMESH3D_H
	#include <amesh3d.hpp>
#endif

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
	const UMesh* m;					///< Data about the original linear mesh. We need this to compute spline reconstruction of the boundary.
	UMesh* mq;						///< Data of the corresponding (straight-faced) quadratic mesh
	Meshmove* mmv;					///< Pointer to parent class for the mesh-movement classes, such RBF, DGM or linear elasticity.
	BoundaryReconstruction* br;		///< Object to reconstruct the boundary using cubic splines.
	std::string brtype;				///< Type of boundary reconstruction, can be FACE or VERTEX
	int degree;						///< Degree of the generated mesh (2 or 3; currently only 2 is supported)
	
	double tol;						///< Tolerance for linear solver used for computing spline coefficients.
	int maxiter;					///< Maximum number of iterations for linear solver used to compute spline coefficients.
	int rbfchoice;					///< Parameters for mesh movement - the type of RBF to be used, if applicable
	amc_real supportradius;			///< Parameters for mesh movement - the support radius to be used, if applicable
	int nummovesteps;				///< Number of steps in which to accomplish the total mesh movement.
	std::string rbfsolver;				///< string describing the method to use for solving the RBF equations

	amc_int nbounpoin;						///< Number if boundary points.
	amc_int ninpoin;						///< Number of interior points.
	amat::Matrix<amc_real> disps;			///< Displacement of midpoint of each face
	amat::Matrix<amc_real> boundisps;		///< Displacement at each boundary point of the quadratic mesh, computed using [disps](@ref disps).
	amat::Matrix<amc_real> bounpoints;
	amat::Matrix<amc_real> inpoints;
	amat::Matrix<amc_int> bflagg;			///< This flag is true if the corresponding mesh node lies on a boundary.
	amat::Matrix<amc_int> toRec;			///< This flag is true if a boundary face is to be reconstructed.
	amat::Matrix<amc_real> allpoint_disps;	///< Initial displacements of all points in the high-order mesh; zero for interior points

public:
	void setup(const UMesh* mesh, UMesh* meshq, std::string br_type, double angle_threshold, double toler, int maxitera, int rbf_choice, amc_real support_radius, int rbf_steps, std::string rbf_solver);

	~CurvedMeshGen();

	void compute_boundary_displacements();

	void generate_curved_mesh();
};

void CurvedMeshGen::setup(const UMesh* mesh, UMesh* meshq, std::string br_type, double angle_threshold, double toler, int maxitera, int rbf_choice, amc_real support_radius, int rbf_steps, std::string rbf_solver)
{
	degree = 2;
	
	m = mesh;
	mq = meshq;

	brtype = br_type;
	if(br_type == "VERTEX")
		br = new VertexCenteredBoundaryReconstruction(m, degree, true, 1.0e2);
	else
		br = new FaceCenteredBoundaryReconstruction(m, degree, true, 1.0e2);

	tol = toler;
	maxiter = maxitera;
	rbfchoice = rbf_choice; supportradius = support_radius;
	nummovesteps = rbf_steps;
	rbfsolver = rbf_solver;
	disps.setup(m->gnface(),m->gndim());
	disps.zeros();
	bflagg.setup(mq->gnpoin(),1);
	allpoint_disps.setup(mq->gnpoin(),mq->gndim());

	// demarcate which faces are to be reconstructed
	/*toRec.setup(m->gnface(),1);
	toRec.zeros();
	for(amc_int iface = 0; iface < m->gnface(); iface++)
		for(int i = 0; i < boundarymarkers.size(); i++)
			for(int j = 0; j < boundarymarkers[i].size(); j++)
				if(m->gbface(iface,m->gnnofa()) == boundarymarkers[i][j])
					toRec(iface) = 1;*/
}

CurvedMeshGen::~CurvedMeshGen()
{
	delete br;
}

/// Computes displacement of midpoint of each face.
/** \note NOTE: currently only for tetrahedral elements!
 */
void CurvedMeshGen::compute_boundary_displacements()
{
	int i, inode, idim;
	amc_real sum;

	br->preprocess();
	br->solve();

	allpoint_disps.zeros();

	amat::Matrix<amc_real>* linearedgepoints;
	amat::Matrix<amc_real>* linearfacepoints;
	if(degree == 2)
	{
		linearedgepoints = new amat::Matrix<amc_real>[m->gnbedge()];
		for(i = 0; i < m->gnbedge(); i++)
			linearedgepoints[i].setup(1,m->gndim());
	}
	else if(degree == 3)
	{
		linearedgepoints = new amat::Matrix<amc_real>[m->gnbedge()];
		for(i = 0; i < m->gnbedge(); i++)
			linearedgepoints[i].setup(2,m->gndim());
		linearfacepoints = new amat::Matrix<amc_real>[m->gnface()];
		for(i = 0; i < m->gnface(); i++)
			linearfacepoints[i].setup(1,m->gndim());
	}

	// get coords of midpoints of each boundary edge and face, and their positions on the reconstructed surface
	
	std::vector<amc_real> recpoint(m->gndim());

	if(degree == 2)
	{
		double uh = 0.5;
		for(amc_int ied = 0; ied < m->gnbedge(); ied++)
		{
			br->getEdgePoint(uh,ied,recpoint);
			
			for(idim = 0; idim < m->gndim(); idim++)
			{
				linearedgepoints[ied](0,idim) = 0;
				
				for(inode = 0; inode < m->gnnoded(); inode++)
					linearedgepoints[ied](0,idim) += m->gcoords(m->gedgepo(ied,inode),idim);
				
				linearedgepoints[ied](0,idim) /= m->gnnoded();
			}

			for(idim = 0; idim < m->gndim(); idim++)
				allpoint_disps(mq->gedgepo(ied,2), idim) = recpoint[idim] - linearedgepoints[ied].get(0,idim);
		}
	}

	if(degree == 3)
	{
		std::vector<amc_real> areacoords(m->gndim());
		if(m->gnnofa() == 3)
		{
			areacoords[0] = areacoords[1] = areacoords[2] = 1.0/3.0;
		}

		std::vector<double> uh(2); uh[0] = 1.0/3.0; uh[1] = 2.0/3.0;

		for(amc_int ied = 0; ied < m->gnbedge(); ied++)
		{
			for(i = 0; i < 2; i++)
			{
				br->getEdgePoint(uh[i],ied,recpoint);
				for(idim = 0; idim < m->gndim(); idim++)
				{
					linearedgepoints[ied](i,idim) = 0;
					for(inode = 0; inode < m->gnnoded(); inode++)
						linearedgepoints[ied](i,idim) += m->gcoords(m->gedgepo(ied,inode),idim);
				}
				for(idim = 0; idim < m->gndim(); idim++)
				{
					linearedgepoints[ied](i,idim) /= m->gnnoded();
					// store displacement in linearedgepoints
					linearedgepoints[ied](i,idim) = recpoint[idim] - linearedgepoints[ied](i,idim);
				}
			
				for(idim = 0; idim < m->gndim(); idim++)
					allpoint_disps(mq->gedgepo(ied,2+i), idim) = linearedgepoints[ied].get(i,idim);
			}
		}
		for(amc_int iface = 0; iface < m->gnface(); iface++)
		{
			br->getFacePoint(areacoords, iface, recpoint);
			for(idim = 0; idim < m->gndim(); idim++)
			{
				sum = 0;
				for(inode = 0; inode < m->gnnofa(); inode++)
					sum += m->gcoords(m->gbface(iface,inode),idim);
				linearfacepoints[iface](0,idim) = recpoint[idim] - sum/m->gnnofa();
				allpoint_disps(mq->gbface(iface,9),idim) = linearfacepoints[iface](0,idim);
			}
		}
	}

	/*std::ofstream fout("testdisps.dat");
	for(i = 0; i < mq->gnpoin(); i++)
	{
		for(idim = 0; idim < m->gndim(); idim++)
			fout << allpoint_disps.get(i,idim) << " ";
		fout << '\n';
	}
	fout.close();*/
	
	/// We do not need the linear mesh once we have the displacements of the faces' midpoints.
	
	delete [] linearedgepoints;
	if(degree == 3)
		delete [] linearfacepoints;
}

/** Uses the previously computed displacements of the face midpoints to curve the mesh.
*/
void CurvedMeshGen::generate_curved_mesh()
{
	/** 
	Note that this function works with the straight quadratic mesh.
	We assume that the face numberings (bface) and edge numberings (intedge) of the linear mesh and the quadratic mesh are the same.
	*/

	amc_int ipoin;
	int idim;

	// check
	/*amat::Matrix<amc_real> diff(m->gnpoin(), m->gndim());
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		for(idim = 0; idim < m->gndim(); idim++)
			diff(ipoin,idim) = mq->gcoords(ipoin,idim) - m->gcoords(ipoin,idim);

	std::ofstream fout("testdiff.dat");
	fout << std::setprecision(14);
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		fout << diff.get(ipoin,0) << " " << diff.get(ipoin,1) << " " << diff.get(ipoin,2) << '\n';
	fout.close();

	std::vector<amc_real> err(3,0);
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		if(m->gbpointsinv(ipoin) >= 0)
			for(idim = 0; idim < 3; idim++)
				err[idim] += diff.get(ipoin,idim)*diff.get(ipoin,idim);
	for(idim = 0; idim < 3; idim++)
		err[idim] = sqrt(err[idim]);
	std::cout << std::setprecision(14);
	std::cout << "acmg3d: Initial error in positions of low-order nodes: " << err[0] << " " << err[1] << " " << err[2] << std::endl;*/


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
	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
		if(bflagg(ipoin))
		{
			for(idim = 0; idim < mq->gndim(); idim++){
				bounpoints(k,idim) = mq->gcoords(ipoin,idim);
				boundisps(k,idim) = allpoint_disps(ipoin,idim);
			}
			k++;
		}
		else
		{
			for(idim = 0; idim < mq->gndim(); idim++)	
				inpoints(l,idim) = mq->gcoords(ipoin,idim);
			l++;
		}
	
	/// We now have all we need to call the mesh-movement functions and generate the curved mesh.
	//Call RBF functions here

	mmv = new RBFmove(&inpoints, &bounpoints, &boundisps, rbfchoice, supportradius, nummovesteps, tol, maxiter, rbfsolver);
	mmv->move();

	bounpoints = mmv->getBoundaryPoints();
	inpoints = mmv->getInteriorPoints();

	/// Finally, we reassemble the coord array for the curved mesh using the mesh-mover's computed values.
	//Get coord array of curved mesh
	amat::Matrix<amc_real> newcoords(mq->gnpoin(),mq->gndim());
	k = 0; l = 0;
	for(ipoin = 0; ipoin < mq->gnpoin(); ipoin++)
	{
		if(bflagg(ipoin)) {
			for(idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = bounpoints(k,idim);
			k++;
		}
		else {
			for(idim = 0; idim < mq->gndim(); idim++)
				newcoords(ipoin,idim) = inpoints(l,idim);
			l++;
		}
	}

	// set it in mesh mq
	mq->setcoords(&newcoords);

	/*for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		for(idim = 0; idim < m->gndim(); idim++)
			diff(ipoin,idim) = mq->gcoords(ipoin,idim) - m->gcoords(ipoin,idim);

	fout.open("testdiff.dat");
	fout << std::setprecision(14);
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		fout << diff.get(ipoin,0) << " " << diff.get(ipoin,1) << " " << diff.get(ipoin,2) << '\n';
	fout.close();

	err.assign(3,0);
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		if(m->gbpointsinv(ipoin) >= 0)
			for(idim = 0; idim < 3; idim++)
				err[idim] += diff.get(ipoin,idim)*diff.get(ipoin,idim);
	for(idim = 0; idim < 3; idim++)
		err[idim] = sqrt(err[idim]);
	std::cout << std::setprecision(14);
	std::cout << "acmg3d: error in positions of low-order nodes: " << err[0] << " " << err[1] << " " << err[2] << std::endl;*/
	
	delete mmv;
}

// ------------ end --------------------
}
