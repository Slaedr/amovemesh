/** @brief Class to convert linear mesh into quadratic mesh using P2 linear elasticity.
 * @author Aditya Kashi
 * @date October 29, 2015
 */

#ifndef __AGEOMETRY_H
	#include <ageometry.hpp>
#endif

#ifndef __ALINALG_H
	#include <alinalg.hpp>
#endif

#ifndef __ALINELAST_P2_STIFFENED_H
	#include <alinelast_p2_stiffened.hpp>
#endif

#define __ALINELAST_CURVEDMESHGEN2D_H 1

using namespace std;
using namespace amat;
using namespace acfd;


/** Class to generate curved mesh from a linear mesh using cubic spline reconstruction and one of the mesh movement techniques. */

class Curvedmeshgen2d
{
	UMesh2d* m;						///< Data about the original linear mesh. We need this to compute spline reconstruction of the boundary.
	UMesh2d* mq;					///< Data of the corresponding (straight-faced) quadratic mesh
	LinElastP2* mmv;					///< Pointer to parent class for the mesh-movement classes, such RBF, DGM or linear elasticity.
	BoundaryReconstruction2d br;	///< Object to reconstruct the boundary using cubic splines.
	
	double tol_g;					///< Tolerance for linear solver of spline solver
	double tol_e;					///< Tolerance for linear elasticity solver
	double maxiter;					///< Maximum number of iterations for linear solvers.
	double young;					///< Young's modulus
	double nu;						///< Poisson's ratio
	double lambda;
	double mu;
	double chi;						///< stiffening exponent for stiffened linear elasticity
	string stiffscheme;				///< stiffening scheme used for stiffened linear elasticity

	int nbounpoin;					///< Number if boundary points.
	int ninpoin;					///< Number of interior points.
	Matrix<double> disps;			///< Displacement of midpoint of each face
	Matrix<double> boundisps;		///< Displacement at each boundary point of the quadratic mesh, computed using [disps](@ref disps).
	Matrix<double> bounpoints;
	Matrix<double> inpoints;
	Matrix<int> bflagg;				///< This flag is true if the corresponding mesh node lies on a boundary.
	Matrix<int> toRec;				///< This flag is true if a boundary face is to be reconstructed.

	/// Stiffening factors
	Matrix<double> stiff;

public:
	void setup(UMesh2d* mesh, UMesh2d* meshq, LinElastP2* mmove, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, double tolg, double tole, double maxitera, double youngsmodulus, double poissonsratio, double xchi, string sscheme); 

	void compute_boundary_displacements();

	void generate_curved_mesh();
};

void Curvedmeshgen2d::setup(UMesh2d* mesh, UMesh2d* meshq, LinElastP2* mmove, int num_parts, vector<vector<int>> boundarymarkers, double angle_threshold, double toler1, double toler2, double maxitera, double youngsmodulus, double poissonsratio, double xchi, string sscheme)
{
	m = mesh;
	mq = meshq;
	mmv = mmove;
	br.setup(m, num_parts, boundarymarkers, angle_threshold);
	tol_g = toler1;
	tol_e = toler2;
	maxiter = maxitera;
	young = youngsmodulus;
	nu = poissonsratio;
	chi = xchi;
	stiffscheme = sscheme;
	lambda = nu*young/((1+nu)*(1-2*nu));
	mu = young/(2*(1+nu));
		
	// stiffening factor
	stiff.setup(m->gnelem(), 1);
	m->compute_jacobians();
	double j0 = 0;
	for(int i = 0; i < m->gnelem(); i++)
		j0 += m->gjacobians(i);
	j0 /= m->gnelem();

	if(stiffscheme == "size")
	{
		for(int iel = 0; iel < m->gnelem(); iel++)
			stiff(iel) = pow(j0/m->gjacobians(iel),chi);
	}
	else
	{
		m->compute_metric_quantities();
		if(stiffscheme == "shape")
		{
			m->linearmetric_shape(&stiff);
			
			// if we have zero metric, we're toast
			//for(int iel = 0; iel < m->gnelem(); iel++)
			//	if(stiff.get(iel) < ZERO_TOL) stiff(iel) = ZERO_TOL;

			for(int iel = 0; iel < m->gnelem(); iel++)
				stiff(iel) = pow(1.0/stiff.get(iel),chi);
		}
		else if(stiffscheme == "shapesize")
		{
			m->linearmetric_shapesize(&stiff);
			for(int iel = 0; iel < m->gnelem(); iel++)
				stiff(iel) = pow(1.0/stiff.get(iel),chi);
		}
	}

	disps.setup(m->gnface(),m->gndim());
	disps.zeros();
	
	// compute boundary boolean vector for the quadratic mesh
	bflagg.setup(mq->gnpoin(),1);
	bflagg.zeros();
	for(int i = 0; i < mq->gnface(); i++)
		for(int j = 0; j < mq->gnnofa(); j++)
			bflagg(mq->gbface(i,j)) = 1;
	
	// demarcate which faces are to be reconstructed
	toRec.setup(m->gnface(),1);
	toRec.zeros();
	for(int iface = 0; iface < m->gnface(); iface++)
		for(int i = 0; i < boundarymarkers.size(); i++)
			for(int j = 0; j < boundarymarkers[i].size(); j++)
				if(m->gbface(iface,m->gnnofa()) == boundarymarkers[i][j])
					toRec(iface) = 1;
}

/** Computes displacement of midpoint of each face. */
void Curvedmeshgen2d::compute_boundary_displacements()
{
	br.preprocess();
	br.detect_corners();
	br.split_parts();
	br.compute_splines(tol_g,maxiter);

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

/** Uses the previously computed displacements of the face midpoints to curve the mesh.
*/
void Curvedmeshgen2d::generate_curved_mesh()
{
	/** 
	Note that this function works with the straight quadratic mesh.
	We assume that the face numberings of the linear mesh and the quadratic mesh are the same.
	*/
	
	/// Get a vector of displacements for each node of the quadratic mesh
	Matrix<double> allpoint_disps(mq->gnpoin(),mq->gndim());
	allpoint_disps.zeros();
	for(int iface = 0; iface < mq->gnface(); iface++)
	{
		for(int idim = 0; idim < mq->gndim(); idim++)
			allpoint_disps(mq->gbface(iface,mq->gnnofa()-1), idim) = disps(iface,idim);
	}

	mmv->setup(mq, mu, lambda, chi, &stiff);
	cout << "Cuvedmeshgen2d: generate_curved_mesh(): Assembling stiffness matrix and load vector." << endl;
	mmv->assembleStiffnessMatrix();
	mmv->assembleLoadVector();
	
	cout << "Cuvedmeshgen2d: generate_curved_mesh(): Applying Dirichlet BCs." << endl;
	mmv->dirichletBC_onAllBface(allpoint_disps, bflagg);

	Matrix<double> alldisps(mq->gnpoin()*2,1);

	SpMatrix A = mmv->stiffnessMatrix();
	Matrix<double> b = mmv->loadVector();

	Matrix<double> xold(2*mq->gnpoin(),1);
	xold.zeros();

	/*ofstream fout("matrix.dat");
	b.fprint(fout);
	fout.close();*/

	alldisps = sparseCG_d(&A, b, xold, tol_e, maxiter);

	/*ofstream ofile("disps.dat");
	for(int i = 0; i < mq->gnpoin(); i++)
	{
		if(dabs(alldisps.get(i,0)) > ZERO_TOL || dabs(alldisps.get(mq->gnpoin()+i,0) > ZERO_TOL))
			ofile << i << " " << alldisps(i,0) << " " << alldisps(mq->gnpoin()+i,0) << '\n';
	}
	ofile.close();*/

	Matrix<double> fpos(mq->gnpoin(),mq->gndim());
	for(int i = 0; i < mq->gnpoin(); i++)
		for(int j = 0; j < mq->gndim(); j++)
			fpos(i,j) = mq->gcoords(i,j);

	cout << "Curvedmeshgen2d: generate_curved_mesh(): Updating position." << endl;
	// update poistions
	for(int j = 0; j < mq->gndim(); j++)
	{
		for(int i = 0; i < mq->gnpoin(); i++)
			fpos(i,j) += alldisps(j*mq->gnpoin() + i);
	}
	//cout << "Curvedmeshgen2d: generate_curved_mesh(): Positions updated. Updating mesh." << endl;
	//update mesh
	mq->setcoords(&fpos);
	
}

// ------------ end --------------------
