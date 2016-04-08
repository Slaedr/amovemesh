/** \brief Generation of 3D curved meshes when an analytical expression for the boundary is available.
 * \author Aditya Kashi
 */

#ifndef __AMESH3D_H
#include <amesh3d.hpp>
#endif

#ifndef __ARBF_H
#include <arbf.hpp>
#endif

using namespace std;
using namespace amat;
using namespace amc;

/// Generation of 3D curved meshes when an analytical expression for the boundary is available
/** RBF mesh movement scheme is used.
 * RBF support radius, tolerance and maximum iterations for solver need to be provided.
 */
class CurvedMeshGeneration
{
	UMesh* m;
	/// interior points to be moved
	Matrix<double> mipoints;
	/// fixed interior points
	Matrix<double> fipoints;
	/// Boundary points
	Matrix<double> bpoints;
	/// Boundary motion
	Matrix<double> bmotion;
	/// Flag indicating whether or not there are internal (non-boundary) fixed points
	bool noInternalFixedPoints;
	/// number of interior points
	int ninpoin;
	/// number of interior points to be moved
	int npomo;
	/// number of  boundary points
	int nbpoin;
	/// npoin-by-1, 1 if corresponding point is a boundary points
	Matrix<int> bflag;
	/// 1 if corresponding point is a high-order point (npoin-by-1)
	Matrix<int> hflag;
	/// 1 if corresponding boundary point is a high-order point (nbpoin-by-1)
	Matrix<int> bhflag;
	Matrix<int> fixedflag;
	/// Contains markers of fixed boundary flags
	Matrix<int> fbf;
	RBFmove rbfm;
	int rbf_c;
	double srad;
	int rbf_nsteps;
	double tol;
	int maxiter;
	/// linear solver to use for RBFmove
	std::string rbfsolver;

public:
	CurvedMeshGeneration(UMesh* mesh, Matrix<int> fixed_boundary_flags, int rbf_steps, int rbf_choice, double support_radius, double rbf_tol, int rbf_maxiter, std::string rbf_solver)
	// Note that the input mesh should be quadratic
	{
		cout << "CurvedMeshGeneration: Pre-processing..." << endl;
		m = mesh;
		fbf = fixed_boundary_flags;
		rbf_c = rbf_choice;
		srad = support_radius;
		rbf_nsteps = rbf_steps;
		tol = rbf_tol;
		maxiter = rbf_maxiter;
		rbfsolver = rbf_solver;
		
		noInternalFixedPoints = false;
		
		cout << "CurvedMeshGeneration: Support radius " << srad << endl;

		Matrix<int> pfixedflag(m->gnpoin(),1);
		pfixedflag.ones();

		// contains 1 if a face is to be curved
		amat::Matrix<int> facesToCurve(m->gnface(),1);
		facesToCurve.ones();
		for(int iface = 0; iface < m->gnface(); iface++)
			for(int ifl = 0; ifl < fbf.msize(); ifl++)
				if(m->gbface(iface,m->gnnofa()) == fbf.get(ifl))
					facesToCurve(iface) = 0;

		// pre-process to get ipoints and bpoints :-
		//cout << "CurvedMeshGeneration: bflag" << endl;
		bflag.setup(m->gnpoin(), 1);
		bflag.zeros();
		// bflag will contain 1 if the corresponding point is a boundary point
		for(int i = 0; i < m->gnface(); i++)
		{
			for(int j = 0; j < m->gnnofa(); j++)
				bflag(m->gbface(i,j)) = 1;
			
			if( facesToCurve.get(i) )
				for(int j = 0; j < m->gnnofa(); j++)
					pfixedflag(m->gbface(i,j)) = 0;			// register this point as belonging to a curved boundary
		}

		//cout << "CurvedMeshGeneration: hflag" << endl;
		hflag.setup(m->gnpoin(), 1);
		hflag.ones();
		for(int i = 0; i < m->gnelem(); i++)
		{
			for(int j = 0; j < m->gnnode()-m->gnfael()-m->gnedel()-1; j++)
				hflag(m->ginpoel(i,j)) = 0;
		}

		//cout << "CurvedMeshGeneration: nbpoin" << endl;
		nbpoin = 0;
		for(int i = 0; i < bflag.rows(); i++)
		{
			if(bflag(i) == 1)
			{
				nbpoin++;
			}
		}

		ninpoin = m->gnpoin() - nbpoin;			// number of interior points
		cout << "CurvedMeshGeneration: Number of internal points = " << ninpoin << ", number of boundary points = " << nbpoin << endl;

		bpoints.setup(nbpoin, m->gndim());
		bmotion.setup(nbpoin, m->gndim());
		fixedflag.setup(nbpoin, 1);
		fixedflag.zeros();
		bhflag.setup(nbpoin, 1);
		bhflag.zeros();

		//cout << "CurvedMeshGeneration: bpoints and hflag" << endl;
		int k = 0;
		for(int i = 0; i < bflag.rows(); i++)
		{
			if(bflag(i) == 1)
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					bpoints(k,idim) = mesh->gcoords(i,idim);
				if(hflag(i)==1) bhflag(k) = 1;					// this point is a high-order point
				if(pfixedflag(i)==1) fixedflag(k) = 1;			// this point belongs to a fixed boundary
				k++;
			}
		}

		npomo = 0;
		for(int i = 0; i < m->gnpoin(); i++)
			if(hflag(i) == 1 && bflag(i)==0) npomo++;		// point is a high order point and not a boundary point

		cout << "CurvedMeshGeneration: Number of movable points = " << npomo << endl;

		mipoints.setup(npomo, m->gndim());					// npomo is the number of interior points to be moved
		if(ninpoin-npomo == 0)
			noInternalFixedPoints = true;
		else
			fipoints.setup(ninpoin-npomo, m->gndim());

		k = 0; int l = 0;
		for(int i = 0; i < m->gnpoin(); i++)
		{
			if(bflag(i)==0 && hflag(i) == 1)				// point is an interior point and a high-order point
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					mipoints(k,idim) = m->gcoords(i,idim);
				k++;
			}
			else if(bflag(i)==0 && hflag(i) == 0 && noInternalFixedPoints==false)			// point is an interior point but not a high order point (fixed)
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					fipoints(l,idim) = m->gcoords(i,idim);
				l++;
			}
		}
	}

	double trueboundary(double x0, double y)
	{
		// TODO: change this block according to the desired boundary motion
		
		// Below, the analytical boundary of the 3D-bump case
		double x = x0 - 0.3*pow(sin(PI*y),4);
		if(x >= 0.3 && x <= 1.2)
			return 0.05*pow(sin(PI*x/0.9-(PI/3.0)),4);
		else return 0.0;
		
		/*if(x0 >= 0.3 && x0 <= 1.2)
			return 0.05*pow((sin(PI*x0/0.9-(PI/3.0))),4);
		else return 0.0;*/
	}

	void compute_boundary_displacement()
	{
		cout << "CurvedMeshGeneration: compute_boundary_displacement(): Calculating bmotion" << endl;
		bmotion.zeros();
		for(int i = 0; i < nbpoin; i++)
		{
			if(bhflag(i) == 1 && fixedflag(i) != 1)		// if point is a high-order point and is not part of a fixed boundary
			{
				// TODO: change this block according to the desired boundary motion

				//bmotion(i,0) = 0;
				//bmotion(i,1) = 0;
				bmotion(i,2) = trueboundary(bpoints(i,0),bpoints(i,1)) - bpoints(i,2);
				//bmotion(i,1) = 0.2;
			}
		}
	}

	/*void compute_boundary_displacement_comp(const UMesh& mesh)
	// gets boundary displacements by comparison with another mesh having same topology.
	{
		cout << "CurvedMeshGeneration: compute_boundary_displacement_comp(): Calculating bmotion" << endl;
		bmotion.zeros();

		// first, get bpoints of second mesh
		int k = 0;
		Matrix<double> bpoints2(nbpoin, m->gndim());
		for(int i = 0; i < bflag.rows(); i++)
		{
			if(bflag(i) == 1)
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					bpoints2(k,idim) = mesh.gcoords(i,idim);
				k++;
			}
		}

		// now assuming that topology of the 2 meshes is exactly same, get boundary displacements

		for(int i = 0; i < nbpoin; i++)
		{
			if(bhflag(i) == 1 && fixedflag(i) != 1)		// if point is a high-order point and is not part of a fixed boundary
			{
				bmotion(i,0) = bpoints2(i,0) - bpoints(i,0);
				bmotion(i,1) = bpoints2(i,1) - bpoints(i,1);
			}
		}

		double norm = 0;
		for(int i = 0; i < nbpoin; i++)
		{
			norm += bmotion(i,0)*bmotion(i,0) + bmotion(i,1)*bmotion(i,1);
		}
		norm = sqrt(norm);
		cout << "Norm of bmotion is " << norm << endl;
		cout << "CurvedMeshGeneration: compute_boundary_displacement_comp(): bmotion calculated." << endl;
	}*/

	void generate()
	{
		/*int ipoin, idim, k = 0;
		amc_real temp;
		for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		{
			if(bflag.get(ipoin) == 1)
			{
				for(idim = 0; idim < m->gndim(); idim++)
				{
					temp = m->gcoords(ipoin,idim);
					m->scoords(ipoin,idim, temp+bmotion.get(k,idim));
				}
				k++;
			}
		}*/

		cout << "CurvedMeshGeneration: generate(): Setting up RBF mesh movement" << endl;
		cout << nbpoin << endl;
		rbfm.setup(&mipoints, &bpoints, &bmotion, rbf_c, srad, rbf_nsteps, tol, maxiter, rbfsolver);
		cout << "CurvedMeshGeneration: generate(): Moving the mesh" << endl;
		rbfm.move();
		bpoints = rbfm.getBoundaryPoints();
		mipoints = rbfm.getInteriorPoints();

		cout << "CurvedMeshGeneration: generate(): Re-assembling coordinate matrix" << endl;
		// re-assemble coords
		int mp = 0, fp = 0, bp = 0;
		Matrix<double> coord(m->gnpoin(), m->gndim());
		for(int i = 0; i < m->gnpoin(); i++)
		{
			if(bflag(i) == 1)
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					coord(i,idim) = bpoints(bp,idim);
				bp++;
			}
			else if(hflag(i) == 1)
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					coord(i,idim) = mipoints(mp,idim);
				mp++;
			}
			else
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					coord(i,idim) = fipoints(fp,idim);
				fp++;
			}
		}
		m->setcoords(&coord);
		cout << "CurvedMeshGeneration: generate(): Done." << endl;
	}
};
