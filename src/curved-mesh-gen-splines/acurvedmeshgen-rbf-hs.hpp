#ifndef __AMESH2DGENERAL_H
#include <amesh2d.hpp>
#endif

#include <arbf.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

class CurvedMeshGeneration
{
	UMesh2d* m;
	Matrix<double> mipoints;		// interior points to be moved
	Matrix<double> fipoints;			// fixed interior points
	Matrix<double> bpoints;
	Matrix<double> bmotion;
	int ninpoin;				// number of interior points
	int npomo;					// number of interior points to be moved
	int nbpoin;					// number of  boundary points
	Matrix<int> bflag;			// 1 if corresponding point is a boundary points (npoin-by-1)
	Matrix<int> hflag;			// 1 if corresponding point is a high-order point (npoin-by-1)
	Matrix<int> bhflag;			// 1 if corresponding boundary point is a high-order point (nbpoin-by-1)
	Matrix<int> fixedflag;
	Matrix<int> fbf;
	RBFmove rbfm;
	int rbf_c;
	double srad;
	int rbf_nsteps;

	HermiteSpline2d h;			// to hold spline data read from file
	Matrix<double> bspmotion;		// to hold displacement for all boundary high-order nodes, and zero for other nodes

public:
	CurvedMeshGeneration(UMesh2d* mesh, Matrix<int> fixed_boundary_flags, string spfile,  int rbf_steps, int rbf_choice, double support_radius)
	// Note that the input mesh should be quadratic
	{
		cout << "CurvedMeshGeneration: Pre-processing..." << endl;
		m = mesh;
		fbf = fixed_boundary_flags;
		rbf_c = rbf_choice;
		srad = support_radius;
		rbf_nsteps = rbf_steps;
		cout << "CurvedMeshGeneration: Support radius " << srad << endl;
		h.readGeomCoeffsFromFile(spfile);
		bspmotion.setup(m->gnpoin(),m->gndim());
		bspmotion.zeros();

		cout << "CurvedMeshGeneration: Calculating boundary displacements..." << endl;
		vector<double> tb(m->gndim());
		for(int iface = 0; iface < m->gnface(); iface++)
		{
			int pno = m->gbface(iface,2);		// get point number of high-order node of this bface
			for(int idim = 0; idim < m->gndim(); idim++)
			{
				tb[idim] = h.spline(iface,idim,0.5);
				bspmotion(pno,idim) = tb[idim] - m->gcoords(pno,idim);
			}
		}

		ofstream fout("ransnaca.dat");
		for(int iface = 0; iface < m->gnface(); iface++)
		{
			bool check = false;
			for(int ibf = 0; ibf < fbf.msize(); ibf++)
				if(m->gbface(iface,m->gnnofa()) == fbf(ibf)) check = true;
			if(check != true)
			{
				for(double u = 0; u <= 1; u+=0.025)
					fout << h.spline(iface,0,u) << " " << h.spline(iface,1,u) << '\n';
			}
		}
		fout.close();

		Matrix<int> pfixedflag(m->gnpoin(),1);
		pfixedflag.zeros();

		// pre-process to get ipoints and bpoints :-
		//cout << "CurvedMeshGeneration: bflag" << endl;
		bflag.setup(m->gnpoin(), 1);
		bflag.zeros();
		// bflag will contain 1 if the corresponding point is a boundary point
		for(int i = 0; i < m->gnface(); i++)
		{
			for(int j = 0; j < m->gnnofa(); j++)
				bflag(m->gbface(i,j)) = 1;
			for(int ifl = 0; ifl < fbf.msize(); ifl++)
			{
				if(m->gbface(i,m->gnnofa()) == fbf(ifl))
				{
					for(int j = 0; j < m->gnnofa(); j++)
						pfixedflag(m->gbface(i,j)) = 1;			// register this point as belonging to a fixed boundary
				}
			}
		}

		//cout << "CurvedMeshGeneration: hflag" << endl;
		hflag.setup(m->gnpoin(), 1);
		hflag.ones();
		for(int i = 0; i < m->gnelem(); i++)
		{
			for(int j = 0; j < m->gnnode()-m->gnfael()-1; j++)		// NOTE: This is assuming there's a cell-centre high-order node as well!!
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

		bpoints.setup(nbpoin, m->gndim());
		bmotion.setup(nbpoin, m->gndim());
		fixedflag.setup(nbpoin, 1);
		fixedflag.zeros();
		bhflag.setup(nbpoin, 1);
		bhflag.zeros();

		//cout << "CurvedMeshGeneration: bpoints, bmotion and hflag" << endl;
		int k = 0;
		for(int i = 0; i < bflag.rows(); i++)
		{
			if(bflag(i) == 1)
			{
				for(int idim = 0; idim < m->gndim(); idim++) {
					bpoints(k,idim) = mesh->gcoords(i,idim);
					bmotion(k,idim) = bspmotion(i,idim);
				}
				if(hflag(i)==1) bhflag(k) = 1;					// this point is a high-order point
				if(pfixedflag(i)==1) fixedflag(k) = 1;			// this point belongs to a fixed boundary
				k++;
			}
		}

		npomo = 0;
		for(int i = 0; i < m->gnpoin(); i++)
			if(hflag(i) == 1 && bflag(i)==0) npomo++;		// point is a high order point and not a boundary point

		mipoints.setup(npomo, m->gndim());					// npomo is the number of interior points to be moved
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
			else if(bflag(i)==0 && hflag(i) == 0)			// pointis an interior point but not a high order point (fixed)
			{
				for(int idim = 0; idim < m->gndim(); idim++)
					fipoints(l,idim) = m->gcoords(i,idim);
				l++;
			}
		}
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
				/*if(!(bpoints(i,0) <= 0.00002721276 && bpoints(i,1) < 0.00092310243 && bpoints(i,1) > -0.00092310243))
					if(bpoints(i,1) <= 0.0)
						bmotion(i,1) = -1*trueboundary(bpoints(i,0)) - bpoints(i,1);
					else
						bmotion(i,1) = trueboundary(bpoints(i,0)) - bpoints(i,1);*/

				bmotion(i,1) = 0.3;
			}
		}
	}

	void compute_boundary_displacement_comp(const UMesh2d& mesh)
	/* gets boundary displacements by comparison with another mesh having same topology. */
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

		/*double norm = 0;
		for(int i = 0; i < nbpoin; i++)
		{
			norm += bmotion(i,0)*bmotion(i,0) + bmotion(i,1)*bmotion(i,1);
		}
		norm = sqrt(norm);
		cout << "Norm of bmotion is " << norm << endl;*/
		cout << "CurvedMeshGeneration: compute_boundary_displacement_comp(): bmotion calculated." << endl;
	}

	void generate()
	{
		cout << "CurvedMeshGeneration: generate(): Setting up RBF mesh movement" << endl;
		cout << nbpoin << endl;
		rbfm.setup(&mipoints, &bpoints, &bmotion, rbf_c, srad, rbf_nsteps, 1e-8, 20000);
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
