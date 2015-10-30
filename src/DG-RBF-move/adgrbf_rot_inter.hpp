/* Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
Aditya Kashi
July 1, 2015
*/
#include <alinalg.hpp>
#include <arotation2db.hpp>
#include "abowyerwatson.hpp"

using namespace std;
using namespace amat;

namespace acfd {

class DGRBFmove
{
	UTriMesh* m;
	Matrix<int> nrot;
	double angle;		// in radians
	MRotation2d rot;
	double xc;
	double yc;

	int ndim;
	int npoin;
	int ninpoin;
	Matrix<double> points;
	int ndgpoin;
	Matrix<double> dgpoints;
	int ndgelem;
	Matrix<int> dginpoel;
	Matrix<int> dgesuel;
	Matrix<double> dgintfac;
	Matrix<double> bmotion;

	Matrix<double>* Mtt;
	Matrix<double>* ralpha;
	Matrix<double>* rs;
	bool isallocMtt;
	bool isallocralpha;
	double (DGRBFmove::*rbf)(double);			// function pointer for radial basis functions
	double srad;			// support radius for RBFs

	double tol;

public:
	Delaunay2D dg;				// object to construct triangulation of boundary points.
	Matrix<double> newcoords;

	DGRBFmove(UTriMesh* mesh, int rbf_type, double support_radius, Matrix<int> n_rot, double angle_rot, double cx, double cy)
	{
		// note that bflags contains as many rows as points in the mesh and contains 1 if the point is a boundary point and 0 otherwise
		// boundary_motion has as many rows as points in the original mesh and contains x and y displacement values for each point (zero for interior points)
		m = mesh;
		angle = angle_rot;
		nrot = n_rot;
		ndim = m->gndim();
		ndgpoin = m->gnbpoin();		// ndgpoin is the number of points in the DG, which equals the number of boundary nodes in the original mesh
		ninpoin = m->gnpoin() - ndgpoin;
		points.setup(ninpoin, ndim+1);	// for each interior point, store coords and containing DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element.
		dgpoints.setup(ndgpoin,ndim);

		int i;
		int k = 0;
		for(i = 0; i < m->gnbpoin(); i++)
			for(int j = 0; j < ndim; j++)
				dgpoints(i,j) = m->gcoords(i,j);

		for( ; i < m->gnpoin(); i++)
			for(int j = 0; j < ndim; j++)
				points(i-m->gnbpoin(),j) = m->gcoords(i,j);

		newcoords.setup(m->gnpoin(),ndim);

		dg.setup(&dgpoints, ndgpoin);

		isallocMtt = false;
		isallocralpha = false;

		switch(rbf_type)
		{
			case(0): rbf = &DGRBFmove::rbf_c0;
			break;
			case(2): rbf = &DGRBFmove::rbf_c2_compact;
			break;
			case(4): rbf = &DGRBFmove::rbf_c4;
			default: rbf = &DGRBFmove::rbf_c2_compact;
		}

		srad = support_radius;

		tol = 1e-20;			// very small number

		cout << "DGRBFmove: Calculating boundary displacements\n";
		rot.setup(m, angle, cx, cy, nrot);
		xc = cx;
		yc = cy;
		bmotion = rot.rhsvect_rotate();

		cout << "DGRBFmove: Setup complete. Angle = " << angle << ", RBF choice = " << rbf_type << endl;
	}

	~DGRBFmove()
	{
		if(isallocMtt == true) delete [] Mtt;
		if(isallocralpha == true) { delete [] ralpha; delete [] rs; }
	}

	// RBFs
	double rbf_c2(double xi)
	{
		return pow(1-xi,4)*(4*xi+1);
	}
	double rbf_c2_compact(double xi)
	{
		double t = xi/srad;
		if(t <= 1.0)
			return pow(1-t,4)*(4*t+1);
		else return 0;
	}
	double rbf_c0(double xi)
	{
		return (1-xi)*(1-xi);
	}
	double rbf_c4(double xi)
	{
		return pow(1-xi,6)*(35*xi*xi + 18*xi + 3);
	}
	double gaussian(double xi)
	{
		return exp(-xi*xi);
	}

	void generateDG()
	{
		// generate the Delaunay graph
		dg.bowyer_watson();
		dg.writeGmsh2("rot-dg.msh");

		ndgelem = dg.elems.size();
		cout << "DGmove: generateDG(): No. of DG elements: " << ndgelem << endl;
		dginpoel.setup(ndgelem,ndim+1);
		dgesuel.setup(ndgelem,ndim+1);

		// populate DG data arrays
		for(int iel = 0; iel < dg.elems.size(); iel++)
		{
			for(int j = 0; j < ndim+1; j++)
			{
				dginpoel(iel,j) = dg.elems[iel].p[j];
				dgesuel(iel,j) = dg.elems[iel].surr[j];
			}
		}

		// for debugging
		/*ofstream dginp("dginpoel.dat");
		dginpoel.fprint(dginp);
		dginp.close();*/

		if(isallocMtt == false)
		{
			Mtt = new Matrix<double>[ndgelem];
			for(int i = 0; i < ndgelem; i++)
				Mtt[i].setup(ndim+1, ndim+1);
			isallocMtt = true;
		}
		if(isallocralpha == false)
		{
			ralpha = new Matrix<double>[ndgelem];
			rs = new Matrix<double>[ndgelem];
			for(int i = 0; i < ndgelem; i++)
			{
				ralpha[i].setup(ndim+1,1);
				rs[i].setup(ndim+1,1);
			}
			isallocralpha = true;
		}
	}

	void set_boundary_rotation()
	{
		if(isallocralpha == false)
		{
			ralpha = new Matrix<double>[ndgelem];
			rs = new Matrix<double>[ndgelem];
			for(int i = 0; i < ndgelem; i++)
			{
				ralpha[i].setup(ndim+1,1);
				rs[i].setup(ndim+1,1);
			}
			isallocralpha = true;
		}

		for(int i = 0; i < ndgelem; i++)
			for(int j = 0; j < ndim+1; j++)
				rs[i](j) = 0;

		for(int i = 0; i < ndgelem; i++)
		{
			for(int j = 0; j < ndim+1; j++)
				for(int k = 0; k < nrot.msize(); k++)
					if(m->gbpoints(dginpoel(i,j),4) == nrot(k))
						rs[i](j) = angle;
		}
	}

	double det3(Matrix<double> A)		// determinant of 3x3 matrix
	{
		double d;
		d = A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) - A(0,1)*(A(1,0)*A(2,2)-A(2,0)*A(1,2)) + A(0,2)*(A(1,0)*A(2,1)-A(2,0)*A(1,1));
		return d;
	}

	Matrix<double> cramer3(Matrix<double> A, Matrix<double> b)
	{
		double ddet = 0;
		Matrix<double> x(b.rows(), b.cols());
		Matrix<double> Anum(A.rows(),A.cols());
		ddet = det3(A);
		if(dabs(ddet) < tol) cout << "! DGRBFmove: cramer3(): Matrix A is singular!!\n";
		for(int j = 0; j < b.cols(); j++)
		{
			for(int i = 0; i < b.rows(); i++)
			{
				Anum = A;
				for(int ii = 0; ii < b.rows(); ii++)
					Anum(ii,i) = b(ii,j);
				x(i,j) = det3(Anum)/ddet;
			}
		}
		return x;
	}

	void calcMtt()
	{
		//check if Mtt is allocated
		if(isallocMtt == false)
		{
			Mtt = new Matrix<double>[ndgelem];
			for(int i = 0; i < ndgelem; i++)
				Mtt[i].setup(ndim+1, ndim+1);
			isallocMtt = true;
		}
		cout << "DGRBFmove: calcrMtt(): Calculating RBFs\n";

		int iel, ipoin, jpoin;
		double msum = 0;
		for(iel = 0; iel < ndgelem; iel++)
		{
			for(int i = 0; i < ndim+1; i++)
			{
				ipoin = dginpoel(iel,i);
				for(int j = 0; j < ndim+1; j++)
				{
					jpoin = dginpoel(iel,j);
					msum = 0.0;
					for(int idim = 0; idim < ndim; idim++)
					 	msum += (dgpoints(ipoin,idim)-dgpoints(jpoin,idim))*(dgpoints(ipoin,idim)-dgpoints(jpoin,idim));
					if(msum < 0) { cout << "! DGRBFmove: calcMtt(): msum is " << msum << endl; msum = 0;}
					Mtt[iel](i,j) = (this->*rbf)(sqrt(msum));
				}
			}
		}
	}

	void calcralpha()
	{
		for(int iel = 0; iel < ndgelem; iel++)
		{
			// first, calculate s
			// this is to be done in calling function

			// calculate alpha
			for(int idim = 0; idim < ndim; idim++)
			{
				ralpha[iel] = cramer3(Mtt[iel],rs[iel]);
			}
		}
	}

	void movemesh()
	{
		generateDG();

		// cycle over interior points
		int contelem;
		cout << "DGRBFmove: movemesh(): Calculating containing elements for each interior point\n";
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			//cout << ipoin << endl;
			// first find containing DG element by "walking-through" the DG
			contelem = dg.find_containing_triangle(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);

			// store DG element in points
			points(ipoin,2) = contelem;
		}

		set_boundary_rotation();

		calcMtt();

		calcralpha();

		ofstream alph("rs-compact.dat");
		for(int iel = 0; iel < ndgelem; iel++)
		{
			for(int inode = 0; inode < ndim+1; inode++)
			{
					alph << rs[iel](inode) << " ";
				alph << '\n';
			}
			alph << '\n';
		}
		alph.close();

		// calculate new positions of interior points by mapping them to deformed DG elements using area coordinates calculated before
		cout << "DGRBFmove: movemesh(): Moving the interior points\n";
		Matrix<double> A(ndim+1,1);
		int elem;
		double* rr = new double[ndim];
		double rrs = 0;
		double sum = 0;
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			rrs = 0;
			elem = points(ipoin,2);
			//cout << "* " << elem << endl;

			// calculate RBFs of point ipoin and store in A
			for(int inode = 0; inode < ndim+1; inode++)
			{
				sum = 0;
				for(int idim = 0; idim < ndim; idim++)
				{
					rr[idim] = dgpoints.get(dginpoel(elem,inode),idim);
					sum += (rr[idim] - points(ipoin,idim)) * (rr[idim] - points(ipoin,idim));
				}

				A(inode) = (this->*rbf)( sqrt(sum) );
			}

			// Calculate angle of ipoin using A and ralpha, and update coordinates of point
			//for(int idim = 0; idim < ndim; idim++)
			//{
				for(int inode = 0; inode < ndim+1; inode++)
					rrs += A(inode)*ralpha[elem](inode);
				points(ipoin,0) = (points(ipoin,0)-xc)*cos(rrs) - (points(ipoin,1)-yc)*sin(rrs) + xc;
				points(ipoin,1) = (points(ipoin,0)-xc)*sin(rrs) + (points(ipoin,1)-yc)*cos(rrs) + yc;
			//}
		}
		delete [] rr;

		// update coordinates in dgpoints using bmotion
		cout << "DGmove: movemesh(): Moving the Delaunay graph\n";
		for(int i = 0; i < ndgpoin; i++)
		{
				for(int j = 0; j < ndim; j++)
					dgpoints(i,j) += bmotion.get(i,j);
		}
	}

	void movedg()
	{
		//moves the DG according to the boundary motion - for visualizing deformed DG
		for(int i = 0; i < dgpoints.rows(); i++)
		{
			dg.nodes[i].x = dgpoints(i,0);
			dg.nodes[i].y = dgpoints(i,1);
		}
	}

	Matrix<double> getcoords()
	{
		// create a coords matrix with same point numbering as initial matrix and return it
		int i;
		for(i = 0; i < ndgpoin; i++)
		{
			for(int idim = 0; idim < ndim; idim++)
				newcoords(i,idim) = dgpoints(i,idim);
		}
		for(i = ndgpoin; i < m->gnpoin(); i++)
		{
			for(int idim = 0; idim < ndim; idim++)
				newcoords(i,idim) = points(i-ndgpoin,idim);
		}
		return newcoords;
	}
};

} // end namespace acfd
