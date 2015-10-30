/* Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
Aditya Kashi
July 1, 2015
*/

#include "abowyerwatson.hpp"

using namespace std;
using namespace amat;

namespace acfd {

class DGRBFmove
{
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
	Matrix<double>* bmotion;

	Matrix<double>* Mtt;
	Matrix<double>* alpha;
	Matrix<double>* s;
	bool isallocMtt;
	bool isallocalpha;
	double (DGRBFmove::*rbf)(double);
	double srad;			// support radius for RBFs

	double tol;

public:
	Delaunay2D dg;

	DGRBFmove() {}
	
	DGRBFmove(Matrix<double>* incoords, Matrix<double>* bouncoords, Matrix<double>* boundary_motion, int rbf_type, double support_radius)
	{
		// note that bflags contains as many rows as points in the mesh and contains 1 if the point is a boundary point and 0 otherwise
		// boundary_motion has as many rows as boundary points in the original mesh and contains x and y displacement values for each point
		bmotion = boundary_motion;
		ndim = incoords->cols();
		ndgpoin = bouncoords->rows();		// ndgpoin is the number of points in the DG, which equals the number of boundary nodes in the original mesh
		ninpoin = incoords->rows();
		points.setup(ninpoin, ndim+1);	// for each interior point, store coords and containing DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.

		for(int i = 0; i < incoords->rows(); i++)
			for(int j = 0; j < ndim; j++)
				points(i,j) = incoords->get(i,j);

		// now get dgpoints using bfac
		dgpoints.setup(ndgpoin,ndim);
		dgpoints = *bouncoords;

		dg.setup(&dgpoints, ndgpoin);

		isallocMtt = false;
		isallocalpha = false;
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
		//cout << "Test: " << rbf(0)<< ' ' << sqrt(0.0) << endl;
	}

	void setup(Matrix<double>* incoords, Matrix<double>* bouncoords, Matrix<double>* boundary_motion, int rbf_type, double support_radius)
	{
		// note that bflags contains as many rows as points in the mesh and contains 1 if the point is a boundary point and 0 otherwise
		// boundary_motion has as many rows as boundary points in the original mesh and contains x and y displacement values for each point
		bmotion = boundary_motion;
		ndim = incoords->cols();
		ndgpoin = bouncoords->rows();		// ndgpoin is the number of points in the DG, which equals the number of boundary nodes in the original mesh
		ninpoin = incoords->rows();
		points.setup(ninpoin, ndim+1);	// for each interior point, store coords and containing DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.

		for(int i = 0; i < incoords->rows(); i++)
			for(int j = 0; j < ndim; j++)
				points(i,j) = incoords->get(i,j);

		// now get dgpoints using bfac
		dgpoints.setup(ndgpoin,ndim);
		dgpoints = *bouncoords;

		dg.setup(&dgpoints, ndgpoin);

		isallocMtt = false;
		isallocalpha = false;
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
		//cout << "Test: " << rbf(0)<< ' ' << sqrt(0.0) << endl;
	}

	~DGRBFmove()
	{
		if(isallocMtt == true) delete [] Mtt;
		if(isallocalpha == true) { delete [] alpha; delete [] s; }
	}

	// RBFs
	double rbf_c2(double xi)
	{
		return pow(1-xi,4)*(4*xi+1);
	}
	double rbf_c2_compact(double xi)
	{
		if(xi/srad <= 1.0)
			return pow(1-xi/srad,4)*(4*xi/srad+1);
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
		if(isallocalpha == false)
		{
			alpha = new Matrix<double>[ndgelem];
			s = new Matrix<double>[ndgelem];
			for(int i = 0; i < ndgelem; i++)
			{
				alpha[i].setup(ndim+1,ndim);
				s[i].setup(ndim+1,ndim);
			}
			isallocalpha = true;
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
		cout << "DGRBFmove: calcMtt(): Calculating RBFs\n";

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

	void calcalpha()
	{
		Matrix<double>* a = new Matrix<double>[ndim];
		// iterate over Delaunay elements
		for(int iel = 0; iel < ndgelem; iel++)
		{
			// first, calculate s
			for(int i = 0; i < ndim+1; i++)
			{
				for(int idim = 0; idim < ndim; idim++)
					s[iel](i,idim) = bmotion->get(dginpoel(iel,i),idim);
			}

			// calculate alpha
			Matrix<double> x0(ndim+1,1); x0.zeros();

			for(int idim = 0; idim < ndim; idim++)
			{
				alpha[iel] = cramer3(Mtt[iel],s[iel]);
			}

			/*for(int i = 0; i < ndim+1; i++)
			{
				for(int idim = 0; idim < ndim; idim++)
					alpha[iel](i,idim) = a[idim](i,0);
			}*/
		}
		delete [] a;
	}

	void movemesh()
	{
		// cycle over interior points
		int contelem;
		cout << "DGmove: movemesh(): Calculating containing elements for each interior point\n";
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			//cout << ipoin << endl;
			// first find containing DG element by "walking-through" the DG
			contelem = dg.find_containing_triangle(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);

			// store DG element in points
			points(ipoin,2) = contelem;
		}

		calcMtt();


		calcalpha();

		ofstream alph("alphas-compact.dat");
		for(int iel = 0; iel < ndgelem; iel++)
		{
			for(int inode = 0; inode < ndim+1; inode++)
			{
				for(int idim = 0; idim < ndim; idim++)
				{
					alph << alpha[iel](inode,idim) << " ";
				}
				alph << '\n';
			}
			alph << '\n';
		}
		alph.close();

		// calculate new positions of interior points by mapping them to deformed DG elements
		cout << "DGmove: movemesh(): Moving the interior points\n";
		Matrix<double> A(ndim+1,1);
		cout << "DGmove: movemesh():  Moving the interior points\n";
		int elem; double* rr = new double[ndim];
		double sum = 0;
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
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

			// Calculate displacement of ipoin using A and alpha, and update coordinates of point
			for(int idim = 0; idim < ndim; idim++)
			{
				rr[idim] = 0.0;
				for(int inode = 0; inode < ndim+1; inode++)
					rr[idim] += A(inode)*alpha[elem](inode, idim);
				points(ipoin,idim) += rr[idim];
			}
		}
		delete [] rr;

		// update coordinates in dgpoints using bmotion
		cout << "DGmove: movemesh(): Moving the Delaunay graph\n";
		for(int i = 0; i < bmotion->rows(); i++)
		{
			for(int j = 0; j < ndim; j++)
				dgpoints(i,j) += bmotion->get(i,j);
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

	Matrix<double> getInteriorPoints()
	{ return points; }

	Matrix<double> getBoundaryPoints()
	{ return dgpoints; }
};

} // end namespace acfd
