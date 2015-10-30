/* Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
Aditya Kashi
July 1, 2015
*/

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

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

	Matrix<double>* Mtt;
	Matrix<double>* alpha;
	Matrix<double>* s;
	bool isallocMtt;
	bool isallocalpha;
	double (DGRBFmove::*rbf)(double);
	double srad;			// support radius for RBFs

	double tol;

	int ndgpoin;
	Matrix<double> dgpoints;
	Matrix<double> testdgpoints;
	int ndgelem;
	Matrix<int> dginpoel;
	Matrix<int> dgesuel;
	Matrix<double> dgintfac;
	bool reflag;

	Matrix<double>* bmotion;
	Matrix<int> bflag;
	Delaunay2D testdg;

	bool check;
	Matrix<double>* bmot;

public:
	Delaunay2D dg;
	Matrix<double> newcoords;

	DGRBFmove(int num_dimn, Matrix<double>* coords, Matrix<int> bflags, Matrix<double>* boundary_motion, int rbf_type, double support_radius)
	{
		// note that bflags contains as many rows as points in the mesh and contains 1 if the point is a boundary point and 0 otherwise
		// boundary_motion has as many rows as points in the original mesh and contains x and y displacement values for each point (zero for interior points)
		bmotion = boundary_motion;
		if(bmotion->rows() != bflags.rows()) cout << "DGMove: !! Error: No. of rows in bflags and bmotion are not equal!\n";
		ndim = coords->cols();
		ndgpoin = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1) ndgpoin++;		// ndgpoin is the number of points in the DG, which equals the number of boundary nodes in the original mesh
		ninpoin = coords->rows() - ndgpoin;
		points.setup(ninpoin, ndim+1);	// for each interior point, store coords, containing DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element.

		int k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 0)
			{
				for(int j = 0; j < ndim; j++)
					points(k,j) = coords->get(i,j);
				k++;
			}

		// now get dgpoints using bfac
		dgpoints.setup(ndgpoin,ndim);
		testdgpoints.setup(ndgpoin,ndim);
		k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1)
			{
				for(int j = 0; j < ndim; j++)
					dgpoints(k,j) = coords->get(i,j);
				k++;
			}

		bflag = bflags;
		newcoords.setup(bflags.rows(),ndim);
		dg.setup(&dgpoints, ndgpoin);
		//testdg = dg;

		// copy bmotion into bmot
		bmot = new Matrix<double>;
		/*bmot->setup(bmotion->rows(),bmotion->cols());
		for(int i = 0; i < bmotion->rows(); i++)
			for(int j = 0; j < bmotion->cols(); j++)
				(*bmot)(i,j) = bmotion->get(i,j);*/

		bmot->setup(ndgpoin, ndim);
		k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1)
			{
				for(int j = 0; j < ndim; j++)
					(*bmot)(k,j) = bmotion->get(i,j);
				k++;
			}

		isallocMtt = false;
		isallocalpha = false;
		reflag = false;
		switch(rbf_type)
		{
			case(0): rbf = &DGRBFmove::rbf_c0_compact;
			break;
			case(2): rbf = &DGRBFmove::rbf_c2_compact;
			break;
			case(4): rbf = &DGRBFmove::rbf_c4_compact;
			case(5): rbf = &DGRBFmove::gaussian_compact;
			default: rbf = &DGRBFmove::rbf_c2_compact;
		}

		srad = support_radius;

		tol = 1e-20;			// very small number
	}

	~DGRBFmove()
	{
		delete bmot;
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
	double rbf_c0_compact(double xi)
	{
		if(xi <= srad)
			return (1-xi/srad)*(1-xi/srad);
		else return 0;
	}
	double rbf_c4_compact(double xi)
	{
		if(xi/srad <= 1.0)
			return pow(1-xi/srad,6)*(35*xi*xi/(srad*srad) + 18*xi/srad + 3);
		else
			return 0;
	}
	double gaussian_compact(double xi)
	{
		if(xi <= srad)
			return exp(-xi*xi/(srad*srad));
		else return 0;
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
					Mtt[iel](i,j) = (this->*rbf)(sqrt(msum));
				}
			}
		}
	}

	void calcalpha()
	{
		//cout << "DGRBFmove: calcalpha(): Calculating coefficients of RBFs\n";
		Matrix<double>* a = new Matrix<double>[ndim];
		for(int iel = 0; iel < ndgelem; iel++)
		{
			// first, calculate s
			for(int i = 0; i < ndim+1; i++)
			{
				for(int idim = 0; idim < ndim; idim++)
					s[iel](i,idim) = bmot->get(dginpoel(iel,i),idim);
			}

			// calculate alpha
			Matrix<double> x0(ndim+1,1); x0.zeros();

			for(int idim = 0; idim < ndim; idim++)
			{
				//a[idim] = gaussseidel(Mtt[iel], s[iel].col(idim), x0, 1e-8, 10000, 'n');
				//a[idim] = gausselim(Mtt[iel], s[iel].col(0));
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
		int num_steps = 1;
		int step = 1;			// just for info of total number of steps
		int numtri = 1;			// number of triangulations needed

		generateDG();

		cout << "DGmove: movemesh():  Calculating containing element for each interior point\n";
		int contelem;
		// cycle over interior points
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			// first find containing DG element by "walking-through" the DG
			contelem = dg.find_containing_triangle(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);

			// store DG element in points
			points(ipoin,2) = contelem;
		}

		// calculate RBFs
		calcMtt();

		while(num_steps > 0)
		{
			//cout << "DGmove: movemesh(): Step " << step << endl;
			cout << "DGmove: movemesh(): number of steps " << num_steps << endl;

			testdg = dg;

			// update coordinates in dgpoints using bmotion
			cout << "DGmove: movemesh():  Moving the Delaunay graph\n";

			//testdgpoints.zeros();
			for(int i = 0; i < ndgpoin; i++)
			{
				for(int j = 0; j < ndim; j++)
					testdgpoints(i,j) = dgpoints(i,j) + bmot->get(i,j);
			}

			// check for validity of Delaunay graph
			movegraph(testdg, testdgpoints);
			testdg.compute_jacobians();
			check = testdg.detect_negative_jacobians();
			cout << "DGmove: movemesh():  Check is " << check << endl;
			if(check == true)								// if there exist elements with negative jacobian
			{
				testdgpoints = dgpoints;					// reset testdgpoints to original values
				testdg = dg;								// reset the test DG as well

				//halve the boundary motion
				for(int i = 0; i < bmot->rows(); i++)
					for(int j = 0; j < bmot->cols(); j++)
						(*bmot)(i,j) = bmot->get(i,j)/2;

				num_steps *= 2;
				reflag = true;
				continue;
			}
			else
			{
				movegraph(dg, testdgpoints);
				num_steps--;
			}

			// calculate displacements of DG nodes and calculate alpha
			calcalpha();

			// calculate new positions of interior points
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

			dgpoints = testdgpoints;

			if(num_steps > 0 && reflag == true)
			{
				delete [] Mtt;
				isallocMtt = false;

				cout << "DGmove: movemesh():  Re-triangulating\n";
				dg.clear();
				dg.setup(&dgpoints,ndgpoin);
				generateDG();
				numtri++;
				int contelem2;
				// re-calculate containing triangles
				for(int ipoin = 0; ipoin < ninpoin; ipoin++)
				{
					// first find containing DG element by "walking-through" the DG
					contelem2 = dg.find_containing_triangle(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);

					// store DG element in points
					points(ipoin,2) = contelem2;
				}
				calcMtt();
				reflag = false;
			}
			step++;
			cout << endl;
		}
		cout << "DGmove: movemesh(): Number of (re)triangulations needed: " << numtri << ".\n";
	}

	void movegraph(Delaunay2D& tdg, Matrix<double>& dgp)
	{
		//moves the DG according to the boundary motion - for visualizing deformed DG
		for(int i = 0; i < dgp.rows(); i++)
		{
			tdg.nodes[i].x = dgp(i,0);
			tdg.nodes[i].y = dgp(i,1);
		}
	}

	Matrix<double> getcoords()
	{
		// create a coords matrix with same point numbering as initial matrix and return it
		int a = 0, b = 0, k = 0;
		for(int i = 0; i < bflag.rows(); i++)
		{
			if(bflag(i) == 0)
			{
				for(int dim = 0; dim < ndim; dim++)
					newcoords(k,dim) = points(a,dim);
				k++;
				a++;
			}
			else
			{
				for(int dim = 0; dim < ndim; dim++)
					newcoords(k,dim) = dgpoints.get(b,dim);
				k++;
				b++;
			}
		}
		return newcoords;
	}
};

} // end namespace acfd
