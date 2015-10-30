/* Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
Aditya Kashi
July 1, 2015

Changelog:
Oct 7, 2015: Changed the way s (prescribed boundary motion at each DG point) is calculated in calcalpha().
				The earlier way was probably wrong - bmotion was used instead of bmotionb.
Oct 7, 2015: Added a similar class, DGRBFrotate, for rotation interpolation.
*/

#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

#ifndef __ABOWYERWATSON_H
#include <abowyerwatson.hpp>
#endif

#ifndef ZERO_TOL
	#define ZERO_TOL 1e-14
#endif

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
	Matrix<double>* bmotion;						// boundary motion of all points (0 for interior points)
	Matrix<double> bmotionb;						// boundary motion of each boundary points
	Matrix<int> bflag;

	Matrix<double>* Mtt;
	Matrix<double>* alpha;
	Matrix<double>* s;
	bool isallocMtt;
	bool isallocalpha;
	double (DGRBFmove::*rbf)(double);
	double srad;			// support radius for RBFs

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
		points.setup(ninpoin, ndim+1);	// for each interior point, store coords and containing DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.

		int k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 0)
			{
				for(int j = 0; j < ndim; j++)
					points(k,j) = coords->get(i,j);
				k++;
			}

		// now get dgpoints and bmotionb using bflags
		dgpoints.setup(ndgpoin,ndim);
		bmotionb.setup(ndgpoin,ndim);
		k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1)
			{
				for(int j = 0; j < ndim; j++){
					dgpoints(k,j) = coords->get(i,j);
					bmotionb(k,j) = bmotion->get(i,j);
				}
				k++;
			}

		bflag = bflags;
		newcoords.setup(bflags.rows(),ndim);
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
		if(dabs(ddet) < ZERO_TOL) cout << "! DGRBFmove: cramer3(): Matrix A is singular!!\n";
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
		for(int iel = 0; iel < ndgelem; iel++)
		{
			// first, calculate s
			for(int i = 0; i < ndim+1; i++)
			{
				for(int idim = 0; idim < ndim; idim++)
					s[iel](i,idim) = bmotionb.get(dginpoel(iel,i),idim);
			}

			// calculate alpha
			Matrix<double> x0(ndim+1,1); x0.zeros();

			for(int idim = 0; idim < ndim; idim++)
			{
				alpha[iel] = cramer3(Mtt[iel],s[iel]);
			}
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

		// calculate new positions of interior points by mapping them to deformed DG elements using area coordinates calculated before
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
		int  k = 0;
		for(int i = 0; i < bmotion->rows(); i++)
		{
			if(bflag(i) == 1)
			{
				for(int j = 0; j < ndim; j++)
					dgpoints(k,j) += bmotion->get(i,j);
				k++;
			}
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

/** Implements rotation interpolation (interpolation of angle of rotation) by the DG-RBF method. */
class DGRBFrotate
{
	int ndim;								///< dimension of geometry
	int rdim;								///< number of components of rotation required (1 for 2d, 3 for 3d)
	int npoin;								///< seems unused...
	int ninpoin;							///< Number of interior points to be moved
	Matrix<double> points;					///< Holds coords of interior points, as well as their containing DG elements.
	int ndgpoin;							///< Number of boundary points, or number of points in the DG
	Matrix<double> dgpoints;				///< Holds coordinates of points in the DG
	int ndgelem;							///< Number of Delaunay elements (elements in the DG)
	Matrix<int> dginpoel;					///< Element connectivity matrix of the Delaunay Graph (DG)
	Matrix<int> dgesuel;					///< Elements surrounding elements for DG
	Matrix<double> dgintfac;				///< Face data structure for DG
	Matrix<double>* bmotion;				///< prescribed boundary motion - angles, in this case, for each point of orginal mesh (zero for interior points).
	Matrix<double> bmotionb;				///< prescribed boundary motion for each DG node.
	Matrix<int> bflag;						///< vector of flags; 1 if corresponding point is a boundary point.

	Matrix<double>* Mtt;					///< RBF LHS matrices for each DG element.
	Matrix<double>* alpha;					///< RBF coefficients for each DG element.
	Matrix<double>* s;						///< Prescribed motion of nodes of each DG element.
	bool isallocMtt;
	bool isallocalpha;
	double (DGRBFrotate::*rbf)(double);		///< Pointer to the specific RBF function to be used
	double srad;							///< support radius for RBFs

	vector<double> rc;						///< coordinates of centre of rotation

public:
	Delaunay2D dg;							///< Object of Delaunay2D class that will create the Delaunay graph (DG) and iterate through it.
	Matrix<double> newcoords;				///< Final coordinates after movement.

	/** Constructor */
	DGRBFrotate(int num_dimn, Matrix<double>* coords, Matrix<int> bflags, Matrix<double>* boundary_angles, vector<double> rce, int rbf_type, double support_radius)
	{
		/** Note that bflags contains as many rows as points in the mesh and contains 1 if the point is a boundary point and 0 otherwise.
			boundary_angles has as many rows as points in the original mesh and contains (in-plane) rotation values for each point (in radians) (zero for interior points).
			It is assumed that boundary_angles has rdim columns (1 for 2d, 3 for 3d).
			rce contains coordinates of centre of rotation.
		*/
		
		bmotion = boundary_angles;
		if(bmotion->rows() != bflags.rows()) cout << "DGMove: !! Error: No. of rows in bflags and bmotion are not equal!\n";
		ndim = coords->cols();
		rc = rce;

		ndgpoin = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1) ndgpoin++;

		ninpoin = coords->rows() - ndgpoin;

		points.setup(ninpoin, ndim+1);	// for each interior point, store coords and containing DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.

		// effective dimension for alpha and s is not necessarily ndim - for 2d, we need only one rotation component.
		if(ndim == 2)
			rdim = 1;
		else
			rdim = ndim;
		
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
		bmotionb.setup(ndgpoin,rdim);
		k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1)
			{
				for(int j = 0; j < ndim; j++) 
				{
					dgpoints(k,j) = coords->get(i,j);
					if(j < rdim)
						bmotionb(k,j) = bmotion->get(i,j);
				}
				k++;
			}
		//bmotionb.mprint();

		bflag = bflags;
		newcoords.setup(bflags.rows(),ndim);
		dg.setup(&dgpoints, ndgpoin);

		isallocMtt = false;
		isallocalpha = false;
		switch(rbf_type)
		{
			case(0): rbf = &DGRBFrotate::rbf_c0;
			break;
			case(2): rbf = &DGRBFrotate::rbf_c2_compact;
			break;
			case(4): rbf = &DGRBFrotate::rbf_c4;
			default: rbf = &DGRBFrotate::rbf_c2_compact;
		}

		srad = support_radius;

	}

	~DGRBFrotate()
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
				alpha[i].setup(ndim+1,rdim);
				s[i].setup(ndim+1,rdim);
			}
			isallocalpha = true;
		}
	}

	/// Computes determinant of 3x3 matrix.
	double det3(Matrix<double> A)
	{
		double d;
		d = A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) - A(0,1)*(A(1,0)*A(2,2)-A(2,0)*A(1,2)) + A(0,2)*(A(1,0)*A(2,1)-A(2,0)*A(1,1));
		return d;
	}

	/// Computes the solution of a 3x3 linear system using det3().
	Matrix<double> cramer3(Matrix<double> A, Matrix<double> b)
	{
		double ddet = 0;
		Matrix<double> x(b.rows(), b.cols());
		Matrix<double> Anum(A.rows(),A.cols());
		ddet = det3(A);
		if(dabs(ddet) < ZERO_TOL) cout << "! DGRBFmove: cramer3(): Matrix A is singular!!\n";
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

	/// Computes LHS matrix for computing RBF coefficients.
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

	/** Computes RBF coefficients uing Mtt.
		Uses bmotionb to get rotation angle of DG nodes.
		Uses cramer3 to solve the linear system associated with each DG element.
	*/
	void calcalpha()
	{
		for(int iel = 0; iel < ndgelem; iel++)
		{
			// first, calculate s
			for(int i = 0; i < ndim+1; i++)
			{
				for(int idim = 0; idim < rdim; idim++)
					s[iel](i,idim) = bmotionb.get(dginpoel(iel,i),idim);
			}

			// calculate alpha
			alpha[iel] = cramer3(Mtt[iel],s[iel]);
		}
	}

	void movemesh()
	{
		// cycle over interior points
		int contelem;
		cout << "DGmove: movemesh(): Calculating containing elements for each interior point\n";
		// This loop can be parallelized
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			// first find containing DG element by "walking-through" the DG
			contelem = dg.find_containing_triangle(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);

			// store DG element in points
			points(ipoin,2) = contelem;
		}

		/// Execute calcMtt() and calcalpha() to get RBF LHS and coefficients.

		calcMtt();

		calcalpha();

		/// Calculate new positions of interior points by mapping them to deformed DG elements using 
		///		RBFs and their coeffecients (alpha).
		Matrix<double> A(ndim+1,1);
		cout << "DGmove: movemesh():  Moving the interior points\n";
		int elem; 
		double* rr = new double[ndim];
		double sum = 0;
		// This loop can be parallelized
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			// get the containing element
			elem = points(ipoin,2);

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
			for(int idim = 0; idim < rdim; idim++)
			{
				// get rotation angle
				rr[idim] = 0.0;
				for(int inode = 0; inode < ndim+1; inode++)
					rr[idim] += A(inode)*alpha[elem](inode, idim);

			}
			
			//calculate new position
			if(ndim == 2) {
				// the formulae below are taken from the paper on DGRBF by Wang, Qin and Zhao.
				//points(ipoin,0) = (points(ipoin,0) - rc[0])*cos(rr[0]) + (points(ipoin,1)-rc[1])*sin(rr[0]) + rc[0];
				//points(ipoin,1) = (points(ipoin,1) - rc[1])*cos(rr[0]) - (points(ipoin,0)-rc[0])*sin(rr[0]) + rc[1];

				// the ones below are derived by me (just rotation in opposite direction to above)
				points(ipoin,0) = (points(ipoin,0) - rc[0])*cos(rr[0]) - (points(ipoin,1)-rc[1])*sin(rr[0]) + rc[0];
				points(ipoin,1) = (points(ipoin,1) - rc[1])*cos(rr[0]) + (points(ipoin,0)-rc[0])*sin(rr[0]) + rc[1];
			}
			else {
				cout << "DGRBFrotate: movemesh(): ! Position update not implemented for 3d!!" << endl;
				// TODO: figure out and implement for 3d
			}
		}
		delete [] rr;
		
		// update coordinates in dgpoints using bmotionb
		cout << "DGRBFrotate: movemesh(): Moving the Delaunay graph\n";
		cout << "DGRBFrotate: movemesh(): Centre of rotation " << rc[0] << " " << rc[1] << endl;
		for(int ipoin = 0; ipoin < ndgpoin; ipoin++)
		{
			if(ndim == 2) {
				dgpoints(ipoin,0) = (dgpoints(ipoin,0) - rc[0])*cos(bmotionb(ipoin)) - (dgpoints(ipoin,1)-rc[1])*sin(bmotionb(ipoin)) + rc[0];
				dgpoints(ipoin,1) = (dgpoints(ipoin,1) - rc[1])*cos(bmotionb(ipoin)) + (dgpoints(ipoin,0)-rc[0])*sin(bmotionb(ipoin)) + rc[1];
			}
			else {
				cout << "DGRBFrotate: movemesh(): ! Position update not implemented for 3d!!" << endl;
				// TODO: figure out and implement for 3d
			}
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

/** Combines rotation interpolation and displacement interpolation to get Wang, Qin and Zhao's DGRBF2 method.
*/
class DGRBF2
{
	int ndim;								///< dimension of geometry
	int rdim;								///< number of components of rotation required (1 for 2d, 3 for 3d)
	int npoin;								///< seems unused...
	int ninpoin;							///< Number of interior points to be moved
	Matrix<double> points;					///< Holds coords of interior points, as well as their containing DG elements.
	int ndgpoin;							///< Number of boundary points, or number of points in the DG
	Matrix<double> dgpoints;				///< Holds coordinates of points in the DG
	int ndgelem;							///< Number of Delaunay elements (elements in the DG)
	Matrix<int> dginpoel;					///< Element connectivity matrix of the Delaunay Graph (DG)
	Matrix<int> dgesuel;					///< Elements surrounding elements for DG
	Matrix<double> dgintfac;				///< Face data structure for DG
	Matrix<double>* bmotiona;				///< prescribed boundary motion - angles, for each point of orginal mesh (zero for interior points).
	Matrix<double>* bmotiond;				///< prescribed boundary motion - displacements
	Matrix<double> bmotionba;				///< prescribed boundary motion (angle) for each DG node.
	Matrix<double> bmotionbd;				///< prescribed boundary motion (displacement) for each DG node.
	Matrix<int> bflag;						///< vector of flags; 1 if corresponding point is a boundary point.

	Matrix<double>* Mtt;					///< RBF LHS matrices for each DG element.
	Matrix<double>* alphaa;					///< RBF angle coefficients for each DG element.
	Matrix<double>* sa;						///< Prescribed angle of nodes of each DG element.
	Matrix<double>* alphad;					///< RBF displacement coefficients for each DG element
	Matrix<double>* sd;						///< Prescribed displacements of nodes for each DG element
	bool isallocMtt;
	bool isallocalpha;
	double (DGRBF2::*rbf)(double);		///< Pointer to the specific RBF function to be used
	double srad;							///< support radius for RBFs

	vector<double> rc;						///< coordinates of centre of rotation
	double tol;								///< for detecting nonsingularity while solving systems by Cramer's rule.

public:
	Delaunay2D dg;
	Matrix<double> newcoords;

	/** Constructor */
	DGRBF2(int num_dimn, Matrix<double>* coords, Matrix<int> bflags, Matrix<double>* boundary_angles, Matrix<double>* boundary_displ, vector<double> rce, int rbf_type, double support_radius)
	{
		/** Note that bflags contains as many rows as points in the mesh and contains 1 if the point is a boundary point and 0 otherwise.
			boundary_angles and boundary_displ have as many rows as points in the original mesh.
			They contain (in-plane) rotation values for each point (in radians) and displacement for each point respectively (zero for interior points).
			It is assumed that boundary_angles has rdim columns (1 for 2d, 3 for 3d). boundary_displ has ndim columns.
			rce contains coordinates of centre of rotation
		*/
		
		bmotiona = boundary_angles;
		bmotiond = boundary_displ;
		if(bmotiond->rows() != bflags.rows()) cout << "DGMove: !! Error: No. of rows in bflags and bmotion are not equal!\n";
		ndim = coords->cols();
		rc = rce;

		ndgpoin = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1) ndgpoin++;

		ninpoin = coords->rows() - ndgpoin;

		points.setup(ninpoin, ndim+1);	// for each interior point, store coords and containing DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.

		// effective dimension for alpha and s is not necessarily ndim - for 2d, we need only one rotation component.
		if(ndim == 2)
			rdim = 1;
		else
			rdim = ndim;
		
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
		bmotionba.setup(ndgpoin,rdim);
		bmotionbd.setup(ndgpoin,ndim);
		k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1)
			{
				for(int j = 0; j < ndim; j++) 
				{
					dgpoints(k,j) = coords->get(i,j);
					if(j < rdim)
						bmotionba(k,j) = bmotiona->get(i,j);
					if(j < ndim)
						bmotionbd(k,j) = bmotiond->get(i,j);
				}
				k++;
			}

		bflag = bflags;
		newcoords.setup(bflags.rows(),ndim);
		dg.setup(&dgpoints, ndgpoin);

		isallocMtt = false;
		isallocalpha = false;
		switch(rbf_type)
		{
			case(0): rbf = &DGRBF2::rbf_c0;
			break;
			case(2): rbf = &DGRBF2::rbf_c2_compact;
			break;
			case(4): rbf = &DGRBF2::rbf_c4;
			default: rbf = &DGRBF2::rbf_c2_compact;
		}

		srad = support_radius;
		tol = 1e-15;

		//cout << "Test: " << rbf(0)<< ' ' << sqrt(0.0) << endl;
	}

	~DGRBF2()
	{
		if(isallocMtt == true) delete [] Mtt;
		if(isallocalpha == true) { delete [] alphaa; delete [] sa; delete [] alphad; delete [] sd; }
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

		if(isallocMtt == false)
		{
			Mtt = new Matrix<double>[ndgelem];
			for(int i = 0; i < ndgelem; i++)
				Mtt[i].setup(ndim+1, ndim+1);
			isallocMtt = true;
		}
		if(isallocalpha == false)
		{
			alphaa = new Matrix<double>[ndgelem];
			sa = new Matrix<double>[ndgelem];
			alphad = new Matrix<double>[ndgelem];
			sd = new Matrix<double>[ndgelem];
			for(int i = 0; i < ndgelem; i++)
			{
				alphaa[i].setup(ndim+1,rdim);
				sa[i].setup(ndim+1,rdim);
				alphad[i].setup(ndim+1,ndim);
				sd[i].setup(ndim+1,ndim);
			}
			isallocalpha = true;
		}
	}

	/// Computes determinant of 3x3 matrix.
	double det3(Matrix<double> A)
	{
		double d;
		d = A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) - A(0,1)*(A(1,0)*A(2,2)-A(2,0)*A(1,2)) + A(0,2)*(A(1,0)*A(2,1)-A(2,0)*A(1,1));
		return d;
	}

	/// Computes the solution of a 3x3 linear system using det3().
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

	/// Computes LHS matrix for computing RBF coefficients.
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

	/** Computes RBF coefficients uing Mtt.
		Uses bmotionb to get motion of DG nodes.
		Uses cramer3 to solve the linear system associated with each DG element.
	*/
	void calcalpha()
	{
		for(int iel = 0; iel < ndgelem; iel++)
		{
			// first, calculate s
			for(int i = 0; i < ndim+1; i++)
			{
				for(int idim = 0; idim < rdim; idim++)
					sa[iel](i,idim) = bmotionba.get(dginpoel(iel,i),idim);
			}

			// calculate alpha
			alphaa[iel] = cramer3(Mtt[iel],sa[iel]);
		}

		for(int iel = 0; iel < ndgelem; iel++)
		{
			// first, calculate s
			for(int i = 0; i < ndim+1; i++)
			{
				for(int idim = 0; idim < ndim; idim++)
					sd[iel](i,idim) = bmotionbd.get(dginpoel(iel,i),idim);
			}

			// calculate alpha
			alphad[iel] = cramer3(Mtt[iel],sd[iel]);
		}
	}

	void movemesh()
	{
		// cycle over interior points
		int contelem;
		cout << "DGmove: movemesh(): Calculating containing elements for each interior point\n";
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			// first find containing DG element by "walking-through" the DG
			contelem = dg.find_containing_triangle(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);

			// store DG element in points
			points(ipoin,2) = contelem;
		}

		/// Execute calcMtt() and calcalpha() to get RBF LHS and coefficients.

		calcMtt();


		calcalpha();

		/// Calculate new positions of interior points by mapping them to deformed DG elements using 
		///		RBFs and their coeffecients (alpha).
		Matrix<double> A(ndim+1,1);
		cout << "DGmove: movemesh():  Moving the interior points\n";
		int elem; 
		double* rra = new double[ndim];
		double* rrd = new double[ndim];
		double sum = 0;
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			elem = points(ipoin,2);

			// calculate RBFs of point ipoin and store in A
			for(int inode = 0; inode < ndim+1; inode++)
			{
				sum = 0;
				for(int idim = 0; idim < ndim; idim++)
				{
					rrd[idim] = dgpoints.get(dginpoel(elem,inode),idim);
					sum += (rrd[idim] - points(ipoin,idim)) * (rrd[idim] - points(ipoin,idim));
				}

				A(inode) = (this->*rbf)( sqrt(sum) );
			}

			// Calculate displacement of ipoin using A and alpha, and update coordinates of point
			for(int idim = 0; idim < rdim; idim++)
			{
				// get rotation angle
				rra[idim] = 0.0;
				for(int inode = 0; inode < ndim+1; inode++)
					rra[idim] += A(inode)*alphaa[elem](inode, idim);
			}

			for(int idim = 0; idim < ndim; idim++)
			{
				//get discplacements
				rrd[idim] = 0.0;
				for(int inode = 0; inode < ndim+1; inode++)
					rrd[idim] += A(inode)*alphad[elem](inode, idim);
			}
			
			//calculate new position based on rotation
			if(ndim == 2) {
				// the formulae below are taken from the paper on DGRBF by Wang, Qin and Zhao.
				//points(ipoin,0) = (points(ipoin,0) - rc[0])*cos(rr[0]) + (points(ipoin,1)-rc[1])*sin(rr[0]) + rc[0];
				//points(ipoin,1) = (points(ipoin,1) - rc[1])*cos(rr[0]) - (points(ipoin,0)-rc[0])*sin(rr[0]) + rc[1];

				// the ones below are derived by me
				points(ipoin,0) = (points(ipoin,0) - rc[0])*cos(rra[0]) - (points(ipoin,1)-rc[1])*sin(rra[0]) + rc[0];
				points(ipoin,1) = (points(ipoin,1) - rc[1])*cos(rra[0]) + (points(ipoin,0)-rc[0])*sin(rra[0]) + rc[1];
			}
			else {
				cout << "DGRBFrotate: movemesh(): ! Position update not implemented for 3d!!" << endl;
				// TODO: figure out and implement for 3d
			}

			// add displacement to rotated position
			for(int idim = 0; idim < ndim; idim++)
			{
				points(ipoin,idim) += rrd[idim];
			}
		}
		delete [] rra;
		delete [] rrd;

		// update coordinates in dgpoints using bmotionb
		cout << "DGmove: movemesh(): Moving the Delaunay graph\n";
		for(int ipoin = 0; ipoin < ndgpoin; ipoin++)
		{
			if(ndim == 2) {
				dgpoints(ipoin,0) = (dgpoints(ipoin,0) - rc[0])*cos(bmotionba(ipoin)) - (dgpoints(ipoin,1)-rc[1])*sin(bmotionba(ipoin)) + rc[0];
				dgpoints(ipoin,1) = (dgpoints(ipoin,1) - rc[1])*cos(bmotionba(ipoin)) + (dgpoints(ipoin,0)-rc[0])*sin(bmotionba(ipoin)) + rc[1];
			}
			else {
				cout << "DGRBFrotate: movemesh(): ! Position update not implemented for 3d!!" << endl;
				// TODO: figure out and implement for 3d
			}
			
			for(int idim = 0; idim < ndim; idim++)
				dgpoints(ipoin,idim) += bmotionbd(ipoin,idim);
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
