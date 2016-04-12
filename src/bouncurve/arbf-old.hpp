/** A mesh-movement method using radial basis function interpolation based on the 2007 paper by de Boer, van der Schoot and Bijl,
   but with one major difference - the interpolation is done using only RBFs; the polynomial part is ignored.
   Aditya Kashi
   August 6, 2015

   Aug 14, 2015: Modified to work with curved mesh generation. Now requires interior points and boundary points as separate inputs.
 */

#ifndef __GLIBCXX_CSTDIO
#include <cstdio>
#endif

#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

using namespace std;
using namespace amat;

namespace acfd {

class RBFmove
{
	Matrix<double> inpoints;
	Matrix<double> bpoints;
	Matrix<double> bmotion;
	Matrix<int> bflag;
	int npoin;			// total number of points
	int ninpoin;		// number of interior points
	int nbpoin;			// number of boundary points
	int ndim;
	double (RBFmove::*rbf)(double);
	double srad;

	int nsteps;
	double tol;
	int maxiter;

	SpMatrix A;					// LHS matrix; contains RBF values relative to pairs of boundary points, and coordinates of boundary points
	Matrix<double>* coeffs;		// contains coefficients of RBFs and linear polynomial for each boundary point
	Matrix<double>* b;			// rhs for each of the dimensions; contains displacements of boundary points

public:

	RBFmove() { }

	RBFmove(Matrix<double>* int_points, Matrix<double>* boun_points, Matrix<double>* boundary_motion, int rbf_ch, double support_radius, int num_steps, double tolerance, int iter)
	// boundary_motion is npoin-by-ndim array - containing displacements corresponding to points in all_points (zero for interior points).
	// boundary_flags is a npoin-by-1 vector, containing 1 if the corresponding row in all_points is a boundary point and zero otherwise.
	{
		cout << "RBFmove: Storing inputs" << endl;
		inpoints = *int_points;
		bpoints = *boun_points;
		npoin = int_points->rows() + boun_points->rows();
		ndim = int_points->cols();
		nbpoin = bpoints.rows();

		int i, k;
		ninpoin = inpoints.rows();
		bmotion = *boundary_motion;
		A.setup(nbpoin,nbpoin);

		nsteps = num_steps;

		//checks

		switch(rbf_ch)
		{
			case(0): rbf = &RBFmove::rbf_c0;
			break;
			case(2): rbf = &RBFmove::rbf_c2_compact;
			break;
			case(4): rbf = &RBFmove::rbf_c4;
			default: rbf = &RBFmove::gaussian;
		}

		b = new Matrix<double>[ndim];		// RHS vectors of linear system for each coordinate direction
		coeffs = new Matrix<double>[ndim];
		for(i = 0; i < ndim; i++)
		{
			b[i].setup(nbpoin,1);
			coeffs[i].setup(nbpoin,1);
		}

		cout << "RBFmove: Preparing RHS vector of linear system to solve" << endl;
		// Fill the first nbpoin entries of b's with the motion for each step
		for(i = 0; i < nbpoin; i++)
		{
			for(int j = 0; j < ndim; j++)
				b[j](i) = bmotion.get(i,j)/nsteps;
		}
		tol = tolerance;
		maxiter = iter;
		srad = support_radius;
	}

	void setup(Matrix<double>* int_points, Matrix<double>* boun_points, Matrix<double>* boundary_motion, int rbf_ch, double support_radius, int num_steps, double tolerance, int iter)
	/// boundary_motion is npoin-by-ndim array - containing displacements corresponding to points in all_points (zero for interior points).
	/// boundary_flags is a npoin-by-1 vector, containing 1 if the corresponding row in all_points is a boundary point and zero otherwise.
	{
		cout << "RBFmove: Storing inputs" << endl;
		inpoints = *int_points;
		bpoints = *boun_points;
		npoin = int_points->rows() + boun_points->rows();
		ndim = int_points->cols();
		nbpoin = bpoints.rows();
		cout << "RBFmove: Number of boundary points " << nbpoin << endl;
		int i, k;
		ninpoin = inpoints.rows();
		bmotion = *boundary_motion;
		A.setup(nbpoin,nbpoin);

		nsteps = num_steps;

		//checks

		switch(rbf_ch)
		{
			case(0): rbf = &RBFmove::rbf_c0;
			break;
			case(2): rbf = &RBFmove::rbf_c2_compact;
			break;
			case(4): rbf = &RBFmove::rbf_c4;
			default: rbf = &RBFmove::gaussian;
		}

		b = new Matrix<double>[ndim];		// RHS vectors of linear system for each coordinate direction
		coeffs = new Matrix<double>[ndim];
		for(i = 0; i < ndim; i++)
		{
			b[i].setup(nbpoin,1);
			coeffs[i].setup(nbpoin,1);
		}

		cout << "RBFmove: Preparing RHS vector of linear system to solve" << endl;
		// Fill the first nbpoin entries of b's with the motion for each step
		for(i = 0; i < nbpoin; i++)
		{
			for(int j = 0; j < ndim; j++)
				b[j](i) = bmotion.get(i,j)/nsteps;
		}
		tol = tolerance;
		maxiter = iter;
		srad = support_radius;
		cout << "RBFmove: RBF to use: " << rbf_ch << endl;
		cout << "RBFmove: Support radius = " << srad << endl;
		cout << "RBFmove: Number of steps = " << nsteps << endl;
	}

	~RBFmove()
	{
		delete [] b;
		delete [] coeffs;
	}

	// RBFs
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

	void assembleLHS()
	{
		cout << "RBFmove:  assembleLHS(): assembling LHS matrix" << endl;
		int i, j;
		double dist;
		double temp;

		SpMatrix* A = &(RBFmove::A);
		Matrix<double>* bpoints = &(RBFmove::bpoints);
		double (RBFmove::*rbfunc)(double) = rbf;
		int nbpoin = RBFmove::nbpoin;
		int ndim = RBFmove::ndim;

		//cout << (this->*rbfunc)(0.0) << endl;

		// set the top nbpoin-by-nbpoin elements of A, ie, M_bb
		//cout << "RBFmove:  assembleLHS(): assembling M_bb" << endl;
		#pragma omp parallel for default(none) private(i,j,dist,temp) shared(A,bpoints,rbfunc,nbpoin,ndim)
		for(i = 0; i < nbpoin; i++)
		{
			A->set_fast(i,i, (this->*rbfunc)(0.0));			// set diagonal element in row i

			//cout << "RBF value = " << (this->*rbf)(0.0) << endl;
			for(j = i+1; j < nbpoin; j++)		// traverse lower triangular matrix
			{
				dist = 0;
				for(int id = 0; id < ndim; id++)
					dist += (bpoints->get(i,id) - bpoints->get(j,id))*(bpoints->get(i,id) - bpoints->get(j,id));
				dist = sqrt(dist);
				temp = (this->*rbfunc)(dist);
				if(dabs(temp) > tol*tol)
				{
					A->set_fast(i,j, temp);
					A->set_fast(j,i, temp);
				}
			}
		}

		/*ofstream outf("A.dat");
		A->fprint(outf);
		outf.close();*/

		/*cout << "RBFmove:  assembleLHS(): assembling P_b" << endl;
		// set P_b and P_b transpose
		for(i = 0; i < nbpoin; i++)
		{
			A.set(i,nbpoin, 1.0);						// the first column of P_b is 1
			A.set(nbpoin,i, 1.0);
			for(j = nbpoin+1; j < nbpoin+ndim+1; j++)
			{
				A.set(i,j, bpoints.get(i,j-nbpoin-1));
				A.set(j,i, bpoints.get(i,j-nbpoin-1));
			}
		}*/
	}

	void move_step()
	{
		// solve for RBF coefficients
		cout << "RBFmove:  move_step(): Solving linear system" << endl;
		Matrix<double> xold(nbpoin,1);
		xold.zeros();
		for(int idim = 0; idim < ndim; idim++)
			coeffs[idim] = sparseCG_d(&A, b[idim], xold, tol, maxiter);
			//coeffs[idim] = sparsePCG(&A, b[idim], xold, "jacobi", tol, maxiter);
			//coeffs[idim] = sparsegaussseidel(&A, b[idim], xold, tol, maxiter);
			//coeffs[idim] = sparseSOR(&A, b[idim], xold, tol, maxiter);

		cout << "RBFmove:  move_step(): Moving interior points" << endl;
		// calculate new positions of interior points
		int i;
		Matrix<double>* co = coeffs;			// first assign local pointers to class variables for OpenMP etc
		Matrix<double>* bp = &bpoints;
		Matrix<double>* ip = &inpoints;
		double (RBFmove::*rbfunc)(double) = rbf;
		int nbpoin = RBFmove::nbpoin;
		int ninpoin = RBFmove::ninpoin;
		int ndim = RBFmove::ndim;

		#pragma omp parallel for default(none) private(i) shared(co, bp, ip, rbfunc, nbpoin, ninpoin, ndim)
		for(i = 0; i < ninpoin; i++)
		{
			double* sum = new double[ndim];		// for storing sum of RBFs corresponding to an interior point
			double* psum = new double[ndim];	// for storing value of linear polynomial corresponding to an interior point

			int j;
			for(j = 0; j < ndim; j++)
			{
				sum[j] = 0;
				psum[j] = 1.0;
			}

			// get RBF part
			for(j = 0; j < nbpoin; j++)
			{
				double dist = 0;
				for(int idim = 0; idim < ndim; idim++)
					dist += (ip->get(i,idim)-bp->get(j,idim))*(ip->get(i,idim)-bp->get(j,idim));
				dist = sqrt(dist);

				for(int idim = 0; idim < ndim; idim++)
					sum[idim] += co[idim].get(j) * (this->*rbfunc)(dist);
			}

			//get polynomial part
			/*for(int idim = 0; idim < ndim; idim++)
				for(int jdim = 0; jdim < ndim; jdim++)
					psum[idim] += co[idim].get(nbpoin+1+jdim) * ip->get(i, jdim);*/

			for(int idim = 0; idim < ndim; idim++)
			{
				(*ip)(i,idim) += sum[idim];
			}

			delete [] sum;
			delete [] psum;

			if(i % 100 == 0)
				printf(" =");
		}
		printf("\n");
	}

	void move()
	{
		int istep;

		for(istep = 0; istep < nsteps; istep++)
		{
			cout << "RBFmove: move(): Step " << istep << endl;
			assembleLHS();

			move_step();

			// move buondary points
			for(int i = 0; i < nbpoin; i++)
				for(int j = 0; j < ndim; j++)
					bpoints(i,j) += b[j](i);
		}
	}

	Matrix<double> getInteriorPoints()
	{
		return inpoints;
	}

	Matrix<double> getBoundaryPoints()
	{
		return bpoints;
	}
};

} // end namespace
