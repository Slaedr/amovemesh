/** A mesh-movement method using radial basis function interpolation based on the 2007 paper by de Boer, van der Schoot and Bijl,
 * but with one major difference - the interpolation is done using only RBFs; the polynomial part is ignored.
 * @author Aditya Kashi
 * @date August 6, 2015
 * 
 * Aug 14, 2015: Modified to work with curved mesh generation. Now requires interior points and boundary points as separate inputs.
 */

#ifndef __GLIBCXX_CSTDIO
#include <cstdio>
#endif

#ifndef __GLIBCXX_STRING
#include <string>
#endif

#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

#define __ARBF_H 1

namespace amc {

class RBFmove
{
	amat::Matrix<double> inpoints;
	amat::Matrix<double> bpoints;
	amat::Matrix<double> bmotion;
	amat::Matrix<int> bflag;
	int npoin;			///< total number of points
	int ninpoin;		///< number of interior points
	int nbpoin;			///< number of boundary points
	int ndim;
	double (RBFmove::*rbf)(double);
	double srad;

	int nsteps;			///< Number of steps in which to carry out the movement. More steps lead to better results upto a certain number of steps.
	double tol;
	int maxiter;

	amat::SpMatrix A;					///< LHS matrix; contains RBF values relative to pairs of boundary points, and coordinates of boundary points
	amat::Matrix<double>* coeffs;		///< contains coefficients of RBFs and linear polynomial for each boundary point
	amat::Matrix<double>* b;			///< rhs for each of the dimensions; contains displacements of boundary points
	bool isalloc;				///< This flag is true if both [b](@ref b) and [coeffs](@ref coeffs) have been allocated
	
	/// string indicating the solver to use - options are 'CG', 'SOR', 'BICGSTAB' or 'LU' (defaults to LU)
	std::string lsolver;

public:

	/// No-arg constructor
	RBFmove();

	RBFmove(amat::Matrix<double>* int_points, amat::Matrix<double>* boun_points, amat::Matrix<double>* boundary_motion, int rbf_ch, double support_radius, int num_steps, 
			double tolerance, int iter, std::string linear_solver);
	///< boundary_motion is nbpoin-by-ndim array - containing displacements corresponding to boundary points.

	/// Sets the data needed
	/** Note that all parameters are deep-copied.
	 * \param int_points is a list of all interior points to be moved
	 * \param boun_points is the array of boundary points
	 * \param boundary_motion is nbpoin-by-ndim array - containing displacements corresponding to boundary points.
	 * \param rbf_ch indicates the RBF to use - 0 : C0, 2 : C2, 4 : C4, default : Gaussian
	 * \param num_steps is the number of steps in which to break up the movement to perform separately (sequentially)
	 * \param linear_solver indicates the linear solver to use to solve the RBF equations - "CG" or "LU"
	 */
	void setup(amat::Matrix<double>* int_points, amat::Matrix<double>* boun_points, amat::Matrix<double>* boundary_motion, int rbf_ch, double support_radius, int num_steps, 
			double tolerance, int iter, std::string linear_solver);

	~RBFmove();

	/// Specific RBFs
	double rbf_c2_compact(double xi);
	double rbf_c0(double xi);
	double rbf_c4(double xi);
	double gaussian(double xi);

	/// Assembles the LHS matrix.
	void assembleLHS();

	/// Executes 1 step of the mesh movement.
	void move_step();

	/// Uses [move_step](@ref move_step) to execute the total number of steps specified.
	/** This function calls all other required functions, so it should be used directly after setup.
	 */
	void move();

	/// Returns new positions of interior points.
	amat::Matrix<double> getInteriorPoints();

	/// Returns new positions of boundary points.
	amat::Matrix<double> getBoundaryPoints();
};
	
RBFmove::RBFmove() {isalloc = false; }

RBFmove::RBFmove(amat::Matrix<double>* int_points, amat::Matrix<double>* boun_points, amat::Matrix<double>* boundary_motion, int rbf_ch, double support_radius, int num_steps, double tolerance, int iter, std::string linear_solver)
// boundary_motion is nbpoin-by-ndim array - containing displacements corresponding to boundary points.
{
	std::cout << "RBFmove: Storing inputs" << std::endl;
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

	switch(rbf_ch)
	{
		case(0): rbf = &RBFmove::rbf_c0;
		break;
		case(2): rbf = &RBFmove::rbf_c2_compact;
		break;
		case(4): rbf = &RBFmove::rbf_c4;
		default: rbf = &RBFmove::rbf_c2_compact;
	}

	b = new amat::Matrix<double>[ndim];		// RHS std::vectors of linear system for each coordinate direction
	coeffs = new amat::Matrix<double>[ndim];
	for(i = 0; i < ndim; i++)
	{
		b[i].setup(nbpoin,1);
		coeffs[i].setup(nbpoin,1);
	}
	isalloc = true;

	std::cout << "RBFmove: Preparing RHS vector of linear system to solve" << std::endl;
	// Fill the first nbpoin entries of b's with the motion for each step
	for(i = 0; i < nbpoin; i++)
	{
		for(int j = 0; j < ndim; j++)
			b[j](i) = bmotion.get(i,j)/nsteps;
	}
	tol = tolerance;
	maxiter = iter;
	srad = support_radius;
	lsolver = linear_solver;
}

void RBFmove::setup(amat::Matrix<double>* int_points, amat::Matrix<double>* boun_points, amat::Matrix<double>* boundary_motion, int rbf_ch, double support_radius, int num_steps, double tolerance, int iter, std::string linear_solver)
// boundary_motion is nbpoin-by-ndim array - containing displacements corresponding to boundary points.
{
	std::cout << "RBFmove: Storing inputs" << std::endl;
	inpoints = *int_points;
	bpoints = *boun_points;
	npoin = int_points->rows() + boun_points->rows();
	ndim = int_points->cols();
	nbpoin = bpoints.rows();
	std::cout << "RBFmove: Number of boundary points " << nbpoin << std::endl;
	int i, k;
	ninpoin = inpoints.rows();
	bmotion = *boundary_motion;
	A.setup(nbpoin,nbpoin);

	nsteps = num_steps;

	switch(rbf_ch)
	{
		case(0): rbf = &RBFmove::rbf_c0;
		break;
		case(2): rbf = &RBFmove::rbf_c2_compact;
		break;
		case(4): rbf = &RBFmove::rbf_c4;
		default: rbf = &RBFmove::gaussian;
	}

	b = new amat::Matrix<double>[ndim];		// RHS std::vectors of linear system for each coordinate direction
	coeffs = new amat::Matrix<double>[ndim];
	for(i = 0; i < ndim; i++)
	{
		b[i].setup(nbpoin,1);
		coeffs[i].setup(nbpoin,1);
	}
	isalloc = true;

	std::cout << "RBFmove: Preparing RHS vector of linear system to solve" << std::endl;
	// Fill the first nbpoin entries of b's with the motion for each step
	for(i = 0; i < nbpoin; i++)
	{
		for(int j = 0; j < ndim; j++)
			b[j](i) = bmotion.get(i,j)/nsteps;
	}
	
	tol = tolerance;
	maxiter = iter;
	srad = support_radius;
	lsolver = linear_solver;
	
	std::cout << "RBFmove: RBF to use: " << rbf_ch << std::endl;
	std::cout << "RBFmove: Support radius = " << srad << std::endl;
	std::cout << "RBFmove: Number of steps = " << nsteps << std::endl;
}

RBFmove::~RBFmove()
{
	if(isalloc) {
		delete [] b;
		delete [] coeffs;
	}
}

// RBFs
double RBFmove::rbf_c2_compact(double xi)
{
	if(xi/srad <= 1.0)
		return pow(1-xi/srad,4)*(4*xi/srad+1);
	else return 0;
}
double RBFmove::rbf_c0(double xi)
{
	return (1-xi)*(1-xi);
}
double RBFmove::rbf_c4(double xi)
{
	return pow(1-xi,6)*(35*xi*xi + 18*xi + 3);
}
double RBFmove::gaussian(double xi)
{
	return exp(-xi*xi);
}

void RBFmove::assembleLHS()
{
	std::cout << "RBFmove:  assembleLHS(): assembling LHS matrix" << std::endl;
	int i, j;
	double dist;
	double temp;

	amat::SpMatrix* A = &(RBFmove::A);
	amat::Matrix<double>* bpoints = &(RBFmove::bpoints);
	double (RBFmove::*rbfunc)(double) = rbf;
	int nbpoin = RBFmove::nbpoin;
	int ndim = RBFmove::ndim;

	// set the top nbpoin-by-nbpoin elements of A, ie, M_bb
	//std::cout << "RBFmove:  assembleLHS(): assembling M_bb" << std::endl;
	//#pragma omp parallel for default(none) private(i,j,dist,temp) shared(A,bpoints,rbfunc,nbpoin,ndim)
	for(i = 0; i < nbpoin; i++)
	{
		A->set(i,i, (this->*rbfunc)(0.0));			// set diagonal element in row i

		//std::cout << "RBF value = " << (this->*rbf)(0.0) << std::endl;
		for(j = i+1; j < nbpoin; j++)		// traverse lower triangular matrix
		{
			dist = 0;
			for(int id = 0; id < ndim; id++)
				dist += (bpoints->get(i,id) - bpoints->get(j,id))*(bpoints->get(i,id) - bpoints->get(j,id));
			dist = sqrt(dist);
			temp = (this->*rbfunc)(dist);
			if(fabs(temp) > tol*tol)
			{
				A->set(i,j, temp);
				A->set(j,i, temp);
			}
		}
	}

	/*std::cout << "RBFmove:  assembleLHS(): assembling P_b" << std::endl;
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

void RBFmove::move_step()
{
	// solve for RBF coefficients
	std::cout << "RBFmove:  move_step(): Solving linear system" << std::endl;
	amat::Matrix<double> xold(nbpoin,1);
	xold.zeros();
		
	if(lsolver == "CG")
		for(int idim = 0; idim < ndim; idim++)
			coeffs[idim] = sparseCG_d(&A, b[idim], xold, tol, maxiter);
	else if(lsolver == "SOR")
		for(int idim = 0; idim < ndim; idim++)
			coeffs[idim] = sparseSOR(&A, b[idim], xold, tol, maxiter);
	else if(lsolver == "BICGSTAB")
		for(int idim = 0; idim < ndim; idim++)
			coeffs[idim] = sparse_bicgstab(&A, b[idim], xold, tol, maxiter);
	else
	{
		amat::Matrix<double> coeffsm(nbpoin,ndim);
		amat::Matrix<double> rhs(nbpoin,ndim);
		amat::Matrix<double> B = A.toDense();
		
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				rhs(i,j) = b[j](i);
		
		coeffsm = gausselim(B, rhs);
		
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				coeffs[j](i) = coeffsm.get(i,j);
	}
	
	
		//coeffs[idim] = sparsePCG(&A, b[idim], xold, "jacobi", tol, maxiter);
		//coeffs[idim] = sparsegaussseidel(&A, b[idim], xold, tol, maxiter);
		//coeffs[idim] = gausselim(B, b[idim]);

	std::cout << "RBFmove:  move_step(): Moving interior points" << std::endl;
	// calculate new positions of interior points
	int i;
	amat::Matrix<double>* co = coeffs;			// first assign local pointers to class variables for OpenMP
	amat::Matrix<double>* bp = &bpoints;
	amat::Matrix<double>* ip = &inpoints;
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
	}
	printf("\n");
}

void RBFmove::move()
{
	int istep;

	for(istep = 0; istep < nsteps; istep++)
	{
		std::cout << "RBFmove: move(): Step " << istep << std::endl;
		assembleLHS();

		move_step();

		// move buondary points
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				bpoints(i,j) += b[j](i);
	}
}

amat::Matrix<double> RBFmove::getInteriorPoints()
{
	return inpoints;
}

amat::Matrix<double> RBFmove::getBoundaryPoints()
{
	return bpoints;
}

} // end namespace
