/* @brief A mesh-movement method using radial basis function interpolation 
 * 
 * based on the 2007 paper by de Boer, van der Schoot and Bijl,
 * but with one major difference - the interpolation is done using only RBFs; the polynomial part is ignored.
 * @author Aditya Kashi
 * @date August 6, 2015
 * 
 * Aug 14, 2015: Modified to work with curved mesh generation. Now requires interior points and boundary points as separate inputs.
 */

#include "arbf.hpp"

namespace amc {
	
RBFmove::RBFmove() {isalloc = false; }

RBFmove::RBFmove(amat::Matrix<double>* int_points, amat::Matrix<double>* boun_points, amat::Matrix<double>* boundary_motion, const int rbf_ch, const double support_radius, 
		const int num_steps, const double tolerance, const int iter, const std::string linear_solver)
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
	
	std::cout << "RBFmove: RBF to use: " << rbf_ch << std::endl;
	std::cout << "RBFmove: Support radius = " << srad << std::endl;
	std::cout << "RBFmove: Number of steps = " << nsteps << std::endl;
}

void RBFmove::setup(amat::Matrix<double>* int_points, amat::Matrix<double>* boun_points, amat::Matrix<double>* boundary_motion, const int rbf_ch, const double support_radius, 
		const int num_steps, const double tolerance, const int iter, const std::string linear_solver)
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
				 std::cout << "RBFmove: setup(): Selected C2 function as RBF" << std::endl;
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
/// Wendland's C2 function
/*double RBFmove::rbf_c2_compact(double xi)
{
	if(xi < srad)
		return pow(1-xi/srad,4)*(4*xi/srad+1);
	else return 0;
}*/
double RBFmove::rbf_c2_compact(double xi)
{
	if(xi < srad)
	{
		double q = xi/srad;
		return (1.0-q)*(1.0-q)*(1.0-q)*(1.0-q)*(4.0*q+1.0);
	}
	else return 0.0;
}

double RBFmove::rbf_c0(double xi)
{
	if(xi < srad)
		return (1-xi/srad)*(1-xi/srad);
	else
		return 0;
}

double RBFmove::rbf_c4(double xi)
{
	if(xi < srad)
		return pow(1-xi/srad,6)*(35*xi*xi/(srad*srad) + 18*xi/srad + 3);
	else
		return 0;
}

double RBFmove::gaussian(double xi)
{
	return exp(-xi*xi);
}

/** Note that an element is only inserted into the sparse LHS matrix if its magnitude is more than tol * tol
 */
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
	
	//#pragma omp parallel for default(none) private(i,j,dist,temp) shared(A,bpoints,rbfunc,nbpoin,ndim)
	for(i = 0; i < nbpoin; i++)
	{
		A->set(i,i, (this->*rbfunc)(0.0));			// set diagonal element in row i
		//A->set(i,i, 1.0);			// set diagonal element in row i

		for(j = i+1; j < nbpoin; j++)		// traverse lower triangular matrix
		{
			dist = 0;
			for(int id = 0; id < ndim; id++)
				dist += (bpoints->get(i,id) - bpoints->get(j,id))*(bpoints->get(i,id) - bpoints->get(j,id));
			dist = sqrt(dist);
			temp = (this->*rbfunc)(dist);
			//std::cout << temp << " ";
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
	
	/*std::ofstream fout("rbfmatrix.dat");
	A->fprint(fout);
	fout.close();*/
}

void RBFmove::move_step()
{
	// solve for RBF coefficients
	std::cout << "RBFmove:  move_step(): Solving linear system" << std::endl;
	amat::Matrix<double> xold(nbpoin,1);
	xold.zeros();
	
	if(lsolver == "CG")
		for(int idim = 0; idim < ndim; idim++)
		{
			coeffs[idim] = sparseCG(&A, b[idim], xold, tol, maxiter);
		}	
	else if(lsolver == "PCG")
		for(int idim = 0; idim < ndim; idim++)
		{
			coeffs[idim] = sparseCG_d(&A, b[idim], xold, tol, maxiter);
		}
	else if(lsolver == "SOR")
		for(int idim = 0; idim < ndim; idim++)
			coeffs[idim] = sparseSOR(&A, b[idim], xold, tol, maxiter);
	else if(lsolver == "BICGSTAB")
		for(int idim = 0; idim < ndim; idim++)
			coeffs[idim] = sparse_bicgstab(&A, b[idim], xold, tol, maxiter);
	else if(lsolver == "DLU")
	{
		amat::Matrix<double> coeffsm(nbpoin,ndim);
		amat::Matrix<double> rhs(nbpoin,ndim);
		amat::Matrix<double> B = A.toDense();
		
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				rhs(i,j) = b[j](i);
		
		gausselim(B, rhs, coeffsm);
		
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				coeffs[j](i) = coeffsm.get(i,j);
	}
#ifdef EIGEN_LIBRARY
	else if(lsolver == "EIGENLU")
	{
		amat::Matrix<double> coeffsm(nbpoin,ndim);
		amat::Matrix<double> rhs(nbpoin,ndim);
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				rhs(i,j) = b[j](i);

		// solve the system by Eigen's SparseLU routine
		coeffsm = gausselim(A,rhs);

		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				coeffs[j](i) = coeffsm.get(i,j);
	}
#ifdef PASTIX_LIBRARY
	else if(lsolver == "PASTIXLDLT")
	{
		amat::Matrix<double> coeffsm(nbpoin,ndim);
		amat::Matrix<double> rhs(nbpoin,ndim);
		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				rhs(i,j) = b[j](i);

		// solve the system by Eigen's PastiX LDL^T solver interface
		pastix_LDLT(A,rhs,coeffsm);

		for(int i = 0; i < nbpoin; i++)
			for(int j = 0; j < ndim; j++)
				coeffs[j](i) = coeffsm.get(i,j);
	}
#endif
#endif
	
	
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
	amc_real sr = RBFmove::srad;

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

			if(dist < srad)	
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
