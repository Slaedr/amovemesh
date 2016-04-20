#include "alinalg.hpp"

#ifdef EIGEN_LIBRARY

#ifndef EIGEN_SPARSE_MODULE_H
#include <Eigen/Sparse>
#endif

#ifndef EIGEN_SVD_MODULE_H
#include <Eigen/SVD>
#endif

#ifdef PASTIX_LIBRARY
	#ifndef EIGEN_PASTIXSUPPORT_MODULE_H
	#include <Eigen/PaStiXSupport>
	#endif
#endif

#endif

namespace amat {

/* Note: Cholesky algorithm only implemented for a row-major matrix */
Matrix<double> cholesky(Matrix<double> A, Matrix<double> b)
{
	Matrix<double> B;

	std::cout << "\ncholesky: Input LHS matrix is " << A.rows() << " x " << A.cols() << std::endl;
	if(A.rows() != b.rows()) { std::cout << "\nInvalid dimensions of A and b!"; return B; }
	int N = A.rows();

	//Part 1: Cholesky decomposition
	B.setup(N,N); B.zeros();

	B(0,0) = sqrt(A(0,0));
	for(int i = 1; i < N; i++)
		B(i,0) = A(i,0)/B(0,0);

	for(int j = 1; j < N; j++)
	{
		double bjk_sum = 0;
		int k = 0;
		do
		{
			bjk_sum += B(j,k)*B(j,k);
			k++;
		}
		while(k <= j-1);

		if(bjk_sum >= A(j,j)) std::cout << "\n! cholesky: Negative argument to sqrt at ("<<j<<","<<j<<")\n";
		B(j,j) = sqrt(A(j,j) - bjk_sum);

		for(int i = j+1; i < N; i++)
		{
			double bsum = 0;
			k=0;
			do
			{	bsum += B(i,k)*B(j,k);
				k++;
			}
			while(k <= j-1);
			B(i,j) = (A(i,j) - bsum)/B(j,j);
		}
	}
	// We now have B, the lower triangular matrix

	// Check if any of the diagonal elements of B are zero
	for(int i = 0; i < N; i++)
		if(abs(B(i,i)) < 1e-10)
		{
			std::cout << "\ncholesky: Element (" << i <<"," << i << ") of lower triangular matrix is near zero!";
			return B;
		}

	// Part 2: forward substitution to obtain intermediate vector y
	Matrix<double> y(N,1);

	y(0,0) = b(0,0)/B(0,0);

	for(int i = 1; i < N; i++)
	{
		double sum = 0;
		int k = 0;
		do
		{	sum += B(i,k)*y(k,0);
			k++;
		} while(k <= i-1);
		y(i,0) = (b(i,0) - sum)/B(i,i);
	}

	//Part 3: back substitution to obtain final solution
	// Note: the final solution is stored in b
	b.zeros();
	//Matrix<double> f(N,1);
	b(N-1,0) = y(N-1,0)/B(N-1,N-1);

	for(int i = N-2; i >= 0; i--)
	{
		double sum = 0;
		int k = i+1;
		do
		{	sum += B(k,i)*b(k,0);
			k++;
		} while(k <= N-1);
		b(i,0) = (y(i,0) - sum)/B(i,i);
	}

	return b;
}

/// Solves Ax = b by Cholesky decomposition
/** The output is finally stored in b.
 * TODO: optimize this function to factorize A in-place
 */
void chol(Matrix<amc_real>& A, Matrix<amc_real>& b)
{
	Matrix<amc_real> B;

	std::cout << "\ncholesky: Input LHS matrix is " << A.rows() << " x " << A.cols() << std::endl;
	if(A.rows() != b.rows()) { std::cout << "\nInvalid dimensions of A and b!"; return; }
	int N = A.rows(), i, j, k;
	amc_real bjk_sum, bsum, sum;

	//Part 1: Cholesky decomposition
	B.setup(N,N); B.zeros();

	B(0,0) = sqrt(A(0,0));
	for(i = 1; i < N; i++)
		B(i,0) = A(i,0)/B(0,0);

	for(j = 1; j < N; j++)
	{
		bjk_sum = 0;
		k = 0;
		do
		{
			bjk_sum += B(j,k)*B(j,k);
			k++;
		}
		while(k <= j-1);

		if(bjk_sum >= A(j,j)) std::cout << "\n! cholesky: Negative argument to sqrt at ("<<j<<","<<j<<")\n";
		B(j,j) = sqrt(A(j,j) - bjk_sum);

		for(i = j+1; i < N; i++)
		{
			bsum = 0;
			k=0;
			do
			{	bsum += B(i,k)*B(j,k);
				k++;
			}
			while(k <= j-1);
			B(i,j) = (A(i,j) - bsum)/B(j,j);
		}
	}
	// We now have B, the lower triangular matrix

#if DEBUG==1
	// Check if any of the diagonal elements of B are zero
	for(i = 0; i < N; i++)
		if(abs(B(i,i)) < ZERO_TOL)
		{
			std::cout << "\ncholesky: Element (" << i <<"," << i << ") of lower triangular matrix is near zero!";
		}
#endif

	// Part 2: forward substitution to obtain intermediate vector y
	Matrix<amc_real> y(N,1);

	y(0,0) = b(0,0)/B(0,0);

	for(i = 1; i < N; i++)
	{
		sum = 0;
		k = 0;
		do
		{	sum += B(i,k)*y(k,0);
			k++;
		} while(k <= i-1);
		y(i,0) = (b(i,0) - sum)/B(i,i);
	}

	//Part 3: back substitution to obtain final solution
	// Note: the final solution is stored in b
	b.zeros();
	b(N-1,0) = y(N-1,0)/B(N-1,N-1);

	for(i = N-2; i >= 0; i--)
	{
		sum = 0;
		int k = i+1;
		do
		{	sum += B(k,i)*b(k,0);
			k++;
		} while(k <= N-1);
		b(i,0) = (y(i,0) - sum)/B(i,i);
	}
}

void gausselim(Matrix<double>& A, Matrix<double>& b, Matrix<double>& x)
{
	//std::cout << "gausselim: Input LHS matrix is " << A.rows() << " x " << A.cols() << std::endl;
	if(A.rows() != b.rows()) { std::cout << "gausselim: Invalid dimensions of A and b!\n"; return; }
	int N = A.rows();

	x.zeros();
	int l;
	int k;
	double ff;

	for(int i = 0; i < N-1; i++)
	{
		double max = dabs(A(i,i));
		int maxr = i;
		for(int j = i+1; j < N; j++)
		{
			if(dabs(A(j,i)) > max)
			{
				max = dabs(A(j,i));
				maxr = j;
			}
		}
		if(max > ZERO_TOL)
		{
			//interchange rows i and maxr 
			for(k = i; k < N; k++)
			{
				double temp = A(i,k);
				A(i,k) = A(maxr,k);
				A(maxr,k) = temp;
			}
			// do the interchange for b as well
			for(k = 0; k < b.cols(); k++)
			{
				double temp = b(i,k);
				b(i,k) = b(maxr,k);
				b(maxr,k) = temp;
			}
		}
		else { std::cout << "! gausselim: Pivot not found!!\n"; return; }

		for(int j = i+1; j < N; j++)
		{
			ff = A(j,i);
			for(l = i; l < N; l++)
				A(j,l) = A(j,l) - ff/A(i,i)*A(i,l);
			for(k = 0; k < b.cols(); k++)
				b(j,k) = b(j,k) - ff/A(i,i)*b(i,k);
		}
	}
	//Thus, A has been transformed to an upper triangular matrix, b has been transformed accordingly.

	//Part 2: back substitution to obtain final solution
	// Note: the solution is stored in x
	double sum;
	for(l = 0; l < b.cols(); l++)
	{
		x(N-1,l) = b(N-1,l)/A(N-1,N-1);

		for(int i = N-2; i >= 0; i--)
		{
			sum = 0;
			k = i+1;
			do
			{	
				sum += A(i,k)*x(k,l);
				k++;
			} while(k <= N-1);
			x(i,l) = (b(i,l) - sum)/A(i,i);
		}
	}
}

#ifdef EIGEN_LIBRARY
Matrix<double> gausselim(const SpMatrix& A, const Matrix<amc_real>& b)
{
	int nr = A.rows(), nc = A.cols(), nrhs = b.cols();
	Matrix<amc_real> x(nr,nrhs);
	if(b.rows() != nr)
	{
		std::cout << "gausselim: ! Dimension mismatch between LHS and RHS!" << std::endl;
		return x;
	}
	
	std::cout << "gausselim: Converting sparse matrix and RHS to Eigen formats..." << std::endl;
	int i,j,k;

	Eigen::Matrix<amc_real, Eigen::Dynamic, Eigen::Dynamic> B(nr,nrhs);
	for(i = 0; i < nr; i++)
		for(j = 0; j < nrhs; j++)
			B(i,j) = b.get(i,j);

	SMatrixCRS<amc_real> lhs;
	A.get_CRS_matrix(lhs);
	
	Eigen::MappedSparseMatrix<amc_real, Eigen::RowMajor> AA(A.rows(), A.cols(), lhs.nnz, lhs.row_ptr, lhs.col_ind, lhs.val);

	std::cout << "gausselim: Solving via Eigen SparseLU..." << std::endl;
	
	Eigen::SparseLU < Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::COLAMDOrdering<int> > eigsolver;
	eigsolver.analyzePattern(AA);
	eigsolver.factorize(AA);

	Eigen::Matrix<amc_real, Eigen::Dynamic, Eigen::Dynamic> xx = eigsolver.solve(B);
	for(i = 0; i < nr; i++)
		for(j = 0; j < nrhs; j++)
			x(i,j) = xx(i,j);

	return x;
}

#ifdef PASTIX_LIBRARY
void pastix_LDLT(const SpMatrix& A, const Matrix<amc_real>& b, Matrix<amc_real>& x)
{
	int nr = A.rows(), nc = A.cols(), nrhs = b.cols();
	if(b.rows() != nr)
	{
		std::cout << "pastix_LDLT: ! Dimension mismatch between LHS and RHS!" << std::endl;
		return;
	}
	
	std::cout << "pastix_LDLT: Converting sparse matrix and RHS to Eigen formats..." << std::endl;
	int i,j,k;

	Eigen::Matrix<amc_real, Eigen::Dynamic, Eigen::Dynamic> B(nr,nrhs);
	for(i = 0; i < nr; i++)
		for(j = 0; j < nrhs; j++)
			B(i,j) = b.get(i,j);

	SMatrixCRS<amc_real> lhs;
	A.get_CRS_matrix(lhs);
	
	Eigen::MappedSparseMatrix<amc_real, Eigen::RowMajor> AA(A.rows(), A.cols(), lhs.nnz, lhs.row_ptr, lhs.col_ind, lhs.val);

	std::cout << "pastix_LDLT: Solving via Eigen PastixLDLT..." << std::endl;
	
	Eigen::PastixLDLT < Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Upper > eigsolver;
	eigsolver.analyzePattern(AA);
	eigsolver.factorize(AA);

	Eigen::Matrix<amc_real, Eigen::Dynamic, Eigen::Dynamic> xx = eigsolver.solve(B);
	for(i = 0; i < nr; i++)
		for(j = 0; j < nrhs; j++)
			x(i,j) = xx(i,j);
}
#endif

#endif

/* Re-stores matrix data in the form needed by SuperLU, and calls the SuperLU routine to solve.
 */
/*void superLU_solve(const SpMatrix* aa, const Matrix<double>* b, Matrix<double>* ans)
{
	if(aa->rows() != b->rows())
	{
		std::cout << "superLU_solve(): Input sizes do not match!!" << std::endl;
		return;
	}
	
	// get LHS matrix in CRS form and store in lhs
	SMatrixCRS<double> lhs;
	aa->getCRSMatrix(lhs);

	SuperMatrix A, B, L, U;
	double *a, *rhs;
	int *col_ind, *row_ptr, *perm_r, *perm_c; 
	a = doubleMalloc(lhs.nnz);
	col_ind = intMalloc(lhs.nnz);
	row_ptr = intMalloc(aa->rows()+1);
	rhs = doubleMalloc(b->rows()*b->cols());
	perm_r = intMalloc(aa->rows());
	perm_c = intMalloc(aa->cols());

	int i,j,k,info,nprocs;
	nprocs = std::atoi(std::getenv("OMP_NUM_THREADS"));

	for(i = 0; i < lhs.nnz; i++)
	{
		a[i] = lhs.val[i];
		col_ind[i] = lhs.col_ind[i];
	}
	for(i = 0; i < aa->rows()+1; i++)
		row_ptr[i] = lhs.row_ptr[i];

	// store the RHS into rhs in column-major format
	k = 0;
	for(j = 0; j < b->cols(); j++)
		for(i = 0; i < b->rows(); i++)
		{
			rhs[k] = b->get(i,j);
			k++;
		}

	dCreate_CompCol_Matrix(&A, aa->rows(), aa->cols(), lhs.nnz, a, col_ind, row_ptr, SLU_NR, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&B, b->rows(), b->cols(), rhs, b->rows(), SLU_DN, SLU_D, SLU_GE);
	
	pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);

	// copy the solution into ans
	double* ansptr = (double*)((DNformat*)B.Store)->nzval;
	k=0;
	for(j = 0; j < b->cols(); j++)
		for(i = 0; i < b->rows(); i++)
		{
			(*ans)(i,j) = ansptr[i + j*b->rows()];
			k++;
		}

	SUPERLU_FREE(rhs);
	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_SCP(&L);
	Destroy_CompCol_NCP(&U);
}*/

//-------------------- Iterative Methods ----------------------------------------//

Matrix<double> pointjacobi(Matrix<double> A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check)
{
	std::cout << "\npointjacobi: Input LHS matrix is " << A.rows() << " x " << A.cols() << std::endl;
	if(A.rows() != b.rows()) { std::cout << "! pointjacobi: Invalid dimensions of A and b!\n"; return b; }
	if(xold.rows() != b.rows()) { std::cout << "! pointjacobi: Invalid dimensions on xold !\n"; return b; }
	int N = A.rows();

	if(check == 'y')
	{
		for(int i = 0; i < A.rows(); i++)
		{
			double sum = 0;
			for(int j = 0; j < A.cols(); j++)
			{
				if(i != j) sum += dabs(A(i,j));
			}
			if(dabs(A(i,i)) <= sum) std::cout << "* pointjacobi(): LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
		}
	}

	Matrix<double> x(b.rows(),1);

	Matrix<double> M(N,N);		// diagonal matrix
	M.zeros();

	//populate diagonal matrix M
	for(int i = 0; i < N; i++)
		M(i,i) = A(i,i);

	A = A - M;

	//invert M
	for(int i = 0; i < N; i++)
		M(i,i) = 1.0/M(i,i);

	//x.zeros();
	//std::cout << "pointjacobi: Initial error = " << (xold-x).dabsmax() << std::endl;

	x = xold;
	int c = 0;

	do
	{
		xold = x;
		x = M * (b - (A*xold));
		c++;
		if(c > maxiter) { std::cout << "pointjacobi: Max iterations exceeded!\n"; break; }
	} while((x-xold).dabsmax() >= tol);
	
	std::cout << "pointjacobi: No. of iterations = " << c << std::endl;

	return x;
}

Matrix<double> gaussseidel(Matrix<double> A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check)
{
	//std::cout << "\ngaussseidel: Input LHS matrix is " << A.rows() << " x " << A.cols() << std::endl;
	if(A.rows() != b.rows()) { std::cout << "! gaussseidel: Invalid dimensions of A and b!\n"; return b; }
	if(xold.rows() != b.rows()) { std::cout << "! gaussseidel: Invalid dimensions on xold !\n"; return b; }
	int N = A.rows();

	if(check == 'y')
	{
		for(int i = 0; i < A.rows(); i++)
		{
			double sum = 0;
			for(int j = 0; j < A.cols(); j++)
			{
				if(i != j) sum += dabs(A(i,j));
			}
			if(dabs(A(i,i)) <= sum) std::cout << "* gaussseidel: LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
		}
	}

	Matrix<double> x(b.rows(),1);

	Matrix<double> M(N,N);		// diagonal matrix
	M.zeros();

	//populate diagonal matrix M
	for(int i = 0; i < N; i++)
		M(i,i) = A(i,i);

	A = A - M;						// diagonal entries of A are now zero

	//invert M
	for(int i = 0; i < N; i++)
		M(i,i) = 1.0/M(i,i);

	//x.zeros();
	//std::cout << "gaussseidel: Initial error = " << (xold-x).dabsmax() << std::endl;

	x = xold;
	int c = 0;
	double initres;
	bool first = true;

	Matrix<double> Axold(N,1);
	Matrix<double> inter(N,1);
	Matrix<double> diff(N,1);	// diff = x - xold
	double error = 1.0;
	do
	{
		xold = x;
		Axold.zeros();
		for(int i = 0; i < N; i++)
		{
			for(int k = 0; k < i; k++)
				Axold(i,0) += A(i,k)*x(k,0);
			for(int k = i; k < N; k++)
				Axold(i,0) += A(i,k)*xold(k,0);
			inter(i,0) = b(i,0) - Axold(i,0);
			x(i,0) = M(i,i) * inter(i,0);
		}
		// NOTE: The above results in FORWARD Gauss-Seidel

		c++;
		if(c > maxiter) { std::cout << "gaussseidel: Max iterations exceeded!\n"; break; }
		for(int i = 0; i < N; i++)
			diff(i,0) = x(i,0) - xold(i,0);

		error = dabs(diff(0,0));
		for(int i = 1; i < N; i++)
			if(dabs(diff(i,0)) > error) error = dabs(diff(i,0));
		//std::cout << "gaussseidel: error = " << error << std::endl;
		if(first == true)
		{	initres = error;
			//std::cout << "gaussseidel: Initial residue = " << initres << std::endl;
			first = false;
		}

		//if(c%20 == 0) std::cout << "gaussseidel(): Step " << c << ", Relative error = " << error/initres << std::endl;

	} while(error/initres >= tol);
	//std::cout << "gaussseidel: No. of iterations = " << c << std::endl;

	return x;
}

Matrix<double> sparsegaussseidel(SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check)
{
	std::cout << "sparsegaussseidel(): Input LHS matrix is " << A->rows() << " x " << A->cols() << std::endl;
	std::cout << "sparsegaussseidel(): b is " << b.rows() << ", and xold is " << xold.rows() << std::endl;
	if(A->rows() != b.rows()) { std::cout << "! gaussseidel: Invalid dimensions of A and b!\n"; return b; }
	if(xold.rows() != b.rows()) { std::cout << "! gaussseidel: Invalid dimensions on xold !\n"; return b; }
	int N = A->rows();

	if(check == 'y')
	{
		for(int i = 0; i < A->rows(); i++)
		{
			double sum = 0;
			for(int j = 0; j < A->cols(); j++)
			{
				if(i != j) sum += dabs(A->get(i,j));
			}
			if(dabs(A->get(i,i)) <= sum) std::cout << "* gaussseidel: LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
		}
	}

	Matrix<double> x(b.rows(),1);

	Matrix<double> M(N,1);		// vector of diagonal elements of A
	M.zeros();

	// diagonal matrix M
	//std::cout << "sparsegaussseidel(): Getting diagonal of sparse matrix\n";
	A->get_diagonal(&M);

	//invert M
	for(int i = 0; i < N; i++)
		M(i) = 1.0/M(i);

	//x.zeros();
	//std::cout << "gaussseidel: Initial error = " << (xold-x).dabsmax() << std::endl;

	for(int i = 0; i < x.rows(); i++)
		x(i) = xold(i);
	int c = 0;
	double initres;
	bool first = true;

	Matrix<double> Axold(N,1);
	Matrix<double> inter(N,1);
	Matrix<double> diff(N,1);	// diff = x - xold
	double error = 1.0;
	//std::cout << "sparsegaussseidel(): Starting iterations\n";
	do
	{
		xold = x;
		Axold.zeros();
		int i;

		//#pragma omp parallel for default(none) private(i) shared(A,b,Axold,x,xold,inter,M,N)
		for(i = 0; i < N; i++)
		{
			//std::cout << "  Calling sparse getelem_multiply_parts\n";
			Axold(i) = A->getelem_multiply_parts(i, &x, &xold, i, xold.get(i));
			//std::cout << "  Setting inter\n";
			inter(i,0) = b(i,0) - Axold(i,0);
			x(i,0) = M(i) * inter(i,0);
		}
		// NOTE: The above results in FORWARD Gauss-Seidel

		//if(c > maxiter) { std::cout << "gaussseidel: Max iterations exceeded!\n"; break; }
		for(int i = 0; i < N; i++)
			diff(i,0) = x(i,0) - xold(i,0);

		error = dabs(diff(0,0));
		for(int i = 1; i < N; i++)
			if(dabs(diff(i,0)) > error) error = dabs(diff(i,0));
		//std::cout << "gaussseidel: error = " << error << std::endl;
		if(first == true)
		{	initres = error;
			if(dabs(initres) < tol*tol)
			{
				std::cout << "sparsegaussseidel(): Initial residue = " << initres << std::endl;
				break;
			}
			first = false;
		}

		if(c%10 == 0 || c == 1) std::cout << "gaussseidel(): Step " << c << ", Relative error = " << error/initres << std::endl;
		c++;

	} while(error/initres >= tol && c <= maxiter);
	std::cout << "gaussseidel: No. of iterations = " << c << ", final error " << error/initres << std::endl;

	return x;
}

Matrix<double> sparseSOR(SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, double w, char check)
// CAUTION: does not work in parallel due to some reason
{
	std::cout << "sparseSOR(): Input LHS matrix is " << A->rows() << " x " << A->cols() << std::endl;
	if(A->rows() != b.rows()) { std::cout << "! sparseSOR(): Invalid dimensions of A and b!\n"; return b; }
	if(xold.rows() != b.rows()) { std::cout << "! sparseSOR(): Invalid dimensions on xold !\n"; return b; }
	int N = A->rows();

	if(check == 'y')
	{
		for(int i = 0; i < A->rows(); i++)
		{
			double sum = 0;
			for(int j = 0; j < A->cols(); j++)
			{
				if(i != j) sum += dabs(A->get(i,j));
			}
			if(dabs(A->get(i,i)) <= sum) std::cout << "* sparseSOR(): LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
		}
	}

	Matrix<double> x(b.rows(),1);

	Matrix<double> M(N,1);		// vector of diagonal elements of A
	M.zeros();

	// diagonal matrix M
	//std::cout << "sparsegaussseidel(): Getting diagonal of sparse matrix\n";
	A->get_diagonal(&M);

	//invert M
	for(int i = 0; i < N; i++)
		M(i) = 1.0/M(i);

	//x.zeros();
	//std::cout << "gaussseidel: Initial error = " << (xold-x).dabsmax() << std::endl;

	x = xold;
	int c = 0;
	double initres;
	bool first = true;

	Matrix<double> Axold(N,1);
	Matrix<double> inter(N,1);
	Matrix<double> diff(N,1);	// diff = x - xold
	double error = 1.0;
	//std::cout << "sparsegaussseidel(): Starting iterations\n";
	do
	{
		xold = x;
		Axold.zeros();
		int i;

		//#pragma omp parallel for default(none) private(i) shared(A,b,Axold,x,xold,inter,M,N,w)
		for(i = 0; i < N; i++)
		{
			//std::cout << "  Calling sparse getelem_multiply_parts\n";
			Axold(i) = A->getelem_multiply_parts(i, &x, &xold, i, xold.get(i));
			//std::cout << "  Setting inter\n";
			inter(i,0) = w*(b(i,0) - Axold(i,0));
			x(i,0) = (1-w)*xold(i,0) + M(i) * inter(i,0);
		}
		// NOTE: The above results in FORWARD Gauss-Seidel

		//if(c > maxiter) { std::cout << "gaussseidel: Max iterations exceeded!\n"; break; }
		for(int i = 0; i < N; i++)
			diff(i,0) = x(i,0) - xold(i,0);

		error = dabs(diff(0,0));
		for(int i = 1; i < N; i++)
			if(dabs(diff(i,0)) > error) error = dabs(diff(i,0));
		//std::cout << "gaussseidel: error = " << error << std::endl;
		if(first == true)
		{	initres = error;
			if(dabs(initres) < tol*tol)
			{
				std::cout << "sparseSOR(): Initial residue = " << initres << std::endl;
				break;
			}
			first = false;
		}

		if(c%10 == 0 || c == 1) std::cout << "sparseSOR(): Step " << c << ", Relative error = " << error/initres << std::endl;
		c++;

	} while(error/initres >= tol && c <= maxiter);
	std::cout << "sparseSOR(): No. of iterations = " << c << ", final error " << error/initres << std::endl;

	return x;
}

/* Calculates solution of Ax=b where A is a SPD matrix in sparse format.
 * No preconditioning is used.
 * NOTE: The parallel version is actually slower, due to some reason.
 */
Matrix<double> sparseCG(const SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter)
{
	std::cout << "sparseCG(): Solving " << A->rows() << "x" << A->cols() << " system by conjugate gradient method with no preconditioner\n";
	if(A->rows() != b.rows() || A->rows() != xold.rows()) std::cout << "sparseCG(): ! Mismatch in number of rows!!" << std::endl;

	Matrix<double> x(A->rows(),1);		// solution vector
	Matrix<double> rold(A->rows(),1);		// initial residual = b - A*xold
	Matrix<double> r(A->rows(),1);			// residual = b - A*x
	Matrix<double> p(A->rows(),1);
	Matrix<double> pold(A->rows(),1);
	Matrix<double> temp(A->rows(),1);
	Matrix<double> diff(A->rows(),1);
	double temp1, temp2;
	double theta;
	double beta;
	double error = 1.0;
	double initres;
	double normalizer = b.l2norm();

	A->multiply(xold, &temp);		// temp := A*xold
	rold = b - temp;
	error  = rold.l2norm();		// initial residue
	if(error < tol)
	{
		std::cout << "sparseCG(): Initial residual is very small. Nothing to do." << std::endl;
		return xold;
	}

	pold = rold;

	int steps = 0;

	std::cout << "SparseCG: Starting loop" << std::endl;
	do
	{
		if(steps % 10 == 0 || steps == 1)
			std::cout << "sparseCG(): Iteration " << steps << ", relative residual = " << error/normalizer << std::endl;
		int i;

		temp1 = rold.dot_product(rold);

		A->multiply(pold, &temp);

		// p^T A p
		temp2 = pold.dot_product(temp);
		
		if(temp2 <= ZERO_TOL) { 
			std::cout << "sparseCG: Matrix A may not be positive-definite!! temp2 is " << temp2 << "\n";
			std::cout << "sparseCG: Magnitude of pold in this iteration is " << pold.l2norm() << std::endl;
		}
		theta = temp1/temp2;

		//#pragma omp parallel for default(none) private(i) shared(x,r,xold,rold,pold,temp,theta)
		for(i = 0; i < x.rows(); i++)
		{
			x(i) = xold.get(i) + pold.get(i)*theta;
			rold(i) = rold.get(i) - temp.get(i)*theta;
		}

		beta = rold.dot_product(rold) / temp1;
		
		for(i = 0; i < x.rows(); i++)
			pold(i) = rold.get(i) + pold.get(i)*beta;

		//calculate ||b - A*x||. 'error' is a misnomer - it's the residual norm.
		error = rold.l2norm();

		// set old variables
		xold = x;

		if(steps > maxiter)
		{
			std::cout << "! sparseCG(): Max iterations reached!\n";
			break;
		}
		steps++;
	} while(error/normalizer > tol);

	std::cout << "sparseCG(): Done. Number of iterations: " << steps << "; final residual " << error/normalizer << ".\n";
	return x;
}

/* Calculates solution of Ax=b where A is a SPD matrix in sparse format.
 * The preconditioner is a diagonal matrix.
 * NOTE: The parallel version is actually slower, due to some reason.
 */
Matrix<double> sparseCG_d(const SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter)
{
	std::cout << "sparseCG_d(): Solving " << A->rows() << "x" << A->cols() << " system by conjugate gradient method with diagonal preconditioner\n";
	if(A->rows() != b.rows() || A->rows() != xold.rows()) std::cout << "sparseCG_d(): ! Mismatch in number of rows!!" << std::endl;

	Matrix<double> x(A->rows(),1);		// solution vector
	Matrix<double> M(A->rows(), 1);		// diagonal preconditioner, or soon, inverse of preconditioner
	Matrix<double> rold(A->rows(),1);		// initial residual = b - A*xold
	Matrix<double> r(A->rows(),1);			// residual = b - A*x
	Matrix<double> z(A->rows(),1);
	Matrix<double> zold(A->rows(),1);
	Matrix<double> p(A->rows(),1);
	Matrix<double> pold(A->rows(),1);
	Matrix<double> temp(A->rows(),1);
	Matrix<double> diff(A->rows(),1);
	double temp1, temp2;
	double theta;
	double beta;
	double error = 1.0;
	double initres;
	double normalizer = b.l2norm();

	M.zeros();
	A->get_diagonal(&M);
	for(int i = 0; i < A->rows(); i++)
	{
		M(i) = 1.0/M(i);
	}

	A->multiply(xold, &temp);		// temp := A*xold
	rold = b - temp;
	error  = rold.l2norm();		// initial residue
	if(error < tol)
	{
		std::cout << "sparseCG_d(): Initial residual is very small. Nothing to do." << std::endl;
		return xold;
	}

	for(int i = 0; i < A->rows(); i++)
		//zold(i) = M(i)*rold(i);				// zold = M*rold
		zold(i) = rold(i);

	pold = zold;

	int steps = 0;

	std::cout << "SparseCG_d: Initial residual norm = " << normalizer << std::endl;
	do
	{
		if(steps % 10 == 0 || steps == 1)
			std::cout << "sparseCG_d(): Iteration " << steps << ", relative residual = " << error/normalizer << std::endl;
		int i;

		temp1 = rold.dot_product(zold);

		A->multiply(pold, &temp);
		//temp.mprint();

		temp2 = pold.dot_product(temp);
		if(temp2 <= ZERO_TOL) { 
			std::cout << "sparseCG_d: Matrix A may not be positive-definite!! temp2 is " << temp2 << "\n";
			std::cout << "sparseCG_D: Magnitude of pold in this iteration is " << pold.l2norm() << std::endl;
		}
		theta = temp1/temp2;

		//#pragma omp parallel for default(none) private(i) shared(x,r,xold,rold,pold,temp,theta)
		for(i = 0; i < x.rows(); i++)
		{
			x(i) = xold.get(i) + pold.get(i)*theta;
			rold(i) = rold.get(i) - temp.get(i)*theta;
		}

		if(steps > 2)
		{
			//#pragma omp parallel for default(none) private(i) shared(zold,M,rold,A)
			for(i = 0; i < A->rows(); i++)
				zold(i) = M(i)*rold(i);
		}
		else
		{
			//#pragma omp parallel for default(none) private(i) shared(zold,M,rold,A)
			for(i = 0; i < A->rows(); i++)
				zold(i) = rold(i);
		}

		beta = rold.dot_product(zold) / temp1;

		//#pragma omp parallel for default(none) private(i) shared(zold,x,p,pold,beta)
		for(i = 0; i < x.rows(); i++)
			pold(i) = zold.get(i) + pold.get(i)*beta;

		//calculate ||b - A*x||. 'error' is a misnomer - it's the residual norm.
		error = rold.l2norm();

		// set old variables
		xold = x;
		/*rold = r;
		zold = z;
		pold = p;*/

		if(steps > maxiter)
		{
			std::cout << "! sparseCG_d(): Max iterations reached!\n";
			break;
		}
		steps++;
	} while(error/normalizer > tol);

	std::cout << "sparseCG_d(): Done. Number of iterations: " << steps << "; final residual " << error/normalizer << ".\n";
	return x;
}

// Does not currently have a prototype in the header file
void precon_jacobi(SpMatrix* A, const Matrix<double>& r, Matrix<double>& z)
// Multiplies r by the Jacobi preconditioner matrix of A, and stores the result in z
{
	Matrix<double> diag(A->rows(), 1);
	A->get_diagonal(&diag);
	
	for(int i = 0; i < A->rows(); i++)
		diag(i) = 1.0/diag(i);
	
	for(int i = 0; i < A->rows(); i++)
		z(i) = diag(i)*r.get(i);
}

// Does not currently have a prototype in the header file
void precon_lusgs(SpMatrix* A, const Matrix<double>& r, Matrix<double>& z)
// Multiplies r by the LU-SGS preconditioner matrix of A, and stores the result in z
{
	SpMatrix L;
	SpMatrix U;
	Matrix<double> D(A->rows(), 1);
	Matrix<double> Dinv(A->rows(), 1);
	Matrix<double> z_initial(z.rows(), 1);
	for(int i = 0; i < z.rows(); i++)
		z_initial(i) = z.get(i);
	
	A->get_diagonal(&D);
	for(int i = 0; i < A->rows(); i++)
		Dinv(i) = 1.0/D(i);
	
	A->get_lower_triangle(L);
	A->get_upper_triangle(U);
	
	double temp = 0;
	Matrix<double> zold(A->rows(),1);
	
	// solve (D+L)*zold = r by forward substitution
	
}

// Does not currently have a prototype in the header file
/* Calculates solution of Ax=b where A is a SPD matrix in sparse format. The preconditioner is supplied by a function pointer.*/
Matrix<double> sparsePCG(SpMatrix* A, Matrix<double> b, Matrix<double> xold, std::string precon, double tol, int maxiter)
{
	std::cout << "sparsePCG(): Solving " << A->rows() << "x" << A->cols() << " system by conjugate gradient method with diagonal preconditioner\n";
	
	void (*precond)(SpMatrix* lhs, const Matrix<double>& r, Matrix<double>& z);
	//const Matrix<double>& z_initial,
	
	if(precon == "jacobi") precond = &precon_jacobi;
	else if(precon == "lusgs") precond = &precon_lusgs;
	
	// check
	//if(A->rows() != b.rows() || A->rows() != xold.rows()) std::cout << "sparseCG_d(): ! Mismatch in number of rows!!" << std::endl;

	Matrix<double> x(A->rows(),1);		// solution vector
	Matrix<double> M(A->rows(), 1);		// diagonal preconditioner, or soon, inverse of preconditioner
	Matrix<double> rold(A->rows(),1);		// initial residual = b - A*xold
	Matrix<double> r(A->rows(),1);			// residual = b - A*x
	Matrix<double> z(A->rows(),1);
	Matrix<double> zold(A->rows(),1);
	Matrix<double> p(A->rows(),1);
	Matrix<double> pold(A->rows(),1);
	Matrix<double> temp(A->rows(),1);
	Matrix<double> diff(A->rows(),1);
	double temp1, temp2;
	double theta;
	double beta;
	double error = 1.0;
	double initres;

	//std::cout << "sparseCG_d(): Declared everything" << std::endl;

	M.ones();		// disable preconditioner
	//std::cout << "sparseCG_d(): preconditioner enabled" << std::endl;

	A->multiply(xold, &temp);		// temp := A*xold
	rold = b - temp;
	error = initres = rold.l2norm();		// initial residue
	if(error < tol)
	{
		std::cout << "sparsePCG(): Initial residual is very small. Nothing to do." << std::endl;
		//x.zeros();
		return xold;
	}

	for(int i = 0; i < A->rows(); i++)
		//zold(i) = M(i)*rold(i);				// zold = M*rold
		zold(i) = rold(i);

	pold = zold;

	int steps = 0;

	do
	{
		if(steps % 10 == 0 || steps == 1)
			std::cout << "sparsePCG(): Iteration " << steps << ", relative residual = " << error << std::endl;
		int i;

		temp1 = rold.dot_product(zold);

		A->multiply(pold, &temp);
		//temp.mprint();

		temp2 = pold.dot_product(temp);
		if(temp2 <= 0) std::cout << "sparsePCG: ! Matrix A is not positive-definite!! temp2 is " << temp2 << "\n";
		theta = temp1/temp2;

		//#pragma omp parallel for default(none) private(i) shared(x,r,xold,rold,pold,temp,theta)
		for(i = 0; i < x.rows(); i++)
		{
			//std::cout << "Number of threads " << omp_get_num_threads();
			x(i) = xold.get(i) + pold.get(i)*theta;
			rold(i) = rold.get(i) - temp.get(i)*theta;
			//diff(i) = x(i) - xold(i);
		}
		//std::cout << "x:\n"; x.mprint();
		//std::cout << "r:\n"; r.mprint();

		if(steps > 5)
		{
			// calculate zold as (M^-1)*rold
			precond(A, rold, zold);
		}
		else
		{
			//#pragma omp parallel for default(none) private(i) shared(zold,M,rold,A)
			for(i = 0; i < A->rows(); i++)
				zold(i) = rold(i);
		}

		beta = rold.dot_product(zold) / temp1;

		//#pragma omp parallel for default(none) private(i) shared(zold,x,p,pold,beta)
		for(i = 0; i < x.rows(); i++)
			pold(i) = zold.get(i) + pold.get(i)*beta;

		//calculate ||b - A*x||
		error = rold.l2norm();

		// set old variables
		xold = x;
		/*rold = r;
		zold = z;
		pold = p;*/

		if(steps > maxiter)
		{
			std::cout << "! sparsePCG(): Max iterations reached!\n";
			break;
		}
		steps++;
	} while(error > tol);

	std::cout << "sparsePCG(): Done. Number of iterations: " << steps << "; final residual " << error << ".\n";
	return xold;
}


///Solves general linear system Ax=b using stabilized biconjugate gradient method of van der Vorst.
/** Jacobi preconditioning is used.
 * NOTE: The initial guess vector xold is modified by this function.
 */
Matrix<double> sparse_bicgstab(const SpMatrix* A, const Matrix<double>& b, Matrix<double>& xold, double tol, int maxiter)
{
	Matrix<double> x(b.rows(), 1);
	// checks
	if(!(A->rows() == A->cols() && A->rows() == b.rows() && b.rows() == xold.rows())) {
		std::cout << "sparse_bicgstab(): ! Input size error!!" << std::endl;
		return x;
	}
	
	Matrix<double> rold(b.rows(), 1);
	Matrix<double> r(b.rows(), 1);
	Matrix<double> rhat(b.rows(), 1);
	Matrix<double> t(b.rows(), 1), u(b.rows(),1);
	Matrix<double> vold(b.rows(), 1);
	Matrix<double> v(b.rows(), 1);
	Matrix<double> pold(b.rows(), 1);
	Matrix<double> p(b.rows(), 1);
	Matrix<double> s(b.rows(), 1);
	Matrix<double> errdiff(b.rows(), 1);
	Matrix<double> M(A->rows(), 1);		// preconditioner. If K is the preconditioning matrix, this holds either K or K^(-1)
	Matrix<double> y(b.rows(), 1);
	Matrix<double> z(b.rows(), 1);
	double rhoold, rho, wold, w, alpha, beta;
	double resnormrel, resnorm, normalizer;

	M.zeros();
	A->get_diagonal(&M);
	for(int i = 0; i < A->rows(); i++)
	{
		M(i) = 1.0/M(i);
	}
	// M now contains K^(-1), where K is the preconditioning matrix
	A->multiply(xold, &t);		// t = A*xold, initially (t is only a temp variable here)
	rold = b - t;
	resnorm = rold.l2norm();
	
	if(resnorm < A_SMALL_NUMBER) {
		std::cout << "sparse_bicgstab(): Initial residual is very small. Exiting." << std::endl;
		return x;
	}

	normalizer = b.l2norm();
	if(normalizer < ZERO_TOL) {
		std::cout << "sparse_bicgstab(): RHS has zero norm. Exiting." << std::endl;
		return x;
	}
	resnormrel = resnorm/normalizer;
	int steps = 0;
	int i;

	vold.zeros(); pold.zeros();

	rhat = rold;

	rhoold = alpha = wold = 1.0;

	while(resnormrel > tol && steps <= maxiter)
	{
		if(steps % 10 == 0 || steps == 1)
			std::cout << "sparse_bicgstab(): Iteration " << steps << ": rel residual norm = " << resnormrel << std::endl;

		rho = rhat.dot_product(rold);
		beta = rho*alpha/(rhoold*wold);

		for(i = 0; i < b.rows(); i++)
		{
			p(i) = rold.get(i) + beta * (pold.get(i) - wold*vold.get(i));
			y(i) = M.get(i)*p.get(i);
		}
		
		A->multiply(y, &v);			// v = A*y
		alpha = rho / rhat.dot_product(v);

		for(i = 0; i < b.rows(); i++)
		{
			s(i) = rold.get(i) - alpha*v.get(i);
			z(i) = M.get(i)*s.get(i);
		}

		A->multiply(z, &t);			// t = A*K^-1*s

		//--------------- modification begins here (this is for purely left-preconditioning)

		for(i = 0; i < b.rows(); i++)
			u(i) = M.get(i)*t.get(i);		// u = K^(-1)*A*K^(-1)*s = K^(-1)*t

		//w = t.dot_product(s)/t.dot_product(t);
		w = u.dot_product(z) / u.dot_product(u);

		//--------------- modification ends here (done according to Van der Vorst's paper)

		for(i = 0; i < b.rows(); i++)
		{
			x(i) = xold.get(i) + alpha*y.get(i) + w*z.get(i);
		}

		for(i = 0; i < b.rows(); i++)
			rold(i) = s.get(i) - w*t.get(i);

		resnormrel = rold.l2norm()/normalizer;

		for(i = 0; i < b.rows(); i++)
		{
			xold(i) = x.get(i);
			vold(i) = v.get(i);
			pold(i) = p.get(i);
		}
		rhoold = rho;
		wold = w;
		steps++;
	}

	std::cout << "sparse_bicgstab(): Done. Iterations: " << steps << ", final relative residual norm: " << resnormrel << std::endl;

	return x;
}

/// solves the least squares problem (finds the minimum point x) \f$ \min(||Ax - b||_2) \f$ by solving the normal equations
/** The columns of the LHS matrix A are first scaled by the respective column-vector norms, so as to improve the condition number.
 */
void leastSquares_NE(amat::Matrix<amc_real>& A, amat::Matrix<amc_real>& b, amat::Matrix<amc_real>& x)
{
	amc_int m = A.rows(), n = A.cols();
#if DEBUG==1
	if(b.rows() != m || x.rows() != n) 
	{
		std::cout << "leastSquares_NE(): ! Size error at input!" << std::endl;
		return;
	}
#endif
	std::vector<amc_real> scale(n);		// for scaling the column of A
	amc_real csum;
	amc_int i,j,k;
	
	// get norms of column-vectors and scale the columns of A
	for(j = 0; j < n; j++)
	{
		csum = 0;
		for(i = 0; i < m; i++)
			csum += A.get(i,j)*A.get(i,j);
		scale[j] = 1.0/sqrt(csum);
		for(i = 0; i < m; i++)
			A(i,j) *= scale[j];
	}

	amat::Matrix<amc_real> B(n,n);
	amat::Matrix<amc_real> c(n,1);
	for(i = 0; i < n; i++)
	{
		c(i) = 0;
		for(k = 0; k < m; k++)
			c(i) += A.get(k,i)*b.get(k);

		for(j = 0; j < n; j++)
		{
			B(i,j) = 0;
			for(k = 0; k < m; k++)
				B(i,j) += A.get(k,i)*A.get(k,j);
		}
	}
	// we now have B = A^T A and c = A^T b
	// solve by Cholesky
	//chol(B, c);
	gausselim(B, c, x);

	for(i = 0; i < n; i++)
		x(i) *= scale[i];
}

/// Computes the QR decomposition of a dense matrix by Householder algorithm given in Trefethen and Bau, Numerical Linear Algebra.
/** \param A is an m x n matrix, replaced by the upper triangular matrix R at the end of the algorithm.
 * \param v is the set of reflection vectors which can be used to compute Q, the action of Q and the action of Q*.
 * v should contain n vectors: The first vector with size m, second with size m-1,...upto m-n+1.
 */
void qr(amat::Matrix<amc_real>& A, std::vector<amc_real>* v)
{
	int m = A.rows(), n = A.cols(), k, i,j;
	//std::cout << "Problem size " << m << " " << n << std::endl;
	amc_real magx, sgnx0;
		
	amat::Matrix<amc_real> P(m,m);
	std::vector<amc_real> x(m);
	std::vector<amc_real> va(m);

	for(k = 0; k < n; k++)
	{
		// get rows k to m of the kth column of A
		for(i = k; i < m; i++)
			x[i-k] = A.get(i,k);

		if(fabs(x[0]) > ZERO_TOL)
			sgnx0 = x[0]/fabs(x[0]);
		else sgnx0 = 1.0;		// arbitrary 1 or -1
		magx = 0;
		for(i = 0; i < m-k; i++)
			magx += x[i]*x[i];
		magx = sqrt(magx);

		v[k][0] = sgnx0 * magx + x[0];
		magx = v[k][0]*v[k][0];

		for(i = 1; i < m-k; i++)
		{
			v[k][i] = x[i];
			magx += v[k][i]*v[k][i];
		}
		magx = sqrt(magx);
		for(i = 0; i < m-k; i++)
			v[k][i] /= magx;

		// compute v* times A(k:m,k:n) to get a row matrix of size (m-k)x1
		for(j = 0; j < n-k; j++)
		{
			va[j] = 0;
			for(i = 0; i < m-k; i++)
				va[j] += v[k][i]*A.get(k+i,k+j);
		}

		// compute outer product of v[k] and va
		for(i = 0; i < m-k; i++)
			for(j = 0; j < n-k; j++)
				P(i,j) = v[k][i]*va[j];

		// finally update A
		for(i = k; i < m; i++)
			for(j = k; j < n; j++)
				A(i,j) -= 2.0*P(i-k,j-k);

		/*for(i = 0; i < m-k; i++)
			std::cout << x[i] << " ";
		std::cout << std::endl;*/
	}
}

void leastSquares_QR(amat::Matrix<amc_real>& A, amat::Matrix<amc_real>& b, amat::Matrix<amc_real>& x)
{
	amc_int m = A.rows(), n = A.cols();
#ifdef DEBUG
	if(b.rows() != m || x.rows() != n) 
	{
		std::cout << "leastSquares_QR(): ! Size error at input!" << std::endl;
		return;
	}
#endif
	std::vector<amc_real> scale(n);		// for scaling the column of A
	amc_real csum, dp;
	amc_int i,j,k;
	
	// get norms of column-vectors of A and scale them
	for(j = 0; j < n; j++)
	{
		csum = 0;
		for(i = 0; i < m; i++)
			csum += A.get(i,j)*A.get(i,j);
		scale[j] = 1.0/sqrt(csum);
		for(i = 0; i < m; i++)
			A(i,j) *= scale[j];
	}

	std::vector<amc_real>* v = new std::vector<amc_real>[n];
	for(i = 0; i < n; i++)
		v[i].resize(m-i);

	// get QR decomposition (R is stored in A, Q is determined by v)
	qr(A,v);

	solve_QR(v, A, b, x);

	for(i = 0; i < n; i++)
		for(j = 0; j < b.cols(); j++)
			x(i,j) *= scale[i];

	delete[] v;
}

void solve_QR(const std::vector<amc_real>* v, const amat::Matrix<amc_real>& R, amat::Matrix<amc_real>& b, amat::Matrix<amc_real>& x)
{
	amc_int m = R.rows(), n = R.cols();
#if DEBUG==1
	if(b.rows() != m || x.rows() < n) 
	{
		std::cout << "solve_QR(): ! Size error at input!" << std::endl;
		return;
	}
#endif
	amc_real csum, dp;
	amc_int i,j,k;
	
	// compute RHS b := Q^T b (note that b possibly has multiple columns corresponding to multiple least-squares problems with a common LHS)
	for(j = 0; j < b.cols(); j++)
		for(k = 0; k < n; k++)
		{
			dp = 0;
			for(i = k; i < m; i++)
				dp += v[k][i-k]*b(i,j);
			
			for(i = k; i < m; i++)
				b(i,j) -= 2.0*dp*v[k][i-k];
		}

	// back-substitution
	for(j = 0; j < b.cols(); j++)
	{
		x(n-1,j) = b(n-1,j)/R.get(n-1,n-1);
		for(i = n-2; i >= 0; i--)
		{
			csum = 0;
			for(k = i+1; k < n; k++)
				csum += R.get(i,k)*x.get(k,j);
			x(i,j) = (b.get(i,j) - csum)/R.get(i,i);
		}
	}
}

#ifdef EIGEN_LIBRARY
/** b, the RHS, can contain several columns corresponding to different RHS vectors, but
 * \note Eigen3 JacobiSVD::solve() does not necessarily support this!
 */
void leastSquares_SVD(Matrix<amc_real>& A, Matrix<amc_real>& b, Matrix<amc_real>& x)
{
	int m = A.rows(), n = A.cols(), nrhs = b.cols();
	int i,j,k;
#	if DEBUG==1
	if(b.rows() != m || x.rows() != n) 
	{
		std::cout << "leastSquares_NE(): ! Size error at input!" << std::endl;
		return;
	}
#	endif

	// first convert A and b into Eigen3 format
	Eigen::MatrixXd Ae(m,n), be(m,nrhs);
	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
			Ae(i,j) = A.get(i,j);
		for(j = 0; j < nrhs; j++)
			be(i,j) = b.get(i,j);
	}

	Eigen::JacobiSVD<Eigen::MatrixXd> esvd(Ae, Eigen::ComputeThinU|Eigen::ComputeThinV);
	auto so = esvd.solve(be);
	Eigen::MatrixXd sol;
	so.evalTo(sol);
	for(i = 0; i < n; i++)
		for(j = 0; j < nrhs; j++)
			x(i,j) = sol(i,j);
}
#endif

}
