/* A library to solve linear systems.
   Aditya Kashi
   Feb 2015
*/

#ifndef __AMATRIX2_H
#include "amatrix2.hpp"
#endif

#ifndef __ASPARSEMATRIX_H
#include "asparsematrix.hpp"
#endif

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifdef _OPENMP
#ifndef OMP_H
#include <omp.h>
#endif
#endif

#define __ALINALG_H 1

using namespace std;
using namespace amat;

namespace amat
{

typedef MatrixCOO SpMatrix;

/* Note: Cholesky algorithm only implemented for a row-major matrix */
Matrix<double> cholesky(Matrix<double> A, Matrix<double> b)
{
	Matrix<double> B;

	cout << "\ncholesky: Input LHS matrix is " << A.rows() << " x " << A.cols() << endl;
	if(A.rows() != b.rows()) { cout << "\nInvalid dimensions of A and b!"; return B; }
	int N = A.rows();

	//Part 1: Cholesky decomposition
	B.setup(N,N,ROWMAJOR); B.zeros();

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

		if(bjk_sum >= A(j,j)) cout << "\n! cholesky: Negative argument to sqrt at ("<<j<<","<<j<<")\n";
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
			cout << "\ncholesky: Element (" << i <<"," << i << ") of lower triangular matrix is near zero!";
			return B;
		}

	// Part 2: forward substitution to obtain intermediate vector y
	Matrix<double> y(N,1,ROWMAJOR);

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
	//Matrix<double> f(N,1,ROWMAJOR);
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

Matrix<double> gausselim(Matrix<double> A, Matrix<double> b, double tol=1e-12)
{
	//cout << "gausselim: Input LHS matrix is " << A.rows() << " x " << A.cols() << endl;
	if(A.rows() != b.rows()) { cout << "gausselim: Invalid dimensions of A and b!\n"; return A; }
	int N = A.rows();

	Matrix<double> P(N,N,ROWMAJOR);			// Permutation matrix
	P.identity();

	for(int i = 0; i < N-1; i++)
	{
		bool pivot_found = false;
		int k = i;
		if(dabs(A(i,i)) < tol)
		{
			for(int j = i+1; j < N; j++)
			{
				if(dabs(A(j,i)) >= tol)
				{
					pivot_found = true;
					//interchange rows i and j
					for(int k = i; k < N; k++)
					{
						double temp = A(i,k);
						A(i,k) = A(j,k);
						A(j,k) = temp;
					}
					// do the interchange for P and b as well
					for(int k = 0; k < N; k++)
					{
						double temp = P(i,k);
						P(i,k) = P(j,k);
						P(j,k) = temp;
					}
					double temp = b(i,0);
					b(i,0) = b(j,0);
					b(j,0) = temp;
					break;
				}
			}
			if(pivot_found == false) { cout << "! gausselim: Pivot not found!!\n";}
		}

		for(int j = i+1; j < N; j++)
		{
			double ff = A(j,i);
			for(int l = i; l < N; l++)
				A(j,l) = A(j,l) - ff/A(i,i)*A(i,l);
			b(j,0) = b(j,0) - ff/A(i,i)*b(i,0);
		}
	}
	//Thus, A has been transformed to an upper triangular matrix, b has been transformed accordingly, and Permutation matrix P has been created.

	//Part 2: back substitution to obtain final solution
	// Note: the solution is stored in x
	Matrix<double> x(N,1);
	x.zeros();
	//Matrix<double> f(N,1,ROWMAJOR);
	x(N-1,0) = b(N-1,0)/A(N-1,N-1);

	for(int i = N-2; i >= 0; i--)
	{
		double sum = 0;
		int k = i+1;
		do
		{	sum += A(i,k)*x(k,0);		// or A(k,i) ??
			k++;
		} while(k <= N-1);
		x(i,0) = (b(i,0) - sum)/A(i,i);
	}
	return x;
}

Matrix<double> LUdecomp(Matrix<double> A, Matrix<double> b)
{
	//Computes solution of Ax=b by LU decomposition with partial (row) pivoting

	cout << "LUdecomp: Input LHS matrix is " << A.rows() << " x " << A.cols() << endl;
	if(A.rows() != b.rows()) { cout << "LUdecomp: Invalid dimensions of A and b!\n"; return A; }
	int N = A.rows();

	Matrix<double> L(N, N, ROWMAJOR);
	L.identity();
	Matrix<double> P(N, N, ROWMAJOR);
	P.identity();

	for(int k = 0; k <= N-2; k++)
	{
		int i = k;
		for(int j = k; j <= N-1; j++)	// find row i s.t. A(i,k) is max among row k (A(k,k)) to last row (A(N-1,k))
			if(dabs(A(j,k)) > dabs(A(i,k))) i = j;
		if(i != k)						// interchange rows i and k from column k to column N-1 (last column)
		{
			for(int l = k; l = N-1; l++)
			{
				double temp = A(i,l);
				A(i,l) = A(k,l);
				A(k,l) = temp;
			}
			for(int l = 0; l = k-1; l++)
			{
				double temp = L(i,l);
				L(i,l) = L(k,l);
				L(k,l) = temp;
			}
			for(int l = 0; l = N-1; l++)
			{
				double temp = P(i,l);
				P(i,l) = P(k,l);
				P(k,l) = temp;
			}
		}
	}
	//INCOMPLETE
	return A;
}



//-------------------- Iterative Methods ----------------------------------------//

Matrix<double> pointjacobi(Matrix<double> A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check)
{
	cout << "\npointjacobi: Input LHS matrix is " << A.rows() << " x " << A.cols() << endl;
	if(A.rows() != b.rows()) { cout << "! pointjacobi: Invalid dimensions of A and b!\n"; return b; }
	if(xold.rows() != b.rows()) { cout << "! pointjacobi: Invalid dimensions on xold !\n"; return b; }
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
			if(dabs(A(i,i)) <= sum) cout << "* pointjacobi(): LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
		}
	}

	Matrix<double> x(b.rows(),1);

	Matrix<double> M(N,N,ROWMAJOR);		// diagonal matrix
	M.zeros();

	//populate diagonal matrix M
	for(int i = 0; i < N; i++)
		M(i,i) = A(i,i);

	A = A - M;

	//invert M
	for(int i = 0; i < N; i++)
		M(i,i) = 1.0/M(i,i);

	//x.zeros();
	//cout << "pointjacobi: Initial error = " << (xold-x).dabsmax() << endl;

	x = xold;
	int c = 0;

	do
	{
		xold = x;
		x = M * (b - (A*xold));
		c++;
		if(c > maxiter) { cout << "pointjacobi: Max iterations exceeded!\n"; break; }
	} while((x-xold).dabsmax() >= tol);

	/* double error = 1.0;
	do
	{
		xold = x;

		Matrix<double> Axold(N,1);
		Axold.zeros();
		for(int i = 0; i < N; i++)
		{
			for(int k = 0; k < N; k++)
				Axold(i,0) += A(i,k)*xold(k,0);
		}
		Matrix<double> inter(N,1);
		for(int i = 0; i < N; i++)
			inter(i,0) = b(i,0) - Axold(i,0);
		for(int i = 0; i < N; i++)
			x(i,0) = M(i,i) * inter(i,0);

		c++;
		if(c > maxiter) { cout << "pointjacobi: Max iterations exceeded!\n"; break; }
		Matrix<double> diff(N,1);	// diff = x - xold
		for(int i = 0; i < N; i++)
			diff(i,0) = x(i,0) - xold(i,0);

		error = dabs(diff(0,0));
		for(int i = 1; i < N; i++)
			if(dabs(diff(i,0)) > error) error = dabs(diff(i,0));
		//cout << "pointjacobi: error = " << error << endl;

	} while(error >= tol); */
	cout << "pointjacobi: No. of iterations = " << c << endl;

	return x;
}

Matrix<double> gaussseidel(Matrix<double> A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check)
{
	//cout << "\ngaussseidel: Input LHS matrix is " << A.rows() << " x " << A.cols() << endl;
	if(A.rows() != b.rows()) { cout << "! gaussseidel: Invalid dimensions of A and b!\n"; return b; }
	if(xold.rows() != b.rows()) { cout << "! gaussseidel: Invalid dimensions on xold !\n"; return b; }
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
			if(dabs(A(i,i)) <= sum) cout << "* gaussseidel: LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
		}
	}

	Matrix<double> x(b.rows(),1);

	Matrix<double> M(N,N,ROWMAJOR);		// diagonal matrix
	M.zeros();

	//populate diagonal matrix M
	for(int i = 0; i < N; i++)
		M(i,i) = A(i,i);

	A = A - M;						// diagonal entries of A are now zero

	//invert M
	for(int i = 0; i < N; i++)
		M(i,i) = 1.0/M(i,i);

	//x.zeros();
	//cout << "gaussseidel: Initial error = " << (xold-x).dabsmax() << endl;

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
		if(c > maxiter) { cout << "gaussseidel: Max iterations exceeded!\n"; break; }
		for(int i = 0; i < N; i++)
			diff(i,0) = x(i,0) - xold(i,0);

		error = dabs(diff(0,0));
		for(int i = 1; i < N; i++)
			if(dabs(diff(i,0)) > error) error = dabs(diff(i,0));
		//cout << "gaussseidel: error = " << error << endl;
		if(first == true)
		{	initres = error;
			//cout << "gaussseidel: Initial residue = " << initres << endl;
			first = false;
		}

		//if(c%20 == 0) cout << "gaussseidel(): Step " << c << ", Relative error = " << error/initres << endl;

	} while(error/initres >= tol);
	//cout << "gaussseidel: No. of iterations = " << c << endl;

	return x;
}

Matrix<double> sparsegaussseidel(SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check)
{
	cout << "sparsegaussseidel(): Input LHS matrix is " << A->rows() << " x " << A->cols() << endl;
	if(A->rows() != b.rows()) { cout << "! gaussseidel: Invalid dimensions of A and b!\n"; return b; }
	if(xold.rows() != b.rows()) { cout << "! gaussseidel: Invalid dimensions on xold !\n"; return b; }
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
			if(dabs(A->get(i,i)) <= sum) cout << "* gaussseidel: LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
		}
	}

	Matrix<double> x(b.rows(),1);

	Matrix<double> M(N,1);		// vector of diagonal elements of A
	M.zeros();

	// diagonal matrix M
	//cout << "sparsegaussseidel(): Getting diagonal of sparse matrix\n";
	A->get_diagonal(&M);

	//invert M
	for(int i = 0; i < N; i++)
		M(i) = 1.0/M(i);

	//x.zeros();
	//cout << "gaussseidel: Initial error = " << (xold-x).dabsmax() << endl;

	x = xold;
	int c = 0;
	double initres;
	bool first = true;

	Matrix<double> Axold(N,1);
	Matrix<double> inter(N,1);
	Matrix<double> diff(N,1);	// diff = x - xold
	double error = 1.0;
	//cout << "sparsegaussseidel(): Starting iterations\n";
	do
	{
		xold = x;
		Axold.zeros();
		for(int i = 0; i < N; i++)
		{
			/*for(int k = 0; k < i; k++)
				Axold(i,0) += A(i,k)*x(k,0);
			for(int k = i; k < N; k++)
				Axold(i,0) += A(i,k)*xold(k,0);*/
			
			//cout << "  Calling sparse getelem_multiply_parts\n";
			Axold(i) = A->getelem_multiply_parts(i, &x, &xold, i);
			//cout << "  Setting inter\n";
			inter(i,0) = b(i,0) - Axold(i,0);
			x(i,0) = M(i) * inter(i,0);
		}
		// NOTE: The above results in FORWARD Gauss-Seidel

		c++;
		if(c > maxiter) { cout << "gaussseidel: Max iterations exceeded!\n"; break; }
		for(int i = 0; i < N; i++)
			diff(i,0) = x(i,0) - xold(i,0);

		error = dabs(diff(0,0));
		for(int i = 1; i < N; i++)
			if(dabs(diff(i,0)) > error) error = dabs(diff(i,0));
		//cout << "gaussseidel: error = " << error << endl;
		if(first == true)
		{	initres = error;
			//cout << "gaussseidel: Initial residue = " << initres << endl;
			first = false;
		}

		if(c%10 == 0) cout << "gaussseidel(): Step " << c << ", Relative error = " << error/initres << endl;

	} while(error/initres >= tol);
	cout << "gaussseidel: No. of iterations = " << c << endl;

	return x;
}


}
