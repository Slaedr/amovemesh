/* A library of functions to solve linear systems, linear least-squares problems etc.
   Aditya Kashi
   Feb 2015
*/

#ifndef __ALINALG_H

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif

#ifndef __ASPARSEMATRIX_H
#include <asparsematrix.hpp>
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

namespace amat
{

/// Computes solution of Ax = b by Gaussian elimination. Reasonably well-tested.
/** A is mxm, b is mxk where k is the number of systems to be solved with the same LHS.
*/
void gausselim(Matrix<double>& A, Matrix<double>& b, Matrix<double>& x);

#ifdef EIGEN_LIBRARY
/// Uses Eigen3's supernodal sparse LU solver to solve Ax = b
Matrix<double> gausselim(const SpMatrix& A, const Matrix<amc_real>& b);

/// Uses Eigen3's interface to PastiX to solve a SPD system
void pastix_LDLT(const SpMatrix& A, const Matrix<amc_real>& b, Matrix<amc_real>& x);

/// Uses Eigen3's SVD module to solve linear least-squares 
void leastSquares_SVD(Matrix<amc_real>& A, Matrix<amc_real>& b, Matrix<amc_real>& x);
#endif

// Uses the SuperLU direct sparse solver to solve Ax = b and stores the solution in ans
//void superLU_solve(const SpMatrix* A, const Matrix<double>* b, Matrix<double>* ans);

/// Old (somewhat tested) implementation of Cholesky direct solver
Matrix<double> cholesky(Matrix<double> A, Matrix<double> b);

/// Solves Ax = b by Cholesky decomposition
/** The output is finally stored in b.
 */
void chol(Matrix<amc_real>& A, Matrix<amc_real>& b);


/// Solves Ax = b for dense A by point-Jacobi iterations
/** \param check should be set to 'y' for testing diagonal dominance before solving.
 */
Matrix<double> pointjacobi(Matrix<double> A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check='n');

/// Solves Ax = b for dense A by forward Gauss-Seidel iterations
/** \param check should be set to 'y' for testing diagonal dominance before solving.
 */
Matrix<double> gaussseidel(Matrix<double> A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check='n');

/// Solves Ax=b for sparse A by forward Gauss-Seidel iterations.
Matrix<double> sparsegaussseidel(SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, char check='n');

/// Solves Ax=b for sparse A by forward SOR.
Matrix<double> sparseSOR(SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter, double w=1.25, char check='n');


/// Calculates solution of Ax=b where A is a SPD matrix in sparse format. This function is well-tested.
/** The preconditioner is a diagonal matrix.
 * NOTE: The parallel version is actually slower, due to some reason.
 */
Matrix<double> sparseCG_d(SpMatrix* A, Matrix<double> b, Matrix<double> xold, double tol, int maxiter);

/**	Solves general linear system Ax=b using stabilized biconjugate gradient method of van der Vorst.
 * NOTE: The initial guess vector xold is modified by this function.
 */
Matrix<double> sparse_bicgstab(const SpMatrix* A, const Matrix<double>& b, Matrix<double>& xold, double tol, int maxiter);


/// solves the least squares problem (finds the minimum point x) \f$ \min(||Ax - b||_2) \f$ by solving the normal equations
void leastSquares_NE(amat::Matrix<amc_real>& A, amat::Matrix<amc_real>& b, amat::Matrix<amc_real>& x);

/// Computes the QR decomposition of a dense matrix by Householder algorithm given in Trefethen and Bau, Numerical Linear Algebra.
/** \param A is an m x n matrix, replaced by the upper triangular matrix R at the end of the algorithm.
 * \param v is the set of reflection vectors which can be used to compute Q, the action of Q and the action of Q*.
 * v should contain n vectors: The first vector with size m, second with size m-1,...upto m-n+1.
 */
void qr(amat::Matrix<amc_real>& A, std::vector<amc_real>* v);

/// Computes solution to a linear least-squares problem by [QR decomposition](@ref qr)
void leastSquares_QR(amat::Matrix<amc_real>& A, amat::Matrix<amc_real>& b, amat::Matrix<amc_real>& x);

/// Given factors Q and R, solve QRx = b
/** \param v is the set of vectors that determines Q
 * \param R is the mxn upper triangular R
 * \param b is the m-vector RHS
 * \param x is the n-vector solution
 */
void solve_QR(const std::vector<amc_real>* v, const amat::Matrix<amc_real>& R, amat::Matrix<amc_real>& b, amat::Matrix<amc_real>& x);

}	// end namespace

#endif

