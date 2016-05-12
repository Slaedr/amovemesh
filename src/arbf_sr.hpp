/** @file arbf_sr.hpp
 * @brief A mesh-movement method using radial basis function interpolation, but using different support radii for each interpolation center.
 * 
 * based on the 2007 paper by de Boer, van der Schoot and Bijl,
 * but with one major difference - the interpolation is done using only RBFs; the polynomial part is ignored.
 * @author Aditya Kashi
 * @date May 9, 2016
 */

#ifndef __ARBF_SR_H

#ifndef _GLIBCXX_CSTDIO
#include <cstdio>
#endif

#ifndef _GLIBCXX_STRING
#include <string>
#endif

#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

#define __ARBF_SR_H 1

namespace amc {

/// Movement of a point-cloud based on radial basis function interpolation of some other points called interpolation centers
/*! Requires specification of support radii separately at each interpolation center (bpoints)
 */
class RBFmove
{
	amat::Matrix<amc_real>* const inpoints;
	amat::Matrix<amc_real>* const bpoints;
	const amat::Matrix<amc_real>* const bmotion;
	amat::Matrix<int> bflag;
	int npoin;			///< total number of points
	int ninpoin;		///< number of interior points
	int nbpoin;			///< number of boundary points
	int ndim;
	double (RBFmove::*rbf)(double, amc_int);
	const amat::Matrix<amc_real>* const srad;

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
	//RBFmove();

	/// Sets the data needed
	/** Note that all parameters are deep-copied.
	 * \param int_points is a list of all interior points to be moved
	 * \param boun_points is the array of boundary points
	 * \param boundary_motion is nbpoin-by-ndim array - containing displacements corresponding to boundary points.
	 * \param rbf_ch indicates the RBF to use - 0 : C0, 2 : C2, 4 : C4, default : Gaussian
	 * \param support_radius contains the support radius to use for each boundary point
	 * \param num_steps is the number of steps in which to break up the movement to perform separately (sequentially)
	 * \param linear_solver indicates the linear solver to use to solve the RBF equations - "DLU", "CG", "LU"
	 */
	RBFmove(amat::Matrix<double>* int_points, amat::Matrix<double>* boun_points, amat::Matrix<double>* boundary_motion, const int rbf_ch, const amat::Matrix<amc_real>* const support_radius, 
			const int num_steps, const double tolerance, const int iter, const std::string linear_solver);

	~RBFmove();

	// Specific RBFs
	/// Wendland's compact C2 function - most tested
	double rbf_c2_compact(double xi, amc_int ibp);

	double rbf_c0(double xi, amc_int ibp);
	double rbf_c4(double xi, amc_int ibp);
	double gaussian(double xi, amc_int ibp);

	/// Assembles the LHS matrix.
	void assembleLHS();

	/// Executes 1 step of the mesh movement.
	void move_step();

	/// Uses [move_step](@ref move_step) to execute the total number of steps specified.
	/** This function calls all other required functions, so it should be used directly after setup.
	 */
	void move();

	/// Returns new positions of interior points.
	amat::Matrix<double>* getInteriorPoints();

	/// Returns new positions of boundary points.
	amat::Matrix<double>* getBoundaryPoints();
};

} // end namespace
#endif
