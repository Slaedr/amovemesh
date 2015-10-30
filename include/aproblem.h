/** Abstract problem class.
	Aditya Kashi
	October 19, 2015
*/

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif

#ifndef __ASPARSEMATRIX_H
#include <asparsematrix.hpp>
#endif

#define __APROBLEM_H 1

using namespace amat;

/** Abstract class from which to derive any class that requires solution to a linear system.
	Pointers to objects of this class can be passed to linear-solver functions, while the pointers actually point to the objects of classes that inherit from this one.
	Virtual functions of this class (implemented in derived classes) can then be invoked in the linear solver.
*/
class Problem
{
public:
	virtual void lhs_action(Matrix<double>* x, Matrix<double>* ans) = 0;
	virtual Matrix<double> get_rhs() = 0;
};
