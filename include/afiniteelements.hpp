/** \file afiniteelements.hpp
 * \brief Classes for finite element basis functions 
 * \author Aditya Kashi
 * \date Feb 29, 2016
 */

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif

/// Namespace for functionality related to finite element basis functions
namespace afem {

/// Abstract class for getting shape function values on an arbitrary reference element
class ReferenceFiniteElement
{
protected:
	int ndim;				///< Dimension of the element
	int nnode;				///< Number of nodes in the element
	int ngauss;				///< Number of Gauss points
public:
	virtual double getBasis(const int ibasisnum, const std::vector& r) const = 0;
	virtual void getBasisDeriv(const int ibasisnum, const std::vector<double>& r. std::vector<double>& derivs) const = 0;
	virtual double getDetJacobian(const std::vector<double>& r) const = 0;
};

class LagrangeFiniteElement : public ReferenceFiniteElement
{
};
	
}	// end namespace afem
