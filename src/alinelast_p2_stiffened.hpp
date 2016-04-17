/** @brief Funtionality for solving equations of linear elasticity on a straight-edged quadratic mesh using P2 Lagrange finite elements (no non-homogeneous Neumann BCs)
 * We use a 'stiffened' formulation, in which smaller elements are stiffer (undergo less deformation) than larger elements. Other kinds of stiffening are also possible, based on some mesh quality metrics, for instance.
 * @date January 27, 2016
 * @author Adtya Kashi
 */

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif
#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif
#ifndef _GLIBCXX_STRING
#include <string>
#endif
#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif
#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif
#ifdef _OPENMP
#ifndef OMP_H
#include <omp.h>
#endif
#endif

#ifndef __AMATRIX_H
#include <amatrix.hpp>
#endif
#ifndef __AMESH2DGENERAL_H
#include <amesh2d.hpp>
#endif

#define __ALINELAST_P2_STIFFENED_H 1

namespace amc {

/// Class implementing solution of linear elasticity system by P2 Lagrange finite elements
/** Currently only works on 2D triangular meshes.
 */
class LinElastP2
{

	UMesh2d* m;
	int ngeoel;
	int ngeofa;
	int ndofe;					///< Number of DOFs per element
	amat::Matrix<double> geoel;		///< holds 2*area of element, and derivatives of barycentric coordinate functions lambdas
	amat::Matrix<double> geofa;		///< holds normals to and length of boundary faces
	amat::SpMatrix K;					///< global stiffness matrix
	amat::Matrix<double> f;			///< global load vector
	amat::Matrix<double>* stiff;		///< stiffening factor for each element

	double muE;					///< isotropic elasticity constant mu
	double lambdaE;				///< isotropic elasticity constant lambda
	double j0;					///< scaling constant for the jacobian-based stiffening
	double chi;					///< exponent for jacobian-based stiffening
	double cbig;				///< for Dirichlet BCs

public:

	/// Sets inputs and computes derivatives of basis functions and face-normals
	/** \note The computations are only for straight-sided P2 elements!
	 * \param xch Exponent to which the stiffening factor is raised in each element stiffness entry.
	 * \param stiffm An array containing the stiffening factor to be used for each element.
	 */
	void setup(UMesh2d* mesh, double mu, double lambd, double xch, amat::Matrix<double>* stiffm)
	{
		m = mesh;
		ngeoel = 7;
		ngeofa = 3;
		muE = mu; lambdaE = lambd;
		chi = xch;
		geoel.setup(m->gnelem(), ngeoel);
		geofa.setup(m->gnface(), ngeofa);
		K.setup(m->gndim()*m->gnpoin(), m->gndim()*m->gnpoin());
		f.setup(m->gndim()*m->gnpoin(),1);
		stiff = stiffm;
		cbig = 1e30;
		
		//loadvec.setup(m.gnpoin(), 1);

		for(int i = 0; i < m->gnelem(); i++)
		{
			geoel(i,0) = m->gcoords(m->ginpoel(i,0),0)*(m->gcoords(m->ginpoel(i,1),1) - m->gcoords(m->ginpoel(i,2),1)) - m->gcoords(m->ginpoel(i,0),1)*(m->gcoords(m->ginpoel(i,1),0)-m->gcoords(m->ginpoel(i,2),0)) + m->gcoords(m->ginpoel(i,1),0)*m->gcoords(m->ginpoel(i,2),1) - m->gcoords(m->ginpoel(i,2),0)*m->gcoords(m->ginpoel(i,1),1);
			double D = geoel(i,0);

			//a1 = (y2 - y3) :
			geoel(i,1) = (m->gcoords(m->ginpoel(i,1), 1) - m->gcoords(m->ginpoel(i,2), 1));
			//a2 = (y3 - y1):
			geoel(i,2) = (m->gcoords(m->ginpoel(i,2), 1) - m->gcoords(m->ginpoel(i,0), 1));
			//a3 = (y1 - y2) :
			geoel(i,3) = -geoel(i,1)-geoel(i,2);
			//b1 = (x3 - x2) :
			geoel(i,4) = (m->gcoords(m->ginpoel(i,2), 0) - m->gcoords(m->ginpoel(i,1), 0));
			//b2 = (x1 - x3) :
			geoel(i,5) = (m->gcoords(m->ginpoel(i,0), 0) - m->gcoords(m->ginpoel(i,2), 0));
			//b3 = (x2 - x1) :
			geoel(i,6) = -geoel(i,5) - geoel(i,4);
			
			j0 += geoel.get(i,0);
		}
		
		// set j0 as the average jacobian
		j0 /= m->gnelem();
		
		for(int i = 0; i < m->gnface(); i++)
		{
			//n_x = y_2 - y_1
			geofa(i,0) = m->gcoords(m->gbface(i,2), 1) - m->gcoords(m->gbface(i,0), 1);
			//n_y = x_1 - x_2
			geofa(i,1) = m->gcoords(m->gbface(i,0), 0) - m->gcoords(m->gbface(i,2), 0);
			//l = sqrt(n_x^2 + n_y^2)
			geofa(i,2) = sqrt(geofa(i,0)*geofa(i,0) + geofa(i,1)*geofa(i,1));
			//geofa(i,2) = sqrt((m.gcoords(m.gbface(i,1)-1, 1) - m.gcoords(m.gbface(i,0)-1, 1))*(m.gcoords(m.gbface(i,1)-1, 1) - m.gcoords(m.gbface(i,0)-1, 1)) + (m.gcoords(m.gbface(i,0)-1, 0) - m.gcoords(m.gbface(i,1)-1, 0))*(m.gcoords(m.gbface(i,0)-1, 0) - m.gcoords(m.gbface(i,1)-1, 0)));
		}

		std::cout << "LinElastP2: Computed derivatives of basis functions, and normals to and lengths of boundary faces.\n";
	}

	amat::Matrix<double> elementstiffnessK11(int iel)
	{
		amat::Matrix<double> K11(6,6);
		double coeff = 2*muE+lambdaE;
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
			{
				if(i==j)
					K11(i,i) = (coeff*geoel(iel,i+1)*geoel(iel,i+1) + muE*geoel(iel,i+4)*geoel(iel,i+4))/(2*geoel(iel,0)) * stiff->get(iel);	// ai^2 / 2D
				else
					K11(i,j) = -1.0*(coeff*geoel(iel,i+1)*geoel(iel,j+1) + muE*geoel(iel,i+4)*geoel(iel,j+4))/(6*geoel(iel,0)) * stiff->get(iel);		// (c*ai*aj + mu*bi*bj) / 6D
			}
		K11(0,3) = (coeff*2*geoel(iel,1)*geoel(iel,2) + muE*2*geoel(iel,4)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K11(0,4) = 0.0;
		K11(0,5) = (coeff*2*geoel(iel,1)*geoel(iel,3) + muE*2*geoel(iel,4)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K11(1,3) = (coeff*2*geoel(iel,2)*geoel(iel,1) + muE*2*geoel(iel,5)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);
		K11(1,4) = (coeff*2*geoel(iel,2)*geoel(iel,3) + muE*2*geoel(iel,5)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K11(1,5) = 0.0;
		K11(2,3) = 0.0;
		K11(2,4) = (coeff*2*geoel(iel,3)*geoel(iel,2) + muE*2*geoel(iel,6)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K11(2,5) = (coeff*2*geoel(iel,3)*geoel(iel,1) + muE*2*geoel(iel,6)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);

		K11(3,3) = 4/(3*geoel(iel,0))*( coeff*(geoel(iel,2)*geoel(iel,2) + geoel(iel,1)*geoel(iel,1) + geoel(iel,2)*geoel(iel,1)) + muE*(geoel(iel,5)*geoel(iel,5)+geoel(iel,4)*geoel(iel,4) + geoel(iel,5)*geoel(iel,4)) ) * stiff->get(iel);
		K11(4,4) = 4/(3*geoel(iel,0))* ( coeff*(geoel(iel,3)*geoel(iel,3) + geoel(iel,2)*geoel(iel,2) + geoel(iel,3)*geoel(iel,2)) + muE*(geoel(iel,6)*geoel(iel,6) + geoel(iel,5)*geoel(iel,5) + geoel(iel,6)*geoel(iel,5)) ) * stiff->get(iel);
		K11(5,5) = 4/(3*geoel(iel,0))* ( coeff*(geoel(iel,1)*geoel(iel,1) + geoel(iel,3)*geoel(iel,3) + geoel(iel,1)*geoel(iel,3)) + muE*(geoel(iel,4)*geoel(iel,4) + geoel(iel,6)*geoel(iel,6) + geoel(iel,4)*geoel(iel,6)) ) * stiff->get(iel);

		K11(3,4) = 2/(3*geoel(iel,0))* ( coeff*(geoel(iel,2)*geoel(iel,3) + geoel(iel,2)*geoel(iel,2) + 2*geoel(iel,1)*geoel(iel,3) + geoel(iel,1)*geoel(iel,2)) + muE*(geoel(iel,5)*geoel(iel,6) + geoel(iel,5)*geoel(iel,5) + 2*geoel(iel,4)*geoel(iel,6) + geoel(iel,4)*geoel(iel,5)) ) * stiff->get(iel);
		K11(3,5) = 2/(3*geoel(iel,0))* ( coeff*(geoel(iel,2)*geoel(iel,1) + 2*geoel(iel,2)*geoel(iel,3) + geoel(iel,1)*geoel(iel,1) + geoel(iel,1)*geoel(iel,3)) + muE*(geoel(iel,5)*geoel(iel,4) + 2*geoel(iel,5)*geoel(iel,6) + geoel(iel,4)*geoel(iel,4) + geoel(iel,4)*geoel(iel,6)) ) * stiff->get(iel);
		K11(4,5) = 2/(3*geoel(iel,0))* ( coeff*(geoel(iel,3)*geoel(iel,1) + geoel(iel,3)*geoel(iel,3) + 2*geoel(iel,2)*geoel(iel,1) + geoel(iel,2)*geoel(iel,3)) + muE*(geoel(iel,6)*geoel(iel,4) + geoel(iel,6)*geoel(iel,6) + 2*geoel(iel,5)*geoel(iel,4) + geoel(iel,5)*geoel(iel,6)) ) * stiff->get(iel);

		//K11 is symmetric
		for(int i = 0; i < 5; i++)
			for(int j = i+1; j < 6; j++)
				K11(j,i) = K11(i,j);

		return K11;
	}

	amat::Matrix<double> elementstiffnessK22(int iel)
	{
		amat::Matrix<double> K22(6,6);
		double coeff = 2*muE+lambdaE;
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
			{
				if(i==j)
					K22(i,i) = (muE*geoel(iel,i+1)*geoel(iel,i+1) + coeff*geoel(iel,i+4)*geoel(iel,i+4))/(2*geoel(iel,0)) * stiff->get(iel);	// ai^2 / 2D
				else
					K22(i,j) = -1.0*(muE*geoel(iel,i+1)*geoel(iel,j+1) + coeff*geoel(iel,i+4)*geoel(iel,j+4))/(6*geoel(iel,0)) * stiff->get(iel);		// -(mu*ai*aj + c*bi*bj) / 6D
			}
		K22(0,3) = (muE*2*geoel(iel,1)*geoel(iel,2) + coeff*2*geoel(iel,4)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K22(0,4) = 0.0;
		K22(0,5) = (muE*2*geoel(iel,1)*geoel(iel,3) + coeff*2*geoel(iel,4)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K22(1,3) = (muE*2*geoel(iel,2)*geoel(iel,1) + coeff*2*geoel(iel,5)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);
		K22(1,4) = (muE*2*geoel(iel,2)*geoel(iel,3) + coeff*2*geoel(iel,5)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K22(1,5) = 0.0;
		K22(2,3) = 0.0;
		K22(2,4) = (muE*2*geoel(iel,3)*geoel(iel,2) + coeff*2*geoel(iel,6)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K22(2,5) = (muE*2*geoel(iel,3)*geoel(iel,1) + coeff*2*geoel(iel,6)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);

		K22(3,3) = 4/(3*geoel(iel,0))* ( muE*(geoel(iel,2)*geoel(iel,2) + geoel(iel,1)*geoel(iel,1) + geoel(iel,2)*geoel(iel,1)) + coeff*(geoel(iel,5)*geoel(iel,5) + geoel(iel,4)*geoel(iel,4) + geoel(iel,5)*geoel(iel,4)) ) * stiff->get(iel);
		K22(4,4) = 4/(3*geoel(iel,0))* ( muE*(geoel(iel,3)*geoel(iel,3) + geoel(iel,2)*geoel(iel,2) + geoel(iel,3)*geoel(iel,2)) + coeff*(geoel(iel,6)*geoel(iel,6) + geoel(iel,5)*geoel(iel,5) + geoel(iel,6)*geoel(iel,5)) ) * stiff->get(iel);
		K22(5,5) = 4/(3*geoel(iel,0))* ( muE*(geoel(iel,1)*geoel(iel,1) + geoel(iel,3)*geoel(iel,3) + geoel(iel,1)*geoel(iel,3)) + coeff*(geoel(iel,4)*geoel(iel,4) + geoel(iel,6)*geoel(iel,6) + geoel(iel,4)*geoel(iel,6)) ) * stiff->get(iel);

		K22(3,4) = 2/(3*geoel(iel,0))* ( muE*(geoel(iel,2)*geoel(iel,3) + geoel(iel,2)*geoel(iel,2) + 2*geoel(iel,1)*geoel(iel,3) + geoel(iel,1)*geoel(iel,2)) + coeff*(geoel(iel,5)*geoel(iel,6) + geoel(iel,5)*geoel(iel,5) + 2*geoel(iel,4)*geoel(iel,6) + geoel(iel,4)*geoel(iel,5)) ) * stiff->get(iel);
		K22(3,5) = 2/(3*geoel(iel,0))* ( muE*(geoel(iel,2)*geoel(iel,1) + 2*geoel(iel,2)*geoel(iel,3) + geoel(iel,1)*geoel(iel,1) + geoel(iel,1)*geoel(iel,3)) + coeff*(geoel(iel,5)*geoel(iel,4) + 2*geoel(iel,5)*geoel(iel,6) + geoel(iel,4)*geoel(iel,4) + geoel(iel,4)*geoel(iel,6)) ) * stiff->get(iel);
		K22(4,5) = 2/(3*geoel(iel,0))* ( muE*(geoel(iel,3)*geoel(iel,1) + geoel(iel,3)*geoel(iel,3) + 2*geoel(iel,2)*geoel(iel,1) + geoel(iel,2)*geoel(iel,3)) + coeff*(geoel(iel,6)*geoel(iel,4) + geoel(iel,6)*geoel(iel,6) + 2*geoel(iel,5)*geoel(iel,4) + geoel(iel,5)*geoel(iel,6)) ) * stiff->get(iel);

		//K22 is symmetric
		for(int i = 0; i < 5; i++)
			for(int j = i+1; j < 6; j++)
				K22(j,i) = K22(i,j);

		return K22;
	}

	amat::Matrix<double> elementstiffnessK12(int iel)
	{
		amat::Matrix<double> K12(6,6);
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
			{
				if(i==j)
					K12(i,i) = (lambdaE*geoel(iel,i+1)*geoel(iel,i+4) + muE*geoel(iel,i+1)*geoel(iel,i+4))/(2*geoel(iel,0)) * stiff->get(iel);	// ai^2 / 2D
				else
					K12(i,j) = -1.0*(lambdaE*geoel(iel,i+4)*geoel(iel,j+1) + muE*geoel(iel,i+1)*geoel(iel,j+4))/(6*geoel(iel,0)) * stiff->get(iel);		// -(c*bi*aj + mu*ai*bj) / 6D
			}

		K12(0,3) = (muE*2*geoel(iel,4)*geoel(iel,2) + lambdaE*2*geoel(iel,1)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(0,4) = 0.0;
		K12(0,5) = (muE*2*geoel(iel,4)*geoel(iel,3) + lambdaE*2*geoel(iel,1)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(1,3) = (muE*2*geoel(iel,5)*geoel(iel,1) + lambdaE*2*geoel(iel,2)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(1,4) = (muE*2*geoel(iel,5)*geoel(iel,3) + lambdaE*2*geoel(iel,2)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(1,5) = 0.0;
		K12(2,3) = 0.0;
		K12(2,4) = (muE*2*geoel(iel,6)*geoel(iel,2) + lambdaE*2*geoel(iel,3)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(2,5) = (muE*2*geoel(iel,6)*geoel(iel,1) + lambdaE*2*geoel(iel,3)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);

		K12(3,3) = 2/(3*geoel(iel,0))* ( muE*(2*geoel(iel,2)*geoel(iel,5) + geoel(iel,2)*geoel(iel,4) + geoel(iel,1)*geoel(iel,5) + 2*geoel(iel,1)*geoel(iel,4)) + lambdaE*(2*geoel(iel,2)*geoel(iel,5) + geoel(iel,2)*geoel(iel,4) + geoel(iel,1)*geoel(iel,5) + 2*geoel(iel,1)*geoel(iel,4)) ) * stiff->get(iel);
		K12(4,4) = 2/(3*geoel(iel,0))* ( muE*(2*geoel(iel,3)*geoel(iel,6) + geoel(iel,3)*geoel(iel,5) + geoel(iel,2)*geoel(iel,6) + 2*geoel(iel,2)*geoel(iel,5)) + lambdaE*(2*geoel(iel,3)*geoel(iel,6) + geoel(iel,3)*geoel(iel,5) + geoel(iel,2)*geoel(iel,6) + 2*geoel(iel,2)*geoel(iel,5)) ) * stiff->get(iel);
		K12(5,5) = 2/(3*geoel(iel,0))* ( muE*(2*geoel(iel,1)*geoel(iel,4) + geoel(iel,1)*geoel(iel,6) + geoel(iel,3)*geoel(iel,4) + 2*geoel(iel,3)*geoel(iel,6)) + lambdaE*(2*geoel(iel,1)*geoel(iel,4) + geoel(iel,1)*geoel(iel,6) + geoel(iel,3)*geoel(iel,4) + 2*geoel(iel,3)*geoel(iel,6)) ) * stiff->get(iel);

		K12(3,4) = 2/(3*geoel(iel,0))* ( muE*(geoel(iel,5)*geoel(iel,3) + geoel(iel,5)*geoel(iel,2) + 2*geoel(iel,4)*geoel(iel,3) + geoel(iel,4)*geoel(iel,2)) + lambdaE*(geoel(iel,2)*geoel(iel,6) + geoel(iel,2)*geoel(iel,5) + 2*geoel(iel,1)*geoel(iel,6) + geoel(iel,1)*geoel(iel,5)) ) * stiff->get(iel);
		K12(3,5) = 2/(3*geoel(iel,0))* ( muE*(geoel(iel,5)*geoel(iel,1) + 2*geoel(iel,5)*geoel(iel,3) + geoel(iel,4)*geoel(iel,1) + geoel(iel,4)*geoel(iel,3)) + lambdaE*(geoel(iel,2)*geoel(iel,4) + 2*geoel(iel,2)*geoel(iel,6) + geoel(iel,1)*geoel(iel,4) + geoel(iel,1)*geoel(iel,6)) ) * stiff->get(iel);
		K12(4,5) = 2/(3*geoel(iel,0))* ( muE*(geoel(iel,6)*geoel(iel,1) + geoel(iel,6)*geoel(iel,3) + 2*geoel(iel,5)*geoel(iel,1) + geoel(iel,5)*geoel(iel,3)) + lambdaE*(geoel(iel,3)*geoel(iel,4) + geoel(iel,3)*geoel(iel,6) + 2*geoel(iel,2)*geoel(iel,4) + geoel(iel,2)*geoel(iel,6)) ) * stiff->get(iel);

		// K12 is not symmetric
		K12(3,0) = (lambdaE*2*geoel(iel,4)*geoel(iel,2) + muE*2*geoel(iel,1)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(4,0) = 0.0;
		K12(5,0) = (lambdaE*2*geoel(iel,4)*geoel(iel,3) + muE*2*geoel(iel,1)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(3,1) = (lambdaE*2*geoel(iel,5)*geoel(iel,1) + muE*2*geoel(iel,2)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(4,1) = (lambdaE*2*geoel(iel,5)*geoel(iel,3) + muE*2*geoel(iel,2)*geoel(iel,6))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(5,1) = 0.0;
		K12(3,2) = 0.0;
		K12(4,2) = (lambdaE*2*geoel(iel,6)*geoel(iel,2) + muE*2*geoel(iel,3)*geoel(iel,5))/(3*geoel(iel,0)) * stiff->get(iel);
		K12(5,2) = (lambdaE*2*geoel(iel,6)*geoel(iel,1) + muE*2*geoel(iel,3)*geoel(iel,4))/(3*geoel(iel,0)) * stiff->get(iel);

		K12(4,3) = 2/(3*geoel(iel,0))* ( lambdaE*(geoel(iel,5)*geoel(iel,3) + geoel(iel,5)*geoel(iel,2) + 2*geoel(iel,4)*geoel(iel,3) + geoel(iel,4)*geoel(iel,2)) + muE*(geoel(iel,2)*geoel(iel,6) + geoel(iel,2)*geoel(iel,5) + 2*geoel(iel,1)*geoel(iel,6) + geoel(iel,1)*geoel(iel,5)) ) * stiff->get(iel);
		K12(5,3) = 2/(3*geoel(iel,0))* ( lambdaE*(geoel(iel,5)*geoel(iel,1) + 2*geoel(iel,5)*geoel(iel,3) + geoel(iel,4)*geoel(iel,1) + geoel(iel,4)*geoel(iel,3)) + muE*(geoel(iel,2)*geoel(iel,4) + 2*geoel(iel,2)*geoel(iel,6) + geoel(iel,1)*geoel(iel,4) + geoel(iel,1)*geoel(iel,6)) ) * stiff->get(iel);
		K12(5,4) = 2/(3*geoel(iel,0))* ( lambdaE*(geoel(iel,6)*geoel(iel,1) + geoel(iel,6)*geoel(iel,3) + 2*geoel(iel,5)*geoel(iel,1) + geoel(iel,5)*geoel(iel,3)) + muE*(geoel(iel,3)*geoel(iel,4) + geoel(iel,3)*geoel(iel,6) + 2*geoel(iel,2)*geoel(iel,4) + geoel(iel,2)*geoel(iel,6)) ) * stiff->get(iel);

		return K12;
	}

	void assembleStiffnessMatrix()
	{
		amat::SpMatrix K11, K22, K12;
		K11.setup(m->gnpoin(),m->gnpoin()); K22.setup(m->gnpoin(), m->gnpoin()); K12.setup(m->gnpoin(), m->gnpoin());
		//K11.zeros(); K22.zeros(); K12.zeros();
		std::vector<int> ip(m->gnnode());
		for(int iel = 0; iel < m->gnelem(); iel++)
		{
			// get element stiffness matrices
			amat::Matrix<double> K11e = elementstiffnessK11(iel);
			amat::Matrix<double> K22e = elementstiffnessK22(iel);
			amat::Matrix<double> K12e = elementstiffnessK12(iel);

			for(int i = 0; i < m->gnnode(); i++)
				ip[i] = m->ginpoel(iel,i);

			double temp;
			for(int i = 0; i < m->gnnode(); i++)
				for(int j = 0; j < m->gnnode(); j++)
				{
					temp = K11.get(ip[i],ip[j]);
					K11.set(ip[i],ip[j], temp+K11e(i,j));
					temp = K22.get(ip[i],ip[j]);
					K22.set(ip[i],ip[j], temp+K22e(i,j));
					temp = K12.get(ip[i],ip[j]);
					K12.set(ip[i],ip[j], temp+K12e(i,j));
				}
		}
		std::cout << "LinElastP2: assembleStiffnessamat::Matrix(): Assembling final 2N by 2N matrix\n";
		/*K11.trim();
		K12.trim();
		K22.trim();*/

		amat::SpMatrix K21(m->gnpoin(), m->gnpoin());
		K21 = K12.transpose();

		K.combine_sparse_matrices(K11,K12,K21,K22);
	}

	void assembleLoadVector()
	{
		//amat::Matrix<double> load(m->gndim()*m->gnpoin(),1);
		f.zeros();
	}

	amat::SpMatrix stiffnessMatrix()
	{
		return K;
	}

	amat::Matrix<double> loadVector()
	{
		return f;
	}

	/** \brief Applies Dirichlet BC on high-order nodes of all boundary faces.
	 *
	 * bdata is a npoin-by-2 vector, containing x- and y-displacements of all mesh points to be imposed as Dirichlet BCs (zero for interior points). The first npoin entries are x-displacements.
	 * bglags is a npoin x 1 vector that contains 1 for boundary points and 0 for interior points, ie, it's a boundary boolean vector for nodes.
	 * IGNORE:extra is a npoin-by-1 vector storing a flag (0 or 1) for each mesh point; if extra(ipoin) == 1, then ipoin is to be kept fixed
	 */
	void dirichletBC_onAllBface(const amat::Matrix<double>& bdata, const amat::Matrix<int>& bflags)
	{
		double temp; int p;
		for(p = 0; p < m->gnpoin(); p++)
		{	
			if(bflags.get(p) == 1)
			{
				//x disp
				temp = K.get(p,p);
				K.set(p,p, temp*cbig);
				//if(ip == m->gnnofa()-1)		// when we encounter the high-order node, set its displacement
				f(p) = temp*cbig*bdata.get(p,0);
				
				//y disp
				temp = K.get(m->gnpoin()+p,m->gnpoin()+p);
				K.set(m->gnpoin()+p, m->gnpoin()+p, temp*cbig);
				//if(ip == m->gnnofa()-1)
				f(m->gnpoin()+p) = temp*cbig*bdata.get(p,1);
			}
			/*if(m->gbface(iface,m->gnnofa()) == 1)		// if non-homog BC is to be applied
			{
				f(m->gbface(iface,1)) = K(m->gbface(iface,1),m->gbface(iface,1)) * bdata(iface);
				f(m->gnpoin() + m->gbface(iface,1)) = K(m->gnpoin() + m->gbface(iface,1), m->gnpoin() + m->gbface(iface,1)) * bdata(m->gnface()+iface);
			}*/
		}
		/*for(int i = 0; i < m->gnpoin(); i++)
		{
			if(extra(i) == 1)
			{
				//x disp
				K(i,i) *= cbig;
				f(i) = 0;		// homogeneous BC by default
				//y disp
				K(m->gnpoin()+i, m->gnpoin()+i) *= cbig;
				f(m->gnpoin()+i) = 0;
			}
		}*/
	}
};

} // end namespace
