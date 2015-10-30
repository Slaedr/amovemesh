// Funtionality for solving equations of linear elasticity on a straight-edged quadratic mesh using P2 Lagrange finite elements
// (no non-homogeneous Neumann BCs)
// June 12, 2015
// Adtya Kashi

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
#define nthreads_elast 8
#endif
#endif

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif
#ifndef __ASPARSEMATRIX_H
#include <asparsematrix.hpp>
#endif
#ifndef __AMESH2D_H
#include "amesh2b.hpp"
#endif

#ifndef __ALINALG_H
typedef MatrixCOO SpMatrix;
#endif

using namespace std;
using namespace amat;
using namespace acfd;

namespace acfd {

class LinElastP1
{
	UTriMesh* m;
	int ngeoel;
	int ngeofa;
	Matrix<double> geoel;		// holds 2*area of element, and derivatives of barycentric coordinate functions lambdas
	Matrix<double> geofa;		// holds normals to and length of boundary faces
	SpMatrix K;					// global stiffness matrix
	Matrix<double> f;			// global load vector

	double muE;					// isotropic elasticity constants
	double lambdaE;
	double cbig;				// for Dirichlet BCs

public:
	LinElastP1(UTriMesh* mesh, double mu, double lambd)
	{
		m = mesh;
		ngeoel = 7;
		ngeofa = 3;
		muE = mu; lambdaE = lambd;
		geoel.setup(m->gnelem(), ngeoel);
		//geofa.setup(m->gnbpoin(), ngeofa);
		K.setup(m->gndim()*m->gnpoin(), m->gndim()*m->gnpoin());
		f.setup(m->gndim()*m->gnpoin(),1);
		cbig = 1e40;

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
		}

		/*for(int i = 0; i < m->gnface(); i++)
		{
			//n_x = y_2 - y_1
			geofa(i,0) = m->gcoords(m->gbface(i,2), 1) - m->gcoords(m->gbface(i,0), 1);
			//n_y = x_1 - x_2
			geofa(i,1) = m->gcoords(m->gbface(i,0), 0) - m->gcoords(m->gbface(i,2), 0);
			//l = sqrt(n_x^2 + n_y^2)
			geofa(i,2) = sqrt(geofa(i,0)*geofa(i,0) + geofa(i,1)*geofa(i,1));
			//geofa(i,2) = sqrt((m.gcoords(m.gbface(i,1)-1, 1) - m.gcoords(m.gbface(i,0)-1, 1))*(m.gcoords(m.gbface(i,1)-1, 1) - m.gcoords(m.gbface(i,0)-1, 1)) + (m.gcoords(m.gbface(i,0)-1, 0) - m.gcoords(m.gbface(i,1)-1, 0))*(m.gcoords(m.gbface(i,0)-1, 0) - m.gcoords(m.gbface(i,1)-1, 0)));
		}*/
		cout << "LinElastP1: Computed derivatives of basis functions, and normals to and lengths of boundary faces.\n";
	}

	Matrix<double> elementstiffnessK11(int iel)
	{
		Matrix<double> K11(3,3);
		double coeff = 2*muE+lambdaE;
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
			{
				if(i==j)
					K11(i,i) = (coeff*geoel(iel,i+1)*geoel(iel,i+1) + muE*geoel(iel,i+4)*geoel(iel,i+4))/(geoel(iel,0)*geoel(iel,0));	// (c*ai^2 + mu*bi^2) / D^2
				else
					K11(i,j) = (coeff*geoel(iel,i+1)*geoel(iel,j+1) + muE*geoel(iel,i+4)*geoel(iel,j+4))/(geoel(iel,0)*geoel(iel,0));		// (c*ai*aj + mu*bi*bj) / D^2
			}

		return K11;
	}

	Matrix<double> elementstiffnessK22(int iel)
	{
		Matrix<double> K22(6,6);
		double coeff = 2*muE+lambdaE;
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
			{
				if(i==j)
					K22(i,i) = (muE*geoel(iel,i+1)*geoel(iel,i+1) + coeff*geoel(iel,i+4)*geoel(iel,i+4))/(geoel(iel,0)*geoel(iel,0));
				else
					K22(i,j) = (muE*geoel(iel,i+1)*geoel(iel,j+1) + coeff*geoel(iel,i+4)*geoel(iel,j+4))/(geoel(iel,0)*geoel(iel,0));
			}

		return K22;
	}

	Matrix<double> elementstiffnessK12(int iel)
	{
		Matrix<double> K12(3,3);
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
			{
				if(i==j)
					K12(i,i) = (lambdaE*geoel(iel,i+1)*geoel(iel,i+4) + muE*geoel(iel,i+1)*geoel(iel,i+4))/(geoel(iel,0)*geoel(iel,0));
				else
					K12(i,j) = (lambdaE*geoel(iel,i+1)*geoel(iel,j+4) + muE*geoel(iel,i+4)*geoel(iel,j+1))/(geoel(iel,0)*geoel(iel,0));
			}

		return K12;
	}

	void assembleStiffnessMatrix()
	{
		SpMatrix K11(m->gnpoin(),m->gnpoin()), K22(m->gnpoin(),m->gnpoin()), K12(m->gnpoin(),m->gnpoin());

		UTriMesh* m = this->m;
		int iel;

		for(iel = 0; iel < m->gnelem(); iel++)
		{
			// get element stiffness matrices
			Matrix<double> K11e = elementstiffnessK11(iel);
			Matrix<double> K22e = elementstiffnessK22(iel);
			Matrix<double> K12e = elementstiffnessK12(iel);

			vector<int> ip(m->gnnode());
			for(int i = 0; i < m->gnnode(); i++)
				ip[i] = m->ginpoel(iel,i);

			double temp;

			for(int i = 0; i < m->gnnode(); i++)
				for(int j = 0; j < m->gnnode(); j++)
				{
					temp = K11.get(ip[i],ip[j]);
					K11.set(ip[i],ip[j], temp + K11e(i,j) );
					temp = K22.get(ip[i],ip[j]);
					K22.set(ip[i],ip[j], temp + K22e(i,j) );
					temp = K12.get(ip[i],ip[j]);
					K12.set(ip[i],ip[j], temp + K12e(i,j) );
				}
		}

		cout << "LinElastP1: assembleStiffnessMatrix(): Assembling final 2N by 2N matrix\n";
		// construct final 2N by 2N stiffness matrix
		int i;
		//SpMatrix* Kg = &(this->K);

		/*for(int i = 0; i < m->gnpoin(); i++)
			for(int j = 0; j < m->gnpoin(); j++)
			{
				K.set(i,j, K11.get(i,j));
				K.set(m->gnpoin()+i, m->gnpoin()+j, K22.get(i,j));
				K.set(i, m->gnpoin()+j, K12.get(i,j));
				K.set(m->gnpoin()+i, j, K12.get(j,i));
			}*/
		SpMatrix K21(m->gnpoin(), m->gnpoin());
		K21 = K12.transpose();

		/*ofstream ofile("checkmat.dat");
		for(int i = 0; i < m->gnpoin(); i++)
		{
			for(int j = 0; j < m->gnpoin(); j++)
				ofile << K11.get(i,j) << " ";
			ofile << endl;
		}
		ofile << endl << endl;
		for(int i = 0; i < m->gnpoin(); i++)
		{
			for(int j = 0; j < m->gnpoin(); j++)
				ofile << K12.get(i,j) << " ";
			ofile << endl;
		}
		ofile.close();*/

		K.combine_sparse_matrices(K11,K12,K21,K22);
		//K.print_data();
	}

	void assembleLoadVector()
	{
		//Matrix<double> load(m->gndim()*m->gnpoin(),1);
		f.zeros();
	}

	SpMatrix stiffnessMatrix()
	{
		return K;
	}

	Matrix<double> loadVector()
	{
		return f;
	}

	/*void dirichletBC_onAllBface(Matrix<double> bdata, Matrix<int> extra, int flag)
	// Applies Dirichlet BC on all boundary faces; homogeneous unless bface(iface,3) is flag for that face
	// bdata is a 2*nface-by-1 vector, containing x- and y-displacements of middle node for each bface. extra is a npoin-by-1 vector storing a flag (0 or 1) for each mesh point; if extra(ipoin) == 1, then ipoin is to be kept fixed
	{
		for(int iface = 0; iface < m->gnface(); iface++)
		{
			for(int ip = 0; ip < m->gnnofa()-1; ip++)
			{
				int p = m->gbface(iface,ip);
				//x disp
				K(p,p) *= cbig;
				f(p) = 0;		// homogeneous BC by default
				//y disp
				K(m->gnpoin()+p, m->gnpoin()+p) *= cbig;
				f(m->gnpoin()+p) = 0;
			}
			if(m->gbface(iface,3) == 1)		// if non-homog BC is to be applied
			{
				f(m->gbface(iface,1)) = K(m->gbface(iface,1),m->gbface(iface,1)) * bdata(iface);
				f(m->gnpoin() + m->gbface(iface,1)) = K(m->gnpoin() + m->gbface(iface,1), m->gnpoin() + m->gbface(iface,1)) * bdata(m->gnface()+iface);
			}
		}
		cout << "LinElastP2: dirichletBC_onAllBface(): Applying Dirichlet BCs\n";
		for(int i = 0; i < m->gnpoin(); i++)
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
		}
	}*/

	//void dirichletBC_onAllBface_2(Matrix<double> bdata, Matrix<int> extra)
	/* Applied dirichlet BCs for all boundary points. bdata is a npoin x 2 array containing x- and y-displacements for each point.
	 * The first nbpoin entries are considered to correspond to boundary points. */
	/*{
		cout << "LinElastP1: Assigning dirichlet BCs. No. of boun points " << m->gnbpoin() << "\n";
		double temp1, temp2;
		for(int i = 0; i < m->gnbpoin(); i++)
		{
			temp1 = K.get(i,i);
			temp2 = K.get(i+m->gnpoin(), i+m->gnpoin());
			f(i) = cbig*bdata(i,0)*temp1;
			f(i + m->gnpoin()) = cbig*bdata(i,1)*temp2;
			K.set(i,i, temp1*cbig);
			K.set(i + m->gnpoin(), i + m->gnpoin(), temp2 * cbig );
		}
	}*/

	void dirichletBC_onAllBface(Matrix<double> bdata, Matrix<int> extra, int flag_move, int flag_fixed)
	// For old domn mesh format: Applied dirichlet BCs for all boundary points. bdata is a npoin x 2 array containing x- and y-displacements for each point.
	// flag_move is the value of bface(i,3) for which inhomogeneous Dirichlet conditions must be imposed. flag_fixed are faces which are fixed.
	{
		cout << "LinElastP1: Assigning dirichlet BCs. No. of boun faces " << m->gnface() << "\n";
		double temp1, temp2;
		int ip[2];
		for(int i = 0; i < m->gnface(); i++)
		{
			ip[0] = m->gbface(i,0);
			ip[1] = m->gbface(i,1);
			for(int j = 0; j < 2; j++)
			{
				temp1 = K.get(ip[j],ip[j]);
				temp2 = K.get(ip[j]+m->gnpoin(), ip[j]+m->gnpoin());

				//if(m->gbface(i,3) == flag_move || m->gbface(i,3) == flag_fixed)
				//{
					K.set(ip[j],ip[j], temp1*cbig);
					K.set(ip[j] + m->gnpoin(), ip[j] + m->gnpoin(), temp2 * cbig );
				//}

				if(m->gbface(i,3)==flag_move)
				{
					f(ip[j]) = cbig*bdata(ip[j],0)*temp1;
					f(ip[j] + m->gnpoin()) = cbig*bdata(ip[j],1)*temp2;
				}
			}
		}
	}
};

} // end namespace
