#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#ifndef __AMESH_H
#include "amesh.hpp"
#endif

#define __ALAPLACELINFEM_H

using namespace std;
using namespace acfd;
using namespace amat;

namespace acfd{

const int ngeoel = 7;		///< no. of elements in each row of geoeol
const int ngeofa = 5;		///< no. of elements in each row of geofa
const int bcstate = 1;		///< value of bface(i,3) for which neumann BC is to be imposed on bface(i,:)

class LaplaceFEspaceP1
{
	/// Array to hold gradients of basis functions for each element of mesh. 
	/// Note: column 0 holds 'D', twice the measure of the element
	/// column-1 holds d(lambda1)/dx, column-2 holds d(lamda1)/dy, and so on.
	Matrix<double> geoel;
	
	/// Array to hold normals and edge-lengths of each boundary face (edge) of mesh.
	/// First column holds x-component of normal to face, 2nd col holds y component of normal and 
	/// 3rd column holds measure of the boundary facet (in case of 2d mesh, the length of the boundary edge)
	Matrix<double> geofa;		///< 4th and 5th columns of this store x and y components of gradient of unknown
	
	UTriMesh mesh;				// TODO: use pointer to mesh rather than a deep copy
	
public:
	
	LaplaceFEspaceP1(UTriMesh m, double dpx, double dpy)
	{
		mesh = m;
		geoel.setup(m.gnelem(), ngeoel, ROWMAJOR);
		geofa.setup(m.gnface(), ngeofa, ROWMAJOR);

		//stiffmat.setup(m.gnpoin(), m.gnpoin(), ROWMAJOR);
		//loadvec.setup(m.gnpoin(), 1, ROWMAJOR);
		
		for(int i = 0; i < m.gnelem(); i++)
		{
			// geoel(i,0) = D(i) = a1*b2 - a2*b1 :
			//geoel(i,0) = (m.gcoords(m.ginpoel(i,1)-1, 1) - m.gcoords(m.ginpoel(i,2)-1, 1))*(m.gcoords(m.ginpoel(i,0)-1, 0)-1 - m.gcoords(m.ginpoel(i,2)-1, 0)) - (m.gcoords(m.ginpoel(i,2)-1, 1) - m.gcoords(m.ginpoel(i,0)-1, 1))*(m.gcoords(m.ginpoel(i,2)-1, 0) - m.gcoords(m.ginpoel(i,1)-1, 0));
			
			geoel(i,0) = m.gcoords(m.ginpoel(i,0)-1,0)*(m.gcoords(m.ginpoel(i,1)-1,1) - m.gcoords(m.ginpoel(i,2)-1,1)) - m.gcoords(m.ginpoel(i,0)-1,1)*(m.gcoords(m.ginpoel(i,1)-1,0)-m.gcoords(m.ginpoel(i,2)-1,0)) + m.gcoords(m.ginpoel(i,1)-1,0)*m.gcoords(m.ginpoel(i,2)-1,1) - m.gcoords(m.ginpoel(i,2)-1,0)*m.gcoords(m.ginpoel(i,1)-1,1);
			double D = dabs(geoel(i,0));
		
			//a1/D = (y2 - y3)/D :
			geoel(i,1) = (m.gcoords(m.ginpoel(i,1)-1, 1) - m.gcoords(m.ginpoel(i,2)-1, 1)) / D;
			//b1/D = (x3 - x2)/D :
			geoel(i,2) = (m.gcoords(m.ginpoel(i,2)-1, 0) - m.gcoords(m.ginpoel(i,1)-1, 0)) / D;
			//a2/D = (y3 - y1)/D:
			geoel(i,3) = (m.gcoords(m.ginpoel(i,2)-1, 1) - m.gcoords(m.ginpoel(i,0)-1, 1)) / D;
			//b2/D = (x1 - x3)/D :
			geoel(i,4) = (m.gcoords(m.ginpoel(i,0)-1, 0) - m.gcoords(m.ginpoel(i,2)-1, 0)) / D;
			//a3/D = (y1 - y2)/D :
			geoel(i,5) = -geoel(i,1)-geoel(i,3);
			//b3/D = (x2 - x1)/D :
			geoel(i,6) = -geoel(i,2) - geoel(i,4);
		}
		
		for(int i = 0; i < m.gnface(); i++)
		{
			//n_x = y_2 - y_1
			geofa(i,0) = m.gcoords(m.gbface(i,1)-1, 1) - m.gcoords(m.gbface(i,0)-1, 1);
			//n_y = x_1 - x_2
			geofa(i,1) = m.gcoords(m.gbface(i,0)-1, 0) - m.gcoords(m.gbface(i,1)-1, 0);
			//l = sqrt(n_x^2 + n_y^2)
			geofa(i,2) = sqrt(geofa(i,0)*geofa(i,0) + geofa(i,1)*geofa(i,1));
			//geofa(i,2) = sqrt((m.gcoords(m.gbface(i,1)-1, 1) - m.gcoords(m.gbface(i,0)-1, 1))*(m.gcoords(m.gbface(i,1)-1, 1) - m.gcoords(m.gbface(i,0)-1, 1)) + (m.gcoords(m.gbface(i,0)-1, 0) - m.gcoords(m.gbface(i,1)-1, 0))*(m.gcoords(m.gbface(i,0)-1, 0) - m.gcoords(m.gbface(i,1)-1, 0)));
		
			// gradients of unknown phi:
			geofa(i,3) = dpx;
			geofa(i,4) = dpy;
		}
		cout << "LaplaceFEspaceP1: Computed derivatives of basis functions, and normals to and lengths of boundary faces.\n";
	}
	
	double ggeoel(int i, int j) { return geoel(i,j); }
	double ggeofa(int i, int j) { return geofa(i,j); }
	
	// Below, i and j are in {1,2,3}
	double elemstiffmat(int ielem, int i, int j)
	{
		// TODO: figure out a better way to calculate components of element stiffness matrix, perhaps by directly
		//		 using a_i and b_i instead of the derivatives of the lamdas (shape functions)
		double val = (geoel(ielem,(i+1)*2-1)*geoel(ielem,(j+1)*2-1) + geoel(ielem,(i+1)*2)*geoel(ielem,(j+1)*2))*abs(geoel(ielem,0))/2;
		return val;
	}

	// geofa(i,3) and geofa(i,4) are x- and y-components of the gradient of the unknown at each face
	double elemloadvec(int iface)
	{
		// return (dpx,dpy).(n_x, n_y) * length(face)/2 when the boundary state is bcstate (==1)
		//double val = (mesh.gbface(iface, 3) == bcstate) ? (geofa(iface,3)*geofa(iface,0) + geofa(iface,4)*geofa(iface,1)) * geofa(iface,2) / 2 : 0.0;
		double val = (mesh.gbface(iface, 3) == bcstate) ? (geofa(iface,3)*geofa(iface,0) + geofa(iface,4)*geofa(iface,1)) / 2 : 0.0;
		return val;
	}

	Matrix<double> stiffnessmatrix()
	{
		cout << "\nAssembling global stiffness matrix...";
		Matrix<double> stiffmat(mesh.gnpoin(), mesh.gnpoin(), ROWMAJOR);
		stiffmat.zeros();
		//stiffmat.setup(mesh.gnpoin(), mesh.gnpoin(), ROWMAJOR);
		for(int ielem = 0; ielem < mesh.gnelem(); ielem++)
		{
			vector<int> ip(mesh.gnnode());
			for(int i = 0; i < mesh.gnnode(); i++)
				ip[i] = mesh.ginpoel(ielem, i)-1;
			
			for(int  i = 0; i < mesh.gnnode(); i++)
				for(int j = 0; j < mesh.gnnode(); j++)
					stiffmat(ip[i],ip[j]) += elemstiffmat(ielem,i,j);
		}
		cout << "\nStiffness matrix assembled.";
		return stiffmat;
	}

	Matrix<double> loadvector()		// works only for 2d triangular elements
	//void loadvector()
	{
		Matrix<double> loadvec(mesh.gnpoin(),1,ROWMAJOR); 
		loadvec.zeros();
	
		for(int iface = 0; iface < mesh.gnface(); iface++)
		{
			vector<int> ip(2);
			for(int i = 0; i < 2; i++)
				ip[i] = mesh.gbface(iface,i)-1;
		
			for(int i = 0; i < 2; i++)
				loadvec(ip[i],0) += elemloadvec(iface);
		}
		cout << "\nLoad vector assembled.";
		return loadvec;
	}
};

Matrix<double> getvelocity(LaplaceFEspaceP1 p, UTriMesh m, Matrix<double> phi)
{
	Matrix<double> ucell(m.gnelem(),2,ROWMAJOR);
	
	for(int ie = 0; ie < m.gnelem(); ie++)
	{
		double sum0 = 0;
		for(int j = 0; j < m.gnnode(); j++)
			sum0 += phi(m.ginpoel(ie,j)-1,0) * p.ggeoel(ie, 2*(j+1)-1);
		ucell(ie,0) = sum0;
		sum0 = 0;
		for(int j = 0; j < m.gnnode(); j++)
			sum0 += phi(m.ginpoel(ie,j)-1,0) * p.ggeoel(ie, 2*(j+1));
		ucell(ie,1) = sum0;
	}
	
	// 0th col contains u_x, 1st col contains u_y, 2nd col contains total weight for that point
	Matrix<double> upoint(m.gnpoin(),3,ROWMAJOR);
	upoint.zeros();
	
	for(int i = 0; i < m.gnelem(); i++)
	{	
		vector<double> ip(m.gnnode());
		for(int j = 0; j < m.gnnode(); j++)
			ip[j] = m.ginpoel(i,j) - 1;
		
		
		for(int j = 0; j < m.gnnode(); j++)
		{
			upoint(ip[j],0) += ucell(i,0)*p.ggeoel(i,0)/2;
			upoint(ip[j],1) += ucell(i,1)*p.ggeoel(i,0)/2;
			upoint(ip[j],2) += abs(p.ggeoel(i,0))/2;				// weights - area
			//TODO: store areas in a separate matrix and return only x- and y-velocities
		}
	}
	
	// Now divide by weights
	for(int i = 0; i < m.gnpoin(); i++)
	{
		upoint(i,0) = upoint(i,0) / upoint(i,2);
		upoint(i,1) = upoint(i,1) / upoint(i,2);
	}
	return upoint;
}

} // end namespace
