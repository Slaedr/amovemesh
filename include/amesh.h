// Data structure and setup for 2D unstructured triangular mesh.
// Aditya Kashi
// Feb 5, 2015

// Feb 26, 2015: Adding mesh-movement function movemesh2d_bs().

#ifndef __AMATRIX2_H
#include "amatrix2.h"
#endif
#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif
#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif
#ifndef _GLIBCXX_STRING
#include <string>
#endif
#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifndef __ALINALG_H
#include "alinalg.h"
#endif
#ifndef __ADATASTRUCTURES_H
#include "adatastructures.h"
#endif

#define __AMESH_H

using namespace std;
using namespace amat;

namespace acfd {

// This function is a cyclic permutation of consecutive integers from 'start' to 'end' (inclusive). It returns the integer (between 'start' and 'end') that is 'off' integers away from 'n' in the cyclic order.
int perm(int start, int end, int n, int off)
{
	if(n > end) { cout << "Permutation point error!\n"; return 0; }
	
	CircList<int> list(start);
	for(int i = start+1; i <= end; i++)
		list.push(i);
	
	Node<int>* nn = list.find(n);
	Node<int>* cur = nn;
	for(int i = 0; i < off; i++)
		cur = cur->next;
	return cur->data;
}
	
const int n_extra_fields_in_bface = 2;

// Unstructured triangular mesh class
class UTriMesh
{
private:
	int npoin;
	int nelem;
	int nface;
	int ndim;
	int nnode;
	Matrix<double> coords;
	Matrix<int> inpoel;
	Matrix<int> bface;
	bool alloc_jacobians;
	public: Matrix<double> jacobians;
	
public:

	UTriMesh() {alloc_jacobians = false;}
	/* UTriMesh()
	{
		npoin = 0; nelem = 0; nface = 0; ndim = 0; nnode = 0;
		coords.setup(1,1,ROWMAJOR);
		inpoel.setup(1,1,ROWMAJOR);
		bface.setup(1,1,ROWMAJOR);
	} */
	
	UTriMesh(ifstream& infile)
	{
		alloc_jacobians = false;
		
		// Do file handling here to populate npoin and nelem
		cout << "\nUTriMesh: Reading mesh file...";
		char ch = '\0'; int dum = 0; double dummy;
		
		infile >> ch; infile >> ch;
		for(int i = 0; i < 4; i++)		//skip 4 lines
			do
				ch = infile.get();
			while(ch != '\n');
		infile >> ndim;
		infile >> nnode;
		infile >> ch;			//get the newline
		do
			ch = infile.get();
		while(ch != '\n');
		infile >> nelem; infile >> npoin; infile >> nface;
		infile >> dummy; 				// get time
		ch = infile.get();			// clear newline
		
		cout << "\nUTriMesh: Number of elements: " << nelem << ", number of points: " << npoin << ", number of nodes per element: " << nnode << endl;
		cout << "Number of boundary faces: " << nface << ", Number of dimensions: " << ndim;
		
		cout << "\nUTriMesh: Allocating coords..";
		coords.setup(npoin, ndim, ROWMAJOR);
		cout << "\nUTriMesh: Allocating inpoel..\n";
		inpoel.setup(nelem, nnode, ROWMAJOR);
		cout << "UTriMesh: Allocating bface...\n";
		bface.setup(nface, ndim + n_extra_fields_in_bface, ROWMAJOR);
		
		cout << "UTriMesh: Allocation done.";
		
		do
			ch = infile.get();
		while(ch != '\n');
			
		//now populate inpoel
		for(int i = 0; i < nelem; i++)
		{
			infile >> dum;
			for(int j = 0; j < nnode; j++)
				infile >> inpoel(i,j);
			do
				ch = infile.get();
			while(ch != '\n');
		}
		cout << "\nUTriMesh: Populated inpoel.";
		
		ch = infile.get();
		do
			ch = infile.get();
		while(ch != '\n');
		
		// populate coords
		for(int i = 0; i < npoin; i++)
		{
			infile >> dum;
			for(int j = 0; j < ndim; j++)
				infile >> coords(i,j);
		}
		cout << "\nUTriMesh: Populated coords.\n";
		
		ch = infile.get();
		for(int i = 0; i < npoin+2; i++)
		{
			do
				ch = infile.get();
			while(ch != '\n');
		}
		for(int i = 0; i < nface; i++)
		{
			infile >> dum;
			for(int j = 0; j < ndim + n_extra_fields_in_bface; j++)
			{
				infile >> bface(i,j);
			}
			do
				ch = infile.get();
			while(ch!='\n');
		}
		cout << "UTriMesh: Populated bface. Done reading mesh.\n";
	}
	
	UTriMesh(UTriMesh& other)
	{
		npoin = other.npoin;
		nelem = other.nelem;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
	}
	
	UTriMesh& operator=(UTriMesh& other)
	{
		npoin = other.npoin;
		nelem = other.nelem;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		return *this;
	}
	
	// Returns x or y coordinate (depending on dim) of node number pointno. Numbering of points begins
	// from 0. Numbering of dimensions begins from 0. 
	double gcoords(int pointno, int dim)
	{
		return coords.get(pointno,dim);
	}
	
	// Returns global node number of locnode th local node of element number elemno. Numberings for both
	// begins from 0
	int ginpoel(int elemno, int locnode)
	{
		return inpoel.get(elemno, locnode);
	}
	
	int gbface(int faceno, int val)
	{
		return bface.get(faceno, val);
	}
	
	int gnpoin() { return npoin; }
	int gnelem() { return nelem; }
	int gnface() { return nface; }
	int gnnode() { return nnode; }
	int gndim() { return ndim; }
	
	void compute_jacobians()
	{
		if (alloc_jacobians == false)
		{
			jacobians.setup(nelem, 1, ROWMAJOR);
			alloc_jacobians = true;
		}
		
		//TODO: implement computing of jacobian of each element here in a way that's compatible with CFD
		for(int i = 0; i < gnelem(); i++)
		{
			// geoel(i,0) = D(i) = a1*b2 - a2*b1 :
			
			jacobians(i,0) = gcoords(ginpoel(i,0)-1,0)*(gcoords(ginpoel(i,1)-1,1) - gcoords(ginpoel(i,2)-1,1)) - gcoords(ginpoel(i,0)-1,1)*(gcoords(ginpoel(i,1)-1,0)-gcoords(ginpoel(i,2)-1,0)) + gcoords(ginpoel(i,1)-1,0)*gcoords(ginpoel(i,2)-1,1) - gcoords(ginpoel(i,2)-1,0)*gcoords(ginpoel(i,1)-1,1);
		}
	}
	
	void detect_negative_jacobians(ofstream& out)
	{
		for(int i = 0; i < nelem; i++)
		{
			if(jacobians(i,0) < 1e-15) out << i << " " << jacobians(i,0) << '\n';
		}
	}
	
	//------------------------- Mesh movement functions -----------------------------------//
	
	// returns stiffness of the edge between global nodes i and j
	double k(int i, int j)
	{
		 return 1/sqrt((gcoords(i,0)-gcoords(j,0))*(gcoords(i,0)-gcoords(j,0)) + (gcoords(i,1)-gcoords(j,1))*(gcoords(i,1)-gcoords(j,1)));
	}
	
	// Moves the mesh according to prescribed boundary displacements xyb (where row 0 contains x-displacements and
	//	and row 1 contains y displacements of ALL nodes in the mesh. Note that the original mesh is overwritten.
	//	NOTE THAT xyb SHOULD BE A 2xN MATRIX, NOT AN Nx2 MATRIX
	/*void movemesh2d_bullshit(Matrix<double> xyb, double tol)
	{
		Matrix<double> d(2, gnpoin());
		d.zeros();
		
		Matrix<double> xnum(1, npoin); xnum.zeros();
		Matrix<double> denom(1, npoin); xdenom.zeros();
		Matrix<double> ynum(1, npoin); ynum.zeros();
		
		Matrix<double> err(1, gnpoin());
	
		int iter = 0;
		do
		{
			double num = 0;
			double denom = 0;
			for(int i = 0; i < nelem; i++)
			{
				//calculate d here
				double ip[nnode];
				for(int j = 0; j < nnode; j++)
					ip[j] = inpoel(i,j);
				
				// Note: Use of the perm function makes calculation of xnum and ynum uniform for any nnode. But it might be more efficient in terms of runtime to just write out everything instead of looping ove j and using perm.
				for(int j = 0; j < nnode; j++)
				{
					if(perm(0,nnode,j,1) >= nnode || perm(0,nnode,j,2) >= nnode)
						cout << "movemesh2d_bs: perm output is out of bounds!\n";
					xnum(0,ip[j]) += k(ip[j],ip[perm(0,nnode,j,1)])*xyb(0,ip[perm(0,nnode,j,1)]) + k(ip[0],ip[perm(0,nnode,j,2)])*xyb(0,ip[perm(0,nnode,j,2)]);
					ynum(0,ip[j]) += k(ip[j],ip[perm(0,nnode,j,1)])*xyb(1,ip[perm(0,nnode,j,1)]) + k(ip[0],ip[perm(0,nnode,j,2)])*xyb(1,ip[perm(0,nnode,j,2)]);
					denom(0,ip[j]) += k(ip[j], ip[perm(0,nnode,j,1)]) + k(ip[j], ip[perm(0,nnode,j,2)]);
				}
			}
			
			for(int i = 0; i < npoin; i++)
			{
				d(0,i) = xnum(0,i)/denom(0,i);
				d(1,i) = ynum(0,1)/denom(0,i);
				e(0,i) = sqrt((d(0,i)-xyb(0,i))*(d(0,i)-xyb(0,i))+(d(1,i)-xyb(1,i))*(d(1,i)-xyb(1,i)));
			}
		
			if(e.absmax() < tol) 
			{ 	cout << "Mesh-movement iterations converged.\n";
				break;
			}
		
			iter++;
			if(iter > 20)
			{	cout << "Mesh-movement did not converge.\n";
				break;
			}
			
			xyb = d;		// update displacements
		} while(true);
	
	}*/
	
	// Returns an array whose ith element is 1 if the ith node is a boundary node
	Matrix<int> bflags()
	{
		Matrix<int> flags(npoin, 1);
		flags.zeros();
		int ip[2];
		
		for(int b = 0; b < nface; b++)
		{
			ip[0] = bface(b,0)-1;
			ip[1] = bface(b,1)-1;
			flags(ip[0],0) = 1;
			flags(ip[1],0) = 1;
		}
		
		return flags;
	}
	
	// Moves the mesh according to prescribed boundary displacements xb and yb
	// blfag contains a 0-or-1 flag that indicates a boundary point  
	// Note that the original mesh is overwritten.
	// NOTE: MAKE SURE to update geoel, geofa and any other downstream data
	void movemesh2d_bs(Matrix<double> xb, Matrix<double> yb)
	{
		Matrix<double> A(npoin, npoin, ROWMAJOR);
		A.zeros();
		
		Matrix<int> bflag = bflags();
		
		cout << "movemesh2d_bs: Assembling stiffness matrix\n";
		
		for(int iel = 0; iel < nelem; iel++)
		{
			int ip[nnode];
			for(int i = 0; i < nnode; i++)
				ip[i] = inpoel(iel,i)-1;
			
			for(int i = 0; i < nnode; i++)
				for(int j = 0; j < nnode; j++)
				{
					if(i == j)
						if(bflag(ip[i],0)==0)
							A(ip[i],ip[i]) += k(ip[i], ip[perm(0,nnode-1,i,1)]);
						else
							A(ip[i],ip[i]) += k(ip[i], ip[perm(0,nnode-1,i,1)]) + k(ip[i], ip[perm(0,nnode-1,i,2)]);
						// the above is wrong for boundary nodes (some edges are counted twice), as elements are not cyclic around a boundary node.
						// However, since a Dirichlet BC is applied to each boundary node, it does not matter.
					else A(ip[i],ip[j]) -= k(ip[i],ip[j]);
				}
		}
		
		ofstream ofile("stiffmat.dat");
		A.fprint(ofile);
		
		cout << "movemesh2d_bs: Setting up boundary conditions\n";
		
		double cbig = 1e30;			// big number for Dirichlet BCs
		// apply Dirichlet BCs using cbig
		for(int i = 0; i < npoin; i++)
		{
			if(bflag(i,0) == 1)
			{
				xb(i,0) *= (cbig * A(i,i));
				yb(i,0) *= (cbig * A(i,i));
				A(i,i) *= cbig;
			}
		}
		
		Matrix<double> dx(npoin,1);
		Matrix<double> dy(npoin,2);
		
		Matrix<double> initvals(npoin,1);
		initvals.zeros();
		
		//solve for the displacements
		cout << "movemesh2d_bs: Solving mesh-movement equations\n";
// 		dx = cholesky(A,xb);
// 		dy = cholesky(A,yb);
		dx = gausselim(A,xb);
		dy = gausselim(A,yb);
// 		dx = pointjacobi(A, xb, initvals, 1e-6, 200, 'y');
// 		dy = pointjacobi(A, yb, initvals, 1e-6, 200, 'y');
		cout << "movemesh2d_bs: Mesh-movement equations solved.\n";
		
		//update coordinates
		for(int i = 0; i < npoin; i++)
		{
			coords(i,0) += dx(i,0);
			coords(i,1) += dy(i,0);
		}
	}
};

} // end namespace acfd
