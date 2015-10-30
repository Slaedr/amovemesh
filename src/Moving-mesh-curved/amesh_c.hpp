// Data structure and setup for 2D unstructured triangular mesh.
// Aditya Kashi
// Feb 5, 2015

// Feb 26, 2015: Adding mesh-movement function movemesh2d_bs().
// June 2, 2015: Adding refine mesh function

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
#ifdef _OPENMP
#ifndef OMP_H
#include <omp.h>
#endif
#endif

#ifndef __AMATRIX2_H
#include "amatrix2.hpp"
#endif
#ifndef __ALINALG_H
#include "alinalg.hpp"
#endif
#ifndef __ADATASTRUCTURES_H
#include "adatastructures.hpp"
#endif

#define __AMESH2D_H

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
	int nnode;		// number of nodes to an element
	int nfael;		// number of faces to an element (equal to number of edges to an element in 2D)
	int nnofa;		// number of node in a face -- needs to be generalized in case of general grids
	int naface;		// total number of (internal and boundary) faces
	int nbface;		// number of boundary faces as calculated by compute_face_data(), as opposed to nface which is read from file
	int naux;		// auxiliary nodes
	int nauxf;		// number of auxiliary nodes per face
	Matrix<double> coords;
	Matrix<int> inpoel;
	Matrix<int> bface;

	// The following 4 are for moving mesh
	Matrix<int> bflag;		// for the basicspring one - useless
	Matrix<int> bflag2;
	Matrix<double> cosa;
	Matrix<double> sina;
public:
	Matrix<int> esup;		// elements surrounding point
	Matrix<int> esup_p;		// stores positions corresponding to points in esup
	Matrix<int> psup;
	Matrix<int> psup_p;
	Matrix<int> esuel;

	Matrix<int> intfac;		// naface x 4 matrix: holds for each face: (1) index of left element, (2) index of right element, (3) index of starting point of face, and (4) index of ending point of face
	Matrix<double> gallfa;	//naface x 3 matrix: holds for each face: (1) x-component of normal (2) y-component of normal (3) Measure (length) of face. Note that normal points from left element towards right element (refer to intfac)

	bool alloc_jacobians;
	Matrix<double> jacobians;

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
		cout << "UTriMesh: Reading mesh file...";
		char ch = '\0'; int dum = 0; double dummy;

		//infile >> ch;
		infile >> dum;
		infile >> ch;
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

		//cout << "\nUTriMesh: Allocating coords..";
		coords.setup(npoin, ndim, ROWMAJOR);
		//cout << "\nUTriMesh: Allocating inpoel..\n";
		inpoel.setup(nelem, nnode, ROWMAJOR);
		//cout << "UTriMesh: Allocating bface...\n";
		bface.setup(nface, ndim + n_extra_fields_in_bface, ROWMAJOR);

		nfael = 3;	// number of faces per element
		nnofa = ndim;	// number of nodes per face

		//cout << "UTriMesh: Allocation done.";

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

		//Correct inpoel:
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
				inpoel(i,j)--;
		}

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
		//coords.mprint();

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
			if (i==nface-1) break;
			do
				ch = infile.get();
			while(ch!='\n');
		}
		cout << "UTriMesh: Populated bface. Done reading mesh.\n";
		//correct first 2 columns of bface
		for(int i = 0; i < nface; i++)
			for(int j = 0; j < 2; j++)
				bface(i,j)--;

		//------------- Calculate other topological properties -----------------
		cout << "UTriMesh: Calculating and storing topological information...\n";
		//1. Elements surrounding points
		cout << "UTriMesh: Elements surrounding points\n";
		esup_p.setup(npoin+1,1,ROWMAJOR);
		esup_p.zeros();

		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
			{
				esup_p(inpoel(i,j)+1,0) += 1;			// inpoel(i,j) + 1 : the + 1 is there because the storage corresponding to the first node begins at 0, not at 1
			}
		}
		// Now make the members of esup_p cumulative
		for(int i = 1; i < npoin+1; i++)
			esup_p(i,0) += esup_p(i-1,0);
		// Now populate esup
		esup.setup(esup_p(npoin,0),1,ROWMAJOR);
		esup.zeros();
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
			{
				int ipoin = inpoel(i,j);
				esup(esup_p(ipoin,0),0) = i;		// now put that element no. in the space pointed to by esup_p(ipoin)
				esup_p(ipoin,0) += 1;				// an element corresponding to ipoin has been found - increment esup_p for that point
			}
		}
		//But now esup_p holds increased values - each member increased by the number elements surrounding the corresponding point.
		// So correct this.
		for(int i = npoin; i >= 1; i--)
			esup_p(i,0) = esup_p(i-1,0);
		esup_p(0,0) = 0;
		// Elements surrounding points is now done.

		//2. Points surrounding points
		cout << "UTriMesh: Points surrounding points\n";
		psup_p.setup(npoin+1,1,ROWMAJOR);
		psup_p.zeros();
		psup_p(0,0) = 0;
		Matrix<int> lpoin(npoin,1);  // The ith member indicates the global point number of which the ith point is a surrounding point
		for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize this vector to -1
		int istor = 0;

		// first pass: calculate storage needed for psup
		for(int ip = 0; ip < npoin; ip++)
		{
			lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
			// Loop over elements surrounding this point
			for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
			{
				int ielem = esup(ie,0);		// element number
				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip)		// test of this point as already been counted as a surrounding point of ip
					{
						istor++;
						//psup(istor,0) = jpoin;	// ! can't do this yet - psup not allocated!
						lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
					}
				}
			}
			psup_p(ip+1,0) = istor;
		}

		psup.setup(istor,1,ROWMAJOR);
		//cout << "+++ " << istor << endl;

		//second pass: populate psup
		istor = 0;
		for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize lpoin to -1
		for(int ip = 0; ip < npoin; ip++)
		{
			lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
			// Loop over elements surrounding this point
			for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
			{
				int ielem = esup(ie,0);		// element number
				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip)		// test of this point as already been counted as a surrounding point of ip
					{
						psup(istor,0) = jpoin;
						istor++;
						lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
					}
				}
			}
			//psup_p(ip+1,0) = istor;
		}
		//Points surrounding points is now done.

		// 3. Elements surrounding elements
		cout << "UTriMesh: Elements surrounding elements...\n";

		esuel.setup(nelem, nfael, ROWMAJOR);
		for(int ii = 0; ii < nelem; ii++)
			for(int jj = 0; jj < nfael; jj++)
				esuel(ii,jj) = -1;
		Matrix<int> lpofa(nfael, nnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in {0,1}, i in {0,1,2})
		lpofa(0,0) = 1; lpofa(0,1) = 2;
		lpofa(1,0) = 2; lpofa(1,1) = 0;
		lpofa(2,0) = 0; lpofa(2,1) = 1;
		Matrix<int> lhelp(nnofa,1);
		lhelp.zeros();
		lpoin.zeros();

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int ifael = 0; ifael < nfael; ifael++)
			{
				for(int i = 0; i < nnofa; i++)
				{
					lhelp(i,0) = inpoel(ielem, lpofa(ifael,i));	// lhelp stores global node nos. of current face of current element
					lpoin(lhelp(i,0)) = 1;
				}
				int ipoin = lhelp(0);
				for(int istor = esup_p(ipoin); istor < esup_p(ipoin+1); istor++)
				{
					int jelem = esup(istor);
					if(jelem != ielem)
					{
						for(int jfael = 0; jfael < nfael; jfael++)
						{
							//Assume that no. of nodes in face ifael is same as that in face jfael
							int icoun = 0;
							for(int jnofa = 0; jnofa < nnofa; jnofa++)
							{
								int jpoin = inpoel(jelem, lpofa(jfael,jnofa));
								if(lpoin(jpoin)==1) icoun++;
							}
							if(icoun == nnofa)		// nnofa is 2
							{
								esuel(ielem,ifael) = jelem;
								esuel(jelem,jfael) = ielem;
							}
						}
					}
				}
				for(int i = 0; i < nnofa; i++)
					lpoin(lhelp(i)) = 0;
			}
		}
		//cout << "UTriMesh: Elements surrounding elements done.\n";

		//for moving mesh stuff:
		bflag = bflags();
		bflag2 = bflags2();
		cout << "UTriMesh: Done.\n";
	}

	void compute_topological()
	{
		//------------- Calculate other topological properties -----------------
		cout << "UTriMesh: Calculating and storing topological information...\n";
		//1. Elements surrounding points
		cout << "UTriMesh: Elements surrounding points\n";
		esup_p.setup(npoin+1,1,ROWMAJOR);
		esup_p.zeros();

		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
			{
				esup_p(inpoel(i,j)+1,0) += 1;			// inpoel(i,j) + 1 : the + 1 is there because the storage corresponding to the first node begins at 0, not at 1
			}
		}
		// Now make the members of esup_p cumulative
		for(int i = 1; i < npoin+1; i++)
			esup_p(i,0) += esup_p(i-1,0);
		// Now populate esup
		esup.setup(esup_p(npoin,0),1,ROWMAJOR);
		esup.zeros();
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
			{
				int ipoin = inpoel(i,j);
				esup(esup_p(ipoin,0),0) = i;		// now put that element no. in the space pointed to by esup_p(ipoin)
				esup_p(ipoin,0) += 1;				// an element corresponding to ipoin has been found - increment esup_p for that point
			}
		}
		//But now esup_p holds increased values - each member increased by the number elements surrounding the corresponding point.
		// So correct this.
		for(int i = npoin; i >= 1; i--)
			esup_p(i,0) = esup_p(i-1,0);
		esup_p(0,0) = 0;
		// Elements surrounding points is now done.

		//2. Points surrounding points
		cout << "UTriMesh: Points surrounding points\n";
		psup_p.setup(npoin+1,1,ROWMAJOR);
		psup_p.zeros();
		psup_p(0,0) = 0;
		Matrix<int> lpoin(npoin,1);  // The ith member indicates the global point number of which the ith point is a surrounding point
		for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize this vector to -1
		int istor = 0;

		// first pass: calculate storage needed for psup
		for(int ip = 0; ip < npoin; ip++)
		{
			lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
			// Loop over elements surrounding this point
			for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
			{
				int ielem = esup(ie,0);		// element number
				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip)		// test of this point as already been counted as a surrounding point of ip
					{
						istor++;
						//psup(istor,0) = jpoin;	// ! can't do this yet - psup not allocated!
						lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
					}
				}
			}
			psup_p(ip+1,0) = istor;
		}

		psup.setup(istor,1,ROWMAJOR);
		//cout << "+++ " << istor << endl;

		//second pass: populate psup
		istor = 0;
		for(int i = 0; i < npoin; i++) lpoin(i,0) = -1;	// initialize lpoin to -1
		for(int ip = 0; ip < npoin; ip++)
		{
			lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
			// Loop over elements surrounding this point
			for(int ie = esup_p(ip,0); ie <= esup_p(ip+1,0)-1; ie++)
			{
				int ielem = esup(ie,0);		// element number
				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip)		// test of this point as already been counted as a surrounding point of ip
					{
						psup(istor,0) = jpoin;
						istor++;
						lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
					}
				}
			}
			//psup_p(ip+1,0) = istor;
		}
		//Points surrounding points is now done.

		// 3. Elements surrounding elements
		cout << "UTriMesh: Elements surrounding elements...\n";

		esuel.setup(nelem, nfael, ROWMAJOR);
		for(int ii = 0; ii < nelem; ii++)
			for(int jj = 0; jj < nfael; jj++)
				esuel(ii,jj) = -1;
		Matrix<int> lpofa(nfael, nnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in {0,1}, i in {0,1,2})
		lpofa(0,0) = 1; lpofa(0,1) = 2;
		lpofa(1,0) = 2; lpofa(1,1) = 0;
		lpofa(2,0) = 0; lpofa(2,1) = 1;
		Matrix<int> lhelp(nnofa,1);
		lhelp.zeros();
		lpoin.zeros();

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int ifael = 0; ifael < nfael; ifael++)
			{
				for(int i = 0; i < nnofa; i++)
				{
					lhelp(i,0) = inpoel(ielem, lpofa(ifael,i));	// lhelp stores global node nos. of current face of current element
					lpoin(lhelp(i,0)) = 1;
				}
				int ipoin = lhelp(0);
				for(int istor = esup_p(ipoin); istor < esup_p(ipoin+1); istor++)
				{
					int jelem = esup(istor);
					if(jelem != ielem)
					{
						for(int jfael = 0; jfael < nfael; jfael++)
						{
							//Assume that no. of nodes in face ifael is same as that in face jfael
							int icoun = 0;
							for(int jnofa = 0; jnofa < nnofa; jnofa++)
							{
								int jpoin = inpoel(jelem, lpofa(jfael,jnofa));
								if(lpoin(jpoin)==1) icoun++;
							}
							if(icoun == nnofa)		// nnofa is 2
							{
								esuel(ielem,ifael) = jelem;
								esuel(jelem,jfael) = ielem;
							}
						}
					}
				}
				for(int i = 0; i < nnofa; i++)
					lpoin(lhelp(i)) = 0;
			}
		}
		//cout << "UTriMesh: Elements surrounding elements done.\n";

	}

	UTriMesh(UTriMesh& other)
	{
		/*npoin = other.npoin;
		nelem = other.nelem;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface; */

		npoin = other.npoin;
		nelem = other.nelem;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		naface = other.naface;
		nbface = other.nbface;
		nfael = other.nfael;
		nnofa = other.nnofa;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		esup = other.esup;
		esup_p = other.esup_p;
		psup = other.psup;
		psup_p = other.psup_p;
		esuel = other.esuel;
		intfac = other.intfac;
		gallfa = other.gallfa;
		alloc_jacobians = other.alloc_jacobians;
		jacobians = other.jacobians;
	}

	UTriMesh& operator=(UTriMesh& other)
	{
		npoin = other.npoin;
		nelem = other.nelem;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		naface = other.naface;
		nbface = other.nbface;
		nfael = other.nfael;
		nnofa = other.nnofa;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		esup = other.esup;
		esup_p = other.esup_p;
		psup = other.psup;
		psup_p = other.psup_p;
		esuel = other.esuel;
		intfac = other.intfac;
		gallfa = other.gallfa;
		alloc_jacobians = other.alloc_jacobians;
		jacobians = other.jacobians;
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

	int gesup(int i) { return esup.get(i); }
	int gesup_p(int i) { return esup_p.get(i); }
	int gpsup(int i) { return psup.get(i); }
	int gpsup_p(int i) { return psup_p.get(i); }
	int gesuel(int ielem, int jnode) { return esuel.get(ielem, jnode); }		//returns element number at face opposite to node number jnode
	double gjacobians(int ielem) { return jacobians(ielem,0); }

	int gnpoin() { return npoin; }
	int gnelem() { return nelem; }
	int gnface() { return nface; }
	int gnbface() { return nbface; }
	int gnnode() { return nnode; }
	int gndim() { return ndim; }
	int gnaface() {return naface; }
	int gnfael() { return nfael; }
	int gnnofa() { return nnofa; }
	int gnaux() { return naux; }

	void compute_jacobians()
	{
		if (alloc_jacobians == false)
		{
			jacobians.setup(nelem, 1, ROWMAJOR);
			alloc_jacobians = true;
		}

		for(int i = 0; i < gnelem(); i++)
		{
			// geoel(i,0) = D(i) = a1*b2 - a2*b1 :

			jacobians(i,0) = gcoords(ginpoel(i,0),0)*(gcoords(ginpoel(i,1),1) - gcoords(ginpoel(i,2),1)) - gcoords(ginpoel(i,0),1)*(gcoords(ginpoel(i,1),0)-gcoords(ginpoel(i,2),0)) + gcoords(ginpoel(i,1),0)*gcoords(ginpoel(i,2),1) - gcoords(ginpoel(i,2),0)*gcoords(ginpoel(i,1),1);
		}
	}

	void detect_negative_jacobians(ofstream& out)
	{
		for(int i = 0; i < nelem; i++)
		{
			if(jacobians(i,0) <= 1e-15) out << i << " " << jacobians(i,0) << '\n';
		}
	}

	void compute_face_data()
	/* Computes, for each face, the elements on either side, the starting node and the ending node of the face. This is stored in intfac. Also computes unit normals to, and lengths of, each face as well as boundary flags of boundary faces, in gallfa.
	NOTE: After this function, esuel holds (nelem + face no.) for each ghost cell, instead of -1 as before.*/
	{
		nbface = naface = 0;
		// first run: calculate nbface
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				//int in1 = perm(0,nnode-1,in,1);
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je == -1)
				{
					//esuel(ie,in) = nelem+nbface;
					nbface++;
				}
			}
		}
		cout << "UTriMesh: compute_face_data(): Number of boundary faces = " << nbface << endl;
		// calculate number of internal faces
		naface = nbface;
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				//int in1 = perm(0,nnode-1,in,1);
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je > ie && je < nelem) naface++;
			}
		}
		cout << "UTriMesh: compute_face_data(): Number of all faces = " << naface << endl;

		//allocate intfac
		intfac.setup(naface,4,ROWMAJOR);

		//reset face totals
		nbface = naface = 0;

		//second run: populate intfac
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				int in1 = perm(0,nnode-1,in,1);
				int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je == -1)
				{
					esuel(ie,in) = nelem+nbface;
					intfac(nbface,0) = ie;
					intfac(nbface,1) = nelem+nbface;
					intfac(nbface,2) = inpoel(ie,in1);
					intfac(nbface,3) = inpoel(ie,in2);

					nbface++;
				}
			}
		}
		naface = nbface;
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				int in1 = perm(0,nnode-1,in,1);
				int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je > ie && je < nelem)
				{
					intfac(naface,0) = ie;
					intfac(naface,1) = je;
					intfac(naface,2) = inpoel(ie,in1);
					intfac(naface,3) = inpoel(ie,in2);
					naface++;
				}
			}
		}

		//Now compute normals and lengths
		gallfa.setup(naface, 3+n_extra_fields_in_bface, ROWMAJOR);
		for(int i = 0; i < naface; i++)
		{
			gallfa(i,0) = coords(intfac(i,3),1) - coords(intfac(i,2),1);
			gallfa(i,1) = -1.0*(coords(intfac(i,3),0) - coords(intfac(i,2),0));
			gallfa(i,2) = sqrt(pow(gallfa(i,0),2) + pow(gallfa(i,1),2));
			//Normalize the normal vector components
			gallfa(i,0) /= gallfa(i,2);
			gallfa(i,1) /= gallfa(i,2);
		}

		//Populate boundary flags in gallfa
		cout << "UTriMesh: compute_face_data(): Storing boundary flags in gallfa..\n";
		for(int ied = 0; ied < nbface; ied++)
		{
			int p1 = intfac(ied,2);
			int p2 = intfac(ied,3);
			// Assumption: order of nodes of boundary faces is such that normal points outside, when normal is calculated as
			//nx = y2 - y1, ny = -(x2-x1).
			if(nbface != nface) { cout << "UTriMesh: compute_face_data(): Calculation of number of boundary faces is wrong!\n"; break; }
			//cout << "UTriMesh: compute_face_data(): Into loop....\n";
			for(int i = 0; i < nface; i++)
			{
				if(bface.get(i,0) == p1 || bface.get(i,1) == p1)
				{
					if(bface.get(i,1) == p2 || bface.get(i,0) == p2)
					{
						gallfa(ied,3) = bface(i,2);
						gallfa(ied,4) = bface(i,3);
					}
				}
			}
		}

		cout << "UTriMesh: compute_face_data(): Done.\n";
	}

	void compute_lengths_and_normals()		// Also happening in compute_face_data
	{
		for(int i = 0; i < naface; i++)
		{
			gallfa(i,0) = coords(intfac(i,3),1) - coords(intfac(i,2),1);
			gallfa(i,1) = -1.0*(coords(intfac(i,3),0) - coords(intfac(i,2),0));
			gallfa(i,2) = sqrt(pow(gallfa(i,0),2) + pow(gallfa(i,1),2));
			//Normalize the normal vector components
			gallfa(i,0) /= gallfa(i,2);
			gallfa(i,1) /= gallfa(i,2);
		}
	}

	void allocate_edge_angles()		// call only after calling compute_face_data()
	{
		cosa.setup(naface,1);
		sina.setup(naface,1);
	}

	void compute_edge_angles()
	{
		cout << "UTriMesh: Computing edge angles...\n";
		// calculate cos(a) and sin(a) where 'a' is angle made by edge (face in 2D) with the x axis
		for(int ied = 0; ied < naface; ied++)
		{
			cosa(ied) = -gallfa(ied,1)/gallfa(ied,2);
			sina(ied) = gallfa(ied,0)/gallfa(ied,2);
		}
	}

	int gintfac(int face, int i) { return intfac.get(face,i); }
	double ggallfa(int elem, int i) {return gallfa.get(elem,i); }
	void printmeshstats()
	{
		cout << "UTriMesh: Stats: npoin = " << npoin << ", nelem = " << nelem << ", nface = " << nface << ", nbface = " << nbface << endl;
		cout << "UTriMesh: Stats: Dimensions of \n\tcoords " << coords.rows() << "x" << coords.cols() << "\n\tinpoel " << inpoel.rows() << "x" << inpoel.cols() << "\n\tbface " << bface.rows() << "x" << bface.cols() << endl;
	}

	void refine_mesh(int fac=1)		// fac is the number of extra points per face - currently works only for fac=1
	{
		compute_face_data();
		//cout << "UTriMesh: refine_mesh(): nbface = " << nbface << ", nface = " << nface << ", ndim = " << ndim << ", nnofa = " << nnofa << endl;

		int elemfac = 4;						// NOTE: This needs to change for higher refinements. 4 corresponds to fac=1
		int facefac = fac+1;					// factor by which naface will increase
		nauxf = fac;
		naux = fac*naface;						// total number of aux points is fac times number of faces
		int nauxel = nnode + nfael*nauxf;		// number of aux nodes per element
		int npoin2 = npoin + naux;
		Matrix<int> inpoel2(nelem*elemfac, nnode);
		Matrix<double> coords2(npoin+naux,ndim);
		Matrix<int> bface2(nbface*facefac, nnofa + n_extra_fields_in_bface);
		//cout << "UTriMesh: refine_mesh(): bface2 is " << bface2.rows() << "x" << bface2.cols() << endl;
		Matrix<double> apel(nelem,nfael);		// holds global point number of aux node on each face of each element

		// the first npoin nodes are the old ones
		for(int i = 0; i < npoin; i++)
		{
			for(int j = 0; j < ndim; j++)
				coords2(i,j) = coords(i,j);
		}

		cout << "UTriMesh: refine_mesh(): Creating auxiliary nodes on internal faces\n";
		int ip1, ip2, ielem, jelem, ep1, ep2, locside;
		for(int ied = nbface; ied < naface; ied++)
		{
			ip1 = gintfac(ied,2);
			ip2 = gintfac(ied,3);
			ielem = gintfac(ied,0);
			jelem = gintfac(ied,1);

			for(int i = 0; i < ndim; i++)
				coords2(npoin + ied,i) = (coords(ip1,i) + coords(ip2,i))/2.0;

			for(int in=0; in<nnode; in++)		// find local node nos. in ielem whose global nos. are ip1 and ip2
			{
				if(ip1 == ginpoel(ielem,in)) ep1 = in;
				if(ip2 == ginpoel(ielem,in)) ep2 = in;
			}
			//cout << "UTriMesh: refine_mesh(): ep1 = " << ep1 << ", ep2 = " << ep2 << endl;
			// find the local node no. of remaining vertex of ielem
			if((ep1 == 0 && ep2 == 1) || (ep1==1 && ep2==0)) locside = 2;
			if((ep1 == 1 && ep2 == 2) || (ep1==2 && ep2==1)) locside = 0;
			if((ep1 == 2 && ep2 == 0) || (ep1==0 && ep2==2)) locside = 1;
			apel(ielem,locside) = npoin + ied;		// store global node no. of the aux point just created; locside is the local node number opposite to the face on which the aux point is created - so locside is like the local face number of the face on which the aux point was created

			for(int in=0; in<nnode; in++)
			{
				if(ip1 == ginpoel(jelem,in)) ep1 = in;
				if(ip2 == ginpoel(jelem,in)) ep2 = in;
			}
			//cout << "UTriMesh: refine_mesh(): ep1 = " << ep1 << ", ep2 = " << ep2 << endl;
			if((ep1 == 0 && ep2 == 1) || (ep1==1 && ep2==0)) locside = 2;
			if((ep1 == 1 && ep2 == 2) || (ep1==2 && ep2==1)) locside = 0;
			if((ep1 == 2 && ep2 == 0) || (ep1==0 && ep2==2)) locside = 1;
			apel(jelem,locside) = npoin + ied;
		}
		cout << "UTriMesh: refine_mesh(): Creating auxiliary nodes on boundary faces and splitting the boundary faces\n";
		for(int ied = 0; ied < nbface; ied++)
		{
			ip1 = gintfac(ied,2);
			ip2 = gintfac(ied,3);
			ielem = gintfac(ied,0);
			for(int i = 0; i < ndim; i++)
				coords2(npoin + ied,i) = (coords(ip1,i) + coords(ip2,i))/2.0;

			for(int in=0; in<nnode; in++)
			{
				if(ip1 == ginpoel(ielem,in)) ep1 = in;
				if(ip2 == ginpoel(ielem,in)) ep2 = in;
			}
			//cout << "UTriMesh: refine_mesh(): ep1 = " << ep1 << ", ep2 = " << ep2 << endl;
			if((ep1 == 0 && ep2 == 1) || (ep1==1 && ep2==0)) locside = 2;
			if((ep1 == 1 && ep2 == 2) || (ep1==2 && ep2==1)) locside = 0;
			if((ep1 == 2 && ep2 == 0) || (ep1==0 && ep2==2)) locside = 1;
			apel(ielem,locside) = npoin + ied;

			//need to update bface2 as well. First find bface number corresponding to ied
			int faceno = -1;
			for(int ifa = 0; ifa < nbface; ifa++)
				if((bface(ifa,0)==ip1 && bface(ifa,1)==ip2) || (bface(ifa,0)==ip2 && bface(ifa,1)==ip1)) faceno = ifa;
			if(faceno == -1) cout << "UTriMesh: refine_mesh(): !! Error !!\n";

			bface2(nbface+faceno,0) = npoin+ied;						// the second half of the split face is assigned a new bface number: nbface+faceno
			bface2(nbface+faceno,1) = ip2;
			for(int iv = 2; iv < 2+n_extra_fields_in_bface; iv++)
				bface2(nbface+faceno,iv) = bface(faceno,iv);

			bface2(faceno,0) = ip1;									// the first half of the split face is assigned the old bface number
			bface2(faceno,1) = npoin+ied;
			for(int iv = 2; iv < 2+n_extra_fields_in_bface; iv++)
				bface2(faceno,iv) = bface(faceno,iv);
		}

		//Next, update inpoel
		cout << "UTriMesh: refine_mesh(): Updating inpoel\n";
		for(int iel = 0; iel < nelem; iel++)
		{
			// element number of central element equals that of the parent element
			for(int ifa = 0; ifa < nnode; ifa++)
				inpoel2(iel,ifa) = apel(iel,ifa);
			// the other 3 new elements are given element numbers after nelem
			inpoel2(nelem+iel*3, 0) = inpoel(iel,0);
			inpoel2(nelem+iel*3, 1) = apel(iel,2);
			inpoel2(nelem+iel*3, 2) = apel(iel,1);

			inpoel2(nelem+iel*3+1, 0) = inpoel(iel,1);
			inpoel2(nelem+iel*3+1, 1) = apel(iel,0);
			inpoel2(nelem+iel*3+1, 2) = apel(iel,2);

			inpoel2(nelem+iel*3+2, 0) = inpoel(iel,2);
			inpoel2(nelem+iel*3+2, 1) = apel(iel,1);
			inpoel2(nelem+iel*3+2, 2) = apel(iel,0);
		}

		// replace this mesh by refined mesh
		cout << "UTriMesh: refine_mesh(): Computation done, replacing mesh by refined mesh\n";
		coords = coords2;
		inpoel = inpoel2;
		bface = bface2;
		npoin = npoin2;
		nelem = nelem*elemfac;
		nface = nface*facefac;
		cout << "UTriMesh: refine_mesh(): Re-computing topological structures like esup and intfac\n";
		compute_topological();
		compute_face_data();
		cout << "UTriMesh: refine_mesh(): Done.\n";
		//inpoel.mprint();
		//printmeshstats();
	}



	//------------------------- Mesh movement functions -----------------------------------//

	// returns stiffness of the edge between global nodes i and j
	double k(int i, int j)
	{
		 return 1/sqrt((gcoords(i,0)-gcoords(j,0))*(gcoords(i,0)-gcoords(j,0)) + (gcoords(i,1)-gcoords(j,1))*(gcoords(i,1)-gcoords(j,1)));
	}

	// Returns an array whose ith element is 1 if the ith node is a boundary node
	Matrix<int> bflags()
	{
		Matrix<int> flags(npoin, 1);
		flags.zeros();
		int ip[2];

		for(int b = 0; b < nface; b++)
		{
			ip[0] = bface(b,0);
			ip[1] = bface(b,1);
			flags(ip[0],0) = 1;
			flags(ip[1],0) = 1;
		}

		return flags;
	}

	Matrix<int> bflags2()
	{
		Matrix<int> flags(2*npoin, 1);
		flags.zeros();
		int ip[2];

		for(int b = 0; b < nface; b++)
		{
			ip[0] = bface(b,0);
			ip[1] = bface(b,1);
			flags(2*ip[0],0) = 1;
			flags(2*ip[0]+1,0) = 1;
			flags(2*ip[1],0) = 1;
			flags(2*ip[1]+1,0) = 1;
		}

		return flags;
	}

	// Moves the mesh according to prescribed boundary displacements xb and yb
	// blfag contains a 0-or-1 flag that indicates a boundary point
	// Note that the original mesh is overwritten.
	// NOTE: MAKE SURE to update geoel, gallfa and any other downstream data

	void movemesh_basicspring(Matrix<double> xb, Matrix<double> yb, Matrix<int> extraflags, string solver, double tol=1e-6, int maxiter=1000)
	{
		Matrix<double> A(npoin, npoin, ROWMAJOR);
		A.zeros();
		bflag = bflags();

		cout << "movemesh_basicspring: Assembling stiffness matrix\n";
		// Essentially, we're solving a homogeneous discrete elliptic equation with Dirichlet BCs

		for(int ip = 0; ip < npoin; ip++)
		{
			for(int i = psup_p(ip,0); i <= psup_p(ip+1,0)-1; i++)
			{
				int ipoin = psup(i,0);
				A(ip,ip) += k(ip,ipoin);
				A(ip,ipoin) -= k(ip,ipoin);
			}
		}

		//ofstream ofile("stiffmat.dat");
		//A.fprint(ofile);
		//ofile.close();

		cout << "movemesh2d_bs: Setting up boundary conditions\n";

		double cbig = 1e30;			// big number for Dirichlet BCs
		// apply Dirichlet BCs using cbig
		for(int i = 0; i < npoin; i++)
		{
			if(bflag(i,0) == 1 || extraflags(i,0)==1)		// extraflags is just like bflags, except that 1-entries need not be boundary nodes
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
		cout << "movemesh_basicspring: Solving mesh-movement equations\n";
		if(solver == "cholesky")
		{
			dx = cholesky(A,xb);
			dy = cholesky(A,yb);
		}
		else if (solver == "gausselim")
		{
			dx = gausselim(A,xb);
			dy = gausselim(A,yb);
		}
		else if (solver == "pointjacobi")
		{
			dx = pointjacobi(A, xb, initvals, tol, maxiter, 'n');
			dy = pointjacobi(A, yb, initvals, tol, maxiter, 'n');
		}
		else if (solver == "gaussseidel")
		{
			dx = gaussseidel(A, xb, initvals, tol, maxiter, 'n');
			dy = gaussseidel(A, yb, initvals, tol, maxiter, 'n');
		}
		else cout << "movemesh_basicspring: Solver " << solver << " not found.\n";
		cout << "movemesh_basicspring: Mesh-movement equations solved.\n";

		//update coordinates
		for(int i = 0; i < npoin; i++)
		{
			coords(i,0) += dx(i,0);
			coords(i,1) += dy(i,0);
		}

		//update jacobians
		compute_jacobians();

		// Also need to update length of faces, normals to faces in gallfa
		compute_lengths_and_normals();
	}

	void movemesh_lineal(Matrix<double>* xb, Matrix<double>* yb, string solver, double tol=1e-6, int maxiter=1000)
	{
		cout << "UTriMesh: movemesh_farhat(): Starting\n";
		Matrix<double> Kg(2*npoin,2*npoin);		// global stiffness matrix
		Kg.zeros();
		//Matrix<double> kli(4,4);				// lineal stiffness matrix for an edge (face in 2D)
		double k11, k12, k22;
		double len;
		int ipx, ipy, jpx, jpy;

		for(int ied = 0; ied < naface; ied++)
		{
			//gather
			len = gallfa(ied,2);
			k11 = cosa(ied)*cosa(ied)/len;
			k12 = sina(ied)*cosa(ied)/len;
			k22 = sina(ied)*sina(ied)/len;
			/*
			kli(0,0) = k11; kli(0,1) = k12; kli(0,2) = -k11; kli(0,3) = -k12;
			kli(1,0) = k12; kli(1,1) = k22; kli(1,2) = -k12; kli(1,3) = -k22;
			kli(2,0) = -k11; kli(2,1) = -k12; kli(2,2) = k11; kli(2,3) = k12;
			kli(3,0) = -k12; kli(3,1) = -k22; kli(3,2) = k12; kli(3,3) = k22;
			*/
			ipx = 2*intfac(ied,2);
			ipy = 2*intfac(ied,2)+1;
			jpx = 2*intfac(ied,3);
			jpy = 2*intfac(ied,3)+1;

			//scatter
			Kg(ipx,ipx) += k11;			// contribution by x-displacement of ip node to x-force at ip node
			Kg(ipx,ipy) += k12;			// contribution by y-displacement of ip node to x-force at ip node
			Kg(ipx,jpx) += -k11;
			Kg(ipx,jpy) += -k12;

			Kg(ipy,ipx) += k12;
			Kg(ipy,ipy) += k22;
			Kg(ipy,jpx) += -k12;
			Kg(ipy,jpy) += -k22;

			Kg(jpx,ipx) += -k11;
			Kg(jpx,ipy) += -k12;
			Kg(jpx,jpx) += k11;
			Kg(jpx,jpy) += k12;

			Kg(jpy,ipx) += -k12;
			Kg(jpy,ipy) += -k22;
			Kg(jpy,jpx) += k12;
			Kg(jpy,jpy) += k22;
		}
		cout << "UTriMesh: movemesh_farhat(): Done iterating over edges.\n";

		//Set boundary condtions
		Matrix<double> b(2*npoin,1);
		for(int i = 0; i < 2*npoin; i=i+2)
		{
			b(i) = (*xb)(i/2);
			b(i+1) = (*yb)(i/2);
		}
		double cbig = 1e30;			// big number for Dirichlet BCs
		// apply Dirichlet BCs using cbig
		for(int i = 0; i < 2*npoin; i++)
		{
			if(bflag2(i,0) == 1)
			{
				b(i,0) *= (cbig * Kg(i,i));
				Kg(i,i) *= cbig;
			}
		}

		cout << "UTriMesh: movemesh_farhat(): BCs set.\n";
		//Solve
		cout << "movemesh_farhat: Solving mesh-movement equations\n";
		Matrix<double> dr(2*npoin,1);
		Matrix<double> initvals(2*npoin,1);
		initvals.zeros();
		if(solver == "cholesky")
		{
			dr = cholesky(Kg,b);
		}
		else if (solver == "gausselim")
		{
			dr = gausselim(Kg,b);
		}
		else if (solver == "pointjacobi")
		{
			dr = pointjacobi(Kg, b, initvals, tol, maxiter, 'n');
		}
		else if (solver == "gaussseidel")
		{
			dr = gaussseidel(Kg, b, initvals, tol, maxiter, 'n');
		}
		else cout << "movemesh_farhat: Solver " << solver << " not found.\n";
		cout << "movemesh_farhat: Mesh-movement equations solved.\n";

		//update coordinates
		for(int i = 0; i < npoin; i++)
		{
			coords(i,0) += dr(2*i);
			coords(i,1) += dr(2*i+1);
		}

		//update jacobians
		compute_jacobians();

		// Also need to update length of faces, normals to faces in gallfa
		compute_lengths_and_normals();
		compute_edge_angles();
	}

	void movemesh_farhat(Matrix<double>* xb, Matrix<double>* yb, Matrix<int> extraflags, string solver, double tol=1e-6, int maxiter=1000)
	{
		allocate_edge_angles();
		compute_edge_angles();
		compute_jacobians();
		bflag2 = bflags2();
		cout << "UTriMesh: movemesh_farhat(): Starting\n";
		Matrix<double> Kg(2*npoin,2*npoin);		// global stiffness matrix
		Kg.zeros();
		//Matrix<double> kli(4,4);				// lineal stiffness matrix for an edge (face in 2D)
		double k11, k12, k22;
		double len;
		int ipx, ipy, jpx, jpy;

		for(int ied = 0; ied < naface; ied++)
		{
			//gather
			len = gallfa(ied,2);
			k11 = cosa(ied)*cosa(ied)/len;
			k12 = sina(ied)*cosa(ied)/len;
			k22 = sina(ied)*sina(ied)/len;
			/*
			kli(0,0) = k11; kli(0,1) = k12; kli(0,2) = -k11; kli(0,3) = -k12;
			kli(1,0) = k12; kli(1,1) = k22; kli(1,2) = -k12; kli(1,3) = -k22;
			kli(2,0) = -k11; kli(2,1) = -k12; kli(2,2) = k11; kli(2,3) = k12;
			kli(3,0) = -k12; kli(3,1) = -k22; kli(3,2) = k12; kli(3,3) = k22;
			*/
			ipx = 2*intfac(ied,2);
			ipy = 2*intfac(ied,2)+1;
			jpx = 2*intfac(ied,3);
			jpy = 2*intfac(ied,3)+1;

			//scatter
			Kg(ipx,ipx) += k11;			// contribution by x-displacement of ip node to x-force at ip node
			Kg(ipx,ipy) += k12;			// contribution by y-displacement of ip node to x-force at ip node
			Kg(ipx,jpx) += -k11;
			Kg(ipx,jpy) += -k12;

			Kg(ipy,ipx) += k12;
			Kg(ipy,ipy) += k22;
			Kg(ipy,jpx) += -k12;
			Kg(ipy,jpy) += -k22;

			Kg(jpx,ipx) += -k11;
			Kg(jpx,ipy) += -k12;
			Kg(jpx,jpx) += k11;
			Kg(jpx,jpy) += k12;

			Kg(jpy,ipx) += -k12;
			Kg(jpy,ipy) += -k22;
			Kg(jpy,jpx) += k12;
			Kg(jpy,jpy) += k22;
		}
		cout << "UTriMesh: movemesh_farhat(): Done iterating over edges.\n";

		//Now for torsional stiffnesses
		Matrix<double> R(3,6);
		Matrix<double> C(3,3); C.zeros();
		Matrix<double> Kt(6,6);		//torsional element stiffness matrix
		double x12, x23, x31, y12, y23, y31, l12, l23, l31;
		int px[3], py[3];		// nnode = 3
		//iterate over elements
		for(int iel = 0; iel < nelem; iel++)
		{
			x12 = gcoords(ginpoel(iel,1),0) - gcoords(ginpoel(iel,0),0);
			x23 = gcoords(ginpoel(iel,2),0) - gcoords(ginpoel(iel,1),0);
			x31 = gcoords(ginpoel(iel,0),0) - gcoords(ginpoel(iel,2),0);

			y12 = gcoords(ginpoel(iel,1),1) - gcoords(ginpoel(iel,0),1);
			y23 = gcoords(ginpoel(iel,2),1) - gcoords(ginpoel(iel,1),1);
			y31 = gcoords(ginpoel(iel,0),1) - gcoords(ginpoel(iel,2),1);

			l12 = x12*x12 + y12*y12;
			l23 = x23*x23 + y23*y23;
			l31 = x31*x31 + y31*y31;

			R(0,0) = -y31/l31 - y12/l12; R(0,1) = x12/l12 + x31/l31; R(0,2) = y12/l12; R(0,3) = -x12/l12; R(0,4) = y31/l31; R(0,5) = -x31/l31;
			R(1,0) = y12/l12; R(1,1) = -x12/l12; R(1,2) = -y12/l12 - y23/l23; R(1,3) = x23/l23 + x12/l12; R(1,4) = y23/l23; R(1,5) = -x23/l23;
			R(2,0) = y31/l31; R(2,1) = -x31/l31; R(2,2) = y23/l23; R(2,3) = -x23/l23; R(2,4) = -y23/l23 - y31/l31; R(2,5) = x31/l31 + x23/l23;

			C(0,0) = l12*l31/(4.0*jacobians(iel)*jacobians(iel));
			C(1,1) = l12*l23/(4.0*jacobians(iel)*jacobians(iel));
			C(2,2) = l31*l23/(4.0*jacobians(iel)*jacobians(iel));

			Kt = R.trans()*(C*R);

			for(int i = 0; i < nnode; i++)
			{
				px[i] = 2*inpoel(iel,i);
				py[i] = 2*inpoel(iel,i) + 1;
			}

			//Add contributions to global stiffness matrix
			for(int i = 0; i < nnode; i++)
				for(int j = 0; j < nnode; j++)
				{
					Kg(px[i],px[j]) += Kt(2*i,2*j);
					Kg(px[i],py[j]) += Kt(2*i,2*j+1);
					Kg(py[i],px[j]) += Kt(2*i+1,2*j);
					Kg(py[i],py[j]) += Kt(2*i+1,2*j+1);
				}
		}
		cout << "UTriMesh: movemesh_farhat(): Stiffness matrix assembled.\n";

		//Set boundary condtions
		Matrix<double> b(2*npoin,1);
		Matrix<double> extf(2*npoin,1);
		for(int i = 0; i < 2*npoin; i=i+2)		// convert boundary data into the "interleaved" format needed
		{
			b(i) = xb->get(i/2);
			b(i+1) = yb->get(i/2);
			extf(i) = extraflags.get(i/2);
			extf(i+1) = extraflags.get(i/2);
		}
		cout << "UTriMesh: movemesh_farhat(): Processed boundary inputs.\n";
		double cbig = 1e30;			// big number for Dirichlet BCs
		// apply Dirichlet BCs using cbig
		for(int i = 0; i < 2*npoin; i++)
		{
			if(bflag2.get(i,0) == 1 || extf.get(i,0)==1)
			{
				b(i,0) *= (cbig * Kg(i,i));
				Kg(i,i) *= cbig;
			}
		}

		cout << "UTriMesh: movemesh_farhat(): BCs set.\n";
		//Solve
		cout << "movemesh_farhat: Solving mesh-movement equations\n";
		Matrix<double> dr(2*npoin,1);
		Matrix<double> initvals(2*npoin,1);
		initvals.zeros();
		if(solver == "cholesky")
		{
			dr = cholesky(Kg,b);
		}
		else if (solver == "gausselim")
		{
			dr = gausselim(Kg,b);
		}
		else if (solver == "pointjacobi")
		{
			dr = pointjacobi(Kg, b, initvals, tol, maxiter, 'n');
		}
		else if (solver == "gaussseidel")
		{
			dr = gaussseidel(Kg, b, initvals, tol, maxiter, 'n');
		}
		else cout << "movemesh_farhat: Solver " << solver << " not found.\n";
		cout << "movemesh_farhat: Mesh-movement equations solved.\n";

		//update coordinates
		for(int i = 0; i < npoin; i++)
		{
			coords(i,0) += dr(2*i);
			coords(i,1) += dr(2*i+1);
		}

		//update jacobians
		compute_jacobians();

		// Also need to update length of faces, normals to faces in gallfa
		compute_lengths_and_normals();
		compute_edge_angles();
	}
};

} // end namespace acfd
