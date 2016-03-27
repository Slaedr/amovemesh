// Data structure and setup for 2D unstructured triangular mesh.
// Aditya Kashi
// Feb 5, 2015

// Feb 26, 2015: Adding mesh-movement function movemesh2d_bs().

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
#define nthreads_mesh 8
#endif
#endif

#define __AMESH2D_H

#ifndef __AMATRIX2_H
#include "amatrix2.hpp"
#endif
#ifndef __ASPARSEMATRIX_H
#include "asparsematrix.hpp"
#endif
#ifndef __ALINALG_H
#include "alinalg.hpp"
#endif
#ifndef __ADATASTRUCTURES_H
#include "adatastructures.hpp"
#endif



using namespace std;
using namespace amat;

namespace acfd {

const int n_extra_fields_in_bface = 2;

// Unstructured triangular mesh class
class UTriMesh
{
private:
	int npoin;
	int nelem;
	int ndim;
	int nnode;		// number of nodes to an element
	int nfael;		// number of faces to an element (equal to number of edges to an element in 2D)
	int nnofa;		// number of node in a face -- needs to be generalized in case of general grids
	int naface;		// total number of (internal and boundary) faces
	int nbface;		// number of boundary faces as calculated by compute_face_data(), as opposed to nface which is read from file
	int nbpoin;		// number of boundary points
	int nface;		// number of boundary faces - not always used
	int nbdata;		// number of integers associated with each boundary point
	Matrix<double> coords;
	Matrix<int> inpoel;
	Matrix<int> bface;
	Matrix<int> bpoints;

	// The following 4 are for moving mesh
	Matrix<int> bflag;		// for the basicspring mesh mpvement, and for output of interior node data
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

	Matrix<double>* xb;			// boundary displacements for mesh-movement
	Matrix<double>* yb;

public:

	UTriMesh() {}

	UTriMesh(Matrix<double>* co, Matrix<int>* inp, int n_poin, int n_elem, int n_node, int n_nofa)
	{
		alloc_jacobians = false;
		coords = *co;
		inpoel = *inp;
		npoin = n_poin;
		nelem = n_elem;
		nnode = n_node;
		nnofa = n_nofa;
	}

	UTriMesh(ifstream& infile)
	{
		alloc_jacobians = false;

		// Do file handling here to populate npoin and nelem
		cout << "UTriMesh: Reading mesh file...\n";
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
		infile >> nelem; infile >> npoin; infile >> nbpoin;
		infile >> dummy; 				// get time
		ch = infile.get();			// clear newline

		cout << "UTriMesh: Number of elements: " << nelem << ", number of points: " << npoin << ", number of nodes per element: " << nnode << endl;
		cout << "Number of boundary points: " << nbpoin << ", Number of dimensions: " << ndim;

		nbdata = 5;

		//cout << "\nUTriMesh: Allocating coords..";
		coords.setup(npoin, ndim, ROWMAJOR);
		//cout << "\nUTriMesh: Allocating inpoel..\n";
		inpoel.setup(nelem, nnode, ROWMAJOR);
		//cout << "UTriMesh: Allocating bface...\n";
		bpoints.setup(nbpoin, nbdata, ROWMAJOR);

		nfael = 3;	// number of faces per element
		nnofa = 2;	// number of nodes per face

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

		for(int i = 0; i < nbpoin; i++)
		{
			infile >> dum;
			for(int j = 0; j < nbdata; j++)
			{
				infile >> bpoints(i,j);
			}
			if (i==nbpoin-1) break;
			do
				ch = infile.get();
			while(ch!='\n');
		}
		cout << "UTriMesh: Populated bpoints. Done reading mesh.\n";
		//correct first 2 columns of bface
		/*for(int i = 0; i < nface; i++)
			for(int j = 0; j < 2; j++)
				bface(i,j)--;*/

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

		cout << "UTriMesh: Done.\n";
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
		nbpoin = other.nbpoin;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		naface = other.naface;
		nbface = other.nbface;
		nfael = other.nfael;
		nnofa = other.nnofa;
		coords = other.coords;
		inpoel = other.inpoel;
		bpoints = other.bpoints;
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
		nbpoin = other.nbpoin;
		nface = other.nface;
		ndim = other.ndim;
		nnode = other.nnode;
		naface = other.naface;
		nbface = other.nbface;
		nfael = other.nfael;
		nnofa = other.nnofa;
		coords = other.coords;
		inpoel = other.inpoel;
		bpoints = other.bpoints;
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

	void readYibinMesh(string meshn)
	{
		/* Reads Yibin Wang's dat mesh format, which does not contain any boundary information. */

		ifstream fin(meshn);

		fin >> npoin >> nelem;

		// some assumptions for properties not specified in the file
		ndim = 2; nnode = 4; nnofa = 2; nfael = 4;
		int a;

		coords.setup(npoin,ndim);
		inpoel.setup(nelem,nnode);

		for(int ipoin = 0; ipoin < npoin; ipoin++)
		{
			for(int idim = 0; idim < ndim; idim++)
				fin >> coords(ipoin,idim);
		}
		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int inode = 0; inode < nnode; inode++) {
				fin >> a;
				inpoel(ielem,inode) = a - 1;
			}
		}

		// finished reading mesh file
		fin.close();

		// Now populate bface and bpoints
		// NOTE: only works for structured mesh!

		Matrix<int> lpoin(npoin,1);
		lpoin.zeros();

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int inode = 0; inode < nnode; inode++)
				lpoin(inpoel(ielem,inode))++;
		}

		nbpoin = 0;
		nface = 0;
		for(int i = 0; i < npoin; i++)
			if(lpoin(i) <= 3) nbpoin++;

		bpoints.setup(nbpoin,1+2);		// stores global point number and surrounding bfaces for each bpoint

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			int nextnode;
			for(int inode = 0; inode < nnode; inode++)
			{
				nextnode = perm(0,nnode-1,inode,1);
				if(lpoin.get(inpoel(ielem,inode)) <= 3 && lpoin.get(inpoel.get(ielem,nextnode)) <= 3) {
					nface++;
				}
			}
		}

		bface.setup(nface,nnofa+n_extra_fields_in_bface);
		bface.zeros();
		Matrix<int> trav(npoin,1);
		trav.zeros();
		int k = 0, l = 0;

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			int thisgnode, nextnode, nextgnode;
			for(int inode = 0; inode < nnode; inode++)
			{
				thisgnode = inpoel(ielem,inode);
				if(lpoin(thisgnode) <= 3 && trav(thisgnode)==0) {
					bpoints(l,0) = thisgnode;
					bpoints(l,1) = 1;
					trav(thisgnode) = 1;
					l++;
				}
				nextnode = perm(0,nnode-1,inode,1);
				nextgnode = inpoel(ielem,nextnode);
				if(lpoin.get(thisgnode) <= 3 && lpoin.get(nextgnode) <= 3) {
					bface(k,0) = nextgnode;
					bface(k,1) = thisgnode;
					bface(k,2) = 1;
					k++;
				}
			}
		}

		// we now have bpoints and bface for a structured mesh.
		cout << "UTriMesh: readYibinMesh(): " << npoin << " points, " << nelem << " elements, " << nface << " boundary faces, " << nbpoin << " boundary points." << endl;
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

	int gbpoints(int pointno, int val)
	{
		return bpoints.get(pointno, val);
	}

	void setcoords(Matrix<double>* c)
	{ coords = *c; }

	void setinpoel(Matrix<int>* inp)
	{ inpoel = *inp; }

	void setbpoints(Matrix<int>* bf)
	{ bpoints = *bf; }

	Matrix<double>* getcoords()
	{ return &coords; }

	int gesup(int i) { return esup.get(i); }
	int gesup_p(int i) { return esup_p.get(i); }
	int gpsup(int i) { return psup.get(i); }
	int gpsup_p(int i) { return psup_p.get(i); }
	int gesuel(int ielem, int jnode) { return esuel.get(ielem, jnode); }		//returns element number at face opposite to node number jnode
	double gjacobians(int ielem) { return jacobians(ielem,0); }

	int gnpoin() { return npoin; }
	int gnelem() { return nelem; }
	int gnbface() { return nbface; }
	int gnbpoin() { return nbpoin; }
	int gnnode() { return nnode; }
	int gndim() { return ndim; }
	int gnaface() {return naface; }
	int gnfael() { return nfael; }
	int gnnofa() { return nnofa; }
	int gnbdata() {return nbdata; }

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
		bool flagj = false;
		for(int i = 0; i < nelem; i++)
		{
			if(jacobians(i,0) <= 1e-15) {
				out << i << " " << jacobians(i,0) << '\n';
				flagj = true;
			}
		}
		if(flagj == true) cout << "UTriMesh: detect_negative_jacobians(): There exist element(s) with negative jacobian!!\n";
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
		gallfa.setup(naface, 3, ROWMAJOR);
		for(int i = 0; i < naface; i++)
		{
			gallfa(i,0) = coords(intfac(i,3),1) - coords(intfac(i,2),1);
			gallfa(i,1) = -1.0*(coords(intfac(i,3),0) - coords(intfac(i,2),0));
			gallfa(i,2) = sqrt(pow(gallfa(i,0),2) + pow(gallfa(i,1),2));
			//Normalize the normal vector components
			gallfa(i,0) /= gallfa(i,2);
			gallfa(i,1) /= gallfa(i,2);
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

	void writeGmsh2(string mfile)
	{
		ofstream outf(mfile);

		outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
		outf << "$Nodes\n" << npoin << '\n';
		for(int ip = 0; ip < npoin; ip++)
		{
			outf << ip+1 << " " << coords(ip,0) << " " << coords(ip,1) << " " << 0 << '\n';
		}
		outf << "$EndNodes\n$Elements\n" << nelem/*+nface*/ << '\n';
		// boundary faces first
		/*for(int iface = 0; iface < nface; iface++)
		{
			outf << iface+1 << " 1 2 0 1";
			for(int i = 0; i < nnofa; i++)
				outf << " " << bface(iface,i)+1;
			outf << '\n';
		}*/
		for(int iel = 0; iel < nelem; iel++)
		{
			outf << /*nface+*/iel+1 << " 2 2 0 2";
			for(int i = 0; i < nnode; i++)
				outf << " " << inpoel(iel,i)+1;
			outf << '\n';
		}
		outf << "$EndElements\n";

		outf.close();
	}



	//------------------------- Mesh movement functions -----------------------------------//

	void movemesh(Matrix<double> disp)
	// Updates coords of mesh by adding displacements to them; the displacements specified by the 2*npoin-by-1 vector disp, which contains x-displacements in the first npoin entries and y-displacements in the rest.
	// Required for elasticity-based mesh movement.
	{
		for(int ip = 0; ip < npoin; ip++)
		{
			coords(ip,0) += disp(ip);
			coords(ip,1) += disp(npoin+ip);
		}
	}
	
	void set_boundary_motion(Matrix<double>* xx, Matrix<double>* yy)
	{
		xb = xx; yb = yy;
	}

	void setbflags2()
	{
		bflag2.setup(2*npoin,1);
		for(int i = 0; i < 2*nbpoin; i+=2)
		{
			bflag2(i) = 1;
			bflag2(i+1) = 1;
		}
	}

	// returns stiffness of the edge between global nodes i and j
	double k(int i, int j)
	{
		 return 1/sqrt((gcoords(i,0)-gcoords(j,0))*(gcoords(i,0)-gcoords(j,0)) + (gcoords(i,1)-gcoords(j,1))*(gcoords(i,1)-gcoords(j,1)));
	}

	// Moves the mesh according to prescribed boundary displacements xb and yb
	// blfag contains a 0-or-1 flag that indicates a boundary point
	// Note that the original mesh is overwritten.
	// NOTE: MAKE SURE to update geoel, gallfa and any other downstream data

	void movemesh_basicspring(Matrix<double> xb, Matrix<double> yb, string solver, double tol=1e-6, int maxiter=1000)
	{
		Matrix<double> A(npoin, npoin, ROWMAJOR);
		A.zeros();

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

	void movemesh_farhat(double tol, int maxiter)
	/* Moves the mesh based on boundary displacements set by the set_boundary_motion() method. The technique is based on Farhat et al's 2D torsional
	   spring analogy method.
	*/
	{
		cout << "UTriMesh: movemesh_farhat(): Starting\n";
		SpMatrix Kg(2*npoin,2*npoin);		// global stiffness matrix
		//Kg.zeros();
		//Matrix<double> kli(4,4);				// lineal stiffness matrix for an edge (face in 2D)
		double k11, k12, k22;
		double len;
		int ipx, ipy, jpx, jpy;
		int ied;

		// copy class variables to local variables for use in OpenMP
		Matrix<double>* cosa = &(UTriMesh::cosa);
		Matrix<double>* sina = &(UTriMesh::sina);
		Matrix<double>* gallfa = &(UTriMesh::gallfa);
		Matrix<int>* intfac = &(UTriMesh::intfac);
		Matrix<double>* jacobians = &(UTriMesh::jacobians);
		#ifdef _OPENMP
		int naface = UTriMesh::naface;
		int nnode = UTriMesh::nnode;
		int nelem = UTriMesh::nelem;
		#endif

		#pragma omp parallel for default(none) private(ied,len,k11,k12,k22,ipx,ipy,jpx,jpy) shared(naface,gallfa,cosa,sina,intfac,Kg) num_threads(nthreads_mesh)
		for(ied = 0; ied < naface; ied++)
		{
			//gather
			len = gallfa->get(ied,2);
			k11 = cosa->get(ied)*cosa->get(ied)/len;
			k12 = sina->get(ied)*cosa->get(ied)/len;
			k22 = sina->get(ied)*sina->get(ied)/len;

			ipx = 2*intfac->get(ied,2);
			ipy = 2*intfac->get(ied,2)+1;
			jpx = 2*intfac->get(ied,3);
			jpy = 2*intfac->get(ied,3)+1;

			//scatter
			double temp;

			#pragma omp critical (omp_lineal)
			{
				temp = Kg.get(ipx,ipx);
				Kg.set(ipx,ipx, temp+k11);			// contribution by x-displacement of ip node to x-force at ip node
				temp = Kg.get(ipx,ipy);
				Kg.set(ipx,ipy, temp+k12);			// contribution by y-displacement of ip node to x-force at ip node
				temp = Kg.get(ipx,jpx);
				Kg.set(ipx,jpx,temp-k11);
				temp = Kg.get(ipx,jpy);
				Kg.set(ipx,jpy, temp-k12);

				temp = Kg.get(ipy,ipx);
				Kg.set(ipy,ipx, temp+k12);
				temp = Kg.get(ipy,ipy);
				Kg.set(ipy,ipy,temp+k22);
				temp = Kg.get(ipy,jpx);
				Kg.set(ipy,jpx, temp-k12);
				temp = Kg.get(ipy,jpy);
				Kg.set(ipy,jpy, temp-k22);

				temp = Kg.get(jpx,ipx);
				Kg.set(jpx,ipx, temp-k11);
				temp = Kg.get(jpx,ipy);
				Kg.set(jpx,ipy, temp-k12);
				temp = Kg.get(jpx,jpx);
				Kg.set(jpx,jpx, temp+k11);
				temp = Kg.get(jpx,jpy);
				Kg.set(jpx,jpy, temp+k12);

				Kg.set(jpy,ipx, Kg.get(jpy,ipx)-k12);
				Kg.set(jpy,ipy, Kg.get(jpy,ipy)-k22);
				Kg.set(jpy,jpx, Kg.get(jpy,jpx)+k12);
				Kg.set(jpy,jpy, Kg.get(jpy,jpy)+k22);
			}

			//if(ied % 40 == 0)
			//	cout << "=" << flush;
		}
		// end parallel
		cout << endl;
		cout << "UTriMesh: movemesh_farhat(): Done iterating over edges.\n";

		//Now for torsional stiffnesses

		double x12, x23, x31, y12, y23, y31, l12, l23, l31;
		int px[3], py[3];		// nnode = 3

		int iel;
		//iterate over elements
		#pragma omp parallel for default(none) private(iel,px,py,x12,x23,x31,y12,y23,y31,l12,l23,l31) shared(naface,nnode,nelem,jacobians,Kg) num_threads(nthreads_mesh)
		for(iel = 0; iel < nelem; iel++)
		{
			Matrix<double> R(3,6);
			Matrix<double> C(3,3); C.zeros();
			Matrix<double> Kt(6,6);		//torsional element stiffness matrix

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

			C(0,0) = l12*l31/(4.0*jacobians->get(iel)*jacobians->get(iel));
			C(1,1) = l12*l23/(4.0*jacobians->get(iel)*jacobians->get(iel));
			C(2,2) = l31*l23/(4.0*jacobians->get(iel)*jacobians->get(iel));

			Kt = R.trans()*(C*R);

			for(int ii = 0; ii < nnode; ii++)
			{
				px[ii] = 2*ginpoel(iel,ii);
				py[ii] = 2*ginpoel(iel,ii) + 1;
			}

			double*** temps = new double**[nnode];
			for(int i = 0; i < nnode; i++)
			{
				temps[i] = new double*[nnode];
				for(int j = 0; j < nnode; j++)
					temps[i][j] = new double[4];
			}

			for(int i = 0; i < nnode; i++)
			{
				for(int j = 0; j < nnode; j++)
				{
					temps[i][j][0] = Kt.get(2*i,2*j);
					temps[i][j][1] = Kt.get(2*i,2*j+1);
					temps[i][j][2] = Kt.get(2*i+1,2*j);
					temps[i][j][3] = Kt.get(2*i+1,2*j+1);
				}
			}

			//Add contributions to global stiffness matrix
			#pragma omp critical (omp_torsion)
			for(int i = 0; i < nnode; i++)
			{
				for(int j = 0; j < nnode; j++)
				{
					Kg.set(px[i],px[j], Kg.get(px[i],px[j])+ temps[i][j][0] );
					Kg.set(px[i],py[j], Kg.get(px[i],py[j])+ temps[i][j][1] );
					Kg.set(py[i],px[j], Kg.get(py[i],px[j])+ temps[i][j][2] );
					Kg.set(py[i],py[j], Kg.get(py[i],py[j])+ temps[i][j][3] );
				}
			}

			for(int i = 0; i < nnode; i++)
			{
				for(int j = 0; j < nnode; j++)
					delete [] temps[i][j];
				delete [] temps[i];
			}
			delete [] temps;

			//if(iel%40==0)
			//	cout << "=" << flush;
		}
		// end parallel

		cout << endl;
		cout << "UTriMesh: movemesh_farhat(): Stiffness matrix assembled.\n";

		setbflags2();
		cout << "UTriMesh: movemesh_farhat(): bflags2 calculated.\n";

		//Set boundary condtions
		Matrix<double> b(2*npoin,1);
		for(int i = 0; i < 2*npoin; i=i+2)
		{
			b(i) = xb->get(i/2);
			b(i+1) = yb->get(i/2);
		}
		double cbig = 1e30;			// big number for Dirichlet BCs
		// apply Dirichlet BCs using cbig
		for(int i = 0; i < 2*nbpoin; i++)
		{
			//if(bflag2(i,0) == 1)
			{
				b(i,0) *= (cbig * Kg.get(i,i));
				Kg.set(i,i, Kg.get(i,i)*cbig);
			}
		}

		cout << "UTriMesh: movemesh_farhat(): BCs set.\n";
		//Solve
		cout << "movemesh_farhat: Solving mesh-movement equations\n";
		Matrix<double> dr(2*npoin,1);
		Matrix<double> initvals(2*npoin,1);
		initvals.zeros();

		dr = sparseCG_d(&Kg, b, initvals, tol, maxiter);
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

	/*double farhat_assemble(int a, int b)
	{
		double val = 0;
		//cout << "UTriMesh: farhat_assemble(): Starting\n";
		//Matrix<double> Kg(2*npoin,2*npoin);		// global stiffness matrix
		//Kg.zeros();
		//Matrix<double> kli(4,4);				// lineal stiffness matrix for an edge (face in 2D)
		double k11, k12, k22;
		double len;
		int ipx, ipy, jpx, jpy;

		for(int ied = 0; ied < naface; ied++)
		{
			ipx = 2*intfac(ied,2);
			ipy = 2*intfac(ied,2)+1;
			jpx = 2*intfac(ied,3);
			jpy = 2*intfac(ied,3)+1;

			if((ipx != a && ipy != a && jpx != a && jpy != a) || (ipx != b && ipy != b && jpx != b && jpy != b)) continue;

			//gather
			len = gallfa(ied,2);
			k11 = cosa(ied)*cosa(ied)/len;
			k12 = sina(ied)*cosa(ied)/len;
			k22 = sina(ied)*sina(ied)/len;

			//scatter
			if(ipx==a && ipx==b) val += k11;			// contribution by x-displacement of ip node to x-force at ip node
			if(ipx==a && ipy==b) val += k12;			// contribution by y-displacement of ip node to x-force at ip node
			if(ipx==a && jpx==b) val += -k11;
			if(ipx==a && jpy==b) val += -k12;

			if(ipy==a && ipx==b) val += k12;
			if(ipy==a && ipy==b) val += k22;
			if(ipy==a && jpx==b) val += -k12;
			if(ipy==a && jpy==b) val += -k22;

			if(jpx==a && ipx==b) val += -k11;
			if(jpx==a && ipy==b) val += -k12;
			if(jpx==a && jpx==b) val += k11;
			if(jpx==a && jpy==b) val += k12;

			if(jpy==a && ipx==b) val += -k12;
			if(jpy==a && ipy==b) val += -k22;
			if(jpy==a && jpx==b) val += k12;
			if(jpy==a && jpy==b) val += k22;
		}
		//cout << "UTriMesh: farhat_assemble(): Done iterating over edges.\n";

		//Now for torsional stiffnesses
		Matrix<double> R(3,6);
		Matrix<double> C(3,3); C.zeros();
		Matrix<double> Kt(6,6);		//torsional element stiffness matrix
		double x12, x23, x31, y12, y23, y31, l12, l23, l31;
		int px[3], py[3];		// nnode = 3
		bool cha = false, chb=false;
		//iterate over elements
		for(int iel = 0; iel < nelem; iel++)
		{
			for(int i = 0; i < nnode; i++)
			{
				px[i] = 2*inpoel(iel,i);
				py[i] = 2*inpoel(iel,i) + 1;
				if(px[i]==a || py[i] == a) cha = true; if(px[i]==b || py[i]==b) chb = true;
			}
			if(cha==false || chb==false) continue;

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

			//Add contributions to global stiffness matrix
			for(int i = 0; i < nnode; i++)
				for(int j = 0; j < nnode; j++)
				{
					if(px[i]==a && px[j]==b) val += Kt(2*i,2*j);
					if(px[i]==a && py[j]==b) val += Kt(2*i,2*j+1);
					if(py[i]==a && px[j]==b) val += Kt(2*i+1,2*j);
					if(py[i]==a && py[j]==b) val += Kt(2*i+1,2*j+1);
				}
		}
		//cout << "UTriMesh: movemesh_farhat(): Stiffness matrix assembled.\n";

		//BCs
		double cbig = 1e-30;		// big number for Dirichlet BCs
		if(a==b && a/2 < nbpoin) val *= cbig;

		return val;
	}

	double farhatRHS(int a)
	{
		// Note that the x- and y-components in the resulting nodal displacements are interleaved in the vector dr, while the input boundary displacements are not interleaved - xb and yb are taken separately
		double val = 0;
		//Set boundary condtions
		if(a%2 == 0)
			val = (*xb)(a/2);
		else
			val = (*yb)(a/2);

		// apply Dirichlet BCs using cbig
		if(a/2 < nbpoin)
			val *= farhat_assemble(a,a);

		//cout << "UTriMesh: movemesh_farhat(): BCs set.\n";
		return val;
	}

	void movemesh_farhat(string solver, double tol=1e-6, int maxiter=1000)
	{
		//Solve
		cout << "movemesh_farhat: Solving mesh-movement equations\n";
		Matrix<double> dr(2*npoin,1);
		Matrix<double> initvals(2*npoin,1);
		initvals.zeros();
		dr = gaussseidel_ng(2*npoin, 2*npoin, 2*npoin, initvals, tol, maxiter, 'n');

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

	Matrix<double> gaussseidel_ng(int arows, int acols, int brows, Matrix<double> xold, double tol, int maxiter, char check)
	{
		cout << "\ngaussseidel: Input LHS matrix is " << arows << " x " << acols << endl;
		Matrix<double> bb;		// dummy
		if(arows != acols) { cout << "! gaussseidel: Invalid dimensions of A and b!\n"; return bb; }
		if(xold.rows() != brows) { cout << "! gaussseidel: Invalid dimensions on xold !\n"; return bb; }
		int N = arows;

		if(check == 'y')
		{
			for(int i = 0; i < arows; i++)
			{
				double sum = 0;
				for(int j = 0; j < acols; j++)
				{
					if(i != j) sum += dabs(farhat_assemble(i,j));
				}
				if(dabs(farhat_assemble(i,i)) <= sum) cout << "* gaussseidel: LHS Matrix is NOT strictly diagonally dominant in row " << i << "!!\n";
			}
		}

		Matrix<double> x(brows,1);

		Matrix<double> M(N,1);		// diagonal elements of A
		M.zeros();

		//populate diagonal matrix M
		cout << "gaussseidel_ng: populate diagonal matrix M\n";
		for(int i = 0; i < N; i++)
			M(i) = farhat_assemble(i,i);

		//A = A - M;

		//invert M
		cout << "gaussseidel_ng: invert M\n";
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
		cout << "gaussseidel_ng: Start iterations\n";
		do
		{
			xold = x;
			Axold.zeros();
			for(int i = 0; i < N; i++)
			{
				for(int k = 0; k < i; k++)
					Axold(i,0) += farhat_assemble(i,k)*x(k,0);

				Axold(i,0) += 0*xold(i,0);
				for(int k = i+1; k < N; k++)
					Axold(i,0) += farhat_assemble(i,k)*xold(k,0);

				inter(i,0) = farhatRHS(i) - Axold(i,0);
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
				cout << "gausssiedel: Initial residue = " << initres << endl;
				first = false;
			}
			//if(c%10 == 0)
				cout << "gausssiedel(): Step " << c << ", Relative error = " << error/initres << endl;

		} while(error/initres >= tol);
		cout << "gaussseidel: No. of iterations = " << c << endl;

		return x;
	}*/

};

} // end namespace acfd
