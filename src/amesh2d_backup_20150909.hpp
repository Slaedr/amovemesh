// Data structure and setup for 2D unstructured mesh.
// Aditya Kashi
// Aug 12, 2015

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
#ifndef _GLIBCXX_VECTOR
#include <vector>
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

#define __AMESH2DGENERAL_H

using namespace std;
using namespace amat;

namespace acfd {

// This function is a cyclic permutation of consecutive integers from 'start' to 'end' (inclusive). It returns the integer (between 'start' and 'end') that is 'off' integers away from 'n' in the cyclic order.
int perm(int start, int end, int n, int off)
{
	if(n > end) { cout << "Permutation point error!\n"; return 0; }
	if(off == 0) return n;

	CircList<int> list(start);
	for(int i = start+1; i <= end; i++)
		list.push(i);

	Node<int>* nn = list.find(n);
	Node<int>* cur = nn;
	for(int i = 0; i < off; i++)
		cur = cur->next;
	return cur->data;
}

/* Class UMesh2d is a general mesh class for unstructured mesh (with 1 kind of element throughout - hybrid mesh is NOT supporte) */
class UMesh2d
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
	int nbpoin;		// number of boundary points !! not used !!
	int nbtag;		// number of tags for each boundary face
	int ndtag;		// number of tags for each element
	Matrix<double> coords;
	Matrix<int> inpoel;
	Matrix<int> bface;
	Matrix<double> vol_regions;		// to hold volume region markers, if any

	Matrix<int> esup_p;
	Matrix<int> esup;
	Matrix<int> psup_p;
	Matrix<int> psup;
	Matrix<int> esuel;
	Matrix<int> intfac;

public:

	UMesh2d() {}

	UMesh2d(const UMesh2d& other)
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
		nbtag = other.nbtag;
		ndtag = other.ndtag;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		vol_regions = other.vol_regions;
		esup = other.esup;
		esup_p = other.esup_p;
		psup = other.psup;
		psup_p = other.psup_p;
		esuel = other.esuel;
		intfac = other.intfac;
		//gallfa = other.gallfa;
		//alloc_jacobians = other.alloc_jacobians;
		//jacobians = other.jacobians;
	}

	UMesh2d& operator=(const UMesh2d& other)
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
		nbtag = other.nbtag;
		ndtag = other.ndtag;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		vol_regions = other.vol_regions;
		esup = other.esup;
		esup_p = other.esup_p;
		psup = other.psup;
		psup_p = other.psup_p;
		esuel = other.esuel;
		intfac = other.intfac;
		//gallfa = other.gallfa;
		//alloc_jacobians = other.alloc_jacobians;
		//jacobians = other.jacobians;
		return *this;
	}

	double gcoords(int pointno, int dim) const
	{
		return coords.get(pointno,dim);
	}
	int ginpoel(int elemno, int locnode)
	{
		return inpoel.get(elemno, locnode);
	}
	int gbface(int faceno, int val)
	{
		return bface.get(faceno, val);
	}

	void setcoords(Matrix<double>* c)
	{ coords = *c; }

	void setinpoel(Matrix<int>* inp)
	{ inpoel = *inp; }

	void setbface(Matrix<int>* bf)
	{ bface = *bf; }

	Matrix<double>* getcoords()
	{ return &coords; }

	int gesup(int i) { return esup.get(i); }
	int gesup_p(int i) { return esup_p.get(i); }
	int gpsup(int i) { return psup.get(i); }
	int gpsup_p(int i) { return psup_p.get(i); }
	int gesuel(int ielem, int jnode) { return esuel.get(ielem, jnode); }
	int gintfac(int face, int i) { return intfac.get(face,i); }
	//double gjacobians(int ielem) { return jacobians(ielem,0); }

	int gnpoin() { return npoin; }
	int gnelem() { return nelem; }
	int gnface() { return nface; }
	int gnbface() { return nbface; }
	int gnnode() { return nnode; }
	int gndim() { return ndim; }
	int gnaface() {return naface; }
	int gnfael() { return nfael; }
	int gnnofa() { return nnofa; }
	int gnbpoin() { cout << "UTriMesh: ! Invalid access to 'gnbpoin()'!!" << endl; return nbpoin; }
	int gnbtag() {return nbtag; }
	int gndtag() { return ndtag; }

	void readDomn(ifstream& infile)
	{
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
		infile >> nfael;
		infile >> nnofa;
		infile >> ch;			//get the newline
		do
			ch = infile.get();
		while(ch != '\n');
		infile >> nelem; infile >> npoin; infile >> nface;
		infile >> dummy; 				// get time
		ch = infile.get();			// clear newline

		cout << "UTriMesh: Number of elements: " << nelem << ", number of points: " << npoin << ", number of nodes per element: " << nnode << endl;
		cout << "Number of boundary faces: " << nface << ", Number of dimensions: " << ndim;

		nbtag = 3;

		//cout << "\nUTriMesh: Allocating coords..";
		coords.setup(npoin, ndim, ROWMAJOR);
		//cout << "\nUTriMesh: Allocating inpoel..\n";
		inpoel.setup(nelem, nnode, ROWMAJOR);
		//cout << "UTriMesh: Allocating bface...\n";
		bface.setup(nface, ndim + nbtag, ROWMAJOR);

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

		for(int i = 0; i < nface; i++)
		{
			infile >> dum;
			for(int j = 0; j < ndim + nbtag; j++)
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
	}

	void readGmsh2(string mfile, int dimensions)
	// Reads mesh from Gmsh 2 format file
	{
		cout << "UMesh2d: readGmsh2(): Reading mesh file...\n";
		int dum; double dummy; string dums; char ch;
		ndim = dimensions;

		ifstream infile(mfile);
		for(int i = 0; i < 4; i++)		//skip 4 lines
			do
				ch = infile.get();
			while(ch != '\n');

		infile >> npoin;
		cout << "UMesh2d: readGmsh2(): No. of points = " << npoin << endl;
		coords.setup(npoin,ndim);

		// read coords of points
		for(int i = 0; i < npoin; i++)
		{
			infile >> dum;
			for(int j = 0; j < ndim; j++)
				infile >> coords(i,j);
			if(ndim < 3) infile >> dummy;
		}
		infile >> dums;
		infile >> dums;
		//cout << "UMesh2d: readGmsh2(): coords read." << endl;

		int nelm, elmtype, nbtags, ntags;
		ndtag = 0; nbtag = 0;
		infile >> nelm;
		Matrix<int> elms(nelm,20);
		nface = 0; nelem = 0;
		//cout << "UMesh2d: readGmsh2(): Total number of elms is " << nelm << endl;

		for(int i = 0; i < nelm; i++)
		{
			infile >> dum;
			infile >> elmtype;
			// assumption here is that elm-type is same for all faces and same for all elements
			switch(elmtype)
			{
				case(1): // linear edge
					nnofa = 2;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(8): // quadratic edge
					nnofa = 3;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(2): // linear triangles
					nnode = 3;
					nfael = 3;
					nnofa = 2;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(3):	// linear quads
					nnode = 4;
					nfael = 4;
					nnofa = 2;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(9):	// quadratic triangles
					nnode = 6;
					nfael = 3;
					nnofa = 3;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(16):	// quadratic quad (8 nodes)
					nnode = 8;
					nfael = 4;
					nnofa = 3;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(10):	// quadratic quad (9 nodes)
					nnode = 9;
					nfael = 4;
					nnofa = 3;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				default:
					cout << "! UMesh2d: readGmsh2(): Element type not recognized. Setting as linear triangle." << endl;
					nnode = 3;
					nfael = 3;
					nnofa = 2;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
			}
		}
		//cout << "UMesh2d: readGmsh2(): Done reading elms" << endl;

		if(nface > 0)
			bface.setup(nface, nnofa+nbtag);
		else cout << "UMesh2d: readGmsh2(): NOTE: There is no data to populate bface!" << endl;

		inpoel.setup(nelem, nnode);
		vol_regions.setup(nelem, ndtag);

		cout << "UMesh2d: readGmsh2(): Done. No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << "\n, number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << endl;

		// write into inpoel and bface
		// the first nface rows to be read are boundary faces
		for(int i = 0; i < nface; i++)
		{
			for(int j = 0; j < nnofa; j++)
				bface(i,j) = elms(i,j)-1;			// -1 to correct for the fact that our numbering starts from zero
			for(int j = nnofa; j < nnofa+nbtag; j++)
				bface(i,j) = elms(i,j);
		}
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
				inpoel(i,j) = elms(i+nface,j)-1;
			for(int j = 0; j < ndtag; j++)
				vol_regions(i,j) = elms(i+nface,j+nnode);
		}
		infile.close();
	}

	void printmeshstats()
	{
		cout << "UMesh2d: No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << ", number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << endl;
	}

	void writeGmsh2(string mfile)
	{
		cout << "UMesh2d: writeGmsh2(): writing mesh to file " << mfile << endl;
		// decide element type first, based on nfael/nnode and nnofa
		int elm_type = 2;
		int face_type = 1;
		if(nnode == 4)
			elm_type = 3;
		else if(nnode == 6)
			elm_type = 9;
		else if(nnode == 8)
			elm_type = 16;
		else if(nnode==9)
			elm_type = 10;

		if(nnofa == 3)
			face_type = 8;

		ofstream outf(mfile);
		//cout << "nodes\n";
		outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
		outf << "$Nodes\n" << npoin << '\n';
		for(int ip = 0; ip < npoin; ip++)
		{
			outf << ip+1;
			for(int j = 0; j < ndim; j++)
			 	outf << " " << coords(ip,j);
			for(int j = 3-ndim; j > 0; j--)
				outf << " " << 0.0;
			outf << '\n';
		}
		outf << "$EndNodes\n";
		
		//cout << "boundary faces\n";
		outf << "$Elements\n" << nelem+nface << '\n';
		// boundary faces first
		for(int iface = 0; iface < nface; iface++)
		{
			outf << iface+1 << " " << face_type << " " << nbtag;
			for(int i = nnofa; i < nnofa+nbtag; i++)
				outf << " " << bface(iface,i);			// write tags
			for(int i = 0; i < nnofa; i++)
				outf << " " << bface(iface,i)+1;		// write nodes
			outf << '\n';
		}
		//cout << "elements\n";
		for(int iel = 0; iel < nelem; iel++)
		{
			outf << nface+iel+1 << " " << elm_type << " " << ndtag;
			for(int i = 0; i < ndtag; i++)
				outf << " " << vol_regions(iel,i);
			for(int i = 0; i < nnode; i++)
				outf << " " << inpoel(iel,i)+1;
			outf << '\n';
		}
		outf << "$EndElements\n";

		outf.close();
	}

	void compute_topological()
	{
		// NOTE: Currently only works for linear mesh
		
		cout << "UMesh2d: compute_topological(): Calculating and storing topological information...\n";
		//1. Elements surrounding points
		//cout << "UMesh2d: compute_topological(): Elements surrounding points\n";
		esup_p.setup(npoin+1,1,ROWMAJOR);
		esup_p.zeros();

		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
			{
				esup_p(inpoel(i,j)+1,0) += 1;	// inpoel(i,j) + 1 : the + 1 is there because the storage corresponding to the first node begins at 0, not at 1
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

		//2. Points surrounding points - NOTE: does not work for anything but linear triangles!!
		cout << "UMesh2d: compute_topological(): Points surrounding points\n";
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
				
				// find local node number of ip in ielem
				int inode;
				for(int jnode = 0; jnode < nnode; jnode++)
					if(inpoel(ielem,jnode) == ip) inode = jnode;
				
				vector<bool> nbd(nnode);		// contains true if that local node number is connected to inode
				for(int j = 0; j < nnode; j++)
					nbd[j] = false;
				
				if(nnode == 3)
					for(int i = 0; i < nbd.size(); i++)
						nbd[i] = true;
				else if(nnode == 4)
					for(int jnode = 0; jnode < nnode; jnode++)
					{
						if(jnode == perm(0,nnode-1,inode,1) || jnode == perm(0,nnode-1, inode, -1))
							nbd[jnode] = true;
					}
				
				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip && nbd[inode])		// test of this point as already been counted as a surrounding point of ip
					{
						istor++;
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
				
				// find local node number of ip in ielem
				int inode;
				for(int jnode = 0; jnode < nnode; jnode++)
					if(inpoel(ielem,jnode) == ip) inode = jnode;
				
				vector<bool> nbd(nnode);		// contains true if that local node number is connected to inode
				for(int j = 0; j < nnode; j++)
					nbd[j] = false;
				
				if(nnode == 3)
					for(int i = 0; i < nbd.size(); i++)
						nbd[i] = true;
				else if(nnode == 4)
					for(int jnode = 0; jnode < nnode; jnode++)
					{
						if(jnode == perm(0,nnode-1,inode,1) || jnode == perm(0,nnode-1, inode, -1))
							nbd[jnode] = true;
					}
				
				//loop over nodes of the element
				for(int inode = 0; inode < nnode; inode++)
				{
					//Get global index of this node
					int jpoin = inpoel(ielem, inode);
					if(lpoin(jpoin,0) != ip && nbd[inode])		// test of this point as already been counted as a surrounding point of ip
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
		//cout << "UMesh2d: compute_topological(): Elements surrounding elements...\n";

		//Matrix<int> lpoin(npoin,1);
		esuel.setup(nelem, nfael, ROWMAJOR);
		for(int ii = 0; ii < nelem; ii++)
			for(int jj = 0; jj < nfael; jj++)
				esuel(ii,jj) = -1;
		Matrix<int> lpofa(nfael, nnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in [0:nnofa], i in [0:nfael])
		/*lpofa(0,0) = 1; lpofa(0,1) = 2;
		lpofa(1,0) = 2; lpofa(1,1) = 0;
		lpofa(2,0) = 0; lpofa(2,1) = 1;*/
		for(int i = 0; i < nfael; i++)
		{
			for(int j = 0; j < nnofa; j++)
			{
				lpofa(i,j) = perm(0,nnode-1,i,j);
			}
		}
		//lpofa.mprint();
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

		/* Computes, for each face, the elements on either side, the starting node and the ending node of the face. This is stored in intfac. Also computes unit normals to, and lengths of, each face as well as boundary flags of boundary faces, in gallfa.
		The orientation of the face is such that the element with smaller index is always to the left of the face, while the element with greater index is always to the right of the face.
		NOTE: After the following portion, esuel holds (nelem + face no.) for each ghost cell, instead of -1 as before.*/

		//cout << "UMesh2d: compute_topological(): Computing intfac..." << endl;
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
		cout << "UMesh2d: compute_topological(): Number of boundary faces = " << nbface << endl;
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
		cout << "UMesh2d: compute_topological(): Number of all faces = " << naface << endl;

		//allocate intfac
		intfac.setup(naface,nnofa+2,ROWMAJOR);

		//reset face totals
		nbface = naface = 0;

		//second run: populate intfac
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nnode; in++)
			{
				int in1 = perm(0,nnode-1,in,1);
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je == -1)
				{
					esuel(ie,in) = nelem+nbface;
					intfac(nbface,0) = ie;
					intfac(nbface,1) = nelem+nbface;
					intfac(nbface,2) = inpoel(ie,in);
					intfac(nbface,3) = inpoel(ie,in1);

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
				//int in2 = perm(0,nnode-1,in1,1);
				int je = esuel(ie,in);
				if(je > ie && je < nelem)
				{
					intfac(naface,0) = ie;
					intfac(naface,1) = je;
					intfac(naface,2) = inpoel(ie,in);
					intfac(naface,3) = inpoel(ie,in1);
					naface++;
				}
			}
		}
	}

	UMesh2d convertLinearToQuadratic()
	{
		cout << "UMesh2d: convertLinearToQuadratic(): Producing quadratic mesh from liner mesh" << endl;
		UMesh2d q;
		if(nnofa != 2) { cout << "! UMesh2d: convertLinearToQuadratic(): Mesh is not linear!!" << endl; return q;}
		
		if(nnode == 3)			// for simplicial mesh
		{
			cout << "UMesh2d: convertLinearToQuadratic(): Simplicial mesh." << endl;
			int parm = 1;		// 1 extra node per face
			q.ndim = ndim;
			q.npoin = npoin + naface;
			q.nelem = nelem;
			q.nface = nface;
			q.nbface = nbface;
			q.naface = naface;
			q.nnofa = nnofa+parm;
			q.nnode = nnode + nfael*parm;
			q.nfael = nfael;
			q.nbtag = nbtag;
			q.ndtag = ndtag;

			q.coords.setup(q.npoin, q.ndim);
			q.inpoel.setup(q.nelem, q.nnode);
			q.bface.setup(q.nface, q.nnofa+q.nbtag);

			for(int i = 0; i < npoin; i++)
				for(int j = 0; j < ndim; j++)
					q.coords(i,j) = coords(i,j);

			for(int i = 0; i < nelem; i++)
				for(int j = 0; j < nnode; j++)
					q.inpoel(i,j) = inpoel(i,j);

			for(int i = 0; i < nface; i++)
			{
				for(int j = 0; j < nnofa; j++)
					q.bface(i,j) = bface(i,j);
				for(int j = nnofa; j < nnofa+nbtag; j++)
					q.bface(i,j+parm) = bface(i,j);
			}

			q.vol_regions = vol_regions;

			int ied, p1, p2, ielem, jelem, idim, inode, lp1, lp2, ifa;

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over boundary faces..." << endl;
			// iterate over boundary faces
			for(ied = 0; ied < nbface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+ied;

				// find the bface that this face corresponds to
				for(ifa = 0; ifa < nface; ifa++)
				{
					if((p1 == bface(ifa,0) && p2 == bface(ifa,1)) || (p1 == bface(ifa,1) && p2 == bface(ifa,0)))	// face found
					{
						q.bface(ifa,nnofa) = npoin+ied;
					}
				}
			}

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over internal faces..." << endl;
			// iterate over internal faces
			for(ied = nbface; ied < naface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+ied;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(jelem,inode)) lp1 = inode;
					if(p2 == inpoel(jelem,inode)) lp2 = inode;
				}

				// in the right element, the new point is in face ip2
				q.inpoel(jelem, nnode+lp2) = npoin+ied;
			}
			cout << "UMesh2d: convertLinearToQuadratic(): Done." << endl;
			return q;
		}
		else 					// for non-simplicial mesh, add extra points at cell-centres as well
		{
			cout << "UMesh2d: convertLinearToQuadratic(): Non-simplicial mesh." << endl;
			int parm = 1;		// 1 extra node per face
			q.ndim = ndim;
			q.npoin = npoin + naface + nelem;
			q.nelem = nelem;
			q.nface = nface;
			q.nbface = nbface;
			q.naface = naface;
			q.nnofa = nnofa+parm;
			q.nnode = nnode + nfael*parm + 1;
			q.nfael = nfael;
			q.nbtag = nbtag;
			q.ndtag = ndtag;

			q.coords.setup(q.npoin, q.ndim);
			q.inpoel.setup(q.nelem, q.nnode);
			q.bface.setup(q.nface, q.nnofa+q.nbtag);

			for(int i = 0; i < npoin; i++)
				for(int j = 0; j < ndim; j++)
					q.coords(i,j) = coords(i,j);

			for(int i = 0; i < nelem; i++)
				for(int j = 0; j < nnode; j++)
					q.inpoel(i,j) = inpoel(i,j);

			for(int i = 0; i < nface; i++)
			{
				for(int j = 0; j < nnofa; j++)
					q.bface(i,j) = bface(i,j);
				for(int j = nnofa; j < nnofa+nbtag; j++)
					q.bface(i,j+parm) = bface(i,j);
			}

			q.vol_regions = vol_regions;

			int ied, p1, p2, ielem, jelem, idim, inode, lp1, lp2, ifa;
			
			// get cell centres
			for(int iel = 0; iel < nelem; iel++)
			{
				double c_x = 0, c_y = 0;
				
				for(int inode = 0; inode < nnode; inode++)
				{
					c_x += coords(inpoel(iel,inode),0);
					c_y += coords(inpoel(iel,inode),1);
				}
				c_x /= nnode;
				c_y /= nnode;
				q.coords(npoin+iel,0) = c_x;
				q.coords(npoin+iel,1) = c_y;
				q.inpoel(iel,q.nnode-1) = npoin+iel;
			}

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over boundary faces..." << endl;
			// iterate over boundary faces
			for(ied = 0; ied < nbface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+nelem+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+nelem+ied;

				// find the bface that this face corresponds to
				for(ifa = 0; ifa < nface; ifa++)
				{
					if((p1 == bface(ifa,0) && p2 == bface(ifa,1)) || (p1 == bface(ifa,1) && p2 == bface(ifa,0)))	// face found
					{
						q.bface(ifa,nnofa) = npoin+nelem+ied;
					}
				}
			}

			//cout << "UMesh2d: convertLinearToQuadratic(): Iterating over internal faces..." << endl;
			// iterate over internal faces
			for(ied = nbface; ied < naface; ied++)
			{
				ielem = intfac(ied,0);
				jelem = intfac(ied,1);
				p1 = intfac(ied,2);
				p2 = intfac(ied,3);

				for(idim = 0; idim < ndim; idim++)
					q.coords(npoin+nelem+ied,idim) = (coords(p1,idim) + coords(p2,idim))/2.0;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(ielem,inode)) lp1 = inode;
					if(p2 == inpoel(ielem,inode)) lp2 = inode;
				}

				// in the left element, the new point is in face ip1 (ie, the face whose first point is ip1 in CCW order)
				q.inpoel(ielem, nnode+lp1) = npoin+nelem+ied;

				for(inode = 0; inode < nnode; inode++)
				{
					if(p1 == inpoel(jelem,inode)) lp1 = inode;
					if(p2 == inpoel(jelem,inode)) lp2 = inode;
				}

				// in the right element, the new point is in face ip2
				q.inpoel(jelem, nnode+lp2) = npoin+nelem+ied;
			}
			cout << "UMesh2d: convertLinearToQuadratic(): Done." << endl;
			return q;
		}
	}
};

} // end namespace acfd
