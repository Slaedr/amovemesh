/// Data structure and setup for 3D unstructured mesh.
/// Aditya Kashi
/// May 16, 2015 (more like August 20, 2015)

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

#define __AMESH3D_H

using namespace std;
using namespace amat;

namespace acfd {

// This function is a cyclic permutation of consecutive integers from 'start' to 'end' (inclusive). It returns the integer (between 'start' and 'end') that is 'off' integers away from 'n' in the cyclic order.
/*int perm(int start, int end, int n, int off)
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
}*/

class UMesh
{
private:
	int npoin;
	int nelem;
	int nface;
	int nedge;
	int ndim;
	int nnode;		///< number of nodes to an element
	int nfael;		///< number of faces to an element
	int nedel;		///< number of edges in an element
	int nnofa;		///< number of nodes in a face -- needs to be generalized in case of general grids
	int nedfa;		///< number of edges in a face
	int nnoded;		///< number of nodes per edge
	int naface;		///< total number of (internal and boundary) faces
	int nbface;		///< number of boundary faces as calculated by compute_face_data(), as opposed to nface which is read from file
	int nbedge;		///< number of boundary edges
	int nbtag;		///< Number of boundary markers for each boundary face
	int ndtag;		///< Number of region markers for each element
	Matrix<double> coords;
	Matrix<int> inpoel;
	Matrix<int> m_inpoel;		// same as inpoel, except that it might contain different node ordering
	Matrix<int> bface;
	Matrix<int> flag_bpoin;		// a boolean flag for each point. Contains 1 if the corresponding point is a boundary point
	Matrix<double> vol_regions;		// to hold volume region markers, if any
	bool alloc_jacobians;
	Matrix<double> jacobians;

	Matrix<int> lpofa;		///< for each face of an element, it contains local node numbers of the nodes forming that face. Assumed to be same for all elements.
	Matrix<int> esup;		///< elements surrounding points
	Matrix<int> esup_p;		///< array containing index of esup where elements surrounding a certain point start
	vector<int>* psup;		///< points surrounding points
	Matrix<int> esuel;		///< elements surrounding elements
	Matrix<int> intedge;	///< edge data structure
	vector<int>* elsed;		///< elements surrounding edge
	Matrix<int> intfac;		///< face data strcture

public:

	/** No-arg constructor. */
	UMesh() {alloc_jacobians = false; psup = nullptr; elsed = nullptr;}

	UMesh(const UMesh& other)
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
		nedel = other.nedel;
		nedfa = other.nedfa;
		nnoded = other.nnoded;
		nedge = other.nedge;
		nbedge = other.nbedge;
		nbtag = other.nbtag;
		ndtag = other.ndtag;
		coords = other.coords;
		inpoel = other.inpoel;
		bface = other.bface;
		flag_bpoin = other.flag_bpoin;
		vol_regions = other.vol_regions;
		esup = other.esup;
		esup_p = other.esup_p;
		if(other.psup != nullptr)
		{
			psup = new vector<int>[npoin];
			for(int i = 0; i < npoin; i++)
				psup[i] = other.psup[i];
		}
		if(other.elsed != nullptr)
		{
			elsed = new vector<int>[nedge];
			for(int i = 0; i < nedge; i++)
				elsed[i] = other.elsed[i];
		}
		esuel = other.esuel;
		intedge = other.intedge;
		elsed = other.elsed;
		intfac = other.intfac;
		alloc_jacobians = other.alloc_jacobians;
		jacobians = other.jacobians;
	}

	~UMesh()
	{
		if(psup != nullptr)
			delete [] psup;
		if(elsed != nullptr)
			delete [] elsed;
	}

	double gcoords(int pointno, int dim)
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
	int gflag_bpoin(int ipoin) { return flag_bpoin(ipoin); }

	void setcoords(Matrix<double>* c)
	{ coords = *c; }

	void setinpoel(Matrix<int>* inp)
	{ inpoel = *inp; }

	void setbface(Matrix<int>* bf)
	{ bface = *bf; }

	Matrix<double>* getcoords()
	{ return &coords; }

	int glpofa(int iface, int ifnode) { return lpofa(iface, ifnode); }
	int gesup(int i) { return esup.get(i); }
	int gesup_p(int i) { return esup_p.get(i); }
	int gpsup(int i, int j) { return psup[i].at(j); }		// get jth surrounding point of ith node
	int gpsupsize(int i) { return psup[i].size(); }
	int gintedge(int iedge, int ipoin) { return intedge.get(iedge,ipoin); }
	int gelsed(int iedge, int ielem) { return elsed[iedge].at(ielem); }
	int gelsedsize(int iedge) { return elsed[iedge].size(); }
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
	int gnedel() { return nedel; }
	int gnbtag() {return nbtag; }
	int gndtag() { return ndtag; }
	int gnedge() { return nedge; }
	int gnbedge() { return nbedge; }

	void readGmsh2(string mfile, int dimensions)
	///< Reads mesh from Gmsh 2 format file. For quadratic meshes, mapping has to be applied for node-ordering.
	{
		cout << "UMesh3d: readGmsh2(): Reading mesh file...\n";
		int dum; double dummy; string dums; char ch;
		ndim = dimensions;

		ifstream infile(mfile);
		for(int i = 0; i < 4; i++)		//skip 4 lines
			do
				ch = infile.get();
			while(ch != '\n');

		infile >> npoin;
		cout << "UMesh3d: readGmsh2(): No. of points = " << npoin << endl;
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
		//cout << "UMesh3d: readGmsh2(): coords read." << endl;

		int nelm, elmtype, nbtags, ntags, nskipped= 0;
		ndtag = 0; nbtag = 0;
		infile >> nelm;
		Matrix<int> elms(nelm,40);
		nface = 0; nelem = 0;
		cout << "UMesh3d: readGmsh2(): Total number of elms is " << nelm << endl;
		Matrix<int> temper(27,1);		// to temporarily store nodes in Gmsh order

		for(int i = 0; i < nelm; i++)
		{
			infile >> dum;
			infile >> elmtype;
			// assumption here is that elm-type is same for all edges, same for all boundary faces and same for all elements
			switch(elmtype)
			{
				case(2): // linear triangle face
					nnofa = 3;
					nnoded = 2;
					nedfa = 3;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(3): // linear quad (4-node) face
					nnofa = 4;
					nnoded = 2;
					nedfa = 4;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(9): // quadratic triangle face
					nnofa = 6;
					nnoded = 3;
					nedfa = 3;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(10): // quadratic quad (9-node) face
					nnofa = 9;
					nnoded = 3;
					nedfa = 4;
					infile >> nbtags;
					if(nbtags > nbtag) nbtag = nbtags;
					for(int j = 0; j < nbtags; j++)
						infile >> elms(i,j+nnofa);		// get tags
					for(int j = 0; j < nnofa; j++)
						infile >> elms(i,j);			// get node numbers
					nface++;
					break;
				case(4): // linear tet
					nnode = 4;
					nfael = 4;
					nnofa = 3;
					nnoded = 2;
					nedel = 6;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(5):	// linear hex
					nnode = 8;
					nfael = 6;
					nnofa = 4;
					nnoded = 2;
					nedel = 12;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(11):	// quadratic tet
					nnode = 10;
					nfael = 4;
					nnofa = 6;
					nnoded = 3;
					nedel = 6;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> elms(i,j);			// get node numbers
					nelem++;
					break;
				case(12):	// quadratic hex (27 nodes)
					// add mapping from gmsh format to Xiaodong format
					nnode = 27;
					nfael = 6;
					nnofa = 9;
					nnoded = 3;
					nedel = 12;
					infile >> ntags;
					if(ntags > ndtag) ndtag = ntags;
					for(int j = 0; j < ntags; j++)
						infile >> elms(i,j+nnode);		// get tags
					for(int j = 0; j < nnode; j++)
						infile >> temper(j);			// get node numbers
					
					for(int j = 0; j <= 8; j++)		// first 8 nodes are same
						elms(i,j) = temper(j);
					
					//cout << "mapping..." << endl;
					elms(i,9) = temper(11);
					elms(i,10) = temper(13);
					elms(i,11) = temper(9);
					elms(i,12) = temper(10);
					elms(i,13) = temper(12);
					elms(i,14) = temper(14);
					elms(i,15) = temper(15);
					elms(i,16) = temper(16);
					elms(i,17) = temper(18);
					elms(i,18) = temper(19);
					elms(i,19) = temper(17);
					elms(i,20) = temper(20);
					elms(i,21) = temper(21);
					elms(i,22) = temper(23);
					elms(i,23) = temper(24);
					elms(i,24) = temper(22);
					elms(i,25) = temper(25);
					elms(i,26) = temper(26);
					
					nelem++;
					break;
					
				default:
					cout << "! UMesh3d: readGmsh2(): Element type not recognized. Skipping." << endl;
					do
						ch = infile.get();
					while(ch != '\n');
					nskipped++;
					break;
			}
		}
		//cout << "UMesh3d: readGmsh2(): Done reading elms" << endl;

		if(nface > 0)
			bface.setup(nface, nnofa+nbtag);
		else cout << "UMesh3d: readGmsh2(): NOTE: There is no data to populate bface!" << endl;

		inpoel.setup(nelem, nnode);
		vol_regions.setup(nelem, ndtag);

		// write into inpoel and bface
		// the first nface rows to be read are boundary faces
		for(int i = 0; i < nface; i++)
		{
			for(int j = 0; j < nnofa; j++)
				bface(i,j) = elms(i+nskipped,j)-1;			// -1 to correct for the fact that our numbering starts from zero
			for(int j = nnofa; j < nnofa+nbtag; j++)
				bface(i,j) = elms(i+nskipped,j);
		}
		for(int i = 0; i < nelem; i++)
		{
			for(int j = 0; j < nnode; j++)
				inpoel(i,j) = elms(i+nface+nskipped,j)-1;
			for(int j = 0; j < ndtag; j++)
				vol_regions(i,j) = elms(i+nface+nskipped,j+nnode);
		}
		infile.close();

		cout << "UMesh3d: readGmsh2(): Setting flag_bpoin..." << endl;

		// set flag_bpoin
		flag_bpoin.setup(npoin,1);
		flag_bpoin.zeros();
		for(int i = 0; i < nface; i++)
			for(int j = 0; j < nnofa; j++)
				flag_bpoin(bface(i,j)) = 1;

		cout << "UMesh3d: readGmsh2(): Done. No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << ",\n number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << ", number of nodes per edge " << nnoded << "." << endl;
	}

	void printmeshstats()
	{
		cout << "UMesh3d: No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << ", number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << endl;
	}

	void mapinpoelXiaodongToGmsh()
	/** Changes node ordering. Use only for hex mesh!! */
	{
		int temp;

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int inode = 0; inode <= 8; inode++)
				m_inpoel(ielem,inode) = inpoel(ielem,inode);

			m_inpoel(ielem,11) = inpoel(ielem,9);
			m_inpoel(ielem,13) = inpoel(ielem,10);
			m_inpoel(ielem,9) = inpoel(ielem,11);
			m_inpoel(ielem,10) = inpoel(ielem,12);
			m_inpoel(ielem,12) = inpoel(ielem,13);
			m_inpoel(ielem,14) = inpoel(ielem,14);
			m_inpoel(ielem,15) = inpoel(ielem,15);
			m_inpoel(ielem,16) = inpoel(ielem,16);
			m_inpoel(ielem,18) = inpoel(ielem,17);
			m_inpoel(ielem,19) = inpoel(ielem,18);
			m_inpoel(ielem,17) = inpoel(ielem,19);
			m_inpoel(ielem,20) = inpoel(ielem,20);
			m_inpoel(ielem,21) = inpoel(ielem,21);
			m_inpoel(ielem,23) = inpoel(ielem,22);
			m_inpoel(ielem,24) = inpoel(ielem,23);
			m_inpoel(ielem,22) = inpoel(ielem,24);
			m_inpoel(ielem,25) = inpoel(ielem,25);
			m_inpoel(ielem,26) = inpoel(ielem,26);
		}
	}

	/** Writes mesh to Gmsh2 file format. */
	void writeGmsh2(string mfile)
	{
		cout << "UMesh2d: writeGmsh2(): writing mesh to file " << mfile << endl;

		m_inpoel.setup(nelem,nnode);
		mapinpoelXiaodongToGmsh();

		// decide element type first, based on nfael/nnode and nnofa
		int elm_type = 4;
		int face_type = 2;
		if(nnode == 8)
			elm_type = 5;
		else if(nnode == 10)
			elm_type = 11;
		else if(nnode == 27)
			elm_type = 12;

		if(nnofa == 4) face_type = 3;
		else if(nnofa == 6) face_type = 9;
		else if(nnofa == 9) face_type = 10;

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

		// "elements"
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

		// cout << "elements\n";
		for(int iel = 0; iel < nelem; iel++)
		{
			outf << nface+iel+1 << " " << elm_type << " " << ndtag;
			for(int i = 0; i < ndtag; i++)
				outf << " " << vol_regions(iel,i);
			for(int i = 0; i < nnode; i++)
				outf << " " << m_inpoel(iel,i)+1;
			outf << '\n';
		}
		outf << "$EndElements\n";

		outf.close();
	}

	/** Computes various connectivity data structures for the mesh.
		These include
		(1) Elements surrounding points (esup and esup_p)
		(2) Points surrounding points (psup)
		(3) Elements surrounding elements (esuel)
		(4) Elements surrounding edge (elsed)
		(5) Edge data structure (intedge)
		(6) Face data structure (intfac)
	*/
	void compute_topological()
	{
		/// NOTE: Currently only works for linear mesh - and psup works only for tetrahedral or hexahedral linear mesh

		cout << "UMesh2d: compute_topological(): Calculating and storing topological information...\n";
		//1. Elements surrounding points
		//cout << "UMesh2d: compute_topological(): Elements surrounding points\n";
		esup_p.setup(npoin+1,1,ROWMAJOR);
		esup_p.zeros();
		
		if(psup != nullptr)
			delete [] psup;
		psup = new vector<int>[npoin];
		for(int i = 0; i < npoin; i++)
			psup[i].reserve(10);

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

		//2. Points surrounding points - works only for tets and hexes!!
		cout << "UMesh2d: compute_topological(): Points surrounding points\n";
		Matrix<int> lpoin(npoin,1);  // The ith member indicates the global point number of which the ith point is a surrounding point
		for(int i = 0; i < npoin; i++)
			lpoin(i,0) = -1;	// initialize this vector to -1

		for(int ip = 0; ip < npoin; ip++)
		{
			//cout << "=" << flush;
			lpoin(ip) = ip;
			// iterate over elements surrounding point
			for(int i = esup_p(ip); i < esup_p(ip+1); i++)
			{
				int ielem = esup(i);
				int inode;

				// find local node number of ip in ielem -- needed for anything except tetrahedral mesh
				for(int jnode = 0; jnode < nnode; jnode++)
					if(inpoel(ielem,jnode) == ip) inode = jnode;

				vector<bool> nbd(nnode);		// contains true if that local node number is connected to inode
				for(int j = 0; j < nnode; j++)
					nbd[j] = false;

				if(nnode == 4)					// for a tet, all local nodes are connected to any given node
				{
					for(int inb = 0; inb < nbd.size(); inb++)
						nbd[inb] = true;
				}
				if(nnode == 8)					// for a hex, hard-code based on local point numbering
				{
					if(inode==0) { nbd[1] = true; nbd[3] = true; nbd[4] = true; }
					else if(inode==1) { nbd[0] = true; nbd[2] = true; nbd[5] = true; }
					else if(inode==2) { nbd[1] = nbd[3] = nbd[6] = true; }
					else if(inode==3) { nbd[2] = nbd[0] = nbd[7] = true; }
					else if(inode==4) { nbd[0] = nbd[5] = nbd[7] = true; }
					else if(inode==5) { nbd[1] = nbd[4] = nbd[6] = true; }
					else if(inode==6) { nbd[2] = nbd[5] = nbd[7] = true; }
					else if(inode==7) { nbd[3] = nbd[6] = nbd[4] = true; }
					else cout << "! UMesh3d: compute_topological(): Error in psup!" << endl;
				}

				for(int jnode = 0; jnode < nnode; jnode++)
				{
					int jpoin = inpoel(ielem, jnode);
					if(lpoin(jpoin) != ip && nbd[jnode] == true)
					{
						psup[ip].push_back(jpoin);
						lpoin(jpoin) = ip;
					}
				}
			}
			for(int i = 0; i < npoin; i++)
				lpoin(i,0) = -1;
		}
		//Points surrounding points is done.

		// 3. calculate number of edges using psup
		//cout << "UMesh3d: compute_topological(): Getting number of edges" << endl;
		nedge = 0; nbedge = 0;

		for(int ipoin = 0; ipoin < npoin; ipoin++)
		{
			lpoin(ipoin) = 1;
			for(int jp = 0; jp < psup[ipoin].size(); jp++)
			{
				int jpoin = psup[ipoin].at(jp);

				if(lpoin(jpoin) != 1)
					nedge++;
			}
		}

		for(int i = 0; i < npoin; i++)
			lpoin(i) = 0;

		cout << "UMesh3d: compute_topological(): Number of edges = " << nedge << endl;

		elsed = new vector<int>[nedge];
		for(int i = 0; i < nedge; i++)
			elsed[i].reserve(8);
		intedge.setup(nedge, nnoded);

		// 4. get intedge
		cout << "UMesh3d: compute_topological(): Calculating intedge" << endl;

		// first, boundary edges
		nbedge = 0;

		for(int ipoin = 0; ipoin < npoin; ipoin++)
		{
			lpoin(ipoin) = 1;
			for(int jp = 0; jp < psup[ipoin].size(); jp++)
			{
				int jpoin = psup[ipoin].at(jp);
				if(lpoin.get(jpoin) != 1 && flag_bpoin.get(ipoin)==1 && flag_bpoin.get(jpoin)==1)
				{
					intedge(nbedge,0) = ipoin;
					intedge(nbedge,1) = jpoin;
					nbedge++;
				}
			}
		}

		for(int i = 0; i < npoin; i++)
			lpoin(i) = 0;

		//cout << "UMesh3d: compute_topological(): Calculating intedge - interior" << endl;
		nedge = nbedge;
		for(int ipoin = 0; ipoin < npoin; ipoin++)
		{
			lpoin(ipoin) = 1;
			for(int jp = 0; jp < psup[ipoin].size(); jp++)
			{
				int jpoin = psup[ipoin].at(jp);
				if(lpoin(jpoin) != 1 && !(flag_bpoin.get(ipoin)==1 && flag_bpoin.get(jpoin)==1) )
				{
					intedge(nedge,0) = ipoin;
					intedge(nedge,1) = jpoin;
					nedge++;
				}
			}
		}

		// 5. Get elsed (elements surrounding each edge) using esup
		cout << "UMesh3d: compute_topological(): Calculating elsed" << endl;
		Matrix<int> lelem(nelem,1);
		int* ip = new int[nnoded];

		for(int ied = 0; ied < nedge; ied++)
		{
			for(int i = 0; i < nnoded; i++)
				ip[i] = intedge(ied,i);

			lelem.zeros();
			for(int iel = esup_p(ip[0]); iel < esup_p(ip[0]+1); iel++)
			{
				lelem(esup(iel)) = 1;
			}

			for(int jel = esup_p(ip[1]+1)-1; jel >= esup_p(ip[1]); jel--)
			{
				if(lelem(esup(jel)) == 1)
					elsed[ied].push_back(esup(jel));
			}
		}

		delete [] ip;

		// 6. Elements surrounding elements
		cout << "UMesh3d: compute_topological(): Elements surrounding elements...\n";

		esuel.setup(nelem, nfael, ROWMAJOR);
		for(int ii = 0; ii < nelem; ii++)
			for(int jj = 0; jj < nfael; jj++)
				esuel(ii,jj) = -1;

		lpofa.setup(nfael, nnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in [0:nnofa], i in [0:nfael])

		if(nnode == 4)								// if tet
			for(int i = 0; i < nfael; i++)
			{
				for(int j = 0; j < nnofa; j++)
				{
					lpofa(i,j) = perm(0,nnode-1,i,j);
				}
			}

		if(nnode == 8)								// if hex
		{
			lpofa(0,0) = 0; lpofa(0,1) = 3; lpofa(0,2) = 2; lpofa(0,3) = 1;
			lpofa(1,0) = 0; lpofa(1,1) = 1; lpofa(1,2) = 5; lpofa(1,3) = 4;
			lpofa(2,0) = 0; lpofa(2,1) = 4; lpofa(2,2) = 7; lpofa(2,3) = 3;
			lpofa(3,0) = 1; lpofa(3,1) = 2; lpofa(3,2) = 6; lpofa(3,3) = 5;
			lpofa(4,0) = 2; lpofa(4,1) = 3; lpofa(4,2) = 7; lpofa(4,3) = 6;
			lpofa(5,0) = 4; lpofa(5,1) = 5; lpofa(5,2) = 6; lpofa(5,3) = 7;
		}

		Matrix<int> lhelp(nnofa,1);
		lhelp.zeros();
		lpoin.zeros();

		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int ifael = 0; ifael < nfael; ifael++)
			{
				for(int i = 0; i < nnofa; i++)
				{
					lhelp(i) = inpoel(ielem, lpofa(ifael,i));		// lhelp stores global node nos. of current face of current element
					lpoin(lhelp(i)) = 1;
				}
				int ipoin = lhelp(0);								// ipoin is the global node number of first point of this face
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
							if(icoun == nnofa)		// if all nnofa points of face ifael are found in face jfael, they both are the same face
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

		/** Computes, for each face, the elements on either side, the starting node and the ending node of the face. This is stored in intfac.
		The orientation of the face is such that the face points towards the element with larger index.
		NOTE: After the following portion, esuel holds (nelem + face no.) for each ghost cell, instead of -1 as before.*/

		//cout << "UMesh2d: compute_topological(): Computing intfac..." << endl;
		nbface = naface = 0;

		// first run: calculate nbface
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nfael; in++)
			{
				int je = esuel(ie,in);
				if(je == -1)
					nbface++;
			}
		}
		cout << "UMesh2d: compute_topological(): Number of boundary faces = " << nbface << endl;
		// calculate number of internal faces
		naface = nbface;
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nfael; in++)
			{
				int je = esuel(ie,in);
				if(je > ie && je < nelem)
					naface++;
			}
		}
		cout << "UMesh2d: compute_topological(): Number of all faces = " << naface << endl;

		//allocate intfac
		intfac.setup(naface,nnofa+2);

		//reset face totals
		nbface = naface = 0;

		//second run: populate intfac
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nfael; in++)
			{
				int* inp = new int[nnofa];
				for(int inofa = 0; inofa < nnofa; inofa++)
					inp[inofa] = lpofa.get(in,inofa);

				int je = esuel.get(ie,in);
				if(je == -1)
				{
					esuel(ie,in) = nelem+nbface;
					intfac(nbface,0) = ie;
					intfac(nbface,1) = nelem+nbface;
					for(int j = 0; j < nnofa; j++)
					{
						intfac(nbface,2+j) = inpoel(ie,inp[j]);
					}

					nbface++;
				}

				delete [] inp;
			}
		}

		naface = nbface;
		for(int ie = 0; ie < nelem; ie++)
		{
			for(int in = 0; in < nfael; in++)
			{
				int* inp = new int[nnofa];
				for(int inofa = 0; inofa < nnofa; inofa++)
					inp[inofa] = lpofa.get(in,inofa);

				int je = esuel(ie,in);
				if(je > ie && je < nelem)
				{
					intfac(naface,0) = ie;
					intfac(naface,1) = je;
					for(int j = 0; j < nnofa; j++)
						intfac(naface,2+j) = inpoel(ie,inp[j]);
					naface++;
				}

				delete [] inp;
			}
		}

		cout << "UMesh3d: compute_topological(): Done." << endl;
	}

	
	/** Creates a UMesh object by adding one extra node at each edge centre, and also, in case of a hex mesh, one extra node at each face centre and cell centre.
	*/
	UMesh convertLinearToQuadratic()
	{
		UMesh q;

		cout << "UMesh3d: convertLinearToQuadratic(): Producing a quadratic mesh from linear mesh..." << endl;
		if(nnofa != 3 && nnofa != 4) {
			cout << "! UMesh2d: convertLinearToQuadratic(): Mesh is not linear or is not supported!!" << endl;
			return q; }

		q.ndim = ndim;
		q.npoin = npoin+nedge+naface+nelem;
		q.nelem = nelem;
		q.nface = nface; q.naface = naface; q.nbface = nbface;
		q.nedfa = nedfa;
		q.nfael = nfael;
		q.nnoded = nnoded + 1;
		q.nedge = nedge;
		q.nbedge = nbedge;

		if(nnode == 8)												// for hex
		{
			q.nnode = nnode+nedel+nfael+1;
			q.nnofa = nnofa+nedfa+1;
		}
		else if(nnode == 4)											// for tet
		{
			q.nnode = nnode+nedel;
			q.nnofa = nnofa+nedfa;
		}
		else  {
			cout << "! UMesh3d: convertLinearToQuadratic(): nnode is neither 4 nor 8 - mesh is not supported!" << endl;
			return q; }

		q.nbtag = nbtag; q.ndtag = ndtag;

		q.coords.setup(q.npoin, q.ndim);
		q.inpoel.setup(q.nelem, q.nnode);
		q.bface.setup(q.nface, q.nnofa+q.nbtag);
		q.vol_regions.setup(q.nelem, q.ndtag);

		// copy nodes, elements and bfaces

		for(int ipoin = 0; ipoin < npoin; ipoin++)
		{
			for(int idim = 0; idim < ndim; idim++)
				q.coords(ipoin, idim) = coords(ipoin,idim);
		}

		for(int iel = 0; iel < nelem; iel++)
		{
			for(int inode = 0; inode < nnode; inode++)
				q.inpoel(iel,inode) = inpoel(iel,inode);
			for(int itag = 0; itag < ndtag; itag++)
				q.vol_regions(iel,itag) = vol_regions(iel,itag);
		}

		for(int iface = 0; iface < nface; iface++)
		{
			for(int j = 0; j < nnofa; j++)
				q.bface(iface,j) = bface(iface,j);
			for(int j = nnofa; j < nnofa+nbtag; j++)
				q.bface(iface,q.nnofa-nnofa+j) = bface(iface,j);
		}

		double* centre = new double[ndim];

		if(nnode == 8)						// for hex mesh, we need face centres and cell centres
		{
			// get the body-centre nodes
			cout << "UMesh3d: convertLinearToQuadratic(): Adding nodes at cell-centres" << endl;

			for(int iel = 0; iel < nelem; iel++)
			{
				for(int i = 0; i < ndim; i++)
					centre[i] = 0;

				for(int inode = 0; inode < nnode; inode++)
				{
					for(int idim = 0; idim < ndim; idim++)
						centre[idim] += coords(inpoel(iel,inode),idim);
				}

				for(int idim = 0; idim < ndim; idim++)
					centre[idim] /= nnode;

				//add centre node to q.coords and update q.inpoel
				for(int idim = 0; idim < ndim; idim++)
					q.coords(npoin+iel,idim) = centre[idim];

				q.inpoel(iel,q.nnode-1) = npoin+iel;
			}

			// face centre nodes
			cout << "UMesh3d: convertLinearToQuadratic(): Adding nodes at face centres." << endl;

			for(int iface = 0; iface < nbface; iface++)
			{
				for(int i = 0; i < ndim; i++)
					centre[i] = 0;

				for(int ifnode = 2; ifnode < 2+nnofa; ifnode++)
					for(int idim = 0; idim < ndim; idim++)
						centre[idim] += coords(intfac(iface,ifnode), idim);
				for(int idim = 0; idim < ndim; idim++)
					centre[idim] /= nnofa;

				// add node to coords
				for(int idim = 0; idim < ndim; idim++)
					q.coords(npoin+nelem+iface, idim) = centre[idim];

				// search for the bface corresponding to this face
				int rbface = -1; bool finmatch;
				bool* match = new bool[nnofa];

				for(int ibface = 0; ibface < nface; ibface++)
				{
					for(int inofa = 0; inofa < nnofa; inofa++)
						match[inofa] = false;

					finmatch = true;

					for(int inofa = 0; inofa < nnofa; inofa++)
					{
						for(int jnofa = 0; jnofa < nnofa; jnofa++)
						{
							if(intfac(iface,inofa+2) == bface(ibface,jnofa)) {	// if ith node of iface matches any node of ibface, flag true
								match[inofa] = true;
								break; }
						}
					}

					/*for(int inofa = 0; inofa < nnofa; inofa++)
						cout << " " << match[inofa];
					cout << endl;*/

					for(int inofa = 0; inofa < nnofa; inofa++)
					{
						if(match[inofa] == false) {					// if any node of iface did not match some node of ibface, there's no match
							finmatch = false;
							break;
						}
					}

					if(finmatch == true) rbface = ibface;
				}

				delete [] match;

				if(rbface == -1) cout << "! UMesh3d: convertLinearToQuadratic(): Bface corresponding to face " << iface << " not found!" << endl;

				// now add node to bface
				q.bface(rbface, q.nnofa-1) = npoin + nelem + iface;

				// --- add point to inpoel ---

				int nodenum = -1;
				int elem = intfac(iface,0);

				bool ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,0)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 20;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// if found, this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,4)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 21;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,4)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 24;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,5)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 22;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,6)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 23;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,4) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,7)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 25;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}
			}

			// next, internal faces
			for(int iface = nbface; iface < naface; iface++)
			{
				for(int i = 0; i < ndim; i++)
					centre[i] = 0;

				for(int ifnode = 2; ifnode < 2+nnofa; ifnode++)
					for(int idim = 0; idim < ndim; idim++)
						centre[idim] += coords(intfac(iface,ifnode), idim);
				for(int idim = 0; idim < ndim; idim++)
					centre[idim] /= nnofa;

				// add node to coords
				for(int idim = 0; idim < ndim; idim++)
					q.coords(npoin+nelem+iface, idim) = centre[idim];

				// add point to inpoel
				int nodenum = -1;
				int elem = intfac(iface,0);

				bool ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,0)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 20;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					//continue;							// if found, this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,4)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 21;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					//continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,4)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 24;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					//continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,5)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 22;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					//continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,6)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 23;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					//continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,4) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,7)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 25;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					//continue;							// this face is done
				}

				// other element
				nodenum = -1;
				elem = intfac(iface,1);

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,0)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 20;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// if found, this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,4)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 21;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,4)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 24;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,5)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 22;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,6)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 23;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}

				ematch = true;
				for(int inode = 0; inode < nnofa; inode++)
				{
					if(!(intfac(iface,2+inode)==inpoel(elem,4) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,7)))
						{ ematch = false; break; }
				}
				if(ematch==true)
				{
					nodenum = 25;
					q.inpoel(elem,nodenum) = npoin + nelem + iface;
					continue;							// this face is done
				}
			}

			// --- next, add points at edge centres ---
			cout << "UMesh3d: convertLinearToQuadratic(): Adding points at edge centres for hexes" << endl;

			// first, boundary edges
			for(int ied = 0; ied < nbedge; ied++)
			{
				for(int i = 0; i < ndim; i++)
					centre[i] = 0;

				for(int ifnode = 0; ifnode < nnoded; ifnode++)
					for(int idim = 0; idim < ndim; idim++)
						centre[idim] += coords(intedge(ied,ifnode), idim);
				for(int idim = 0; idim < ndim; idim++)
					centre[idim] /= nnoded;

				// add node to coords
				int cono = npoin+nelem+naface+ied;
				for(int idim = 0; idim < ndim; idim++)
					q.coords(cono, idim) = centre[idim];

				// add to elements surrounding edge
				//cout << "add to elements surr edge" << endl;
				for(int ielem = 0; ielem < elsed[ied].size(); ielem++)
				{
					int elem = elsed[ied].at(ielem);

					if((intedge.get(ied,0)==inpoel(elem,0)&&intedge.get(ied,1)==inpoel(elem,1)) || (intedge.get(ied,1)==inpoel(elem,0)&&intedge.get(ied,0)==inpoel(elem,1)))
						q.inpoel(elem,8) = cono;
					else if((intedge.get(ied,0)==inpoel(elem,1)&&intedge.get(ied,1)==inpoel(elem,2)) || (intedge.get(ied,1)==inpoel(elem,1)&&intedge.get(ied,0)==inpoel(elem,2)))
						q.inpoel(elem,9) = cono;
					else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,3)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,3)))
						q.inpoel(elem,10) = cono;
					else if((intedge(ied,0)==inpoel(elem,3)&&intedge(ied,1)==inpoel(elem,0)) || (intedge(ied,1)==inpoel(elem,3)&&intedge(ied,0)==inpoel(elem,0)))
						q.inpoel(elem,11) = cono;
					else if((intedge(ied,0)==inpoel(elem,0)&&intedge(ied,1)==inpoel(elem,4)) || (intedge(ied,1)==inpoel(elem,0)&&intedge(ied,0)==inpoel(elem,4)))
						q.inpoel(elem,12) = cono;
					else if((intedge(ied,0)==inpoel(elem,1)&&intedge(ied,1)==inpoel(elem,5)) || (intedge(ied,1)==inpoel(elem,1)&&intedge(ied,0)==inpoel(elem,5)))
						q.inpoel(elem,13) = cono;
					else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,6)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,6)))
						q.inpoel(elem,14) = cono;
					else if((intedge(ied,0)==inpoel(elem,3)&&intedge(ied,1)==inpoel(elem,7)) || (intedge(ied,1)==inpoel(elem,3)&&intedge(ied,0)==inpoel(elem,7)))
						q.inpoel(elem,15) = cono;
					else if((intedge(ied,0)==inpoel(elem,4)&&intedge(ied,1)==inpoel(elem,5)) || (intedge(ied,1)==inpoel(elem,4)&&intedge(ied,0)==inpoel(elem,5)))
						q.inpoel(elem,16) = cono;
					else if((intedge(ied,0)==inpoel(elem,5)&&intedge(ied,1)==inpoel(elem,6)) || (intedge(ied,1)==inpoel(elem,5)&&intedge(ied,0)==inpoel(elem,6)))
						q.inpoel(elem,17) = cono;
					else if((intedge(ied,0)==inpoel(elem,6)&&intedge(ied,1)==inpoel(elem,7)) || (intedge(ied,1)==inpoel(elem,6)&&intedge(ied,0)==inpoel(elem,7)))
						q.inpoel(elem,18) = cono;
					else if((intedge(ied,0)==inpoel(elem,7)&&intedge(ied,1)==inpoel(elem,4)) || (intedge(ied,1)==inpoel(elem,7)&&intedge(ied,0)==inpoel(elem,4)))
						q.inpoel(elem,19) = cono;
				}
	
				//cout << "find bface" << endl;
				// find bfaces that this edge belongs to
				vector<int> edfa;
				bool bmatch1, bmatch2;
				for(int ibface = 0; ibface < nface; ibface++)
				{
					bmatch1 = bmatch2 = false;
					for(int inode = 0; inode < nnofa; inode++)
						if(intedge(ied,0)==bface(ibface,inode))
							bmatch1 = true;
					for(int inode = 0; inode < nnofa; inode++)
						if(intedge(ied,1)==bface(ibface,inode))
							bmatch2 = true;

					if(bmatch1 && bmatch2) edfa.push_back(ibface);
				}

				//cout << "add new point to bfaces" << endl;
				// add new point to the bfaces that were found
				for(int ibf = 0; ibf < edfa.size(); ibf++)
				{
					int ibface = edfa.at(ibf);

					if((intedge.get(ied,0)==bface.get(ibface,0) && intedge.get(ied,1)==bface.get(ibface,1)) || (intedge.get(ied,1)==bface.get(ibface,0) && intedge.get(ied,0)==bface.get(ibface,1)))
						q.bface(ibface,4) = cono;
					else if((intedge(ied,0)==bface(ibface,1) && intedge(ied,1)==bface(ibface,2)) || (intedge(ied,1)==bface(ibface,1) && intedge(ied,0)==bface(ibface,2)))
						q.bface(ibface,5) = cono;
					else if((intedge(ied,0)==bface(ibface,2) && intedge(ied,1)==bface(ibface,3)) || (intedge(ied,1)==bface(ibface,2) && intedge(ied,0)==bface(ibface,3)))
						q.bface(ibface,6) = cono;
					else if((intedge(ied,0)==bface(ibface,3) && intedge(ied,1)==bface(ibface,0)) || (intedge(ied,1)==bface(ibface,3) && intedge(ied,0)==bface(ibface,0)))
						q.bface(ibface,7) = cono;
				}
			}

			// internal edges
			//cout << "internal edges" << endl;
			for(int ied = nbedge; ied < nedge; ied++)
			{
				for(int i = 0; i < ndim; i++)
					centre[i] = 0;

				for(int ifnode = 0; ifnode < nnoded; ifnode++)
					for(int idim = 0; idim < ndim; idim++)
						centre[idim] += coords(intedge(ied,ifnode), idim);
				for(int idim = 0; idim < ndim; idim++)
					centre[idim] /= nnoded;

				// add node to coords
				int cono = npoin+nelem+naface+ied;
				for(int idim = 0; idim < ndim; idim++)
					q.coords(cono, idim) = centre[idim];

				// add to elements surrounding edge
				for(int ielem = 0; ielem < elsed[ied].size(); ielem++)
				{
					int elem = elsed[ied].at(ielem);

					if((intedge(ied,0)==inpoel(elem,0)&&intedge(ied,1)==inpoel(elem,1)) || (intedge(ied,1)==inpoel(elem,0)&&intedge(ied,0)==inpoel(elem,1)))
						q.inpoel(elem,8) = cono;
					else if((intedge(ied,0)==inpoel(elem,1)&&intedge(ied,1)==inpoel(elem,2)) || (intedge(ied,1)==inpoel(elem,1)&&intedge(ied,0)==inpoel(elem,2)))
						q.inpoel(elem,9) = cono;
					else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,3)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,3)))
						q.inpoel(elem,10) = cono;
					else if((intedge(ied,0)==inpoel(elem,3)&&intedge(ied,1)==inpoel(elem,0)) || (intedge(ied,1)==inpoel(elem,3)&&intedge(ied,0)==inpoel(elem,0)))
						q.inpoel(elem,11) = cono;
					else if((intedge(ied,0)==inpoel(elem,0)&&intedge(ied,1)==inpoel(elem,4)) || (intedge(ied,1)==inpoel(elem,0)&&intedge(ied,0)==inpoel(elem,4)))
						q.inpoel(elem,12) = cono;
					else if((intedge(ied,0)==inpoel(elem,1)&&intedge(ied,1)==inpoel(elem,5)) || (intedge(ied,1)==inpoel(elem,1)&&intedge(ied,0)==inpoel(elem,5)))
						q.inpoel(elem,13) = cono;
					else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,6)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,6)))
						q.inpoel(elem,14) = cono;
					else if((intedge(ied,0)==inpoel(elem,3)&&intedge(ied,1)==inpoel(elem,7)) || (intedge(ied,1)==inpoel(elem,3)&&intedge(ied,0)==inpoel(elem,7)))
						q.inpoel(elem,15) = cono;
					else if((intedge(ied,0)==inpoel(elem,4)&&intedge(ied,1)==inpoel(elem,5)) || (intedge(ied,1)==inpoel(elem,4)&&intedge(ied,0)==inpoel(elem,5)))
						q.inpoel(elem,16) = cono;
					else if((intedge(ied,0)==inpoel(elem,5)&&intedge(ied,1)==inpoel(elem,6)) || (intedge(ied,1)==inpoel(elem,5)&&intedge(ied,0)==inpoel(elem,6)))
						q.inpoel(elem,17) = cono;
					else if((intedge(ied,0)==inpoel(elem,6)&&intedge(ied,1)==inpoel(elem,7)) || (intedge(ied,1)==inpoel(elem,6)&&intedge(ied,0)==inpoel(elem,7)))
						q.inpoel(elem,18) = cono;
					else if((intedge(ied,0)==inpoel(elem,7)&&intedge(ied,1)==inpoel(elem,4)) || (intedge(ied,1)==inpoel(elem,7)&&intedge(ied,0)==inpoel(elem,4)))
						q.inpoel(elem,19) = cono;
				}
			}
			delete [] centre;
			return q;
		}		// end if for hex

		// Tets

		// first, boundary edges
		for(int ied = 0; ied < nbedge; ied++)
		{
			for(int i = 0; i < ndim; i++)
				centre[i] = 0;

			for(int ifnode = 0; ifnode < nnoded; ifnode++)
				for(int idim = 0; idim < ndim; idim++)
					centre[idim] += coords(intedge(ied,ifnode), idim);
			for(int idim = 0; idim < ndim; idim++)
				centre[idim] /= nnoded;

			// add node to coords
			int cono = npoin+ied;
			for(int idim = 0; idim < ndim; idim++)
				q.coords(cono, idim) = centre[idim];

			// add to elements surrounding edge NOTE: ordering of nodes is taken from Gmsh docs
			for(int ielem = 0; ielem < elsed[ied].size(); ielem++)
			{
				int elem = elsed[ied].at(ielem);

				if((intedge(ied,0)==inpoel(elem,0)&&intedge(ied,1)==inpoel(elem,1)) || (intedge(ied,1)==inpoel(elem,0)&&intedge(ied,0)==inpoel(elem,1)))
					q.inpoel(elem,4) = cono;
				else if((intedge(ied,0)==inpoel(elem,1)&&intedge(ied,1)==inpoel(elem,2)) || (intedge(ied,1)==inpoel(elem,1)&&intedge(ied,0)==inpoel(elem,2)))
					q.inpoel(elem,5) = cono;
				else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,3)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,3)))
					q.inpoel(elem,8) = cono;
				else if((intedge(ied,0)==inpoel(elem,3)&&intedge(ied,1)==inpoel(elem,0)) || (intedge(ied,1)==inpoel(elem,3)&&intedge(ied,0)==inpoel(elem,0)))
					q.inpoel(elem,7) = cono;
				else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,0)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,0)))
					q.inpoel(elem,6) = cono;
				else if((intedge(ied,0)==inpoel(elem,1)&&intedge(ied,1)==inpoel(elem,3)) || (intedge(ied,1)==inpoel(elem,1)&&intedge(ied,0)==inpoel(elem,3)))
					q.inpoel(elem,9) = cono;
			}

			// find bfaces that this edge belongs to
			vector<int> edfa;
			bool bmatch1, bmatch2;
			for(int ibface = 0; ibface < nface; ibface++)
			{
				bmatch1 = bmatch2 = false;

				for(int inode = 0; inode < nnofa; inode++)
					if(intedge(ied,0)==bface(ibface,inode))
						bmatch1 = true;
				for(int inode = 0; inode < nnofa; inode++)
					if(intedge(ied,1)==bface(ibface,inode))
						bmatch2 = true;

				if(bmatch1 && bmatch2) edfa.push_back(ibface);
			}

			// add new point to found bfaces
			for(int ibf = 0; ibf < edfa.size(); ibf++)
			{
				int ibface = edfa.at(ibf);

				if((intedge(ied,0)==bface(ibface,0) && intedge(ied,1)==bface(ibface,1)) || (intedge(ied,1)==bface(ibface,0) && intedge(ied,0)==bface(ibface,1)))
					q.bface(ibface,3) = cono;
				else if((intedge(ied,0)==bface(ibface,1) && intedge(ied,1)==bface(ibface,2)) || (intedge(ied,1)==bface(ibface,1) && intedge(ied,0)==bface(ibface,2)))
					q.bface(ibface,4) = cono;
				else if((intedge(ied,0)==bface(ibface,2) && intedge(ied,1)==bface(ibface,0)) || (intedge(ied,1)==bface(ibface,2) && intedge(ied,0)==bface(ibface,0)))
					q.bface(ibface,5) = cono;
			}
		}

		// internal edges
		for(int ied = nbedge; ied < nedge; ied++)
		{
			for(int i = 0; i < ndim; i++)
				centre[i] = 0;

			for(int ifnode = 0; ifnode < nnoded; ifnode++)
				for(int idim = 0; idim < ndim; idim++)
					centre[idim] += coords(intedge(ied,ifnode), idim);
			for(int idim = 0; idim < ndim; idim++)
				centre[idim] /= nnoded;

			// add node to coords
			int cono = npoin+ied;
			for(int idim = 0; idim < ndim; idim++)
				q.coords(cono, idim) = centre[idim];

			// add to elements surrounding edge NOTE: ordering of nodes is taken from Gmsh docs
			for(int ielem = 0; ielem < elsed[ied].size(); ielem++)
			{
				int elem = elsed[ied].at(ielem);

				if((intedge(ied,0)==inpoel(elem,0)&&intedge(ied,1)==inpoel(elem,1)) || (intedge(ied,1)==inpoel(elem,0)&&intedge(ied,0)==inpoel(elem,1)))
					q.inpoel(elem,4) = cono;
				else if((intedge(ied,0)==inpoel(elem,1)&&intedge(ied,1)==inpoel(elem,2)) || (intedge(ied,1)==inpoel(elem,1)&&intedge(ied,0)==inpoel(elem,2)))
					q.inpoel(elem,5) = cono;
				else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,3)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,3)))
					q.inpoel(elem,8) = cono;
				else if((intedge(ied,0)==inpoel(elem,3)&&intedge(ied,1)==inpoel(elem,0)) || (intedge(ied,1)==inpoel(elem,3)&&intedge(ied,0)==inpoel(elem,0)))
					q.inpoel(elem,7) = cono;
				else if((intedge(ied,0)==inpoel(elem,2)&&intedge(ied,1)==inpoel(elem,0)) || (intedge(ied,1)==inpoel(elem,2)&&intedge(ied,0)==inpoel(elem,0)))
					q.inpoel(elem,6) = cono;
				else if((intedge(ied,0)==inpoel(elem,1)&&intedge(ied,1)==inpoel(elem,3)) || (intedge(ied,1)==inpoel(elem,1)&&intedge(ied,0)==inpoel(elem,3)))
					q.inpoel(elem,9) = cono;
			}
		}

		delete [] centre;
		return q;
	}
};

}	// end namespace acfd
