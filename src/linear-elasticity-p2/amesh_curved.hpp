// Data structure and setup for 2D unstructured triangular mesh.
// Aditya Kashi
// Feb 5, 2015

// Feb 26, 2015: Adding mesh-movement function movemesh2d_bs().
// June 2, 2015: Modifying for high-order mesh

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
#include <amatrix2.hpp>
#endif
#ifndef __AMESH2D_H
#include <amesh_c.hpp>
#endif
#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif
#ifndef __ADATASTRUCTURES_H
#include <adatastructures.hpp>
#endif

#define __AMESH2DCURVED_H

using namespace std;
using namespace amat;

namespace acfd {

// Unstructured triangular mesh class
class UTriMeshCurved
{
private:
	int npoin;		// number of real points
	int nelem;
	int nface;		// number of boundary faces
	int ndim;
	int nnode;		// number of nodes to an element
	int nfael;		// number of faces to an element (equal to number of edges to an element in 2D)
	int nnofa;		// number of node in a face -- needs to be generalized in case of general grids
	int naface;		// total number of (internal and boundary) faces
	int nbface;		// number of boundary faces as calculated by compute_face_data(), as opposed to nface which is read from file
	int naux;		// number of auxiliary (high-order) points
	Matrix<double> coords;
	Matrix<int> inpoel;
	Matrix<int> bface;

	// The following 4 are for moving mesh
	Matrix<int> bflag;		// for the basicspring one - useless
	Matrix<int> bflag2;
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

	UTriMeshCurved() {alloc_jacobians = false;}

	void readGmsh(ifstream& mfile)
	{

	}

	void readDomn(ifstream& infile)
	{
		// Do file handling here to populate npoin and nelem
		cout << "UTriMeshCurved: Reading mesh file...";
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

		cout << "\nUTriMeshCurved: Number of elements: " << nelem << ", number of points: " << npoin << ", number of nodes per element: " << nnode << endl;
		cout << "Number of boundary faces: " << nface << ", Number of dimensions: " << ndim;
		nfael = 3;	// number of faces per element
		nnofa = 3;	// number of nodes per face

		//cout << "\nUTriMeshCurved: Allocating coords..";
		coords.setup(npoin, ndim, ROWMAJOR);
		//cout << "\nUTriMeshCurved: Allocating inpoel..\n";
		inpoel.setup(nelem, nnode, ROWMAJOR);
		//cout << "UTriMeshCurved: Allocating bface...\n";
		bface.setup(nface, nnofa + n_extra_fields_in_bface, ROWMAJOR);

		//cout << "UTriMeshCurved: Allocation done.";

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
		cout << "\nUTriMeshCurved: Populated inpoel.";

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
		cout << "\nUTriMeshCurved: Populated coords.\n";
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
			for(int j = 0; j < nnofa + n_extra_fields_in_bface; j++)
			{
				infile >> bface(i,j);
			}
			if (i==nface-1) break;
			do
				ch = infile.get();
			while(ch!='\n');
		}
		cout << "UTriMeshCurved: Populated bface. Done reading mesh.\n";
		//correct first nnofa columns of bface
		for(int i = 0; i < nface; i++)
			for(int j = 0; j < nnofa; j++)
				bface(i,j)--;
	}

	void generateFromRefinedLinearMesh(UTriMesh* m)
	{
		cout << "UTriMeshCurved: Setting up..\n";
		npoin = m->gnpoin();
		nelem = m->gnelem()/4;
		nface = m->gnface()/2;
		nnode = m->gnnode()*2;
		nnofa = m->gnnofa()+1;		// 3 nodes in a face now
		ndim = m->gndim();
		nfael = m->gnfael();
		naux = m->gnaux();
		coords.setup(npoin,ndim);
		inpoel.setup(nelem,nnode);
		bface.setup(nface,nnofa+n_extra_fields_in_bface);

		cout << "UTriMeshCurved: Copying coords..\n";
		for(int ip = 0; ip < npoin; ip++)
		{
			for(int j = 0; j < ndim; j++)
				coords(ip,j) = m->gcoords(ip,j);
		}

		cout << "UTriMeshCurved: Merging elements to create 2nd order elements\n";

		int el[nfael];		// holds element nos of surrounding elements
		for(int ielem = 0; ielem < nelem; ielem++)
		{
			for(int j = 0; j < nfael; j++)
				el[j] = m->gesuel(ielem,j);
			for(int i = 0; i < m->gnnode(); i++)				// populate auxiliary nodes first, in positions 4,5,6 (ie, 3,4,5 in 0-numering)
				inpoel(ielem,i+3) = m->ginpoel(ielem,i);

			//new node 0 will be opposite to old node 1:
			for(int inode = 0; inode < m->gnnode(); inode++)
			{
				//select the node not shared by old ielem element
				if(m->ginpoel(el[1],inode) != inpoel(ielem,3) && m->ginpoel(el[1],inode) != inpoel(ielem,5))
				{
					inpoel(ielem,0) = m->ginpoel(el[1],inode);
					break;
				}
			}
			//new node 1 will be opposite to old node 2
			for(int inode = 0; inode < m->gnnode(); inode++)
			{
				//select the node not shared by old ielem element
				if(m->ginpoel(el[2],inode) != inpoel(ielem,3) && m->ginpoel(el[2],inode) != inpoel(ielem,4))
				{
					inpoel(ielem,1) = m->ginpoel(el[2],inode);
					break;
				}
			}
			//new node 2 will be opposite to old node 0
			for(int inode = 0; inode < m->gnnode(); inode++)
			{
				//select the node not shared by old ielem element
				if(m->ginpoel(el[0],inode) != inpoel(ielem,5) && m->ginpoel(el[0],inode) != inpoel(ielem,4))
				{
					inpoel(ielem,2) = m->ginpoel(el[0],inode);
					break;
				}
			}
		}

		cout << "UTriMeshCurved: Populating bface\n";
		for(int iface = 0; iface < nface; iface++)		// cycle over first nface elements of m->bface
		{
			bface(iface,0) = m->gbface(iface,0);
			bface(iface,1) = m->gbface(iface,1);
			bface(iface,2) = m->gbface(nface+iface,1);
			if(m->gbface(iface,1) != m->gbface(nface+iface,0)) cout << "UTriMeshCurved: !! Error: Boundary node mismatch!\n";

			//now copy over boundary flags from one of the two old boundary faces
			for(int i = 2; i < 2+n_extra_fields_in_bface; i++)
			{
				if(m->gbface(iface,i) != m->gbface(nface+iface,i)) cout << "UTriMeshCurved: ! Warning: Boundary flags on high order boundary face are not unique!!\n";
			}
			for(int i = 3; i < 3+n_extra_fields_in_bface; i++)
			{
				bface(iface,i) = m->gbface(iface,i-1);
			}
		}
		cout << "UTriMeshCurved: No. of nodes: " << npoin << ", number of elements: " << nelem << ", number of boundary faces: " << nface << endl;
		//inpoel.mprint();
		//bface.mprint();
	}

	UTriMeshCurved(UTriMeshCurved& other)
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
	}

	UTriMeshCurved& operator=(UTriMeshCurved& other)
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

	void writeQuadraticGmsh2(string mfile)
	{
		ofstream outf(mfile);

		outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
		outf << "$Nodes\n" << npoin << '\n';
		for(int ip = 0; ip < npoin; ip++)
		{
			outf << ip+1 << " " << coords(ip,0) << " " << coords(ip,1) << " " << 0 << '\n';
		}
		outf << "$Elements\n" << nelem+nface << '\n';
		// boundary faces first
		for(int iface = 0; iface < nface; iface++)
		{
			outf << iface+1 << " 8 2 0 1";
			for(int i = 0; i < nnofa; i++)
				outf << " " << bface(iface,i)+1;
			outf << '\n';
		}
		for(int iel = 0; iel < nelem; iel++)
		{
			outf << nface+iel+1 << " 9 2 0 2";
			for(int i = 0; i < nnode; i++)
				outf << " " << inpoel(iel,i)+1;
			outf << '\n';
		}
		outf << "$EndElements\n";

		outf.close();
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

	void movemesh(Matrix<double> disp)
	// Updates coords of mesh by adding displacements to them; the displacements specified by the 2*npoin-by-1 vector disp, which contains x-displacements in the first npoin entries and y-displacements in the rest.
	{
		for(int ip = 0; ip < npoin; ip++)
		{
			coords(ip,0) += disp(ip);
			coords(ip,1) += disp(npoin+ip);
		}
	}

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
		cout << "UTriMesh: compute_face_data(): Storing boundary flags in gallfa...\n";
		for(int ied = 0; ied < nbface; ied++)
		{
			int p1 = intfac(ied,2);
			int p2 = intfac(ied,3);
			// Assumption: order of nodes of boundary faces is such that normal points outside, when normal is calculated as
			//nx = y2 - y1, ny = -(x2-x1).
			if(nbface != nface) { cout << "UTriMeshCurved: Calculation of number of boundary faces is wrong!\n"; break; }
			for(int i = 0; i < nface; i++)
			{
				if(bface(i,0) == p1 || bface(i,1) == p1)
				{
					if(bface(i,1) == p2 || bface(i,0) == p2)
					{
						gallfa(ied,3) = bface(i,2);
						gallfa(ied,4) = bface(i,3);
					}
				}
			}
		}

		cout << "UTriMeshCurved: compute_face_data(): Done.\n";
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

	int gintfac(int face, int i) { return intfac.get(face,i); }
	double ggallfa(int elem, int i) {return gallfa.get(elem,i); }
};

} // end namespace acfd
