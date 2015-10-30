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
#include "amesh_c.hpp"
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

	UTriMeshCurved() {alloc_jacobians = false;}
	/* UTriMesh()
	{
		npoin = 0; nelem = 0; nface = 0; ndim = 0; nnode = 0;
		coords.setup(1,1,ROWMAJOR);
		inpoel.setup(1,1,ROWMAJOR);
		bface.setup(1,1,ROWMAJOR);
	} */

	UTriMeshCurved(UTriMesh* m)
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
			if(nbface != nface) { cout << "UTriMesh: Calculation of number of boundary faces is wrong!\n"; break; }
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

	void refine_mesh()
	{
		// introduce one node on each edge
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

	void movemesh_farhat(Matrix<double>* xb, Matrix<double>* yb, string solver, double tol=1e-6, int maxiter=1000)
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
};

} // end namespace acfd
