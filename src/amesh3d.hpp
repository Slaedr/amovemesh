
/// @author Aditya Kashi
/// @date August 20, 2015

#ifndef __AMESH3D_H

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

#ifndef __ADATASTRUCTURES_H
#include "adatastructures.hpp"
#endif

#define __AMESH3D_H

/**
 * \namespace amc
 * \brief Namespace for all identifiers in AMoCurve, except for arrays/matrices and linear algebra.
 * 
 * Includes classes for meshes, mesh movement, geometry reconstruction and more.
 */

namespace amc {

/// Data structure and setup for 3D unstructured mesh
class UMesh
{
private:
	amc_int npoin;
	amc_int nelem;
	amc_int nface;
	amc_int nedge;		///< Number of edges in the mesh
	int ndim;
	int nnode;			///< number of nodes to an element
	int nfael;			///< number of faces to an element
	int nedel;			///< number of edges in an element
	int nnofa;			///< number of nodes in a face -- needs to be generalized in case of general grids
	int nedfa;			///< number of edges in a face
	int nnoded;			///< number of nodes per edge
	amc_int naface;		///< total number of (internal and boundary) faces
	amc_int nbface;		///< number of boundary faces as calculated by compute_face_data(), as opposed to nface which is read from file
	amc_int nbedge;		///< number of boundary edges
	int nbtag;			///< Number of boundary markers for each boundary face
	int ndtag;			///< Number of region markers for each element
	amc_int nbpoin;		///< Number of boundary points
	amat::Matrix<amc_real> coords;
	amat::Matrix<amc_int> inpoel;
	amat::Matrix<amc_int> m_inpoel;			// same as inpoel, except that it might contain different node ordering
	amat::Matrix<amc_int> bface;
	amat::Matrix<amc_int> flag_bpoin;		///< a boolean flag for each point. Contains 1 if the corresponding point is a boundary point
	amat::Matrix<amc_int> vol_regions;		///< to hold volume region markers, if any
	bool alloc_jacobians;
	amat::Matrix<amc_real> jacobians;

	amat::Matrix<amc_real> lpofa;		///< for each face of an element, it contains local node numbers of the nodes forming that face. Assumed to be same for all elements.
	amat::Matrix<amc_int> esup;			///< elements surrounding points
	amat::Matrix<int> esup_p;			///< array containing index of esup where elements surrounding a certain point start
	std::vector<int>* psup;				///< points surrounding points
	amat::Matrix<amc_int> esuel;		///< elements surrounding elements
	amat::Matrix<amc_int> edgepo;		///< edge data structure. Stores the indices of the points making up the edge.
	std::vector<int>* elsed;			///< elements surrounding edge; note that ordering of edges is same as in [edgepo](@ref edgepo)
	amat::Matrix<amc_int> intfac;		///< face data strcture
	
	amat::Matrix<amc_int> bpoints;		///< an ordering of the boundary points, containing corresponding node numbers in coords
	amat::Matrix<amc_int> bpointsinv;	///< inverse of bpoints; contains bpoint number for each node which is a boundary node, and -1 for interior nodes
	amat::Matrix<amc_int> bfacebp;		///< stores [bpoint](@ref bpoints) numbers of nodes of each [boundary face](@ref bface)
	amat::Matrix<amc_int> bfsubp;		///< boundary faces surrounding boundary point
	amat::Matrix<amc_int> bfsubp_p;		///< contains pointers into [bfsubp](@ref bfsubp) for each boundary point
	amat::Matrix<amc_int> bfsubf;		///< boundary faces surrounding boundary face
	amat::Matrix<amc_int> bpsubp;		///< boundary points surrounding boundary point
	amat::Matrix<amc_int> bpsubp_p;		///< contains pointers into [bpsubp](@ref bpsubp) for each boundary point
	amat::Matrix<amc_int> intbedge;		/**< boundary edge - boundary face connectivity. Ordering of edges is same as that of [edgepo](@ref edgepo). 
											NOTE: The edge may not always point from smaller index cell to larger index cell! 
											The edge direction may not be consistent. */
	amat::Matrix<amc_int> bfaceHostCell;		///< Stores the host element number for each bface, computed (inefficiently) by findBfaceHostCell()

public:

	/** No-arg constructor. */
	UMesh();

	UMesh(const UMesh& other);

	~UMesh();

	/// Access point coordinates
	amc_real gcoords(amc_int pointno, int dim) const
	{
		return coords.get(pointno,dim);
	}
	int ginpoel(amc_int elemno, int locnode) const
	{
		return inpoel.get(elemno, locnode);
	}
	int gbface(amc_int faceno, int val) const
	{
		return bface.get(faceno, val);
	}
	int gflag_bpoin(amc_int ipoin) const { return flag_bpoin.get(ipoin); }

	amc_int gbpoints(amc_int ipoin) const { return bpoints.get(ipoin); }
	amc_int gbpointsinv(amc_int ipoin) const { return bpointsinv.get(ipoin); }
	//amc_int gbedge(amc_int ibedge, int index) const { return bedge.get(ibedge, index); }

	/// set coordinates of a certain point; 'set' counterpart of the 'get' function [gcoords](@ref gcoords).
	void scoords(const amc_int pointno, const int dim, const amc_real value)
	{
		coords(pointno,dim) = value;
	}

	void setcoords(amat::Matrix<amc_real>* c)
	{ coords = *c; }

	void setinpoel(amat::Matrix<amc_int>* inp)
	{ inpoel = *inp; }

	void setbface(amat::Matrix<amc_int>* bf)
	{ bface = *bf; }

	const amat::Matrix<amc_real>* getcoords() const
	{ return &coords; }

	int glpofa(amc_int iface, int ifnode) const { return lpofa.get(iface, ifnode); }
	amc_int gesup(amc_int i) const { return esup.get(i); }
	amc_int gesup_p(amc_int i) const { return esup_p.get(i); }
	amc_int gpsup(amc_int i, int j) const { return psup[i].at(j); }		// get jth surrounding point of ith node
	amc_int gpsupsize(amc_int i) const { return psup[i].size(); }
	amc_int gedgepo(amc_int iedge, int ipoin) const { return edgepo.get(iedge,ipoin); }
	amc_int gelsed(amc_int iedge, int ielem) const { return elsed[iedge].at(ielem); }
	amc_int gelsedsize(amc_int iedge) const { return elsed[iedge].size(); }
	amc_int gesuel(amc_int ielem, int jnode) const { return esuel.get(ielem, jnode); }
	amc_int gintfac(amc_int face, int i) const { return intfac.get(face,i); }
	amc_int gbfsubp_p(amc_int i) const { return bfsubp_p.get(i); }
	amc_int gbfsubp(amc_int i) const { return bfsubp.get(i); }
	amc_int gbfsubf(amc_int iface, int isurr) const { return bfsubf.get(iface,isurr); }
	amc_int gbpsubp_p(amc_int i) const { return bpsubp_p.get(i); }
	amc_int gbpsubp(amc_int i) const { return bpsubp.get(i); }
	amc_int gintbedge(amc_int ied, int j) const { return intbedge(ied,j); }
	amc_real gjacobians(amc_int ielem) const { return jacobians.get(ielem,0); }

	amc_int gnpoin() const { return npoin; }
	amc_int gnelem() const { return nelem; }
	amc_int gnface() const { return nface; }
	amc_int gnbface() const { return nbface; }
	amc_int gnbpoin() const { return nbpoin; }
	int gnnode() const { return nnode; }
	int gndim() const { return ndim; }
	amc_int gnaface() const {return naface; }
	int gnfael() const { return nfael; }
	int gnnofa() const { return nnofa; }
	int gnnoded() const { return nnoded; }
	int gnedel() const { return nedel; }
	int gnbtag() const {return nbtag; }
	int gndtag() const { return ndtag; }
	amc_int gnedge() const { return nedge; }
	amc_int gnbedge() const { return nbedge; }

	/// Reads mesh from Gmsh 2 format file. 
	/** For hexahedral quadratic meshes, mapping has to be applied for node-ordering, as the inpoel node-ordering here is same as rDGFlo's.
	 */
	void readGmsh2(std::string mfile, int dimensions);

	/// Reads a special rDGFLO domn format containing only interconnectivity matrix and point coordinates
	/** Computes boundary face data (bface) using intfac, for which esup, esuel and intfac are computed.
	 * An example of the format is given below.
	 *
	 *     nelem   npoin   nnode nedel nfael nnofa nedfa nnoded       time           itime  (except for time and itime, everything else is required)
	 *     175053   31630     4     6     4     3     3      2         0.0000E+00         0
	 *     intmat: nodal points corresponding to each element 
	 *     1        59      1067      1126      3928  1
	 *     2        59      1126      1185      3928  1
	 *     .
	 *     .
	 *     coordinates of points (necessary 1 line between last row of interconnectivity matrix and point coordinates)
	 *     .
	 *     .
	 */
	void readDomn(std::string mfile);

	void printmeshstats();

	/// Changes node ordering. Use only for hex mesh!!
	void mapinpoelRDGFloToGmsh();

	/// lpofal contains, on output, local node numbers corresponding to nodes of the local faces of an element
	void findLocalFaceLocalPointConnectivityLinearElements(amat::Matrix<int>& lpofal);
	
	/// Finds the host element for each boundary face and stores in bfaceHostCell
	void findBfaceHostCell();

	/// Writes mesh to Gmsh2 file format
	void writeGmsh2(std::string mfile);

	/// Writes mesh to regular rDGFlo format domn file
	/** \param[in] mfile is the name of the file to be created
	 * \param[in] farfieldnumber contains the physical face numbers for farfield boundaries.
	 * \param[in] symmetrynumber contains physical face numbers of symmetry boundaries
	 * \param[in] wallnumber contains physical face numbers of the adiabatic walls
	 *
	 * The format is given as a Fortran snippet below.
	 *
	 * 		write(10,*)'npoin,ntetr,npyra,npris,nhexa,ntria,nquad,time'
	 * 		write(10,'(7i9,e12.5e2)') npoin,0,0,0,nhexa,0,nquad,0.0  
	 * 		
	 * 		write(10,*)'element connectivity'
	 * 		do i=1,nhexa
	 * 			write(10,'(i9,27i8)') i,hex(i,1:27)
	 * 		enddo
	 * 		
	 * 		write(10,*)' face connectivity'
	 * 		do i=1,nquad
	 * 			write(10,'(15i8)') i,quad(i,10),quad(i,12),0,0,quad(i,11),bfact(i,1:9) ! face-num, boun-marker, farfield(1 or 0), 0,0, host-element, bface points
	 * 		enddo
	 *
	 * 		write(10,*)' coordinates of the nodes'
	 * 		do i=1,nnode
	 * 			write(10,'(i8,3e20.10e2)') i,(xyz(i,j),j=1,3)
	 * 		end do
	 */
	void writeDomn(std::string mfile, std::vector<int> farfieldnumber, std::vector<int> symmetrynumber, std::vector<int> wallnumber);
	
	/// Computes jacobians for linear elements
	/** Currently only for tetrahedra
	 */
	void compute_jacobians();

	/** \brief Computes various connectivity data structures for the mesh.
	 *
	 * These include
	 * - Elements surrounding points (esup and esup_p)
	 * - Points surrounding points (psup)
	 * - Elements surrounding elements (esuel)
	 * - Elements surrounding edge (elsed)
	 * - Edge data structure (edgepo)
	 * - Face data structure (intfac)
	 * 
	 * \note NOTE: Currently only works for linear mesh - and psup works only for tetrahedral or hexahedral linear mesh
	 */
	void compute_topological();

	/// Computes topological properties of the boundary (surface) mesh
	/** Currently only works for linear mesh.
	 * - Boundary points (bpoints and bpointsinv)
	 * - Boundary faces surrounding boudnary point (bfsubp and bfsubp_p)
	 * - Boundary points surrounding boundary point
	 * - Boundary faces surrounding boundary face (bfsubf)
	 * - Boundary edge connectivity (intbedge) giving the two boundary faces adjoining each edge and global node numbers of nodes making up each edge
	 */
	void compute_boundary_topological();
	
	/// Creates a UMesh object by adding one extra node at each edge centre, and also, in case of a hex mesh, one extra node at each face centre and cell centre.
	/** \note Only works for tetrahedral and hexahedral elements.
	 */
	UMesh convertLinearToQuadratic();
};


}	// end namespace
#endif
