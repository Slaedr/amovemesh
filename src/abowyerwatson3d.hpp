/**
@file abowyerwatson3d.hpp
@brief This file contains a class that implements the Bowyer-Watson algorithm for Delaunay tesselation in 3D.
@author Aditya Kashi
@date November 3, 2015
*/

/**
 * \class Delaunay3d
 * \brief A class for Delaunay tetrahedralization of a given set of points based on Bowyer-Watson algorithm.
 *
 * It has also been referenced from Wikipedia's Bowyer-Watson algorithm page.
 *
 * Notes:
 * \todo Change points to a pointer to const Matrix rather than a matrix itself, in the interest of efficiency.
 *
 * Currently uses a std::std::vector to store elements, faces etc. 
 * \todo Change the deletion and insertion of faces and elements such that we're not moving entire lists (easier said than done).
 *
 * It might be better to use a std::list or std::forward_list instead.
 * \note A graph data structure should also be seriously considered for storing elements, faces, bad elements and the void polygon.
 */

#ifndef __ABOWYERWATSON3D_H

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
#ifndef __AMATRIX_H
#include <amatrix.hpp>
#endif

#define __ABOWYERWATSON3D_H 1

namespace amc {

typedef std::vector<double> Point;

/// Structure for triangular face.
struct Face
{
	int p[3];			///< indices of endpoints of the face
	int elem[2];		///< elements on either side of the face
	//int lfel[2];		///< the local face numbers of this face corresponding to the two elements in elem
};

/// Class representing a tetrahedron.
class Tet
{
public:	
	int p[4];			///< Indices of vertices.
	Point centre;		///< Coords of circumcenter of the tet.
	int surr[4];		///< Indices of surrounding tets. Note that the neighbor corresponding to surr[3] is opposite the vertex p[3].
	double D;			///< 6*volume of tet.
	double radius;		///< Square of radius of circumcircle of tet.

	Tet() {
		centre.resize(3,0.0);
	}

	Tet(const Tet& other) {
		for(int i = 0; i < 4; i++)
		{
			p[i] = other.p[i];
			surr[i] = other.surr[i];
		}
		centre = other.centre;
		D = other.D;
		radius = other.radius;
	}
};

/// Intended to encapsulate data required by 'walk-through' algorithms.
/**  Not needed for mesh generation, but needed, for instance, in independent application of the walk-through subroutine find_containing_tet_and_barycentric_coords().
*/
struct Walkdata
{
	int elem;
	double areacoords[4];
};

/// "Delaunay kernel" for generating the tetrahedral Delaunay tessellation of the convex hull of a set of points.
class Delaunay3d
{
	int cap;
	int badcap;
	double tol;
	int nnode;
	int ndim;

	/// Value used to provide tolerance for the Delaunay criterion.
	/** This value is multiplied by machine epsilon (approx. 2e-16), and the resulting value is used as a tolerance for the Delaunay criterion in \ref bowyer_watson .
	 */
	double zero_scale;

	/// Stores the face-point relationship of a tetrahedron.
	/** Row i contains the local node numbers of the face opposite to node i. The local node numbers are ordered so that the face points outwards. */
	amat::Matrix<int> lpofa;
	void setlpofa();

public:
	amat::Matrix<double> points;
	std::vector<Point> nodes;		///< List of nodes in the Delaunay graph.
	std::vector<Tet> elems;			///< List of all elements ([tetrahedra](@ref Tet)) in the Delaunay graph.
	std::vector<int> badelems;		///< Collection of 'bad elements', that are to be removed while adding a point; its are integers that index members index [elems](@ref elems).
	std::vector<Face> faces;		///< List of all [faces](@ref Face) in the Delaunay graph.
	std::vector<int> voidpoly;		///< Collection of faces that bounds the void obtained after removing bad elements while adding a point; its members are integers that index [faces](@ref faces).
	amat::Matrix<double> jacobians;

	int npoints;

	Delaunay3d();
	Delaunay3d(amat::Matrix<double>* _points);
	Delaunay3d(const Delaunay3d& other);
	Delaunay3d& operator=(const Delaunay3d& other);
	void setup(amat::Matrix<double>* _points);
	
	double l2norm(const std::vector<double>& a) const;
	
	double dot(const std::vector<double>& a, const std::vector<double>& b) const;
	
	/// Computes the jacobian of a tet formed from a point and a face of a tetrahedron.
	double det4(const int ielem, const int i, const std::vector<double>& r) const;
	
	void compute_jacobian(Tet& elem);
	
	void cross_product3(std::vector<double>& c, const std::vector<double>& a, const std::vector<double>& b) const;
	
	/// Computes circumcentre and square of circumradius of a tet
	void compute_circumsphere(Tet& elem);

	/// Computes circumcenter and square of circumradius according to Dr Luo's method
	void compute_circumsphere_contra(Tet& elem);
	
	/// Returs the Jacobian (6 * volume) of a tetrahedron formed by 4 points taken as arguments
	double tetvol(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) const;
	
	/// Unused
	int find_containing_tet_old(const std::vector<double>& r, const int startelement) const;

	/// Locates the Delaunay graph (DG) element which contains the input point (improved).
	/** \param xx is the point which needs to be located in the DG.
	 * \param startelement is the index of the element from which to start the search.
	 *
	 * For each DG element encountered, the 4 tetrahedra formed by the point and each of the 4 faces of the current element are considered.
	 * If the ratios of the Jacobians of all 4 of these new tetrahedra to the Jacobian of the current DG tetrahedron are positive, the current element is the containing element.
	 * The minimum of the 4 ratios is found. 
	 * If the minimum is negative, then the new element to be checked is taken as the DG element adjacent to the face (of the current element) corresponding to the minimum value.
	 */
	int find_containing_tet(const std::vector<double>& r, const int startelement) const;
	
	/// Returns the local face number (number of node opposite to the face) of elem's face given by the second argument.
	/**  If no faces of elem match, it returns -1.
	 */
	int check_face_tet(const Tet& elem, const Face& face) const;

	/// Computes the Delaunay triangulation (tetrahedralization, in this case).
	/** We first scale the point set by the largest coordinate magnitudes in each direction, so as to non-dimensionalize it.
	 */
	void bowyer_watson();
	
	void clear();					///< Reset the Delaunay3d object, except for input data
	
	/// Writes the Delaunay graph to a Gmsh file.
	void writeGmsh2(const std::string mfile) const;
	
	/// Finds the DG element containing a given point and return the area coordinates in that element
	Walkdata find_containing_tet_and_barycentric_coords(const std::vector<double>& rr, const int startelement) const;
	
	/// Computes the jacobian of all elements in the triangulation using cross products
	void compute_jacobians();
	
	bool detect_negative_jacobians() const;

	void write_jacobians(const std::string fname) const;

	/// Carry out checks on the mesh produced
	void check();
};

}	// end namespace amc

#endif
