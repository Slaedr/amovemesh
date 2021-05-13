/**
 * @file abowyerwatson3d.hpp
 * @brief This file contains a class that implements the Bowyer-Watson algorithm
 *        for Delaunay tesselation in 3D.
 * @author Aditya Kashi
 * @date November 3, 2015
*/

/**
 * \class Delaunay3d
 * \brief A class for Delaunay tetrahedralization of a given set of points based on
 *        the Bowyer-Watson algorithm.
 *
 * It has also been referenced from Wikipedia's Bowyer-Watson algorithm page.
 *
 * Notes:
 * \todo Change points to a pointer to const Matrix rather than a matrix itself,
 * in the interest of efficiency.
 *
 * Currently uses a std::vector to store elements, faces etc.
 * \todo Change the deletion and insertion of faces and elements such that
 * we're not moving entire lists (easier said than done).
 *
 * It might be better to use a std::list or std::forward_list instead.
 * \note A graph data structure should also be seriously considered for storing
 * elements, faces, bad elements and the void polygon.
 */

#ifndef AMC_BOWYERWATSON3D_H
#define AMC_BOWYERWATSON3D_H 1

#include <string>
#include <vector>

#include "amatrix.hpp"

namespace amc {

typedef std::array<double,3> Point;

/// Structure for triangular face.
struct Face
{
	int p[3];			///< indices of endpoints of the face
	int elem[2];		///< elements on either side of the face
};

/// Class representing a tetrahedron.
struct Tet
{
	/// Indices of vertices.
	std::array<int,4> p;

	/// Coords of circumcenter of the tet.
	Point centre{};

	/** Indices of surrounding tets.
	 *
	 * Note that the neighbor corresponding to surr[i] is opposite the vertex p[i].
	 */
	std::array<int,4> surr;

	/// 6 times the volume of the tet.
	double D;

	/// Square of radius of circumcircle of tet.
	double radius;
};

/// Intended to encapsulate data required by 'walk-through' algorithms.
/**  Not needed for mesh generation, but needed, for instance,
 * in independent application of the walk-through subroutine
 * find_containing_tet_and_barycentric_coords().
 */
struct Walkdata
{
	int elem;
	double areacoords[4];
};

/** "Delaunay kernel" for generating the tetrahedral Delaunay tessellation of
 * the convex hull of a set of points.
 */
class Delaunay3d
{
public:
	static constexpr int ndim = 3;
	static constexpr int nnode = 4;

	Delaunay3d();
	Delaunay3d(amat::Matrix<double>* _points);
	Delaunay3d(const Delaunay3d& other);
	Delaunay3d& operator=(const Delaunay3d& other);
	void setup(amat::Matrix<double>* _points);

	/// Computes the Delaunay triangulation (tetrahedralization, in this case).
	/** We first scale the point set by the largest coordinate magnitudes in
	 * each direction, so as to non-dimensionalize it.
	 */
	void bowyer_watson();

	/// Carry out checks on the mesh produced
	void check();

	bool detect_negative_jacobians() const;

	void write_jacobians(const std::string fname) const;

	/// Writes the Delaunay graph to a Gmsh file.
	void writeGmsh2(const std::string mfile) const;

private:
	/// Stores the face-point relationship of a tetrahedron.
	/** Row i contains the local node numbers of the face opposite to node i.
	 * The local node numbers are ordered so that the face points outwards.
	 */
	amat::Matrix<int> lpofa;

	void setlpofa();

	static constexpr int cap = 50;
	static constexpr int badcap = 50;
	static constexpr double tol = 1e-10;

	amat::Matrix<double> points;

	/// List of nodes in the Delaunay graph.
	std::vector<Point> nodes;

	/// List of all elements ([tetrahedra](@ref Tet)) in the Delaunay graph.
	std::vector<Tet> elems;

	/// Collection of 'bad elements', that are to be removed while adding a point
	/** Its are integers that index members index [elems](@ref elems).
	 */
	std::vector<int> badelems;

	/// List of all [faces](@ref Face) in the Delaunay graph.
	std::vector<Face> faces;

	/** Collection of faces that bounds the void obtained after removing
	 * bad elements while adding a point.
	 *
	 * Its members are integers that index [faces](@ref faces).
	 */
	std::vector<int> voidpoly;

	amat::Matrix<double> jacobians;

	int npoints;

	/// Returns the dot product of two vectors.
	double dot(const Point& a, const Point& b) const;

	/// Computes the jacobian of a tet formed from a point and a face of a tetrahedron.
	double det4(const int ielem, const int i_face, const Point& r) const;

	void compute_jacobian(Tet& elem);
	
	void cross_product3(Point& c, const Point& a, const Point& b) const;
	
	/// Computes circumcentre and square of circumradius of a tet
	void compute_circumsphere(Tet& elem);

	/// Computes circumcenter and square of circumradius according to Dr Luo's method
	void compute_circumsphere_contra(Tet& elem);
	
	/// Returs the Jacobian (6 * volume) of a tetrahedron formed by 4 points taken as arguments
	double tetvol(const Point& a, const Point& b, const Point& c, const Point& d) const;

	/// Locates the Delaunay graph (DG) element which contains the input point (improved).
	/** \param xx  The point which needs to be located in the DG.
	 * \param startelement  The index of the element from which to start the search.
	 *
	 * For each DG element encountered, the 4 tetrahedra formed by the point and
	 * each of the 4 faces of the current element are considered.
	 * If the ratios of the Jacobians of all 4 of these new tetrahedra to
	 * the Jacobian of the current DG tetrahedron are positive,
	 * the current element is the containing element. The minimum of the 4 ratios is found.
	 * If the minimum is negative, then the new element to be checked is taken as
	 * the DG element adjacent to the face (of the current element) corresponding to
	 * the minimum value.
	 */
	int find_containing_tet(const Point& r, const int startelement) const;
	
	/** \brief Returns the local face number (number of node opposite to the face)
	 *  of elem's face given by the second argument.
	 *
	 *  If no faces of elem match, it returns -1.
	 */
	int check_face_tet(const Tet& elem, const Face& face) const;

	/// Reset the Delaunay3d object, except for input data
	void clear();

	/** \brief Finds the DG element containing a given point.
	 *
	 * \param rr  The coords of the point to be located in the mesh.
	 * \param startelement  The mesh element ID to start the search from.
	 * \return  The area coordinates of the point in the found element.
	 */
	Walkdata find_containing_tet_and_barycentric_coords(const Point& rr,
														const int startelement) const;
	
	/// Computes the jacobian of all elements in the triangulation using cross products
	void compute_jacobians();
};

}	// end namespace amc

#endif
