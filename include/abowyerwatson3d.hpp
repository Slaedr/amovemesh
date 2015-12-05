/**
@file abowyerwatson3d.hpp
@brief This file contains a class that implements the Bowyer-Watson algorithm for Delaunay tesselation in 3D.

This file is part of KaMoCurve.
@author Aditya Kashi
@date November 3, 2015
*/

/**
\class Delaunay3d
\brief A class for Delaunay tetrahedralization of a given set of points based loosely on Bowyer-Watson algorithm.

It has also been referenced from Wikipedia's Bowyer-Watson algorithm page.

Notes:
  Currently uses a std::vector to store elements, faces etc. 
  It might be better to use a std::list or std::forward_list instead.
  \todo A graph data structure should also be seriously considered for storing elements, faces, bad elements and the void polygon.
*/

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
#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif

#define __ABOWYERWATSON3D_H 1

#ifndef DEBUGBW
#define DEBUGBW 1
#endif

using namespace std;
using namespace amat;

typedef vector<double> Point;

/// Structure for triangular face.
struct Face
{
	int p[3];			///< indices of endpoints of the face
	int elem[2];		///< elements on either side of the face
};

/// Class representing a tetrahedron.
class Tet
{
public:	
	int p[4];			///< Indices of vertices.
	Point centre;		///< Coords of circumcenter of the tet.
	int surr[4];		///< Indices of surrounding tets. Note that the neighbor corresponding to surr[3] is opposite the vertex p[3].
	double D;			///< 6*volume of tet.
	double radius;		///< Radius of circumcircle of tet.

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

/// Intended to encapsulate data required by 'walk-through' algorothms.
/**  Not needed for mesh generation, but needed, for instance, in independent application of the walk-through subroutine find_containing_tet_and_barycentric_coords().
*/
struct Walkdata
{
	int elem;
	double areacoords[4];
};

class Delaunay3d
{
	int cap;
	int badcap;
	double tol;
	int nnode;
	int ndim;

	/// Stores the face-point relationship of a tetrahedron.
	/** Row i contains the local node numbers of the face opposite to node i. The local node numbers are ordered so that the face points outwards. */
	Matrix<int> lpofa;

public:
	Matrix<double> points;
	std::vector<Point> nodes;		///< List of nodes in the Delaunay graph.
	std::vector<Tet> elems;			///< List of all elements ([tetrahedra](@ref Tet)) in the Delaunay graph.
	std::vector<int> badelems;		///< Collection of 'bad elements', that are to be removed while adding a point; its are integers that index members index [elems](@ref elems).
	std::vector<Face> faces;		///< List of all [faces](@ref Face) in the Delaunay graph.
	std::vector<int> voidpoly;		///< Collection of faces that bounds the void obtained after removing bad elements while adding a point; its members are integers that index [faces](@ref faces).
	Matrix<double> jacobians;

	int npoints;

	Delaunay3d();
	Delaunay3d(Matrix<double>* _points, int num_points);
	Delaunay3d(const Delaunay3d& other);
	Delaunay3d& operator=(const Delaunay3d& other);
	void setup(Matrix<double>* _points, int num_points);
	
	double l2norm(const vector<double>& a);
	
	double dot(const vector<double>& a, const vector<double>& b);
	
	/// Computes the jacobian of a tet formed from a point and a face of a tetrahedron.
	double det4(int ielem, int i, const vector<double>& r) const;
	
	void compute_jacobian(Tet& elem);
	
	void cross_product3(vector<double>& c, const vector<double>& a, const vector<double>& b);
	
	void compute_circumsphere(Tet& elem);

	void compute_circumsphere_contra(Tet& elem);
	
	/// Returs the Jacobian (6 * volume) of a tetrahedron formed by 4 points taken as arguments
	double tetvol(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) const;
	
	int find_containing_tet(const vector<double>& r, int startelement) const;
	
	int check_face_tet(const Tet& elem, const Face& face) const;

	void bowyer_watson();
	
	void clear();					///< Reset the Delaunay3d object, except for input data
	
	/// Writes the Delaunay graph to a Gmsh file.
	void writeGmsh2(string mfile);
	
	/// Finds the DG element containing a given point.
	Walkdata find_containing_tet_and_barycentric_coords(const vector<double>& rr, int startelement) const;
	
	void compute_jacobians();
	
	bool detect_negative_jacobians();
};

Delaunay3d::Delaunay3d() {}

Delaunay3d::Delaunay3d(Matrix<double>* _points, int num_points)
{
	nnode = 4;
	ndim = 3;
	cap = 1000;
	badcap = 50;
	tol = 1e-10;
	elems.reserve(cap);
	faces.reserve(cap);
	points = *_points;
	npoints = num_points;
	nodes.reserve(num_points+3);

	lpofa.setup(nnode,nnode-1);
	lpofa(0,0) = 1; lpofa(0,1) = 2; lpofa(0,2) = 3;
	lpofa(1,0) = 2; lpofa(1,1) = 0; lpofa(1,2) = 3;
	lpofa(2,0) = 3; lpofa(2,1) = 0; lpofa(2,2) = 1;
	lpofa(3,0) = 0; lpofa(3,1) = 2; lpofa(3,2) = 1;

	cout << setprecision(12);
}

Delaunay3d::Delaunay3d(const Delaunay3d& other)
{
	cap = other.cap;
	badcap = other.badcap;
	tol = other.tol;
	nnode = other.nnode;
	ndim = other.ndim;
	points = other.points;
	nodes = other.nodes;
	elems = other.elems;
	faces = other.faces;
	jacobians = other.jacobians;
	npoints = other.npoints;
}

Delaunay3d& Delaunay3d::operator=(const Delaunay3d& other)
{
	cap = other.cap;
	badcap = other.badcap;
	tol = other.tol;
	nnode = other.nnode;
	ndim = other.ndim;
	points = other.points;
	nodes = other.nodes;
	elems = other.elems;
	faces = other.faces;
	jacobians = other.jacobians;
	npoints = other.npoints;
	return *this;
}

void Delaunay3d::setup(Matrix<double>* _points, int num_points)
{
	nnode = 4;
	ndim = 3;
	cap = 1000;
	badcap = 50;
	tol = 1e-10;
	elems.reserve(cap);
	faces.reserve(cap);
	points = *_points;
	npoints = num_points;
	nodes.reserve(num_points+3);

	lpofa.setup(nnode,nnode-1);
	lpofa(0,0) = 1; lpofa(0,1) = 2; lpofa(0,2) = 3;
	lpofa(1,0) = 2; lpofa(1,1) = 0; lpofa(1,2) = 3;
	lpofa(2,0) = 3; lpofa(2,1) = 0; lpofa(2,2) = 1;
	lpofa(3,0) = 0; lpofa(3,1) = 2; lpofa(3,2) = 1;

	cout << setprecision(12);
}

double Delaunay3d::tetvol(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) const
{
	/*Matrix<double> vold(ndim,ndim);
	vold(0,0) = a[0]-b[0]; vold(0,1) = b[0]-c[0]; vold(0,2) = c[0]-d[0];
	vold(1,0) = a[1]-b[1]; vold(1,1) = b[1]-c[1]; vold(1,2) = c[1]-d[1];
	vold(2,0) = a[2]-b[2]; vold(2,1) = b[2]-c[2]; vold(2,2) = c[2]-d[2];
	return determinant(vold);*/

	double x1, y1, z1, x21, y21, z21, x31, y31, z31, x41, y41, z41;
	x1 = a[0]; y1 = a[1]; z1 = a[2];

	x21 = b[0] - x1;
	y21 = b[1] - y1;
	z21 = b[2] - z1;

	x31 = c[0] - x1;
	y31 = c[1] - y1;
	z31 = c[2] - z1;

	x41 = d[0] - x1;
	y41 = d[1] - y1;
	z41 = d[2] - z1;

	return x21*(y31*z41-z31*y41) + x31*(z21*y41-y21*z41) + x41*(y21*z31-z21*y31);
}

/// Computes 6*volume, ie the Jacobian of any tetrahedron.
void Delaunay3d::compute_jacobian(Tet& elem)
{	
	/*double ret;
	ret = nodes[elem.p[0]][0] * ( nodes[elem.p[1]][1]*(nodes[elem.p[2]][2]-nodes[elem.p[3]][2]) -nodes[elem.p[1]][2]*(nodes[elem.p[2]][1]-nodes[elem.p[3]][1]) + nodes[elem.p[2]][1]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][1] );
	ret -= nodes[elem.p[0]][1] * (nodes[elem.p[1]][0]*(nodes[elem.p[2]][2]-nodes[elem.p[3]][2]) -nodes[elem.p[1]][2]*(nodes[elem.p[2]][0]-nodes[elem.p[3]][0]) + nodes[elem.p[2]][0]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][0] );
	ret += nodes[elem.p[0]][2] * (nodes[elem.p[1]][0]*(nodes[elem.p[2]][1]-nodes[elem.p[3]][1]) -nodes[elem.p[1]][1]*(nodes[elem.p[2]][0]-nodes[elem.p[3]][0]) + nodes[elem.p[2]][0]*nodes[elem.p[3]][1] - nodes[elem.p[2]][1]*nodes[elem.p[3]][0] );
	ret -= nodes[elem.p[1]][0]*( nodes[elem.p[2]][1]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][1] ) -nodes[elem.p[1]][1]*( nodes[elem.p[2]][0]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][0] ) +nodes[elem.p[1]][2]*( nodes[elem.p[2]][0]*nodes[elem.p[3]][1] - nodes[elem.p[2]][1]*nodes[elem.p[3]][0] );
	
	elem.D = ret;*/

	elem.D = tetvol(nodes[elem.p[0]],nodes[elem.p[1]],nodes[elem.p[2]],nodes[elem.p[3]]);
}

/// Computes cross product c of two 3-vectors a and b.
inline void Delaunay3d::cross_product3(vector<double>& c, const vector<double>& a, const vector<double>& b)
{
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
}

/// Returns the dot product of two vectors.
inline double Delaunay3d::dot(const vector<double>& a, const vector<double>& b)
{
	double val = 0;
	for(int i = 0; i < ndim; i++)
		val += a[i]*b[i];
	return val;
}

/// Returns the \f$ l^2 \f$ norm of a vector.
inline double Delaunay3d::l2norm(const vector<double>& a)
{
	// a.size() should ideally equal ndim
	double norm = 0;
	for(int i = 0; i < a.size(); i++)
		norm += a[i]*a[i];
	norm = sqrt(norm);
	return norm;
}

/// Cmoputes circumcentre and circumradius of a tet
void Delaunay3d::compute_circumsphere(Tet& elem)
{
	// normals to the perp bisectors of 3 non-coplanar edges of the tet
	Matrix<double> normal(ndim,ndim), rhs(ndim,1);
	// jth component of midpoint of an edge
	double midpointj;

	int i,j,k;
	for(i = 0; i < ndim; i++)		// loop over three non-coplanar edges of the tet
	{
		for(j = 0; j < ndim; j++)
			normal(i,j) = nodes[elem.p[i+1]][j] - nodes[elem.p[0]][j];
		
		rhs(i) = 0;
		for(j = 0; j < ndim; j++) 
		{
			midpointj = (nodes[elem.p[i+1]][j]+nodes[elem.p[0]][j])/2.0;
			rhs(i) += normal(i,j)*midpointj;
		}
	}

	// solve the linear system by Cramer's rule to get coords of the circumcentre
	vector<double> rc(ndim);
	Matrix<double> numer(ndim,ndim);
	double detnum;
	//double detdenom = determinant(normal);
	double detdenom = 0;
	detdenom = normal(0,0)*( normal(1,1)*normal(2,2)-normal(1,2)*normal(2,1) ) -normal(0,1)*( normal(1,0)*normal(2,2)-normal(1,2)*normal(2,0)) +normal(0,2)*( normal(1,0)*normal(2,1)-normal(1,1)*normal(2,0));
	if(fabs(detdenom) <= ZERO_TOL) cout << "Delaunay3d: compute_circumsphere(): ! System is inconsistent!! Jacobian of elem is " << elem.D << endl;
	
	for(k = 0; k < ndim; k++)
	{
		for(j = 0; j < ndim; j++)
		{
			if(j == k) continue;
			for(i = 0; i < ndim; i++)
				numer(i,j) = normal.get(i,j);
		}
		for(i = 0; i < ndim; i++)
			numer(i,k) = rhs.get(i);
		
		detnum = numer(0,0)*( numer(1,1)*numer(2,2)-numer(1,2)*numer(2,1) ) -numer(0,1)*( numer(1,0)*numer(2,2)-numer(1,2)*numer(2,0)) +numer(0,2)*( numer(1,0)*numer(2,1)-numer(1,1)*numer(2,0));
		//detnum = determinant(numer);
		rc[k] = detnum/detdenom;
	}

	elem.centre = rc;
	elem.radius = 0.0;
	
	// now calculate radius
	for(j = 0; j < ndim; j++)
		elem.radius += (rc[j]-nodes[elem.p[0]][j])*(rc[j]-nodes[elem.p[0]][j]);

	// check
	Matrix<double> rhstest(ndim,1); rhstest.zeros();
	for(i = 0; i < ndim; i++)
		for(j = 0; j < ndim; j++)
			rhstest(i) += normal(i,j)*rc[j];
	
	rhstest.mprint();
	rhs.mprint();

	cout << "Circumsphere data : centre " << elem.centre[0] << "," << elem.centre[1] << "," << elem.centre[2] << ", radius^2 " << elem.radius << endl;
	cout << "Delaunay3d: compute_circumsphere(): Element jacobian = " << elem.D << endl;
}

void Delaunay3d::compute_circumsphere_contra(Tet& elem)
{
	vector<double> g1(ndim), g2(ndim), nrm(ndim),  g2co(ndim), cold(ndim), rdiff(ndim), rave(ndim);
	double tp;
	int k;
	
	for(k = 0; k < ndim; k++)
	{
		g1[k] = nodes[elem.p[1]][k] - nodes[elem.p[0]][k];
		g2[k] = nodes[elem.p[2]][k] - nodes[elem.p[0]][k];
	}
	
	// normal vector to "base" face 0-1-2
	nrm[0] = g1[1]*g2[2] - g1[2]*g2[1];
	nrm[1] = g1[2]*g2[0] - g1[0]*g2[2];
	nrm[2] = g1[0]*g2[1] - g1[1]*g2[0];
	
	// normalize the normal vector
	double normnrm = 0;
	for(k = 0; k < ndim; k++)
		normnrm += nrm[k]*nrm[k];
	
	#if DEBUGW==1
	if(normnrm < ZERO_TOL) cout << "Delaunay3d: compute_circumsphere_contra(): ! Element is degenerate!!" << endl;
	#endif
	
	normnrm = 1.0/sqrt(normnrm);
	for(k = 0; k < ndim; k++)
		nrm[k] *= normnrm;
	
	// "contravariant" vector g2co = nrm _cross_ g1
	g2co[0] = nrm[1]*g1[2] - nrm[2]*g1[1];
	g2co[1] = nrm[2]*g1[0] - nrm[0]*g1[2];
	g2co[2] = nrm[0]*g1[1] - nrm[1]*g1[0];

	// compute "factor t for base", ie, tp
	double tpd = 0;
	tp = 0;
	for(k = 0; k < ndim; k++)
	{
		tpd += g2co[k]*g2[k];
		tp += (g2[k]-g1[k])*g2[k];
	}
	tp = 0.5 * tp / tpd;

	// circumcentre of base face
	for(k = 0; k < ndim; k++)
		cold[k] = nodes[elem.p[0]][k] + 0.5*g1[k] + tp*g2co[k];
	
	// now average and difference of node 3 (the final node, not part of the base face) and node 0
	for(k = 0; k < ndim; k++)
	{
		rdiff[k] = nodes[elem.p[3]][k] - nodes[elem.p[0]][k];
		rave[k] = 0.5*( nodes[elem.p[3]][k] + nodes[elem.p[0]][k] );
	}

	// final t parameter
	tpd = 0; tp = 0;
	for(k = 0; k < ndim; k++)
	{
		tpd += nrm[k]*rdiff[k];
		tp += (rave[k]-cold[k])*rdiff[k];
	}
	tp = tp/tpd;

	// finally, the circumcentre
	for(k = 0; k < ndim; k++)
		elem.centre[k] = cold[k] + tp*nrm[k];
	
	// radius
	elem.radius = 0;
	for(k = 0; k < ndim; k++)
		elem.radius += (elem.centre[k] - nodes[elem.p[0]][k])*(elem.centre[k] - nodes[elem.p[0]][k]);
	
	//cout << "Circumsphere data : centre " << elem.centre[0] << "," << elem.centre[1] << "," << elem.centre[2] << ", radius^2 " << elem.radius << endl;
	//cout << "Delaunay3d: compute_circumsphere_contra(): Element jacobian = " << elem.D << endl;
}

///	Calculates the jacobian of the tetrahedron formed by point r and a face of tetrahedron ielem.
/** The face is selected by i between 0 and 3. Face i is the face opposite to local node i of the tetrahedron.
*/
double Delaunay3d::det4(int ielem, int i, const vector<double>& r) const
{
	#if DEBUGW==1
	if(i > 3) {
		std::cout << "Delaunay3D: det4(): ! Second argument is greater than 3!" << std::endl;
		return 0;
	}
	#endif

	Tet elem = elems[ielem];
	/*double ret = 0;
	switch(i)
	{
		case(0):
		ret = r[0] * ( nodes[elem.p[1]][1]*(nodes[elem.p[2]][2]-nodes[elem.p[3]][2]) -nodes[elem.p[1]][2]*(nodes[elem.p[2]][1]-nodes[elem.p[3]][1]) + nodes[elem.p[2]][1]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][1] );
		ret -= r[1] * (nodes[elem.p[1]][0]*(nodes[elem.p[2]][2]-nodes[elem.p[3]][2]) -nodes[elem.p[1]][2]*(nodes[elem.p[2]][0]-nodes[elem.p[3]][0]) + nodes[elem.p[2]][0]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][0] );
		ret += r[2] * (nodes[elem.p[1]][0]*(nodes[elem.p[2]][1]-nodes[elem.p[3]][1]) -nodes[elem.p[1]][1]*(nodes[elem.p[2]][0]-nodes[elem.p[3]][0]) + nodes[elem.p[2]][0]*nodes[elem.p[3]][1] - nodes[elem.p[2]][1]*nodes[elem.p[3]][0] );
		ret -= nodes[elem.p[1]][0]*( nodes[elem.p[2]][1]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][1] ) -nodes[elem.p[1]][1]*( nodes[elem.p[2]][0]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][0] ) +nodes[elem.p[1]][2]*( nodes[elem.p[2]][0]*nodes[elem.p[3]][1] - nodes[elem.p[2]][1]*nodes[elem.p[3]][0] );
		break;

		case(1):
		ret = nodes[elem.p[0]][0] * ( r[1]*(nodes[elem.p[2]][2]-nodes[elem.p[3]][2]) -r[2]*(nodes[elem.p[2]][1]-nodes[elem.p[3]][1]) + nodes[elem.p[2]][1]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][1] );
		ret -= nodes[elem.p[0]][1] * (r[0]*(nodes[elem.p[2]][2]-nodes[elem.p[3]][2]) -r[2]*(nodes[elem.p[2]][0]-nodes[elem.p[3]][0]) + nodes[elem.p[2]][0]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][0] );
		ret += nodes[elem.p[0]][2] * (r[0]*(nodes[elem.p[2]][1]-nodes[elem.p[3]][1]) -r[1]*(nodes[elem.p[2]][0]-nodes[elem.p[3]][0]) + nodes[elem.p[2]][0]*nodes[elem.p[3]][1] - nodes[elem.p[2]][1]*nodes[elem.p[3]][0] );
		ret -= r[0]*( nodes[elem.p[2]][1]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][1] ) -r[1]*( nodes[elem.p[2]][0]*nodes[elem.p[3]][2] - nodes[elem.p[2]][2]*nodes[elem.p[3]][0] ) +r[2]*( nodes[elem.p[2]][0]*nodes[elem.p[3]][1] - nodes[elem.p[2]][1]*nodes[elem.p[3]][0] );
		break;

		case(2):
		ret = nodes[elem.p[0]][0] * ( nodes[elem.p[1]][1]*(r[2]-nodes[elem.p[3]][2]) -nodes[elem.p[1]][2]*(r[1]-nodes[elem.p[3]][1]) + r[1]*nodes[elem.p[3]][2] - r[2]*nodes[elem.p[3]][1] );
		ret -= nodes[elem.p[0]][1] * (nodes[elem.p[1]][0]*(r[2]-nodes[elem.p[3]][2]) -nodes[elem.p[1]][2]*(r[0]-nodes[elem.p[3]][0]) + r[0]*nodes[elem.p[3]][2] - r[2]*nodes[elem.p[3]][0] );
		ret += nodes[elem.p[0]][2] * (nodes[elem.p[1]][0]*(r[1]-nodes[elem.p[3]][1]) -nodes[elem.p[1]][1]*(r[0]-nodes[elem.p[3]][0]) + r[0]*nodes[elem.p[3]][1] - r[1]*nodes[elem.p[3]][0] );
		ret -= nodes[elem.p[1]][0]*( r[1]*nodes[elem.p[3]][2] - r[2]*nodes[elem.p[3]][1] ) -nodes[elem.p[1]][1]*( r[0]*nodes[elem.p[3]][2] - r[2]*nodes[elem.p[3]][0] ) +nodes[elem.p[1]][2]*( r[0]*nodes[elem.p[3]][1] - r[1]*nodes[elem.p[3]][0] );
		break;

		case(3):
		ret = nodes[elem.p[0]][0] * ( nodes[elem.p[1]][1]*(nodes[elem.p[2]][2]-r[2]) -nodes[elem.p[1]][2]*(nodes[elem.p[2]][1]-r[1]) + nodes[elem.p[2]][1]*r[2] - nodes[elem.p[2]][2]*r[1] );
		ret -= nodes[elem.p[0]][1] * (nodes[elem.p[1]][0]*(nodes[elem.p[2]][2]-r[2]) -nodes[elem.p[1]][2]*(nodes[elem.p[2]][0]-r[0]) + nodes[elem.p[2]][0]*r[2] - nodes[elem.p[2]][2]*r[0] );
		ret += nodes[elem.p[0]][2] * (nodes[elem.p[1]][0]*(nodes[elem.p[2]][1]-r[1]) -nodes[elem.p[1]][1]*(nodes[elem.p[2]][0]-r[0]) + nodes[elem.p[2]][0]*r[1] - nodes[elem.p[2]][1]*r[0] );
		ret -= nodes[elem.p[1]][0] *( nodes[elem.p[2]][1]*r[2] - nodes[elem.p[2]][2]*r[1] ) -nodes[elem.p[1]][1]*( nodes[elem.p[2]][0]*r[2] - nodes[elem.p[2]][2]*r[0] ) + nodes[elem.p[1]][2]*( nodes[elem.p[2]][0]*r[1] - nodes[elem.p[2]][1]*r[0] );
		break;

		default:
		cout << "Delaunay3D: det4(): ! Invalid argument i! Should be between 0 and 3 inclusive." << endl;
		ret = -1;
	}
	return ret;*/

	vector<double> ta[4];
	for(int j = 0; j < nnode; j++)
		ta[j] = nodes[elem.p[j]];
	ta[i] = r;
	return tetvol(ta[0],ta[1],ta[2],ta[3]);
}

/// Locates the Delaunay graph (DG) element which contains the input point.
/** \param xx is the point which needs to be located in the DG.
* \param startelement is the index of the element from which to start the search.
*
* For each DG element encountered, the 4 tetrahedra formed by the point and each of the 4 faces of the current element are considered. 
* If the ratios of the Jacobians of all 4 of these new tetrahedra to the Jacobian of the current DG tetrahedron are positive, the current element is the containing element. If one of these is negative, the new element to be checked is taken as the DG element adjacent to the face (of the current element) corresponding to the negative value.
*/
int Delaunay3d::find_containing_tet(const vector<double>& xx, int startelement) const
{
	if(xx.size() < ndim) {
		std::cout << "Delaunau3D: find_containing_triangle(): ! Input vector is not long enough!\n";
		return -1;
	}
	int ielem = startelement;
	Tet super;
	double l;
	bool found;
	
	while(1)
	{
		found = true;

		if(ielem < 0 || ielem >= elems.size()) { cout << "Delaunay3d:   !! Reached an element index that is out of bounds!! Index is " << ielem << "\n"; return ielem; }
		super = elems[ielem];

		for(int inode = 0; inode < nnode; inode++)
		{
			// get jacobian
			l = det4(ielem,inode,xx);
			#if DEBUGW==1
			if(dabs(l) < ZERO_TOL) cout << "Delaunay3d: find_containing_tet(): ! Degenerate case (type 1) " << inode << "!!\n";
			#endif
			if(l/super.D < 0)
			{
				ielem = super.surr[inode];
				found = false;
				break;
			}
		}

		// if all 4 area-ratios are positive, we've found our element
		if(found) break;
	}
	cout << "Delaunay2D:   Containing triangle found as " << ielem << endl;
	return ielem;
}

/// Returns the local face number (number of node opposite to the face) of elem's face that is the same as the second argument.
/**  If no faces of elem match, it returns -1.
*/
inline int Delaunay3d::check_face_tet(const Tet& elem, const Face& face) const
{
	// iterate over each face of elem
	for(int i = 0; i < 4; i++)
	{
		// check all 6 permutations by which elem's ith face could be same as the input face.
		if(face.p[0]==elem.p[(0+i)%4] && face.p[1]==elem.p[(1+i)%4] && face.p[2]==elem.p[(2+i)%4]) return (3+i)%4;
		if(face.p[0]==elem.p[(0+i)%4] && face.p[1]==elem.p[(2+i)%4] && face.p[2]==elem.p[(1+i)%4]) return (3+i)%4;
		if(face.p[0]==elem.p[(1+i)%4] && face.p[1]==elem.p[(0+i)%4] && face.p[2]==elem.p[(2+i)%4]) return (3+i)%4;
		if(face.p[0]==elem.p[(1+i)%4] && face.p[1]==elem.p[(2+i)%4] && face.p[2]==elem.p[(0+i)%4]) return (3+i)%4;
		if(face.p[0]==elem.p[(2+i)%4] && face.p[1]==elem.p[(0+i)%4] && face.p[2]==elem.p[(1+i)%4]) return (3+i)%4;
		if(face.p[0]==elem.p[(2+i)%4] && face.p[1]==elem.p[(1+i)%4] && face.p[2]==elem.p[(0+i)%4]) return (3+i)%4;
	}
	return -1;
}

/// Computes the Delaunay triangulation (tetrahedralization, in this case).
/** Make sure 'points' has space for three more points when passing to this sub. 'N' is the actual number of real points.
*/
void Delaunay3d::bowyer_watson()
{
	// add super triangle
	//	find minimum and maximum x and y of the point set
	vector<double> rmin(ndim,0), rmax(ndim,0);
	for(int i = 1; i < npoints; i++)
	{
		for(int idim = 0; idim < ndim; idim++) {
			if(points(i,idim) > rmax[idim]) rmax[idim] = points(i,idim);
			if(points(i,idim) < rmin[idim]) rmin[idim] = points(i,idim);
		}
	}
	
	cout << "Bounds of the point set are ";
	for(int idim = 0; idim < ndim; idim++)
		cout << " " << rmin[idim] << " " << rmax[idim] << endl;
	cout << "**\n";

	// factor by which to scale rdelt, for providing a factor of safety.
	double factor = 10.0;
	
	vector<double> rdelt(ndim);		// stores some extra length to have factor of safety in deciding the first 4 points
	for(int idim = 0; idim < ndim; idim++){
		rdelt[idim] = A_SMALL_NUMBER + factor*(rmax[idim]-rmin[idim]);
		rmax[idim] += rdelt[idim];
		rmin[idim] -= rdelt[idim];
	}

	cout << "Bounds of the extended cuboid are ";
	for(int idim = 0; idim < ndim; idim++)
		cout << " " << rmin[idim] << " " << rmax[idim] << endl;
	cout << "**\n";

	// now rmax and rmin denote the expanse of a cuboid that contains all the points.

	Tet super;
	vector<Point> pp(nnode);
	for(int i = 0; i < nnode; i++)
		pp[i].resize(ndim);

	pp[0][0] = rmin[0] - (rmax[2]-rmin[2])*0.5;
	pp[0][1] = rmin[1];
	pp[0][2] = (rmax[2]+rmin[2])*0.5 + rmax[1] - rmin[1];

	pp[1][0] = rmax[0] + rmax[1] - rmin[1];
	pp[1][1] = rmin[1];
	pp[1][2] = rmax[0]-rmin[0] + 2.0*(rmax[1]-rmin[1]) + rmax[2];

	pp[2][0] = rmax[0] + rmax[1] - rmin[1];
	pp[2][1] = rmin[1];
	pp[2][2] = rmin[2] - (rmax[0]-rmin[0]);

	pp[3][0] = rmin[0] - 0.5*(rmax[2]-rmin[2]);
	pp[3][1] = rmax[0]-rmin[0] + 0.5*(rmax[2]-rmin[2]) + rmax[2];
	pp[3][2] = rmin[0] + rmin[2] - rmax[0];

	cout << "Delaunay3d: bowyer_watson(): Coordinates of vertices of the super triangle are:\n";
	for(int inode = 0; inode < nnode; inode++)
	{
		for(int idim = 0; idim < ndim; idim++)
			cout << pp[inode][idim] << "\t";
		cout << endl;
	}
	
	for(int i = 0; i < nnode; i++)
	{
		nodes.push_back(pp[i]);
		super.p[i] = i;
	}
	for(int inode = 0; inode < nnode; inode++)
		super.surr[inode] = -1;

	compute_jacobian(super);
	compute_circumsphere(super);

	elems.push_back(super);			// add super to elems list

	// set up initial face list
	//cout << "Delaunay3d: set up initial face list\n";
	Face f[ndim+1];
	f[0].p[0] = super.p[1]; f[0].p[1] = super.p[0]; f[0].p[2] = super.p[2];
	f[1].p[0] = super.p[2]; f[1].p[1] = super.p[0]; f[1].p[2] = super.p[3];
	f[2].p[0] = super.p[3]; f[2].p[1] = super.p[1]; f[2].p[2] = super.p[2];
	f[3].p[0] = super.p[3]; f[3].p[1] = super.p[0]; f[3].p[2] = super.p[1];
	for(int i = 0; i < nnode; i++)
	{
		f[i].elem[0] = 0;
		f[i].elem[1] = -1;
		faces.push_back(f[i]);
	}

	int newpoinnum, contelem;
	vector<double> newpoin(ndim);

	// iterate through points
	cout << "Delaunay3d: Starting iteration over points\n";
	for(int ipoin = 0; ipoin < npoints; ipoin++)
	{
		cout << "New point coords: " << points.get(ipoin,0) << " " << points.get(ipoin,1) << " " << points.get(ipoin,2) << endl;
		for(int idim = 0; idim < ndim; idim++)
			newpoin[idim] = points.get(ipoin,idim);
		nodes.push_back(newpoin);
		newpoinnum = nodes.size()-1;

		/// First, find the element containing the new point
		contelem = find_containing_tet(newpoin,elems.size()-1);

		/// Second, search among neighbors for other triangles whose circumcircles contain this point
		int curelem;
		vector<int> stk;					// stack to hold the indices of tets to be checked
		double dist;						// square of distance from point to circumcentre
		stk.push_back(contelem);			// add the containing element to the list of bad elements
		vector<int> flags(elems.size());	// stores 1 at an index if the corresponding element has been checked for the Delaunay criterion
		for(int i = 0; i < elems.size(); i++) 
			flags[i] = 0;

		while(stk.empty() == false)
		{
			curelem = stk.back();			// access last element in stack of elements to be checked

			// if this element has already been checked, just remove it and continue
			if(flags[curelem] == 1)
			{
				stk.pop_back();
				continue;
			}

			flags[curelem] = 1;				// curelem will now be checked

			//calculate distance between circumcentre and the point
			dist = 0;
			for(int idim = 0; idim < ndim; idim++)
				dist += (newpoin[idim] - elems[curelem].centre[idim])*(newpoin[idim] - elems[curelem].centre[idim]);

			#if DEBUGBW==1
			if(dabs(dist - elems[curelem].radius) < ZERO_TOL) cout << "Delaunay2D: Degenerate case (type 2)!!\n";
			#endif
			
			// FOR DEBUG
			cout << "Delaunay3d: bowyer_watson(): Dist^2 and radius^2 are " << dist << ", " << elems[curelem].radius << endl;

			if(dist <= elems[curelem].radius)		// if point lies inside circumcircle, ie, Delaunay criterion is violated
			{
				badelems.push_back(curelem);
				stk.pop_back();
				for(int j = 0; j < ndim+1; j++)
					if(elems[curelem].surr[j] >= 0 && flags[elems[curelem].surr[j]] == 0)		// add surrounding elements (only which are not checked) to "to be checked" stack
						stk.push_back(elems[curelem].surr[j]);
			}
			else
			{
				stk.pop_back();
			}
		}

		// Output badelems for debug purpose
		cout << "Delaunay3d:  Badelems: ";
		for(int i = 0; i < badelems.size(); i++)
			cout << badelems[i] << " ";
		cout << endl;

		/// Third, we store the faces that will be obtained after removal of bad elements
		flags.assign(faces.size(),-1);
		cout << "** Flags : (" << faces.size() << ", " << flags.size() << ")\n";

		for(int ifa = 0; ifa < faces.size(); ifa++)
		{
			for(int itri = 0; itri < badelems.size(); itri++)
			{
				if(faces[ifa].elem[0] == badelems[itri])		//this face belongs to at least one bad element
				{
					if(flags[ifa] != 1)							// this face does not belong to two bad elements
						flags[ifa] = 0;
					else										// this face belongs to two bad elements
						flags[ifa] = -10;
				}
				else if(faces[ifa].elem[1] == badelems[itri])	// this face belongs to at least one bad element
				{
					if(flags[ifa] != 0)							// this face is NOT shared by two bad elements
						flags[ifa] = 1;
					else										// this face is shared by two bad elements
						flags[ifa] = -10;
				}
			}
			if(flags[ifa] == 0)
			{
				faces[ifa].elem[0] = -2;		// -2 will be replaced by the index of a newly-formed element later
				voidpoly.push_back(ifa);
			}
			if(flags[ifa] == 1)
			{
				faces[ifa].elem[1] = -2;
				voidpoly.push_back(ifa);
			}
		}
		
		/** Delete faces which are between two bad elements.
		*  NOTE: This is one place that is ineffecient because of use of array stacks (std::vectors) as it needs deletion of arbitrary members.
		*/
		for(int ifa = 0; ifa < faces.size(); ifa++)
		{
			if(flags[ifa] == -10)				// if face belongs to two bad elements, delete face
			{
				faces.erase(faces.begin()+ifa);
				flags.erase(flags.begin()+ifa);
				// now adjust voidpoly for the deleted face
				for(int i = 0; i < voidpoly.size(); i++)
				{
					if(voidpoly[i] > ifa) voidpoly[i]--;
				}
				ifa--;
			}
		}

		/** Fourth, delete bad elements.
		*  This is another place where array stacks (vectors) of elems, badelems etc make the program slower.
		*/
		//cout << "Delaunay2D:  Fourth, delete bad elements\n";
		for(int ibe = 0; ibe < badelems.size(); ibe++)
		{
			elems.erase(elems.begin()+badelems[ibe]);

			//scan badelems for elements with indices greater than the one just deleted
			for(int i = ibe+1; i < badelems.size(); i++)
				if(badelems[i] > badelems[ibe]) badelems[i]--;

			//scan surrounding elements in elems -- probably not efficient
			for(int i = 0; i < elems.size(); i++)
			{
				for(int j = 0; j < ndim+1; j++)
				{
					if(elems[i].surr[j] == badelems[ibe]) elems[i].surr[j] = -3;
					if(elems[i].surr[j] > badelems[ibe]) elems[i].surr[j]--;
				}
			}

			// adjust face data as well
			for(int i = 0; i < faces.size(); i++)
			{
				if(faces[i].elem[0] > badelems[ibe]) faces[i].elem[0]--;
				if(faces[i].elem[1] > badelems[ibe]) faces[i].elem[1]--;
			}
		}

		// output voidpoly for debug purposes
		cout << "Delaunay3d:  voidpoly:";
		for(int i = 0; i < voidpoly.size(); i++)
			cout << " " << voidpoly[i];
		cout << endl;

		/// Fifth, add new elements; these are formed by the faces in voidpoly and the new point. Also correspondingly update 'faces'.
		//cout << "Delaunay2D:  Fifth, add new elements; Also correspondingly update 'faces'\n";
		vector<int> newfaces;				// new faces formed from new elements created
		int temp;
		for(int ifa = 0; ifa < voidpoly.size(); ifa++)		// for each face in void polygon
		{
			Tet nw;
			nw.p[0] = newpoinnum;
			if(faces[voidpoly[ifa]].elem[0] == -2)	
			{	// if the new element is to the left of this face, the orientation of the remaining points of the tet is same as the orientation of the face
				nw.p[1] = faces[voidpoly[ifa]].p[0];
				nw.p[2] = faces[voidpoly[ifa]].p[1];
				nw.p[3] = faces[voidpoly[ifa]].p[2];
				// we also need to reverse the face's orientation to make it point to the new element
				temp = faces[voidpoly[ifa]].elem[1];
				faces[voidpoly[ifa]].elem[1] = elems.size();		// this new element has not been pushed into elems yet, so we need elems.size() - 1 + 1
				faces[voidpoly[ifa]].elem[0] = temp;
				faces[voidpoly[ifa]].p[0] = nw.p[1];
				faces[voidpoly[ifa]].p[1] = nw.p[3];
				faces[voidpoly[ifa]].p[2] = nw.p[2];
			}
			else if(faces[voidpoly[ifa]].elem[1] == -2)	// if the new element is to the right of the face
			{
				nw.p[1] = faces[voidpoly[ifa]].p[1];
				nw.p[2] = faces[voidpoly[ifa]].p[0];
				nw.p[3] = faces[voidpoly[ifa]].p[2];
				faces[voidpoly[ifa]].elem[1] = elems.size();
			}
			else cout << "Delaunay2D: !! Error while creating new element - face " << voidpoly[ifa] << " in voidpoly does not have -2 as either left elem or right elem.\n";

			compute_jacobian(nw);
			compute_circumsphere_contra(nw);

			if(nw.D < ZERO_TOL)
				cout << "Delaunay3d: bowyer_watson(): New elem is degenerate or inverted! Points are " << nw.p[0] << " " << nw.p[1] << " " << nw.p[2] << " " << nw.p[3] << endl;

			cout << "Delaunay3d:  New element created." << endl;
			// Push new element into the elements' list
			elems.push_back(nw);

			Face fc;

			vector<bool> val(ndim+1, false);	
			//^ flag for each face face of the new element - but we actually only need it for those faces that include the new point 0, ie faces 1, 2 and 3.
			//^ val[i] contains true if a newface corresponding to the ith local face of nw has been found.

			int localface, jface;
			for(int jfa = 0; jfa < newfaces.size(); jfa++)
			{
				//Instead of separately comparing each face of nw, use the function to test whether a face is part of a tet.
				localface = check_face_tet(nw, faces[newfaces[jfa]]);
				
				// If this newface does not match any face of the new element, go to the next newface.
				if(localface == -1) continue;

				// Face formed by new point and p[1]
				else if(localface==0) {
					cout << "Delaunay3d: bowyer_watson(): ! Error: A member of newfaces is not actually a new face!!" << endl;
					continue;
				}

				// if this newface has -4 as the right element, then replace this by the new element
				//if(faces[newfaces[jfa]].elem[1] == -4)		// Do we really need this check? I don't think so.
				//{
				
				faces[newfaces[jfa]].elem[1] = elems.size()-1;
				elems.back().surr[localface] = faces[newfaces[jfa]].elem[0];

				// Also, find which local face of newfaces[jfa]'s left element is the same as newfaces[jfa].
				// Set the new element as a surrounding element of of the left element of newfaces[jfa].
				jface = check_face_tet(elems[faces[newfaces[jfa]].elem[0]], faces[newfaces[jfa]]);
				if(jface == -1) {cout << "Delaunay3d: bowyer_watson(): ! Error while setting surrounding element of a new element!" << endl; }
				elems[faces[newfaces[jfa]].elem[0]].surr[jface] = elems.size()-1;
				
				//}

				val[localface] = true;
				
				// if faces corresponding to all 3 new faces of the new element have been found, quit the jfa loop
				bool brk = true;
				for(int i = 1; i < ndim+1; i++)
					if(!val[i]) brk = false;
				if(brk) break;
			}

			for(int iface = 1; iface < ndim+1; iface++)
			{
				// if face corresponding to local node iface does not exist in newfaces, add the face.
				if(val[iface] == false)
				{
					for(int j = 0; j < nnode-1; j++)	
						fc.p[j] = nw.p[lpofa.get(iface,j)];

					fc.elem[0] = elems.size()-1;
					fc.elem[1] = -4;			// the element to the right of the new face is undefined as of now.
					
					// Add this new face to list of faces and the current list of newfaces.
					faces.push_back(fc);
					newfaces.push_back(faces.size()-1);
				}
			}

			// Surrounding element of this new element - across pre-existing face
			//elems.back().surr[0] = (faces[voidpoly[ifa]].elem[0] != elems.size()-1) ? faces[voidpoly[ifa]].elem[1] : faces[voidpoly[ifa]].elem[0];
			elems.back().surr[0] = faces[voidpoly[ifa]].elem[0];

			// Now to set the new element as a surrounding element of the element neighboring this void face
			int nbor = elems.back().surr[0];
			
			if(nbor >= 0)
			{
				localface = check_face_tet(elems[nbor], faces[voidpoly[ifa]]);
				if(localface < 0) cout << "Delaunay3d: bowyer_watson(): ! Error in locating void face on pre-existing element!" << endl;
				elems[nbor].surr[localface] = elems.size()-1;
			}
		}

		//cout << "Delaunay2D:  Added point " << ipoin << ".\n";
		//empty badelems and voidpoly
		badelems.clear();
		voidpoly.clear();
	} 
	
	// end iteration over points

	cout << "Delaunay3d: bowyer_watson(): Number of elements before removing super points is " << elems.size() << endl;
	// Remove super triangle
	//cout << "Delaunay2D:  Remove super triangle\n";
	for(int ielem = 0; ielem < elems.size(); ielem++)
	{
		//vector<bool> val(nnode,false);
		bool finval = false;
		for(int i = 0; i < nnode; i++)
		{
			if(elems[ielem].p[i] == 0 || elems[ielem].p[i] == 1 || elems[ielem].p[i] == 2 || elems[ielem].p[i] == 3) {
				finval = true;
				break;
			}
		}

		// if finval is true, it means ielem contains one of the first 4 nodes
		if(finval)
		{
			elems.erase(elems.begin()+ielem);
			// re-adjust surr[] of each element. This is yet another place where we would benefit from a graph data structure.
			for(int iel = 0; iel < elems.size(); iel++)
			{
				for(int j = 0; j < nnode; j++)
				{
					if(elems[iel].surr[j] == ielem) elems[iel].surr[j] = -1;
					if(elems[iel].surr[j] > ielem) elems[iel].surr[j]--;
				}
			}
			/// \note Perhaps we should readjust the two elems of each face while deleting elements containing the super nodes.

			ielem--;
		}
	}
	
	// remove super nodes
	nodes.erase(nodes.begin(),nodes.begin()+nnode);
	
	for(int ielem = 0; ielem < elems.size(); ielem++)
	{
		for(int i = 0; i < nnode; i++)
		{
			elems[ielem].p[i] = elems[ielem].p[i]-nnode;
		}
	}
	// remove super faces - not needed
	cout << "Delaunay3d: Triangulation done.\n";

	// print surrounding elements
	/*for(int i = 0; i < elems.size(); i++)
	{
		cout << "Element " << i << ": ";
		for(int j = 0; j < 3; j++)
		{
			cout << elems[i].surr[j] << " ";
		}
		cout << endl;
	}*/
}

void Delaunay3d::clear()					// reset the Delaunay2D object, except for input data
{
	nodes.clear();
	elems.clear();
	faces.clear();
	badelems.clear();
	voidpoly.clear();
}

void Delaunay3d::writeGmsh2(string mfile)
{
	ofstream outf(mfile);

	outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
	outf << "$Nodes\n" << nodes.size() << '\n';
	for(int ip = 0; ip < nodes.size(); ip++)
	{
		outf << ip+1 << " " << nodes[ip][0] << " " << nodes[ip][1] << " " << nodes[ip][2] << '\n';
	}
	outf << "$Elements\n" << elems.size() << '\n';
	for(int iel = 0; iel < elems.size(); iel++)
	{
		outf << iel+1 << " 4 2 0 2";
		for(int i = 0; i < nnode; i++)
			outf << " " << elems[iel].p[i]+1;
		outf << '\n';
	}
	outf << "$EndElements\n";

	outf.close();
}

Walkdata Delaunay3d::find_containing_tet_and_barycentric_coords(const vector<double>& xx, int startelement) const
/* Note that the local node numbering is not assumed to be consistent. So checking the sign of the area of the triangle formed by the new point and an edge is not enough.
   Rather, we compare the corresponding area-ratio. If the sign of the area of the triangle created by the new point changes because of opposite orientation, so does the area of the triangle being checked. */
{
	Walkdata dat;
	if(xx.size() < 3) {
		std::cout << "Delaunau3D: find_containing_triangle(): ! Input vector is not long enough!\n";
		return dat;
	}
	int ielem = startelement;
	Tet super;
	vector<double> l(ndim+1);
	bool found;
	cout << "Delaunay3D:   Finding containing tet..." << endl;
	
	while(1)
	{
		found = true;

		if(ielem < 0 || ielem >= elems.size()) { cout << "Delaunay3d:   !! Reached an element index that is out of bounds!! Index is " << ielem << "\n"; }
		super = elems[ielem];

		for(int inode = 0; inode < nnode; inode++)
		{
			// get jacobian
			l[inode] = det4(ielem,inode,xx);
			#if DEBUGW==1
			if(dabs(l[inode]) < ZERO_TOL) cout << "Delaunay3D: find_containing_tet(): ! Degenerate case (type 1) for l " << inode << "!!\n";
			#endif
			if(l[inode]/super.D < 0)
			{
				ielem = super.surr[inode];
				found = false;
				break;
			}
		}

		// if all 4 area-ratios are positive, we've found our element
		if(found) break;
	}
	//cout << "Delaunay2D:   Containing triangle found: " << ielem << ".\n";
	dat.elem = ielem; 
	for(int inode = 0; inode < nnode; inode++)
		dat.areacoords[inode] = l[inode]/super.D;
	return dat;
}

void Delaunay3d::compute_jacobians()
{
	jacobians.setup(elems.size(),1);
	for(int i = 0; i < elems.size(); i++)
	{
		compute_jacobian(elems[i]);
		jacobians(i) = elems[i].D;
	}
}

bool Delaunay3d::detect_negative_jacobians()
{
	bool flagj = false;
	for(int i = 0; i < elems.size(); i++)
	{
		if(jacobians(i,0) < ZERO_TOL) {
			//out << i << " " << jacobians(i,0) << '\n';
			flagj = true;
		}
	}
	if(flagj == true) cout << "Delaunay2D: detect_negative_jacobians(): There exist element(s) with negative jacobian!!\n";
	return flagj;
}
