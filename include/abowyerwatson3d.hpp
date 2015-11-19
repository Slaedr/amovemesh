/**
A class for Delaunay tetrahedralization of a given set of points based loosely on Bowyer-Watson algorithm; also referenced from Wikipedia's Bowyer-Watson algorithm page.
Aditya Kashi
November 3, 2015

Notes:
  Currently uses a std::vector to store elements, faces etc. 
  It might be better to use a std::list or std::forward_list instead.
  A graph data structure should also be seriously considered for storing elements, faces, bad elements and the void polygon.
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
#define DEBUGBW 0
#endif

using namespace std;
using namespace amat;

typedef vector<double> Point;

/** Structure for triangular face. */
struct Face
{
	int p[3];			// indices of endpoints of the face
	int elem[2];		// elements on either side of the face
};

/** Class representing a tetrahedron. */
class Tet
{
public:	
	int p[4];			///< Indices of vertices.
	Point centre;		///< Coords of circumcenter of the tet.
	int surr[4];		///< Indices of surrounding tets. Note that the neighbor corresponding to surr[3] is opposite the vertex p[3].
	double D;			///< 2*volume of tet.
	double radius;		///< Radius of circumcircle of tet.

	Tet() {
		centre.reserve(3,0.0);
	}

	Tet(const Tet& other) {
		for(int i = 0; i < 4; i++)
		{
			p[i] = other.p[i];
			surr[i] = other.surr[i];
		}
		centre = other.centre;
		D = other.D
		radius = other.radius;
	}
};

/** Intended to encapsulate data required by 'walk-through' algorothms.
*  Not needed for mesh generation, but in independent application of the walk-through subroutine find_containing_triangle_and_area_coords().
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
public:
	Matrix<double> points;
	std::vector<Point> nodes;
	std::vector<Tet> elems;
	std::vector<int> badelems;		// collection of 'bad elements', that are to be removed while adding a point; members index 'elems'
	std::vector<Face> faces;
	std::vector<int> voidpoly;		// collection of faces that bounds the void obtained after removing bad elements while adding a point; members index 'faces'
	Matrix<double> jacobians;

	int npoints;

	Delaunay3d();
	Delaunay3d(Matrix<double>* _points, int num_points);
	Delaunay3d(const Delaunay3d& other);
	Delaunay3d& operator=(const Delaunay3d& other);
	void setup(Matrix<double>* _points, int num_points);
	
	void compute_jacobian(Tet& elem);
	void compute_circumsphere(Tet& elem);
	double det4(int ielem, int i, vector<double> r) const;
	int find_containing_triangle(const vector<double>& r, int startelement) const;
	int check_face_tet(const Tet& elem, const Face& face) const;

	void bowyer_watson();
	
	void clear();					///< Reset the Delaunay3d object, except for input data
	void writeGmsh2(string mfile);
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
}

/** Computes 2*volume of any tetrahedron. */
void Delaunay3d::compute_jacobian(Tet& elem)
{	
	double ret;
	ret = elem.p[0](0) * ( elem.p[1][1]*(elem.p[2][2]-elem.p[3][2]) -elem.p[1][2]*(elem.p[2][1]-elem.p[3][1]) + elem.p[2][1]*elem.p[3][2] - elem.p[2][2]*elem.p[3][1] );
	ret -= elem.p[0][1] * (elem.p[1][0]*(elem.p[2][2]-elem.p[3][2]) -elem.p[1][2]*(elem.p[2][0]-elem.p[3][0]) + elem.p[2][0]*elem.p[3][2] - elem.p[2][2]*elem.p[3][0] );
	ret += elem.p[0][2] * (elem.p[1][0]*(elem.p[2][1]-elem.p[3][1]) -elem.p[1][1]*(elem.p[2][0]-elem.p[3][0]) + elem.p[2][0]*elem.p[3][1] - elem.p[2][1]*elem.p[3][0] );
	ret -= elem.p[1][0]*( elem.p[2][1]*elem.p[3][2] - elem.p[2](2)*elem.p[3][1] ) -elem.p[1][1]*( elem.p[2][0]*elem.p[3][2] - elem.p[2][2]*elem.p[3][0] ) +elem.p[1][2]*( elem.p[2][0]*elem.p[3][1] - elem.p[2][1]*elem.p[3][0] );
	elem.D = ret;
}

/** Computes cross product c of two 3-vectors a and b.*/
inline void Delaunay3d::cross_product3(vector<double>& c, const vector<double>& a, const vector<double>& b)
{
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
}

inline double Delaunay3d::dot(const vector<double>& a, const vector<double>& b)
{
	double val = 0;
	for(int i = 0; i < ndim; i++)
		val += a[i]*b[i];
	return val;
}

inline double Delaunay3d::l2norm(const vector<double>& a)
{
	// a.size() should ideally equal ndim
	double norm = 0;
	for(int i = 0; i < a.size(); i++)
		norm += a[i]*a[i];
	norm = sqrt(norm);
	return norm;
}

/** Computes circumcentre and circumradius of a [tetrahedron](@ref Tet).
* NOTE: The tetrahedron's jacobian (2*volume) should be stored in elem.D beforehand.
*/
void Delaunay3d::compute_circumcircle(Tet& elem)
{
	vector<double> a(ndim), b(ndim), c(ndim), n1(idim), n2(idim), n3(idim), fin(idim);
	for(int i = 0; i < ndim; i++)
	{
		a[idim] = elem.p[0](idim)-elem.p[3](idim);
		b[idim] = elem.p[1](idim)-elem.p[3](idim);
		c[idim] = elem.p[2](idim)-elem.p[3](idim);
	}

	cross_product(n1,b,c);
	cross_product(n2,c,a);
	cross_product(n3,a,b);
	for(int i = 0; i < ndim; i++)
	{
		fin[idim] = dot(a,a)*n1[idim] + dot(b,b)*n2[idim] + dot(c,c)*n3[idim];
		fin[idim] /= 6*elem.D;
		elem.centre(idim) = fin[idim];
	}
	elem.radius = l2norm(fin);
}

/**	Calculates the jacobian of the tetrahedron formed by point r and a face of tetrahedron ielem. The face is selected by i between 0 and 3.
	Face i is the face opposite to local node i of the tetrahedron.
*/
double Delaunay3d::det4(int ielem, int i, vector<double> r) const
{
	#if DEBUGW==1
	if(i > 3) {
		std::cout << "Delaunay3D: det4(): ! Second argument is greater than 3!" << std::endl;
		return 0;
	}
	#endif

	double ret = 0;
	Tet elem = elems[ielem];
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
	return ret;
}

int Delaunay3d::find_containing_tet(const vector<double>& xx, int startelement) const
{
	if(xx.size() < 3) {
		std::cout << "Delaunau3D: find_containing_triangle(): ! Input vector is not long enough!\n";
		return -1;
	}
	int ielem = startelement;
	double l;
	bool found;
	cout << "Delaunay3D:   Finding containing triangle..." << endl;
	
	while(1)
	{
		found = true;

		if(ielem < 0 || ielem >= elems.size()) { cout << "Delaunay3d:   !! Reached an element index that is out of bounds!! Index is " << ielem << "\n"; }
		Tet super = elems[ielem];

		for(int inode = 0; inode < nnode; inode++)
		{
			// get jacobian
			l = det4(ielem,inode,xx);
			#if DEBUGW==1
			if(dabs(l) < ZERO_TOL) cout << "Delaunay3D: find_containing_tet(): ! Degenerate case (type 1) for l " << inode << "!!\n";
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
	//cout << "Delaunay2D:   Containing triangle found.\n";
	return ielem;
}

/** Returns the local face number (number of node opposite to the face) of elem's face that is the same as the second argument.
*  If no faces of elem match, it returns -1.
*/
inline int check_face_tet(const Tet& elem, const Face& face) const
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

/** Computes the Delaunay triangulation (tetrahedralization, in this case).
* Make sure 'points' has space for three more points when passing to this sub. 'N' is the actual number of real points.
*/
void Delaunay3d::bowyer_watson()
{
	// add super triangle
	//	find minimum and maximum x and y of the point set
	vector<double> rmin(ndim), rmax(ndim);
	for(int i = 1; i < npoints; i++)
	{
		for(int idim = 0; idim < ndim; idim++) {
			if(points(i,idim) > rmax[idim]) rmax[idim] = points(i,idim);
			if(points(i,idim) < rmin[idim]) rmin[idim] = points(i,idim);
		}
	}
	//cout << "Delaunay3D: bowyer_watson(): xmax, xmin, ymax, ymin " << xmax << " " << xmin << " " << ymax << " " << ymin << '\n';
	vector<double> rdelt(ndim);		// stores some extra length to have factor of safety in deciding the first 4 points
	for(int idim = 0; idim < ndim; idim++)
		rdelt[idim] = 1.0+(rmax[idim]-rmin[idim]);
	Tet super;
	vector<Point> pp(nnode);
	
	pp[0](0) = rmin[0]-rdelt[0]; pp[0](1) = rmax[1] + rdelt[1]; pp[0](2) = rmin[2] - rdelt[2];

	pp[1](0) = rmin[0] + 2*(rmax[0]-rmin[0]) + 2*rdelt[0]; pp[1](1) = pp[0](1); pp[1](2) = pp[0](1);

	pp[2](0) = pp[0](0); pp[2](1) = pp[0](1); pp[2](2) = rmin[2] + 2*(rmax[2] - rmin[2]) + 2.0*rdelt[2];

	pp[3](0) = pp[0](0); pp[3](1) = rmax[1] - 2.0*(rmax[1]-rmin[1]) - rdelt[1]; pp[3](2) = pp[0](2);

	cout << "Delaunay3d: bowyer_watson(): Coordinates of vertices of the super triangle are:\n";
	for(int inode = 0; inode < nnode; inode++)
	{
		for(int idim = 0; idim < ndim; idim++)
			cout << pp[inode](idim) << "\t";
		cout << endl;
	}
	
	for(int i = 0; i < nnode; i++)
	{
		nodes.push_back(pp[i]);
		super.p[i] = i;
	}
	for(int inode = 0; inode < nnode; inode++)
		super.surr[inode] = -1;

	super.D = compute_jacobian(super);
	super.compute_circumcircle();

	elems.push_back(super);			// add super to elems list

	// set up initial face list
	//cout << "Delaunay3d: set up initial face list\n";
	Face f[ndim+1];
	f[0].p[0] = super.p[2]; f[0].p[1] = super.p[2]; f[0].p[2] = super.p[1];
	f[1].p[0] = super.p[2]; f[1].p[1] = super.p[0]; f[1].p[2] = super.p[3];
	f[2].p[0] = super.p[3]; f[2].p[1] = super.p[1]; f[2].p[2] = super.p[2];
	f[3].p[0] = super.p[3]; f[3].p[1] = super.p[0]; f[3].p[2] = super.p[1];
	for(int i = 0; i < 3; i++)
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
		for(int idim = 0; idim < ndim; idim++)
			newpoin[idim] = points.get(ipoin,idim);
		nodes.push_back(newpoin);
		newpoinnum = nodes.size()-1;

		/// First, find the element containing the new point
		contelem = find_containing_tet(newpoin,0);
		//cout << "Delaunay2D:  Containing element is " << contelem << endl;

		/// Second, search among neighbors for other triangles whose circumcircles contain this point
		//cout << "Delaunay2D:  Second, search among neighbors for other triangles whose circumcircles contain this point\n";
		int curelem;
		vector<int> stk;					// stack to hold the indices of tets to be checked
		double dist;						// square of distance from point to circumcentre
		stk.push_back(contelem);			// add the containing element to the list of bad elements
		vector<int> flags(elems.size());	// stores 1 at an index if the corresponding element has been checked for the Delaunay criterion
		for(int i = 0; i < elems.size(); i++) flags[i] = 0;

		while(stk.empty() == false)
		{
			curelem = stk.back();			// access last element in stack of elements to be checked
			flags[curelem] = 1;				// curelem will now be checked

			//calculate distance between circumcentre and the point
			dist = 0;
			for(int idim = 0; idim < ndim; idim++)
				dist += (newpoin[idim] - elems[curelem].centre[idim])*(newpoin[idim] - elems[curelem].centre[idim]);

			#if DEBUGBW==1
			if(dabs(dist - elems[curelem].radius*elems[curelem].radius) < ZERO_TOL) cout << "Delaunay2D: Degenerate case (type 2)!!\n";
			#endif
			if(dist < elems[curelem].radius*elems[curelem].radius)		// if point lies inside circumcircle, ie, Delaunay criterion is violated
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
		/*cout << "Delaunay2D:  Badelems: ";
		for(int i = 0; i < badelems.size(); i++)
			cout << badelems[i] << " ";
		cout << endl;*/

		/// Third, store the faces that will be obtained after removal of bad triangles
		//cout << "Delaunay2D:  Third, store the faces that will be obtained after removal of bad triangles\n";
		flags.assign(faces.size(),-1);
		/*cout << "**Flags : (" << faces.size() << ", " << flags.size() << ")";
		for(int i = 0; i < flags.size(); i++)
			cout << " " << flags[i];
		cout << endl;
		cout << "** Details of face 0: " << faces[0].p[0] << " " << faces[0].p[1] << ", " << faces[0].elem[0] << " " << faces[0].elem[1] << endl;*/

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
		/*cout << "Delaunay2D:  voidpoly:";
		for(int i = 0; i < voidpoly.size(); i++)
			cout << " " << voidpoly[i];
		cout << endl;*/

		/// Fifth, add new elements; these are formed by the faces in voidpoly and the new point. Also correspondingly update 'faces'.
		//cout << "Delaunay2D:  Fifth, add new elements; Also correspondingly update 'faces'\n";
		vector<int> newfaces;				// new faces formed from new elements created
		for(int ifa = 0; ifa < voidpoly.size(); ifa++)		// for each face in void polygon
		{
			Tet nw;
			nw.p[0] = newpoinnum;
			if(faces[voidpoly[ifa]].elem[0] == -2)	
			{	// if the new element is to the left of this face, the orientation of the remaining points of the tet is same as the orientation of the face
				nw.p[1] = faces[voidpoly[ifa]].p[0];
				nw.p[2] = faces[voidpoly[ifa]].p[1];
				nw.p[3] = faces[voidpoly[ifa]].p[2];
				faces[voidpoly[ifa]].elem[0] = elems.size();		// this new element has not been pushed into elems yet, so we need elems.size() - 1 + 1
			}
			else if(faces[voidpoly[ifa]].elem[1] == -2)	// if the new element is to the right of the face
			{
				nw.p[1] = faces[voidpoly[ifa]].p[1];
				nw.p[2] = faces[voidpoly[ifa]].p[0];
				nw.p[3] = faces[voidpoly[ifa]].p[2];
				faces[voidpoly[ifa]].elem[1] = elems.size();
			}
			else cout << "Delaunay2D: !! Error while creating new element - face " << voidpoly[ifa] << " in voidpoly does not have -2 as either left elem or right elem.\n";

			nw.D = compute_jacobian(nw);
			compute_circumcircle(nw);

			//cout << "Delaunay3d:  New element circumcentre and radius: " << nw.centre.x << ", " << nw.centre.y << "; " << nw.radius << endl;
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
				if(faces[newfaces[jfa]].elem[1] == -4)		// Do we really need this check? I don't think so.
				{
					faces[newfaces[jfa]].elem[1] = elems.size()-1;
					elems.back().surr[localface] = faces[newfaces[jfa]].elem[0];

					// Also, find which local face of newfaces[jfa]'s left element is the same as newfaces[jfa].
					// Set the new element as a surrounding element of of the left element of newfaces[jfa].
					jface = check_face_tet(elems[faces[newfaces[jfa]].elem[0]], faces[newfaces[jfa]]);
					if(jface == -1) {cout << "Delaunay3d: bowyer_watson(): ! Error while setting surrounding element of a new element!" << endl; }
					elems[faces[newfaces[jfa]].elem[0]].surr[jface] = elems.size()-1;
				}

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
					fc.p[0] = nw.p[0];
					fc.p[1] = nw.p[3];
					fc.p[2] = nw.p[2];
					fc.elem[0] = elems.size()-1;
					fc.elem[1] = -4;			// the element to the right of the new face is undefined as of now.
					// Add this new faces to list of faces and the current list of newfaces.
					faces.push_back(fc);
					newfaces.push_back(faces.size()-1);
				}
			}

			// Surrounding element of this new element - across pre-existing face
			elems.back().surr[0] = (faces[voidpoly[ifa]].elem[0] != elems.size()-1) ? faces[voidpoly[ifa]].elem[0] : faces[voidpoly[ifa]].elem[1];

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
	} // end iteration over points

	// Remove super triangle
	//cout << "Delaunay2D:  Remove super triangle\n";
	for(int ielem = 0; ielem < elems.size(); ielem++)
	{
		//vector<bool> val(nnode,false);
		bool finval = false;
		for(int i = 0; i < nnode; i++)
		{
			if(elem[ielem].p[nnode] == 0 || elem[ielem].p[nnode] == 1 || elem[ielem].p[nnode] == 2 || elem[ielem].p[nnode] == 3) {
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
				for(int j = 0; j < 3; j++)
				{
					if(elems[iel].surr[j] == ielem) elems[iel].surr[j] = -1;
					if(elems[iel].surr[j] > ielem) elems[iel].surr[j]--;
				}
			}
			//NOTE: Perhaps we should readjust the two elems of each face.

			ielem--;
		}
	}
	// remove super nodes
	nodes.erase(nodes.begin(),nodes.begin()+3);
	for(int ielem = 0; ielem < elems.size(); ielem++)
	{
		for(int i = 0; i < nnode; i++)
		{
			elems[ielem].p[i] = elems[ielem].p[i]-3;
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

void Delaunay3D::clear()					// reset the Delaunay2D object, except for input data
{
	nodes.clear();
	elems.clear();
	faces.clear();
	badelems.clear();
	voidpoly.clear();
}

void Delaunay3D::writeGmsh2(string mfile)
{
	ofstream outf(mfile);

	outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
	outf << "$Nodes\n" << npoints << '\n';
	for(int ip = 0; ip < npoints; ip++)
	{
		outf << ip+1 << " " << nodes[ip][0] << " " << nodes[ip][1] << " " << nodes[ip][2] << '\n';
	}
	outf << "$Elements\n" << elems.size() << '\n';
	for(int iel = 0; iel < elems.size(); iel++)
	{
		outf << iel+1 << " 3 2 0 2";
		for(int i = 0; i < nnode; i++)
			outf << " " << elems[iel].p[i]+1;
		outf << '\n';
	}
	outf << "$EndElements\n";

	outf.close();
	//cout << "Delaunay2D: Number of faces finally = " << faces.size() << endl;
}

Walkdata Delaunay3d::find_containing_tet_and_barycentric_coords(vector<double> xx, int startelement) const
/* Note that the local node numbering is not assumed to be consistent. So checking the sign of the area of the triangle formed by the new point and an edge is not enough.
   Rather, we compare the corresponding area-ratio. If the sign of the area of the triangle created by the new point changes because of opposite orientation, so does the area of the triangle being checked. */
{
	Walkdata dat;
	if(xx.size() < 3) {
		std::cout << "Delaunau3D: find_containing_triangle(): ! Input vector is not long enough!\n";
		return -1;
	}
	int ielem = startelement;
	vector<double> l(ndim+1);
	bool found;
	cout << "Delaunay3D:   Finding containing tet..." << endl;
	
	while(1)
	{
		found = true;

		if(ielem < 0 || ielem >= elems.size()) { cout << "Delaunay3d:   !! Reached an element index that is out of bounds!! Index is " << ielem << "\n"; }
		Tet super = elems[ielem];

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

void Delaunay3D::compute_jacobians()
{
	jacobians.setup(elems.size(),1);
	for(int i = 0; i < elems.size(); i++)
	{
		elems[i].D = compute_jacobian(elems[i]);
		jacobians(i) = elems[i].D;
	}
}

bool Delaunay3D::detect_negative_jacobians()
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
