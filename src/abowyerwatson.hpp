/** \file abowyerwatson.hpp
 * \brief Delaunay triangulation of a given set of points, based on Bowyer-Watson algorithm; also referenced from Wikipedia's Bowyer-Watson algorithm page.
 * \author Aditya Kashi
 * \date June 22, 2015
 * 
 * Notes: \note
 * Currently uses a std::vector to store elements. It is probably better to use a std::list or std::forward_list instead.
 * 
 * Changelog:
 * June 26, 2015: Orientation is now preserved.
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
#ifndef __AMATRIX_H
#include <amatrix.hpp>
#endif

#define __ABOWYERWATSON_H 1

namespace amc {

struct Point2
{
	double x;
	double y;
};

/// Holds indices of points that make up a face and the left and right elements
struct Face
{
	int p[2];			///< indices of endpoints of the face
	int elem[2];		///< elements on either side of the face
};

/// Contains data characterizing a 2D Delaunay element
struct Triangle
{
	int p[3];			///< indices of vertices of triangle
	Point2 centre;		///< coords of circumcentre of triangle
	int surr[3];		///< indices of surrounding triangles
	double D;			///< 2*area of triangle
	double radius;		///< radius of circumcircle of triangle
};

/// Holds an element index and the barycentric coordinates of point w.r.t that element
/** Not needed for mesh generation, but in independent application of the walk-through subroutine find_containing_triangle_and_area_coords()
 */
struct Walkdata
{
	int elem;
	double areacoords[3];
};

/// Delaunay triangulation of a set of 2D points
class Delaunay2D
{
	int cap;
	int badcap;
	double tol;
	int nnode;
public:
	amat::Matrix<double> points;
	std::vector<Point2> nodes;
	std::vector<Triangle> elems;
	std::vector<int> badelems;		// collection of 'bad elements', that are to be removed while adding a point; members index 'elems'
	std::vector<Face> faces;
	std::vector<int> voidpoly;		// collection of faces that bounds the void obtained after removing bad elements while adding a point; members index 'faces'
	amat::Matrix<double> jacobians;

	int npoints;

	Delaunay2D() {}

	Delaunay2D(amat::Matrix<double>* _points)
	{
		nnode = 3;
		cap = 1000;
		badcap = 50;
		tol = A_SMALL_NUMBER;
		elems.reserve(cap);
		faces.reserve(cap);
		//nelem = 0;
		points = *_points;
		npoints = _points->rows();
		nodes.reserve(npoints+3);
	}

	Delaunay2D(const Delaunay2D& other)
	{
		cap = other.cap;
		badcap = other.badcap;
		tol = other.tol;
		nnode = other.nnode;
		points = other.points;
		nodes = other.nodes;
		elems = other.elems;
		faces = other.faces;
		jacobians = other.jacobians;
		npoints = other.npoints;
	}

	Delaunay2D& operator=(const Delaunay2D& other)
	{
		cap = other.cap;
		badcap = other.badcap;
		tol = other.tol;
		nnode = other.nnode;
		points = other.points;
		nodes = other.nodes;
		elems = other.elems;
		faces = other.faces;
		jacobians = other.jacobians;
		npoints = other.npoints;
		return *this;
	}

	void setup(const amat::Matrix<double>* _points)
	{
		nnode = 3;
		cap = 1000;
		badcap = 50;
		tol = A_SMALL_NUMBER;
		elems.reserve(cap);
		faces.reserve(cap);
		points = *_points;
		npoints = _points->rows();
		nodes.reserve(npoints+3);
	}

	/// Walk through the mesh from element to element, until the element containing a given point is found.
	/** NOTE: It will be better to find the minimum of the area coordinates and move to the next element according to that, as done in the 3D code.
	 */
	int find_containing_triangle(double xx, double yy, int startelement)
	{
		int ielem = startelement;
		double l1, l2, l3;
		//std::cout << "Delaunay2D:   Finding containing triangle...\n";
		while(1)
		{
			if(ielem < 0 || ielem >= elems.size()) { std::cout << "Delaunay2D:   !! Reached an element index that is out of bounds!! Index is " << ielem << "\n"; }
			Triangle super = elems[ielem];

			l3 = xx*(nodes[super.p[0]].y - nodes[super.p[1]].y) - yy*(nodes[super.p[0]].x - nodes[super.p[1]].x) + nodes[super.p[0]].x*nodes[super.p[1]].y - nodes[super.p[1]].x*nodes[super.p[0]].y;
			#if DEBUG==1
			if(fabs(l3) < tol) std::cout << "Delaunay2D:   Degenerate case (type 1) l3!!\n";
			#endif
			if(l3/super.D < 0)
			{
				ielem = super.surr[2];
				continue;
			}

			l1 = xx*(nodes[super.p[1]].y - nodes[super.p[2]].y) - yy*(nodes[super.p[1]].x - nodes[super.p[2]].x) + nodes[super.p[1]].x*nodes[super.p[2]].y - nodes[super.p[2]].x*nodes[super.p[1]].y;
			#if DEBUG==1
			if(fabs(l1) < tol) std::cout << "Delaunay2D:   Degenerate case (type 1) l1!!\n";
			#endif
			if(l1/super.D < 0)
			{
				ielem = super.surr[0];
				continue;
			}

			l2 = xx*(nodes[super.p[2]].y - nodes[super.p[0]].y) - yy*(nodes[super.p[2]].x - nodes[super.p[0]].x) + nodes[super.p[2]].x*nodes[super.p[0]].y - nodes[super.p[0]].x*nodes[super.p[2]].y;
			#if DEBUG==1
			if(fabs(l2) < tol) std::cout << "Delaunay2D:   Degenerate case (type 1) l2!!\n";
			#endif
			if(l2/super.D < 0)
			{
				ielem = super.surr[1];
				continue;
			}
			break;		// if all three area-ratios are positive, we've found our element
		}
		//std::cout << "Delaunay2D:   Containing triangle found.\n";
		return ielem;
	}

	/// Implements the Bowyer-Watson algorithm
	void bowyer_watson()
	/* Make sure 'points' has space for three more points when passing to this sub. 'N' is the actual number of real points. */
	{
		// add super triangle
		//find minimum and maximum x and y of the point set
		double xmin = points(0,0), xmax = points(0,0), ymin = points(0,1), ymax = points(0,1);
		for(int i = 1; i < npoints; i++)
		{
			if(points(i,0) > xmax) xmax = points(i,0);
			if(points(i,0) < xmin) xmin = points(i,0);
			if(points(i,1) > ymax) ymax = points(i,1);
			if(points(i,1) < ymin) ymin = points(i,1);
		}
		//std::cout << "Delaunay2D: bowyer_watson(): xmax, xmin, ymax, ymin " << xmax << " " << xmin << " " << ymax << " " << ymin << '\n';

		double xdelt = 1+(xmax-xmin), ydelt = 1+(ymax-ymin);
		double delt = (xdelt >= ydelt) ? xdelt : ydelt;
		Triangle super;
		Point2 pp[3];
		pp[0].x = xmin-xdelt; pp[0].y = ymax+ydelt;
		pp[1].x = xmin-xdelt; pp[1].y = ymax - 2*(ymax-ymin) - 2*ydelt;
		pp[2].x = xmin + 2*(xmax-xmin) + 2*xdelt; pp[2].y = ymax + ydelt;
		for(int i = 0; i < 3; i++)
		{
			nodes.push_back(pp[i]);
			super.p[i] = i;
		}
		super.surr[0] = super.surr[1] = super.surr[2] = -1;

		super.D = nodes[super.p[0]].x*(nodes[super.p[1]].y - nodes[super.p[2]].y) - nodes[super.p[0]].y*(nodes[super.p[1]].x - nodes[super.p[2]].x) + nodes[super.p[1]].x*nodes[super.p[2]].y - nodes[super.p[2]].x*nodes[super.p[1]].y;

		double Dd = 2*(nodes[super.p[0]].x*(nodes[super.p[1]].y-nodes[super.p[2]].y) + nodes[super.p[1]].x*(nodes[super.p[2]].y-nodes[super.p[0]].y) + nodes[super.p[2]].x*(nodes[super.p[0]].y-nodes[super.p[1]].y));
		super.centre.x = ((nodes[super.p[0]].x*nodes[super.p[0]].x+nodes[super.p[0]].y*nodes[super.p[0]].y)*(nodes[super.p[1]].y-nodes[super.p[2]].y) + (nodes[super.p[1]].x*nodes[super.p[1]].x+nodes[super.p[1]].y*nodes[super.p[1]].y)*(nodes[super.p[2]].y-nodes[super.p[0]].y) + (nodes[super.p[2]].x*nodes[super.p[2]].x+nodes[super.p[2]].y*nodes[super.p[2]].y)*(nodes[super.p[0]].y-nodes[super.p[1]].y)) / Dd;
		super.centre.y = ((nodes[super.p[0]].x*nodes[super.p[0]].x+nodes[super.p[0]].y*nodes[super.p[0]].y)*(nodes[super.p[2]].x-nodes[super.p[1]].x) + (nodes[super.p[1]].x*nodes[super.p[1]].x+nodes[super.p[1]].y*nodes[super.p[1]].y)*(nodes[super.p[0]].x-nodes[super.p[2]].x) + (nodes[super.p[2]].x*nodes[super.p[2]].x+nodes[super.p[2]].y*nodes[super.p[2]].y)*(nodes[super.p[1]].x-nodes[super.p[0]].x)) / Dd;

		double a2 = (nodes[super.p[0]].x-nodes[super.p[1]].x)*(nodes[super.p[0]].x-nodes[super.p[1]].x) + (nodes[super.p[0]].y-nodes[super.p[1]].y)*(nodes[super.p[0]].y-nodes[super.p[1]].y);
		double b2 = (nodes[super.p[1]].x-nodes[super.p[2]].x)*(nodes[super.p[1]].x-nodes[super.p[2]].x) + (nodes[super.p[1]].y-nodes[super.p[2]].y)*(nodes[super.p[1]].y-nodes[super.p[2]].y);
		double c2 = (nodes[super.p[2]].x-nodes[super.p[0]].x)*(nodes[super.p[2]].x-nodes[super.p[0]].x) + (nodes[super.p[2]].y-nodes[super.p[0]].y)*(nodes[super.p[2]].y-nodes[super.p[0]].y);
		super.radius = sqrt(a2*b2*c2)/(2*super.D);

		elems.push_back(super);			// add super to elems list

		// set up initial face list
		//std::cout << "Delaunay2D: set up initial face list\n";
		Face f[3];
		f[0].p[0] = super.p[1]; f[0].p[1] = super.p[2];
		f[1].p[0] = super.p[2]; f[1].p[1] = super.p[0];
		f[2].p[0] = super.p[0]; f[2].p[1] = super.p[1];
		for(int i = 0; i < 3; i++)
		{
			f[i].elem[0] = 0;
			f[i].elem[1] = -1;
			faces.push_back(f[i]);
		}

		// iterate through points
		std::cout << "Delaunay2D: Starting iteration over points\n";
		for(int ipoin = 0; ipoin < npoints; ipoin++)
		{
			Point2 newpoin; newpoin.x = points(ipoin,0); newpoin.y = points(ipoin,1);
			nodes.push_back(newpoin);
			int newpoinnum = nodes.size()-1;

			/// First, find the element containing the new point
			int contelem = find_containing_triangle(points(ipoin,0),points(ipoin,1), elems.size()-1);

#if DEBUG==1
			if(ipoin % 20 == 0)
				std::cout << "Delaunay2D:  Point " << ipoin << ", Containing element is " << contelem << std::endl;
#endif

			/// Second, search among neighbors for other triangles whose circumcircles contain this point
			int curelem;
			std::vector<int> stk;					// stack to hold the indices of triangles to be checked
			double dist;						// square of distance from point to circumcentre
			stk.push_back(contelem);			// add the containing element to the list of bad elements
			std::vector<int> flags(elems.size());	// stores 1 at an index if the corresponding element has been checked for the Delaunay criterion
			for(int i = 0; i < elems.size(); i++) flags[i] = 0;

			while(stk.empty() == false)
			{
				curelem = stk.back();			// access last element in stack of elements to be checked
				flags[curelem] = 1;				// curelem will now be checked

				//calculate distance between circumcentre and the point
				dist = (points(ipoin,0) - elems[curelem].centre.x)*(points(ipoin,0) - elems[curelem].centre.x) + (points(ipoin,1) - elems[curelem].centre.y)*(points(ipoin,1) - elems[curelem].centre.y);

				#if DEBUG==1
				if(fabs(dist - elems[curelem].radius*elems[curelem].radius) < tol) std::cout << "Delaunay2D: Degenerate case (type 2)!!\n";
				#endif
				if(dist < elems[curelem].radius*elems[curelem].radius)		// if point lies inside circumcircle, ie, Delaunay criterion is violated
				{
					badelems.push_back(curelem);
					stk.pop_back();
					for(int j = 0; j < 3; j++)
						if(elems[curelem].surr[j] >= 0 && flags[elems[curelem].surr[j]] == 0)		// add surrounding elements (only which are not checked) to "to be checked" stack
							stk.push_back(elems[curelem].surr[j]);
				}
				else
				{
					stk.pop_back();
				}
			}

			/// Third, store the faces that will be obtained after removal of bad triangles
			flags.assign(faces.size(),-1);

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
						if(flags[ifa] != 0)							// this face does NOT shared by two bad elements
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
			// delete faces that were shared by 2 bad triangles
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

			// Fourth, delete bad elements
			for(int ibe = 0; ibe < badelems.size(); ibe++)
			{
				elems.erase(elems.begin()+badelems[ibe]);

				//scan badelems for elements with indices greater than the one just deleted
				for(int i = ibe+1; i < badelems.size(); i++)
					if(badelems[i] > badelems[ibe]) badelems[i]--;

				//scan surrounding elements in elems -- probably not efficient
				for(int i = 0; i < elems.size(); i++)
				{
					for(int j = 0; j < 3; j++)
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

			// Fifth, add new elements; these are formed by the faces in voidpoly and the new point. Also correspondingly update 'faces'.
			std::vector<int> newfaces;				// new faces formed from new elements created
			for(int ifa = 0; ifa < voidpoly.size(); ifa++)		// for each face in void polygon
			{
				Triangle nw;
				nw.p[0] = newpoinnum;
				if(faces[voidpoly[ifa]].elem[0] == -2)		// if the new element is to the left of this face
				{
					nw.p[1] = faces[voidpoly[ifa]].p[0];
					nw.p[2] = faces[voidpoly[ifa]].p[1];
					faces[voidpoly[ifa]].elem[0] = elems.size();		// this new element has not been pushed into elems yet, so we need elems.size() - 1 + 1
				}
				else if(faces[voidpoly[ifa]].elem[1] == -2)	// if the new element is to the right of the face
				{
					nw.p[1] = faces[voidpoly[ifa]].p[1];
					nw.p[2] = faces[voidpoly[ifa]].p[0];
					faces[voidpoly[ifa]].elem[1] = elems.size();
				}
				else std::cout << "Delaunay2D: !! Error while creating new element - face " << voidpoly[ifa] << " in voidpoly does not have -2 as either left elem or right elem.\n";

				nw.D = nodes[nw.p[0]].x*(nodes[nw.p[1]].y - nodes[nw.p[2]].y) - nodes[nw.p[0]].y*(nodes[nw.p[1]].x - nodes[nw.p[2]].x) + nodes[nw.p[1]].x*nodes[nw.p[2]].y - nodes[nw.p[2]].x*nodes[nw.p[1]].y;

				double Dd = 2*(nodes[nw.p[0]].x*(nodes[nw.p[1]].y-nodes[nw.p[2]].y) + nodes[nw.p[1]].x*(nodes[nw.p[2]].y-nodes[nw.p[0]].y) + nodes[nw.p[2]].x*(nodes[nw.p[0]].y-nodes[nw.p[1]].y));
				nw.centre.x = ((nodes[nw.p[0]].x*nodes[nw.p[0]].x+nodes[nw.p[0]].y*nodes[nw.p[0]].y)*(nodes[nw.p[1]].y-nodes[nw.p[2]].y) + (nodes[nw.p[1]].x*nodes[nw.p[1]].x+nodes[nw.p[1]].y*nodes[nw.p[1]].y)*(nodes[nw.p[2]].y-nodes[nw.p[0]].y) + (nodes[nw.p[2]].x*nodes[nw.p[2]].x+nodes[nw.p[2]].y*nodes[nw.p[2]].y)*(nodes[nw.p[0]].y-nodes[nw.p[1]].y)) / Dd;
				nw.centre.y = ((nodes[nw.p[0]].x*nodes[nw.p[0]].x+nodes[nw.p[0]].y*nodes[nw.p[0]].y)*(nodes[nw.p[2]].x-nodes[nw.p[1]].x) + (nodes[nw.p[1]].x*nodes[nw.p[1]].x+nodes[nw.p[1]].y*nodes[nw.p[1]].y)*(nodes[nw.p[0]].x-nodes[nw.p[2]].x) + (nodes[nw.p[2]].x*nodes[nw.p[2]].x+nodes[nw.p[2]].y*nodes[nw.p[2]].y)*(nodes[nw.p[1]].x-nodes[nw.p[0]].x)) / Dd;

				double a2 = (nodes[nw.p[0]].x-nodes[nw.p[1]].x)*(nodes[nw.p[0]].x-nodes[nw.p[1]].x) + (nodes[nw.p[0]].y-nodes[nw.p[1]].y)*(nodes[nw.p[0]].y-nodes[nw.p[1]].y);
				double b2 = (nodes[nw.p[1]].x-nodes[nw.p[2]].x)*(nodes[nw.p[1]].x-nodes[nw.p[2]].x) + (nodes[nw.p[1]].y-nodes[nw.p[2]].y)*(nodes[nw.p[1]].y-nodes[nw.p[2]].y);
				double c2 = (nodes[nw.p[2]].x-nodes[nw.p[0]].x)*(nodes[nw.p[2]].x-nodes[nw.p[0]].x) + (nodes[nw.p[2]].y-nodes[nw.p[0]].y)*(nodes[nw.p[2]].y-nodes[nw.p[0]].y);
				nw.radius = sqrt(a2*b2*c2)/(2*nw.D);
				//std::cout << "Delaunay2D:  New element circumcentre and radius: " << nw.centre.x << ", " << nw.centre.y << "; " << nw.radius << std::endl;

				elems.push_back(nw);

				Face f1, f2; bool val1 = false, val2 = false;
				for(int jfa = 0; jfa < newfaces.size(); jfa++)
				{
					// Face formed by new point and p[1]
					if((faces[newfaces[jfa]].p[0] == nw.p[1] && faces[newfaces[jfa]].p[1] == nw.p[0]) || (faces[newfaces[jfa]].p[0] == nw.p[0] && faces[newfaces[jfa]].p[1] == nw.p[1])) // if the face formed by nw.p[0] and nw.p[1] exists in newfaces
					{
						if(faces[newfaces[jfa]].elem[1] == -4)
						{
							faces[newfaces[jfa]].elem[1] = elems.size()-1;
							elems.back().surr[2] = faces[newfaces[jfa]].elem[0];

							if(elems[faces[newfaces[jfa]].elem[0]].surr[0] == -2)
								if( (elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[0] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[1]) || (elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[1] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[0]) )
									elems[faces[newfaces[jfa]].elem[0]].surr[0] = elems.size()-1;		// set the new element as an element surrounding the neighboring element

							if(elems[faces[newfaces[jfa]].elem[0]].surr[1] == -2)
								if( (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[0] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[1]) || (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[1] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[0]) )
									elems[faces[newfaces[jfa]].elem[0]].surr[1] = elems.size()-1;		// set the new element as an element surrounding the neighboring element

							if(elems[faces[newfaces[jfa]].elem[0]].surr[2] == -2)
								if( (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[0] && elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[1]) || (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[1] && elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[0]) )
									elems[faces[newfaces[jfa]].elem[0]].surr[2] = elems.size()-1;		// set the new element as an element surrounding the neighboring element
						}
						val1 = true;
					}

					// Face formed by new point and p[2]
					if((faces[newfaces[jfa]].p[0] == nw.p[2] && faces[newfaces[jfa]].p[1] == nw.p[0]) || (faces[newfaces[jfa]].p[0] == nw.p[0] && faces[newfaces[jfa]].p[1] == nw.p[2])) 
						// if the face formed by nw.p[0] and nw.p[2] exists in newfaces
					{
						if(faces[newfaces[jfa]].elem[1] == -4)
						{
							faces[newfaces[jfa]].elem[1] = elems.size()-1;
							elems.back().surr[1] = faces[newfaces[jfa]].elem[0];

							if(elems[faces[newfaces[jfa]].elem[0]].surr[0] == -2)
								if( (elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[0] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[1]) || (elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[1] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[0]) )
									elems[faces[newfaces[jfa]].elem[0]].surr[0] = elems.size()-1;		// set the new element as an element surrounding the neighboring element

							if(elems[faces[newfaces[jfa]].elem[0]].surr[1] == -2)
								if( (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[0] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[1]) || (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[1] && elems[faces[newfaces[jfa]].elem[0]].p[2]==faces[newfaces[jfa]].p[0]) )
									elems[faces[newfaces[jfa]].elem[0]].surr[1] = elems.size()-1;		// set the new element as an element surrounding the neighboring element

							if(elems[faces[newfaces[jfa]].elem[0]].surr[2] == -2)
								if( (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[0] && elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[1]) || (elems[faces[newfaces[jfa]].elem[0]].p[0]==faces[newfaces[jfa]].p[1] && elems[faces[newfaces[jfa]].elem[0]].p[1]==faces[newfaces[jfa]].p[0]) )
									elems[faces[newfaces[jfa]].elem[0]].surr[2] = elems.size()-1;		// set the new element as an element surrounding the neighboring element
						}
						else std::cout << "Delaunay2D: !! Error!\n";
						val2 = true;
					}

					if(val1==true && val2==true) break;
				}

				if(val1==false)			// face formed by nw.p[0] and nw.p[1] does not exist in newfaces (or newfaces is empty)
				{
					f1.p[0] = nw.p[0]; f1.p[1] = nw.p[1];
					f1.elem[0] = elems.size()-1; f1.elem[1] = -4;
					faces.push_back(f1);
					val1 = true;
					newfaces.push_back(faces.size()-1);
					elems.back().surr[2] = -2;
				}
				if(val2==false)			// face formed by nw.p[0] and nw.p[2] does not exist in newfaces
				{
					f2.p[0] = nw.p[2]; f2.p[1] = nw.p[0];
					f2.elem[0] = elems.size()-1; f2.elem[1] = -4;
					faces.push_back(f2);
					val2 = true;
					newfaces.push_back(faces.size()-1);
					elems.back().surr[1] = -2;
				}

				// Surrounding element of this new element - across pre-existing face
				elems.back().surr[0] = (faces[voidpoly[ifa]].elem[0] != elems.size()-1) ? faces[voidpoly[ifa]].elem[0] : faces[voidpoly[ifa]].elem[1];

				// Now to set the new element as a surrounding element of the element neighboring this void face
				int nbor = elems.back().surr[0];
				if(nbor >= 0)
				{
					if( (elems[nbor].p[0] == faces[voidpoly[ifa]].p[0] && elems[nbor].p[1] == faces[voidpoly[ifa]].p[1]) || (elems[nbor].p[0] == faces[voidpoly[ifa]].p[1] && elems[nbor].p[1] == faces[voidpoly[ifa]].p[0]) )	// if the neighbor's face 2 matches the current void face
					{
						//if(elems[nbor].surr[2] == -2)
							elems[nbor].surr[2] = elems.size()-1;
					}
					else if((elems[nbor].p[1] == faces[voidpoly[ifa]].p[0] && elems[nbor].p[2] == faces[voidpoly[ifa]].p[1]) || (elems[nbor].p[1] == faces[voidpoly[ifa]].p[1] && elems[nbor].p[2] == faces[voidpoly[ifa]].p[0]))	// if the neighbor's face 0 matches the current void face
					{
						//if(elems[nbor].surr[0] == -2)
							elems[nbor].surr[0] = elems.size()-1;
					}
					else if((elems[nbor].p[2] == faces[voidpoly[ifa]].p[0] && elems[nbor].p[0] == faces[voidpoly[ifa]].p[1]) || (elems[nbor].p[2] == faces[voidpoly[ifa]].p[1] && elems[nbor].p[0] == faces[voidpoly[ifa]].p[0]))	// if the neighbor's face 1 matches the current void face
					{
						//if(elems[nbor].surr[1] == -2)
							elems[nbor].surr[1] = elems.size()-1;
					}
					else std::cout << "Delaunay2D:  !! No match for -3 face!\n";
				}
			}

			//empty badelems and voidpoly
			badelems.clear();
			voidpoly.clear();
		}
		// end iteration over points

		// Remove super triangle
		for(int ielem = 0; ielem < elems.size(); ielem++)
		{
			for(int i = 0; i < 3; i++)
				if( (nodes[elems[ielem].p[i]].x == nodes[0].x && nodes[elems[ielem].p[i]].y == nodes[0].y) || (nodes[elems[ielem].p[i]].x == nodes[1].x && nodes[elems[ielem].p[i]].y == nodes[1].y) || (nodes[elems[ielem].p[i]].x == nodes[2].x && nodes[elems[ielem].p[i]].y == nodes[2].y) )
				{
					elems.erase(elems.begin()+ielem);
					// re-adjust surr[] of each element
					for(int iel = 0; iel < elems.size(); iel++)
					{
						for(int j = 0; j < 3; j++)
						{
							if(elems[iel].surr[j] == ielem) elems[iel].surr[j] = -1;
							if(elems[iel].surr[j] > ielem) elems[iel].surr[j]--;
						}
					}
					//we could also readjust elem of each face here, but it's not needed

					ielem--;
					break;		// once the element has been deleted, we don't want to check other points
				}
		}
		// remove super nodes
		nodes.erase(nodes.begin(),nodes.begin()+3);
		for(int ielem = 0; ielem < elems.size(); ielem++)
		{
			for(int i = 0; i < 3; i++)
			{
				elems[ielem].p[i] = elems[ielem].p[i]-3;
			}
		}
		// remove super faces - not needed
		std::cout << "Delaunay2D: Triangulation done.\n";
	}

	void clear()					// reset the Delaunay2D object, except for input data
	{
		nodes.clear();
		elems.clear();
		faces.clear();
		badelems.clear();
		voidpoly.clear();
	}

	/// Write the Delaunay graph in Gmsh format
	void writeGmsh2(std::string mfile)
	{
		std::ofstream outf(mfile);

		outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
		outf << "$Nodes\n" << npoints << '\n';
		outf << std::setprecision(MESHDATA_DOUBLE_PRECISION);
		for(int ip = 0; ip < npoints; ip++)
		{
			outf << ip+1 << " " << nodes[ip].x << " " << nodes[ip].y << " " << 0 << '\n';
		}
		outf << "$Elements\n" << elems.size() /*+nface*/ << '\n';
		for(int iel = 0; iel < elems.size(); iel++)
		{
			outf << iel+1 << " 2 2 0 2";
			for(int i = 0; i < nnode; i++)
				outf << " " << elems[iel].p[i]+1;
			outf << '\n';
		}
		outf << "$EndElements\n";

		outf.close();
		//std::cout << "Delaunay2D: Number of faces finally = " << faces.size() << std::endl;
	}

	/// Locates a given point with coordinates xx and yy in the Delaunay graph and returns its area coordinates in the containing element
	Walkdata find_containing_triangle_and_area_coords(double xx, double yy, int startelement)
	{
		Walkdata dat;
		int ielem = startelement;
		double l1, l2, l3;
		Triangle super;
		while(1)
		{
			if(ielem < 0 || ielem >= elems.size()) { 
				std::cout << "Delaunay2D: find_containing_triangle..(): !!Reached an element index that is out of bounds!! Index is " << ielem << "\n";
				break;
			}
			//super = elems[ielem];
			super = elems.at(ielem);			// at() function checks index out-of-bounds, unlike overloaded [] operator (so it's also a bit slower).

			l3 = xx*(nodes[super.p[0]].y - nodes[super.p[1]].y) - yy*(nodes[super.p[0]].x - nodes[super.p[1]].x) + nodes[super.p[0]].x*nodes[super.p[1]].y - nodes[super.p[1]].x*nodes[super.p[0]].y;
			#if DEBUG==1
			if(fabs(l3) < tol) std::cout << "Delaunay2D:   Degenerate case (type 1) l3!!\n";
			#endif
			if(l3/super.D < 0)
			{
				ielem = super.surr[2];
				continue;
			}

			l1 = xx*(nodes[super.p[1]].y - nodes[super.p[2]].y) - yy*(nodes[super.p[1]].x - nodes[super.p[2]].x) + nodes[super.p[1]].x*nodes[super.p[2]].y - nodes[super.p[2]].x*nodes[super.p[1]].y;
			#if DEBUG==1
			if(fabs(l1) < tol) std::cout << "Delaunay2D:   Degenerate case (type 1) l1!!\n";
			#endif
			if(l1/super.D < 0)
			{
				ielem = super.surr[0];
				continue;
			}

			l2 = xx*(nodes[super.p[2]].y - nodes[super.p[0]].y) - yy*(nodes[super.p[2]].x - nodes[super.p[0]].x) + nodes[super.p[2]].x*nodes[super.p[0]].y - nodes[super.p[0]].x*nodes[super.p[2]].y;
			#if DEBUG==1
			if(fabs(l2) < tol) std::cout << "Delaunay2D:   Degenerate case (type 1) l2!!\n";
			#endif
			if(l2/super.D < 0)
			{
				ielem = super.surr[1];
				continue;
			}
			break;		// if all three area-ratios are positive, we've found our element
		}
		dat.elem = ielem; dat.areacoords[0] = l1/super.D; dat.areacoords[1] = l2/super.D; dat.areacoords[2] = l3/super.D;
		return dat;
	}

	void compute_jacobians()
	{
		jacobians.setup(elems.size(),1);
		for(int i = 0; i < elems.size(); i++)
		{
			jacobians(i,0) = nodes[elems[i].p[0]].x * (nodes[elems[i].p[1]].y - nodes[elems[i].p[2]].y) - nodes[elems[i].p[0]].y*(nodes[elems[i].p[1]].x-nodes[elems[i].p[2]].x) + nodes[elems[i].p[1]].x*nodes[elems[i].p[2]].y - nodes[elems[i].p[2]].x*nodes[elems[i].p[1]].y;
		}
	}

	bool detect_negative_jacobians()
	{
		bool flagj = false;
		amc_int numneg = 0;
		for(int i = 0; i < elems.size(); i++)
		{
			if(jacobians(i,0) < ZERO_TOL) {
				//out << i << " " << jacobians(i,0) << '\n';
				flagj = true;
				numneg++;
				std::cout << i << '\n';
			}
		}
		if(flagj == true) std::cout << "Delaunay2D: detect_negative_jacobians(): There exist " << numneg << " element(s) with negative jacobian!!\n";
		else std::cout << "Delaunay2d: detect_negative_jacobians(): DG is fine." << std::endl;
		return flagj;
	}
};

} // end namespace amc
