/** A class I wrote while writing code for Delaunay triangulation, that seems useless now.
*/

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef DEBUGW
#define DEBUGW 1
#endif

/** Class to store an n-dimensional point. */
class Point
{
	int dim;		///< Dimension of the space this point is part of.
	double* x;		///< vector of coordinates.
public:
	Point(int dimspace);
	Point(const Point& other);
	~Point();
	void setup(int dimspace);
	Point& operator=(const Point& other);
	double& operator()(int idim);
	double get(int idim) const;
};

Point::Point(int dimspace)
{
	dim = dimspace;
	x = new double[dim];
}

Point::Point(const Point& other)
{
	dim = other.dim;
	x = new double[dim];
	for(int i = 0; i < dim; i++
		x[i] = other.x[i];
}

Point::~Point()
{
	delete [] x;
}

void Point::setup(int dimspace)
{
	dim = dimspace;
	x = new double[dim];
}

Point& Point::operator=(const Point& other)
{
	dim = other.dim;
	x = new double[dim];
	for(int i = 0; i < dim; i++
		x[i] = other.x[i];
	return *this;
}

inline double& Point::operator()(int idim)
{
	#if DEBUGW==1
	if(idim >= dim) {
		std::cout << "! Point: (): Dimension " << idim << " does not exist!" << std::endl;
		return x[0];
	}
	#endif
	return x[idim];
}

inline double Point::get(int idim) const
{
	#if DEBUGW==1
	if(idim >= dim) {
		std::cout << "! Point: (): Dimension " << idim << " does not exist!" << std::endl;
		return x[0];
	}
	#endif
	return x[idim];
}
// --------------------- class Point over -------------------------------------
