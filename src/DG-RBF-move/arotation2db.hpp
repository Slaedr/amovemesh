#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif
#ifndef __AMESH2D_H
#include <amesh2b.hpp>
#endif

#define __AROTATION2D_H 1

using namespace amat;
using namespace acfd;
using namespace std;


class MRotation2d
{
	double theta;
	double cx;
	double cy;
	Matrix<int> nrot;				// a vector containing the integer flags of the boundary-sets (bodies) to be rotated

	UTriMesh* m;

public:
	MRotation2d() {}

	MRotation2d(UTriMesh* mesh, double angle, double centre_x, double centre_y, Matrix<int> nums_to_rotate)
	{
		cout << "MRotation2d: Initializing..\n";
		m = mesh;
		theta = angle; cx = centre_x; cy = centre_y;
		nrot = nums_to_rotate;
	}

	void setup(UTriMesh* mesh, double angle, double centre_x, double centre_y, Matrix<int> nums_to_rotate)
	{
		cout << "MRotation2d: Initializing..\n";
		m = mesh;
		theta = angle; cx = centre_x; cy = centre_y;
		nrot = nums_to_rotate;
	}

	double drx(double x, double y)
	{
		return x*cos(theta) - y*sin(theta) - x;
	}

	double dry(double x, double y)
	{
		return x*sin(theta) + y*cos(theta) - y;
	}

	//NOTE: In the following I assume that bface(i,3) contains 0 if that bface is part of the surface to be rotated, and
	//		it is 1 if bface 'i' is part of the fixed surface.
	Matrix<double> rhsvect_rotate()
	{
		cout << "MRotation2d: Calculating RHS vector\n";
		Matrix<double> rhs(m->gnpoin(), m->gndim());
		rhs.zeros();

		for(int b = 0; b < m->gnbpoin(); b++)
		{
			for(int i = 0; i < nrot.msize(); i++)
			{
				if(m->gbpoints(b,4)==nrot(i))
				{
					rhs(b,0) = drx(m->gcoords(b,0)-cx, m->gcoords(b, 1)-cy);
					rhs(b,1) = dry(m->gcoords(b,0)-cx, m->gcoords(b, 1)-cy);
				}
			}
		}
		cout << "MRotation2d: RHS vector calculated.\n";
		return rhs;
	}
};
