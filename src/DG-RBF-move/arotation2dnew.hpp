#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif
#ifndef __AMESH2DGENERAL_H
#include <amesh2d.hpp>
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

	UMesh2d* m;

public:
	MRotation2d() {}

	MRotation2d(UMesh2d* mesh, double angle, double centre_x, double centre_y, Matrix<int> nums_to_rotate)
	{
		cout << "MRotation2d: Initializing..\n";
		m = mesh;
		theta = angle; cx = centre_x; cy = centre_y;
		nrot = nums_to_rotate;
	}

	void setup(UMesh2d* mesh, double angle, double centre_x, double centre_y, Matrix<int> nums_to_rotate)
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
				if(m->gbface(m->gbpointsb(b,2),m->gnnofa())==nrot(i))			// check one of the faces surrounding boundary point b
				{
					rhs(b,0) = drx(m->gcoords(b,0)-cx, m->gcoords(b, 1)-cy);
					rhs(b,1) = dry(m->gcoords(b,0)-cx, m->gcoords(b, 1)-cy);
				}
			}
		}
		cout << "MRotation2d: RHS vector calculated.\n";
		return rhs;
	}

	/** Return a vector of rotation angles for each point in the mesh. 
		Zero for all points except those in nrot.
	*/
	Matrix<double> rhsvect_angles()
	{
		cout << "MRotation2d: Calculating RHS vector" << endl;
		Matrix<double> rhs(m->gnpoin(),1);
		rhs.zeros();
		for(int i = 0; i < m->gnface(); i++)
		{
			for(int im = 0; im < nrot.msize(); im++)
				if(m->gbface(i,m->gnnofa()) == nrot(im))
				{
					for(int j = 0; j < m->gnnofa(); j++)
						rhs(m->gbface(i,j)) = theta;
				}
		}
		return rhs;
	}
};
