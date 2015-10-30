#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif
#ifndef __AMESH2DGENRERAL_H
#include <amesh2d.hpp>
#endif
#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

#define __AROTATION2D_H 1

using namespace acfd;
using namespace std;


class MRotation2d
{
	double theta;
	double cx;
	double cy;
	Matrix<int> flag;

	UMesh2d* m;

public:
	MRotation2d(UMesh2d* mesh, double t, double x, double y, Matrix<int> flag_to_rotate)
	{
		cout << "MRotation2d: Initializing..\n";
		m = mesh;
		theta = t; cx = x; cy = y;
		flag = flag_to_rotate;
	}

	double drx(double x, double y)
	{
		return x*cos(theta) - y*sin(theta) - x;
	}

	double dry(double x, double y)
	{
		return x*sin(theta) + y*cos(theta) - y;
	}

	//NOTE: In the following I assume that bface(i,3) contains 'flag' if that bface is part of the surface to be rotated, and
	//		it is 1 if bface 'i' is part of the fixed surface.
	Matrix<double> rhsvect_rotate()
	{
		cout << "MRotation2d: Calculating RHS vector\n";
		Matrix<double> rhs(m->gnpoin(), 2, ROWMAJOR);		//NOTE: Column major!!
		rhs.zeros();
		int ip[2];

		for(int b = 0; b < m->gnface(); b++)
		{
			for(int i = 0; i < m->gnnofa(); i++)
			{
				for(int im = 0; im < flag.msize(); im++)
					if(m->gbface(b,m->gnnofa()) == flag(im))
					{
						ip[i] = m->gbface(b,i);
						rhs(ip[i],0) = drx(m->gcoords(ip[i],0)-cx, m->gcoords(ip[i], 1)-cy);
						rhs(ip[i],1) = dry(m->gcoords(ip[i],0)-cx, m->gcoords(ip[i], 1)-cy);
					}
			}
		}
		cout << "MRotation2d: RHS vector calculated.\n";
		return rhs;
	}
	
	/** Return a vector of rotation angles for each point in the mesh. 
		Zero for all points except those in flag.
	*/
	Matrix<double> rhsvect_angles()
	{
		cout << "MRotation2d: Calculating RHS vector" << endl;
		Matrix<double> rhs(m->gnpoin(),1);
		rhs.zeros();
		for(int i = 0; i < m->gnface(); i++)
		{
			for(int im = 0; im < flag.msize(); im++)
				if(m->gbface(i,m->gnnofa()) == flag(im))
				{
					for(int j = 0; j < m->gnnofa(); j++)
						rhs(m->gbface(i,j)) = theta;
				}
		}
		cout << "MRotation2d: RHS vector calculated.\n";
		return rhs;
	}
};
