#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif
#ifndef __AMESH2D_H
#include <amesh2.hpp>
#endif
#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

#define __AMOVEMESH2D_H 1

using namespace acfd;
using namespace std;

class TorsionSpring2D
{
	UTriMesh* m;
	
public:

	TorsionSpring2D(UTriMesh* mesh)
	{
		m = mesh;
	}
};
		

class MRotation2d
{
	double theta;
	double cx;
	double cy;
	
	UTriMesh* m;
	
public:
	MRotation2d(UTriMesh* mesh, double t, double x, double y)
	{
		cout << "MRotation2d: Initializing..\n";
		m = mesh;
		theta = t; cx = x; cy = y;
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
		Matrix<double> rhs(m->gnpoin(), 2, ROWMAJOR);		//NOTE: Column major!!
		rhs.zeros();
		int ip[2];
		
		for(int b = 0; b < m->gnface(); b++)
		{
			for(int i = 0; i < 2; i++)
			{
				ip[i] = m->gbface(b,i);
				rhs(ip[i],0) = (1 - m->gbface(b,3)) * drx(m->gcoords(ip[i],0), m->gcoords(ip[i], 1));
				rhs(ip[i],1) = (1 - m->gbface(b,3)) * dry(m->gcoords(ip[i],0), m->gcoords(ip[i], 1));
			}
		}
		cout << "MRotation2d: RHS vector calculated.\n";
		return rhs;
	}
};
