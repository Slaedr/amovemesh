/** @class MRotation2d
 * Determines boundary displacements corresponding to rotation of a mesh boundary by a given angle.
 */

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

namespace amc {

class MRotation2d
{
	double theta;
	double cx;
	double cy;
	amat::Matrix<int> flag;

	UMesh2d* m;

public:
	/** Sets the required variables.
	 * \param t angle of rotation in radians
	 * \param x x-coordinate of rotation center
	 * \param y y-coordinate of rotation center
	 * \param flag_to_rotate contains boundary markers of the faces to be rotated
	 */
	MRotation2d(UMesh2d* mesh, double t, double x, double y, amat::Matrix<int> flag_to_rotate)
	{
		std::cout << "MRotation2d: Initializing..\n";
		m = mesh;
		theta = t; cx = x; cy = y;
		flag = flag_to_rotate;
	}

	/// Returns x-displacement of a point with coords (x,y)
	double drx(double x, double y)
	{
		return x*cos(theta) - y*sin(theta) - x;
	}

	/// Returns y-displacement of a point with coords (x,y)
	double dry(double x, double y)
	{
		return x*sin(theta) + y*cos(theta) - y;
	}

	/// In this, I assume that the first boundary tag of a bface contains 'flag' if that bface is part of the surface to be rotated
	amat::Matrix<double> rhsvect_rotate()
	{
		std::cout << "MRotation2d: Calculating boundary displacements\n";
		amat::Matrix<double> rhs(m->gnpoin(), 2);
		rhs.zeros();
		int ip[2];

		for(int b = 0; b < m->gnface(); b++)
		{
			for(int im = 0; im < flag.msize(); im++)
				if(m->gbface(b,m->gnnofa()) == flag(im))
				{
					for(int i = 0; i < m->gnnofa(); i++)
					{
						ip[i] = m->gbface(b,i);
						rhs(ip[i],0) = drx(m->gcoords(ip[i],0)-cx, m->gcoords(ip[i],1)-cy);
						rhs(ip[i],1) = dry(m->gcoords(ip[i],0)-cx, m->gcoords(ip[i],1)-cy);
					}
					break;
				}
		}
		std::cout << "MRotation2d: boundary displacements calculated.\n";
		return rhs;
	}
	
	/// Return a vector of rotation angles for each point in the mesh. 
	/** Zero for all points except those in flag.
	 */
	amat::Matrix<double> rhsvect_angles()
	{
		std::cout << "MRotation2d: Calculating RHS vector" << std::endl;
		amat::Matrix<double> rhs(m->gnpoin(),1);
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
		std::cout << "MRotation2d: RHS vector calculated.\n";
		return rhs;
	}
};

}		// end namespace
