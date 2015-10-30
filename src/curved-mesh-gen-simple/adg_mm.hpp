/* Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
Aditya Kashi
July 1, 2015
*/

#include "abowyerwatson.hpp"

using namespace std;
using namespace amat;

namespace acfd {

class DGmove
{
	int ndim;
	int npoin;
	int ninpoin;
	Matrix<double> points;
	int ndgpoin;
	Matrix<double> dgpoints;
	int ndgelem;
	Matrix<int> dginpoel;
	Matrix<int> dgesuel;
	Matrix<double> dgintfac;
	Matrix<double>* bmotion;
	Matrix<int> bflag;

public:
	Delaunay2D dg;
	Matrix<double> newcoords;

	DGmove(int num_dimn, Matrix<double>* coords, Matrix<int> bflags, Matrix<double>* boundary_motion)
	{
		// note that bflags contains as many rows as points in the mesh and contains 1 if the point is a boundary point and 0 otherwise
		// boundary_motion has as many rows as points in the original mesh and contains x and y displacement values for each point (zero for interior points)
		bmotion = boundary_motion;
		if(bmotion->rows() != bflags.rows()) cout << "DGMove: !! Error: No. of rows in bflags and bmotion are not equal!\n";
		ndim = coords->cols();
		ndgpoin = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1) ndgpoin++;		// ndgpoin is the number of points in the DG, which equals the number of boundary nodes in the original mesh
		ninpoin = coords->rows() - ndgpoin;
		points.setup(ninpoin, ndim+1+ndim+1);	// for each interior point, store coords, containing DG element, and area coordinates w.r.t. that DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.

		int k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 0)
			{
				for(int j = 0; j < ndim; j++)
					points(k,j) = coords->get(i,j);
				k++;
			}

		// now get dgpoints using bfac
		dgpoints.setup(ndgpoin,ndim);
		k = 0;
		for(int i = 0; i < bflags.rows(); i++)
			if(bflags(i) == 1)
			{
				for(int j = 0; j < ndim; j++)
					dgpoints(k,j) = coords->get(i,j);
				k++;
			}

		bflag = bflags;
		newcoords.setup(bflags.rows(),ndim);
		dg.setup(&dgpoints, ndgpoin);
	}

	void generateDG()
	{
		// generate the Delaunay graph
		dg.bowyer_watson();

		ndgelem = dg.elems.size();
		cout << "DGmove: generateDG(): No. of DG elements: " << ndgelem << endl;
		dginpoel.setup(ndgelem,ndim+1);
		dgesuel.setup(ndgelem,ndim+1);

		// populate DG data arrays
		for(int iel = 0; iel < dg.elems.size(); iel++)
		{
			for(int j = 0; j < ndim+1; j++)
			{
				dginpoel(iel,j) = dg.elems[iel].p[j];
				dgesuel(iel,j) = dg.elems[iel].surr[j];
			}
		}

		// for debugging
		cout << "DGmove: generateDG(): Checking\n";
		for(int iel = 0; iel < dg.elems.size(); iel++)
		{
			for(int j = 0; j < ndim+1; j++)
			{
				if(dgesuel(iel,j) == 3429)
					cout << "DGmove: generateDG(): esuel contains out-of-bounds entry at (" << iel << ", " << j << "), with entry " << dgesuel(iel,j) << endl;
			}
		}
	}

	void movemesh()
	{
		// cycle over interior points
		Walkdata dat;
		cout << "DGmove: movemesh(): Calculating containing elements and area coordinates for each interior point\n";
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			//cout << ipoin << endl;
			// first find containing DG element by "walking-through" the DG
			dat = dg.find_containing_triangle_and_area_coords(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);

			// store DG element and area coords in points
			points(ipoin,2) = dat.elem;
			for(int i = 0; i < ndim+1; i++)
				points(ipoin,i+3) = dat.areacoords[i];
		}

		// update coordinates in dgpoints using bmotion
		cout << "DGmove: movemesh(): Moving the Delaunay graph\n";
		int  k = 0;
		for(int i = 0; i < bmotion->rows(); i++)
		{
			if(bflag(i) == 1)
			{
				for(int j = 0; j < ndim; j++)
					dgpoints(k,j) += bmotion->get(i,j);
				k++;
			}
		}

		// calculate new positions of interior points by mapping them to deformed DG elements using area coordinates calculated before
		cout << "DGmove: movemesh(): Moving the interior points\n";
		int elem;
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			elem = points(ipoin,2);
			for(int dim = 0; dim < ndim; dim++)
			{
				points(ipoin,dim) = 0.0;
				for(int i = 0; i < ndim+1; i++)
					points(ipoin,dim) += points(ipoin,i+3)*dgpoints(dginpoel(elem,i),dim);
			}
		}
	}

	void movedg()
	{
		//moves the DG according to the boundary motion - for visualizing deformed DG
		for(int i = 0; i < dgpoints.rows(); i++)
		{
			dg.nodes[i].x = dgpoints(i,0);
			dg.nodes[i].y = dgpoints(i,1);
		}
	}

	Matrix<double> getcoords()
	{
		// create a coords matrix with same point numbering as initial matrix and return it
		int a = 0, b = 0, k = 0;
		for(int i = 0; i < bflag.rows(); i++)
		{
			if(bflag(i) == 0)
			{
				for(int dim = 0; dim < ndim; dim++)
					newcoords(k,dim) = points(a,dim);
				k++;
				a++;
			}
			else
			{
				for(int dim = 0; dim < ndim; dim++)
					newcoords(k,dim) = dgpoints.get(b,dim);
				k++;
				b++;
			}
		}
		return newcoords;
	}
};

} // end namespace acfd
