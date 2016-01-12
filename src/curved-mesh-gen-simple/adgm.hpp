/* Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
Aditya Kashi
July 1, 2015
*/

#include <abowyerwatson.hpp>

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

	DGmove() {}

	DGmove(int num_dimn, Matrix<double>* incoords, Matrix<double>* bouncoords, Matrix<double>* boundary_motion)
	{
		// boundary_motion has as many rows as bouncoords and contains x and y displacement values for each boun point
		bmotion = boundary_motion;
		if(bmotion->rows() != bouncoords->rows() || bmotion->cols() != bouncoords->cols())
			cout << "! DGmove: Dimensions of boundary point coordinate array and boundary displacement array do not match!!" << endl;
		ndim = incoords->cols();
		ndgpoin = bouncoords->rows();
		ninpoin = incoords->rows();
		points.setup(ninpoin, ndim+1+ndim+1);	// for each interior point, store coords, containing DG element, and area coordinates w.r.t. that DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.
		dgpoints = *bouncoords;					// copy bouncoords

		for(int i = 0; i < incoords->rows(); i++)
			for(int j = 0; j < incoords->cols(); j++)
				points(i,j) = incoords->get(i,j);

		dg.setup(&dgpoints, ndgpoin);
	}

	void setup(int num_dimn, Matrix<double>* incoords, Matrix<double>* bouncoords, Matrix<double>* boundary_motion)
	{
		// boundary_motion has as many rows as bouncoords and contains x and y displacement values for each boun point
		bmotion = boundary_motion;
		if(bmotion->rows() != bouncoords->rows() || bmotion->cols() != bouncoords->cols())
			cout << "! DGmove: Dimensions of boundary point coordinate array and boundary displacement array do not match!!" << endl;
		ndim = incoords->cols();
		ndgpoin = bouncoords->rows();
		ninpoin = incoords->rows();
		points.setup(ninpoin, ndim+1+ndim+1);	// for each interior point, store coords, containing DG element, and area coordinates w.r.t. that DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.
		dgpoints = *bouncoords;					// copy bouncoords

		for(int i = 0; i < incoords->rows(); i++)
			for(int j = 0; j < incoords->cols(); j++)
				points(i,j) = incoords->get(i,j);

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
		/*cout << "DGmove: generateDG(): Checking\n";
		for(int iel = 0; iel < dg.elems.size(); iel++)
		{
			for(int j = 0; j < ndim+1; j++)
			{
				if(dgesuel(iel,j) >= ndgpoin)
					cout << "DGmove: generateDG(): esuel contains out-of-bounds entry at (" << iel << ", " << j << "), with entry " << dgesuel(iel,j) << endl;
			}
		}*/
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
			//cout << "DGmove: movemesh(): Point " << ipoin << endl;
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
			for(int j = 0; j < ndim; j++)
				dgpoints(i,j) += bmotion->get(i,j);
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
		cout << "DGmove: movedg(): Checking jacobians of the DG\n";
		dg.compute_jacobians();
		dg.detect_negative_jacobians();
	}

	Matrix<double> getInteriorPoints()
	{ return points; }

	Matrix<double> getBoundaryPoints()
	{ return dgpoints; }

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
