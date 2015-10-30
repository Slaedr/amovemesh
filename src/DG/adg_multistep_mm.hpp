/* Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
Aditya Kashi
July 1, 2015
*/

#include "abowyerwatson.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

namespace acfd {

class DGmove
{
	int ndim;
	int npoin;
	int ninpoin;
	Matrix<double> points;
	int ndgpoin;
	Matrix<double> dgpoints;
	Matrix<double> testdgpoints;
	int ndgelem;
	Matrix<int> dginpoel;
	Matrix<int> dgesuel;
	Matrix<double> dgintfac;
	Matrix<double>* bmotion;
	Matrix<int> bflag;
	Delaunay2D testdg;

	Walkdata dat;
	bool check;
	Matrix<double>* bmot;

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
		testdgpoints.setup(ndgpoin,ndim);
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
		//testdg = dg;

		// copy bmotion into bmot
		bmot = new Matrix<double>;
		bmot->setup(bmotion->rows(),bmotion->cols());
		for(int i = 0; i < bmotion->rows(); i++)
			for(int j = 0; j < bmotion->cols(); j++)
				(*bmot)(i,j) = bmotion->get(i,j);
	}

	~DGmove()
	{
		delete bmot;
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
	}

	void movemesh()
	{
		int num_steps = 1;
		int step = 1;			// just for info of total number of steps

		generateDG();
		testdg = dg;

		cout << "DGmove: movemesh():  Calculating containing elements and area coordinates for each interior point\n";
		// cycle over interior points
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			// first find containing DG element by "walking-through" the DG
			dat = dg.find_containing_triangle_and_area_coords(points(ipoin,0), points(ipoin,1), ninpoin/2);

			// store DG element and area coords in points
			points(ipoin,2) = dat.elem;
			for(int i = 0; i < ndim+1; i++)
				points(ipoin,i+3) = dat.areacoords[i];
		}


		while(num_steps > 0)
		{
			cout << "DGmove: movemesh(): Step " << step << endl;

			// update coordinates in dgpoints using bmotion
			cout << "DGmove: movemesh():  Moving the Delaunay graph\n";
			int  k = 0;
			//testdgpoints.zeros();
			for(int i = 0; i < bmot->rows(); i++)
			{
				if(bflag(i) == 1)
				{
					for(int j = 0; j < ndim; j++)
						testdgpoints(k,j) = dgpoints(k,j) + bmot->get(i,j);
					k++;
				}
			}

			// check for validity of Delaunay graph
			movegraph(testdg, testdgpoints);
			testdg.compute_jacobians();
			check = testdg.detect_negative_jacobians();
			cout << "DGmove: movemesh():  Check is " << check << endl;
			if(check == true)
			{
				testdgpoints = dgpoints;					// reset testdgpoints to original values
				testdg = dg;								// reset the DG as well

				//halve the boundary motion
				for(int i = 0; i < bmot->rows(); i++)
					for(int j = 0; j < bmot->cols(); j++)
						(*bmot)(i,j) = bmot->get(i,j)/2;

				num_steps++;
				continue;
			}
			else
			{
				dgpoints = testdgpoints;
				movegraph(dg, dgpoints);
				num_steps--;
			}

			// calculate new positions of interior points by mapping them to deformed DG elements using area coordinates calculated before
			cout << "DGmove: movemesh():  Moving the interior points\n";
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
			if(num_steps > 0)
			{
				cout << "DGmove: movemesh():  Re-triangulating\n";
				dg.clear();
				dg.setup(&dgpoints,ndgpoin);
				generateDG();
				// re-calculate area-coords
				for(int ipoin = 0; ipoin < ninpoin; ipoin++)
				{
					// first find containing DG element by "walking-through" the DG
					dat = dg.find_containing_triangle_and_area_coords(points(ipoin,0), points(ipoin,1), ninpoin/2);

					// store DG element and area coords in points
					points(ipoin,2) = dat.elem;
					for(int i = 0; i < ndim+1; i++)
						points(ipoin,i+3) = dat.areacoords[i];
				}
			}
			step++;
		}
	}

	void movegraph(Delaunay2D& tdg, Matrix<double>& dgp)
	{
		//moves the DG according to the boundary motion - for visualizing deformed DG
		for(int i = 0; i < dgp.rows(); i++)
		{
			tdg.nodes[i].x = dgp(i,0);
			tdg.nodes[i].y = dgp(i,1);
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
