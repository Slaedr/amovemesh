/** @file adgm.hpp
 * @brief Mesh movement using Delaunay graph (DG) mapping technique of Liu, Qin and Xia.
 * @author Aditya Kashi
 * @date July 1, 2015
 */

#ifndef __ABOWYERWATSON_H
#include <abowyerwatson.hpp>
#endif

// for returning the Delaunay graph as a mesh
#ifndef __AMESH2DGENERAL_H
#include <amesh2d.hpp>
#endif

#define __ADGM_H 1

namespace amc {

/// Class to carry out mesh-movement using Delaunay graph mapping.
class DGmove
{
	int ndim;
	int npoin;
	int ninpoin;
	amat::Matrix<double> points;
	int ndgpoin;
	amat::Matrix<double> dgpoints;
	int ndgelem;
	amat::Matrix<int> dginpoel;
	amat::Matrix<int> dgesuel;
	amat::Matrix<double> dgintfac;
	amat::Matrix<double>* bmotion;
	amat::Matrix<int> bflag;

public:
	Delaunay2D dg;
	amat::Matrix<double> newcoords;

	DGmove() {}

	DGmove(int num_dimn, amat::Matrix<double>* incoords, amat::Matrix<double>* bouncoords, amat::Matrix<double>* boundary_motion)
	{
		// boundary_motion has as many rows as bouncoords and contains x and y displacement values for each boun point
		bmotion = boundary_motion;
		if(bmotion->rows() != bouncoords->rows() || bmotion->cols() != bouncoords->cols())
			std::cout << "! DGmove: Dimensions of boundary point coordinate array and boundary displacement array do not match!!" << std::endl;
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

	void setup(int num_dimn, amat::Matrix<double>* incoords, amat::Matrix<double>* bouncoords, amat::Matrix<double>* boundary_motion)
	{
		// boundary_motion has as many rows as bouncoords and contains x and y displacement values for each boun point
		bmotion = boundary_motion;
		if(bmotion->rows() != bouncoords->rows() || bmotion->cols() != bouncoords->cols())
			std::cout << "! DGmove: Dimensions of boundary point coordinate array and boundary displacement array do not match!!" << std::endl;
		ndim = incoords->cols();
		ndgpoin = bouncoords->rows();
		ninpoin = incoords->rows();
		points.setup(ninpoin, ndim+1+ndim+1);	// for each interior point, store coords, containing DG element, and area coordinates w.r.t. that DG element
		//columns 0 and 1 contain x- and y-coords, column 2 contains index of containing element, and columns 3,4,5 contain area coordinates.
		
		dgpoints = *bouncoords;					// deep copy bouncoords

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
		std::cout << "DGmove: generateDG(): No. of DG elements: " << ndgelem << std::endl;
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

		/*std::cout << "DGmove: generateDG(): Checking DG.\n";
		dg.compute_jacobians();
		dg.detect_negative_jacobians();*/

		// for debugging
		/*std::cout << "DGmove: generateDG(): Checking\n";
		for(int iel = 0; iel < dg.elems.size(); iel++)
		{
			for(int j = 0; j < ndim+1; j++)
			{
				if(dgesuel(iel,j) >= ndgpoin)
					std::cout << "DGmove: generateDG(): esuel contains out-of-bounds entry at (" << iel << ", " << j << "), with entry " << dgesuel(iel,j) << std::endl;
			}
		}*/
	}

	void movemesh()
	{
		// cycle over interior points
		Walkdata dat;
		std::cout << "DGmove: movemesh(): Calculating containing elements and area coordinates for each interior point\n";
		for(int ipoin = 0; ipoin < ninpoin; ipoin++)
		{
			//std::cout << ipoin << std::endl;
			// first find containing DG element by "walking-through" the DG
			//std::cout << "DGmove: movemesh(): Point " << ipoin << std::endl;
			dat = dg.find_containing_triangle_and_area_coords(points(ipoin,0), points(ipoin,1), dg.elems.size()/2);
			
			if(dat.elem < 0 || dat.elem >= dg.elems.size())
			{
				std::cout << "DGmove: movemesh(): Error in locating point " << ipoin << "!!" << std::endl;
				return;
			}

			// store DG element and area coords in points
			points(ipoin,2) = dat.elem;
			for(int i = 0; i < ndim+1; i++)
				points(ipoin,i+3) = dat.areacoords[i];
		}

		// update coordinates in dgpoints using bmotion
		std::cout << "DGmove: movemesh(): Moving the Delaunay graph\n";
		int  k = 0;
		for(int i = 0; i < bmotion->rows(); i++)
		{
			for(int j = 0; j < ndim; j++)
				dgpoints(i,j) += bmotion->get(i,j);
		}

		// calculate new positions of interior points by mapping them to deformed DG elements using area coordinates calculated before
		std::cout << "DGmove: movemesh(): Moving the interior points\n";
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
		std::cout << "DGmove: movedg(): Checking jacobians of the DG\n";
		dg.compute_jacobians();
		dg.detect_negative_jacobians();
	}

	/// Combining the 3 steps of movement - generate DG, move mesh and move DG into 1 function.
	void move()
	{
		generateDG();
		movemesh();
		movedg();
	}

	amat::Matrix<double> getInteriorPoints()
	{ return points; }

	amat::Matrix<double> getBoundaryPoints()
	{ return dgpoints; }

	/**amat::Matrix<double> getcoords()
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
	}*/
	
	/// Return the DG as a mesh
	UMesh2d getDelaunayGraph()
	{
		UMesh2d dgmesh;
		dgmesh.setcoords(&dgpoints);
		dgmesh.setinpoel(&dginpoel);
		dgmesh.setnface(0);
		return dgmesh;
	}
};

} // end namespace amc
