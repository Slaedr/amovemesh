/** \file ageometry.cpp
 *
 * Implementation of reconstruction od a C^1 (or C^2) piecewise polynomial (cubic) boundary from a piecewise linear boundary given by a linear mesh.
 * \author Aditya Kashi
 * \date September 4, 2015
 */

#include <ageometry.hpp>

namespace amc {

// --- CSpline implementation ---

void CSpline::setup(UMesh2d* mesh, amat::Matrix<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter)
{
	face_list_available = false;
	
	m = mesh;
	dim = m->gndim();
	rfl = recmarkers;
	isClosed = closed;
	issequenced = sequenced;
	tol = _tol;
	maxiter = _maxiter;
	ndf = 4;

	nseg = 0;
	toRec.setup(m->gnface(),1);
	toRec.zeros();

	// calculate number of spline segments, which equals the number of faces to be reconstructed, and populate toRec
	for(int iface = 0; iface < m->gnface(); iface++)
	{
		for(int im = 0; im < rfl.msize(); im++)
		{
			if(m->gbface(iface,m->gnnofa()) == rfl(im)) { 
				toRec(iface) = 1;
				nseg++;
			}
		}
	}
	
	/*nseg = 0;
	for(int iface = 0; iface < m->gnface(); iface++)
	{
		for(int ibf = 0; ibf < rfl.msize(); ibf++)
			if(m->gbface(iface,m->gnnofa()) == rfl(ibf)) nseg++;
	}*/

	if(isClosed) nspoin = nseg;
		else nspoin = nseg+1;
	
	if(isClosed)
		std::cout << "CSpline: setup(): Closed curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << std::endl;
	else
		std::cout << "CSpline: setup(): Open curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << std::endl;

	// allocate stuff
	
	scf = new amat::Matrix<double>[dim];
	D = new amat::Matrix<double>[dim];
	srhs = new amat::Matrix<double>[dim];
	for(int idim = 0; idim < dim; idim++) {
		scf[idim].setup(nseg,ndf);
		D[idim].setup(nseg+1,1);
		srhs[idim].setup(nseg+1,1);
	}
	
	slhs.setup(nseg+1,nseg+1);
	seq_bface.setup(m->gnface(),1);
	segface.setup(nseg,1);
	seq_spoin.setup(nseg+1, 1);
	for(int i = 0; i < m->gnface(); i++)
		seq_bface(i) = -1;
	for(int i = 0; i < nspoin; i++)
		seq_spoin(i) = -1;
}

void CSpline::setup(UMesh2d* mesh, std::vector<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter)
{
	face_list_available = true;
	
	m = mesh;
	dim = m->gndim();
	facelist = recmarkers;
	isClosed = closed;
	issequenced = sequenced;
	tol = _tol;
	maxiter = _maxiter;
	ndf = 4;

	nseg = 0;
	toRec.setup(m->gnface(),1);
	toRec.zeros();

	// calculate number of spline segments, which equals the number of faces to be reconstructed, and populate toRec
	nseg = facelist.size();
	for(int iseg = 0; iseg < nseg; iseg++)
		toRec(facelist[iseg]) = 1;
	

	if(isClosed) nspoin = nseg;
		else nspoin = nseg+1;
	
	if(isClosed)
		std::cout << "CSpline: setup(): Closed curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << std::endl;
	else
		std::cout << "CSpline: setup(): Open curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << std::endl;

	// allocate stuff
	
	scf = new amat::Matrix<double>[dim];
	D = new amat::Matrix<double>[dim];
	srhs = new amat::Matrix<double>[dim];
	for(int idim = 0; idim < dim; idim++) {
		scf[idim].setup(nseg,ndf);
		D[idim].setup(nseg+1,1);
		srhs[idim].setup(nseg+1,1);
	}
	
	slhs.setup(nseg+1,nseg+1);
	seq_bface.setup(m->gnface(),1);
	segface.setup(nseg,1);
	seq_spoin.setup(nseg+1, 1);
	for(int i = 0; i < m->gnface(); i++)
		seq_bface(i) = -1;
	for(int i = 0; i < nseg+1; i++)
		seq_spoin(i) = -1;
}

CSpline::~CSpline()
{
	delete [] scf;
	delete [] D;
	delete [] srhs;
}

void CSpline::sequence()
{
	if(!issequenced)
	{
		int curpoin = 0, nextpoin = 0, startpoin, startface, curface, nextface;

		// first find a face that belongs to the boundary to be reconstructed
		for(int iface = 0; iface < m->gnface(); iface++)
			if(toRec(iface) == 1) startface = iface;
		
		int nf = 0, np = 0;			// number of faces and number of points respectively

		if(isClosed)
		{
			nextface = startface;
			while(true)
			{
				curface = nextface;
				seq_bface(curface) = nf;
				segface(nf) = curface;
				seq_spoin(np) = m->gbface(curface,0);
				nf++; np++;

				nextpoin = m->gbface(curface,1);

				// search for the other face which contains nextpoin
				nextface = -1;
				for(int iface = 0; iface < m->gnface(); iface++)
					if(toRec(iface) == 1 && m->gbface(iface,0) == nextpoin) nextface = iface;
				if(nextface == -1) std::cout << "UMesh2d: sequence(): ! Point not found!" << std::endl;
				if(nextface == startface) break;
			}
			seq_spoin(np) = seq_spoin(0);		// close the point list
		}
		else		// curve is open
		{
			nextface = startface;

			//first find the last element
			std::cout << "UMesh2d: sequence(): Finding first face of open boundary." << std::endl;
			while(true)
			{
				curface = nextface;
				nextpoin = m->gbface(curface,0);
				
				// search for the other face which contains nextpoin
				nextface = -1;
				for(int iface = 0; iface < m->gnface(); iface++)
					if(toRec(iface) == 1 && m->gbface(iface,1) == nextpoin) nextface = iface;
				
				if(nextface == -1) {
					std::cout << "UMesh2d: sequence(): Reached start of boundary to be reconstructed." << std::endl;
					break;
				}
			}

			nextface = curface;

			// now that the first face has been found, sequence the boundary
			while(true)
			{
				curface = nextface;
				seq_bface(curface) = nf;
				segface(nf) = curface;
				seq_spoin(np) = m->gbface(curface,0);
				nf++; np++;

				nextpoin = m->gbface(curface,1);

				// search for the other face which contains nextpoin
				nextface = -1;
				for(int iface = 0; iface < m->gnface(); iface++)
					if(toRec(iface) == 1 && m->gbface(iface,0) == nextpoin) nextface = iface;
				
				if(nextface == -1) {
					std::cout << "UMesh2d: sequence(): " << nf << " boundary faces and " << np << " boundary points sequenced."  << std::endl;
					if(nf != nseg) std::cout << "UMesh2d: sequence(): ! All faces are not accounted for!" << std::endl;
					break;
				}
			}

			// finally, add last point
			seq_spoin(np) = nextpoin;
		}
	}

	else
	{
		if(!face_list_available) {
			std::cout << "CSpline: sequence(): Boundary part is sequenced, but facelist is not available!!" << std::endl;
			return;
		}
		
		// facelist is already sequenced; just use it.
		for(int iseg = 0; iseg < nseg; iseg++)
		{
			seq_bface(facelist[iseg]) = iseg;
			segface(iseg) = facelist[iseg];
			seq_spoin(iseg) = m->gbface(facelist[iseg],0);
		}
		seq_spoin(nseg) = m->gbface(facelist[nseg-1], 1);
	}
}

void CSpline::compute()
{
	if(isClosed)
	{
		// Assemble RHS std::vector for each dimension
		for(int idim = 0; idim < dim; idim++)
		{
			srhs[idim](0) = 3.0*(m->gcoords(seq_spoin(1),idim) - m->gcoords(seq_spoin(nseg-1),idim));
			for(int i = 1; i < nseg; i++)
			{
				srhs[idim](i) = 3.0*(m->gcoords(seq_spoin(i+1),idim) - m->gcoords(seq_spoin(i-1),idim));
			}
			srhs[idim](nseg) = 3.0*(m->gcoords(seq_spoin(1),idim) - m->gcoords(seq_spoin(nseg-1),idim));
		}

		// Assemble LHS matrix
		slhs.set(0,0, 4.0);
		slhs.set(0,1, 1.0);
		slhs.set(0,nseg-1, 1.0);

		for(int i = 1; i < nseg; i++)
		{
			slhs.set(i,i-1, 1.0);
			slhs.set(i,i, 4.0);
			slhs.set(i,i+1, 1.0);
		}

		slhs.set(nseg,1, 1.0);
		slhs.set(nseg,nseg-1, 1.0);
		slhs.set(nseg,nseg, 4.0);
	}

	else
	{
		// Assemble RHS std::vector for each dimension
		for(int idim = 0; idim < dim; idim++)
		{
			srhs[idim](0) = 3.0*(m->gcoords(seq_spoin(1),idim) - m->gcoords(seq_spoin(0),idim));
			for(int i = 1; i < nseg; i++)
			{
				srhs[idim](i) = 3.0*(m->gcoords(seq_spoin(i+1),idim) - m->gcoords(seq_spoin(i-1),idim));
			}
			srhs[idim](nseg) = 3.0*(m->gcoords(seq_spoin(nseg),idim) - m->gcoords(seq_spoin(nseg-1),idim));
		}

		// Assemble LHS matrix
		slhs.set(0,0, 2.0);
		slhs.set(0,1, 1.0);

		for(int i = 1; i < nseg; i++)
		{
			slhs.set(i,i-1, 1.0);
			slhs.set(i,i, 4.0);
			slhs.set(i,i+1, 1.0);
		}

		slhs.set(nseg,nseg-1, 1.0);
		slhs.set(nseg,nseg, 2.0);
	}

	// solve
	std::cout << "CSpline: compute(): Solving linear systems" << std::endl;
	amat::Matrix<double> d0(nseg+1,1);
	d0.zeros();
	for(int idim = 0; idim < dim; idim++)
		D[idim] = sparseCG_d(&slhs, srhs[idim], d0, tol, maxiter);
	
	// get coeffs
	//std::cout << "CSpline: compute(): Getting spline coefficients" << std::endl;
	double yi, yip;
	for(int idim = 0; idim < dim; idim++)
	{
		for(int i = 0; i < nseg; i++)
		{
			// first get idim-coordinates of the ith and (i+1)th spline points
			yi = m->gcoords(seq_spoin.get(i),idim);
			yip = m->gcoords(seq_spoin.get(i+1),idim);

			scf[idim](i,0) = yi;
			scf[idim](i,1) = D[idim](i);
			scf[idim](i,2) = 3*(yip - yi) - 2*D[idim](i) - D[idim].get(i+1);
			scf[idim](i,3) = 2*(yi - yip) + D[idim](i) + D[idim].get(i+1);
		}
	}
	std::cout << "CSpline: compute(): Done." << std::endl;
}

double CSpline::getspline(int iface, int idim, double t)
{
	return scf[idim].get(seq_bface.get(iface),0) + scf[idim].get(seq_bface.get(iface),1)*t + scf[idim].get(seq_bface.get(iface),2)*t*t + scf[idim].get(seq_bface.get(iface),3)*t*t*t;
}

// --------------------- End of class CSpline --------------------------------------------------------------------------------------------//


void BoundaryReconstruction2d::setup(UMesh2d* mesh, int num_parts, std::vector<std::vector<int>> boundary_markers, double angle_threshold)
{
	m = mesh;
	nparts = num_parts;
	marks = boundary_markers;
	cangle = angle_threshold;

	nparts = marks.size();
	ncorners.resize(nparts,0);
	isClosed.resize(nparts);
	startface.resize(nparts);
	toRec.setup(nparts,m->gnface());

	corners.resize(nparts);
	facepart.setup(m->gnface(),2);
}

BoundaryReconstruction2d::~BoundaryReconstruction2d()
{
	delete [] sparts;
}

void BoundaryReconstruction2d::preprocess()
{
	// populate toRec
	toRec.zeros();
	for(int ipart = 0; ipart < nparts; ipart++)
	{
		for(int iface = 0; iface < m->gnface(); iface++)
		{
			for(int i = 0; i < marks[ipart].size(); i++)
				if(m->gbface(iface,m->gnnofa()) == marks[ipart][i])
					toRec(ipart,iface) = 1;
		}
	}
	std::cout << "BoundaryReconstruction2d: preprocess(): toRec populated." << std::endl;
	for(int ipart = 0; ipart < nparts; ipart++)
	{
		int startfac, nextface, curface, nextpoin;
		bool closed;

		//find some face belonging to this part
		for(int iface = 0; iface < m->gnface(); iface++)
			if(toRec(ipart,iface)) {
				startfac = iface;
				break;
			}

		//next, determine if this part is open or closed
		//std::cout << "BoundaryReconstruction2d: preprocess(): Determining whether part " << ipart << " is open or closed" << std::endl;
		nextface = startfac;
		while(true)
		{
			curface = nextface;
			nextpoin = m->gbface(curface,0);

			// search for the other face which contains nextpoin
			nextface = m->gbpointsb(curface,1);

			if(nextface == startfac) {
				closed = true;
				break;
			}
			
			bool onpart = false;	
			
			if(toRec(ipart,nextface))
				onpart = true;

			if(!onpart) { 
				closed = false;
				break;
			}
		}

		if(closed) {
			isClosed[ipart] = true;
			startface[ipart] = startfac;
		}
		else {
			isClosed[ipart] = false;
			startface[ipart] = curface;
		}
	}
	std::cout << "BoundaryReconstruction2d: preprocess(): Done." << std::endl;
}

/*	Detects corner points by computing dot-products of normals of the two faces associated with a point and comparing it with cos(theta) for each point,
	for a given angle theta.
	Computes the corners data structure.
	Each element of corners is a list of corner points structures for each corner point, so
		corner[ipart][0] refers to the first corner point of part number ipart.
	Each corner point structure contains the global point number of the corner and the two faces associated with that point. So,
		corner[ipart][0][0] refers the global point number of the first corner of part ipart,
		corner[ipart][0][1] contains the left face associated with the corner, and
		corner[ipart][0][2] contians the right face associated with the corner.
*/
void BoundaryReconstruction2d::detect_corners()
{
	int p1, p2, p3;
	double n1x, n1y, n2x, n2y, mag1, mag2, dp;

	std::cout << "BoundaryReconstruction2d: detect_corners(): Searching each part for corners" << std::endl;
	
	// Scan each part. At each boundary point, calculate dot product of the two normals of the two faces containing that point.
	for(int ipart = 0; ipart < nparts; ipart++)
	{
		for(int ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
		{
			// check left and right faces
			if(toRec(ipart,m->gbpointsb(ipoin,1)) && toRec(ipart,m->gbpointsb(ipoin,2)))
			{
				p1 = m->gbface(m->gbpointsb(ipoin,1),0);
				p2 = m->gbpointsb(ipoin,0);
				p3 = m->gbface(m->gbpointsb(ipoin,2),1);
				
				// check dot-product of normals
				n1x = m->gcoords(p2,1) - m->gcoords(p1,1);
				n1y = m->gcoords(p1,0) - m->gcoords(p2,0);
				n2x = m->gcoords(p3,1) - m->gcoords(p2,1);
				n2y = m->gcoords(p2,0) - m->gcoords(p3,0);

				mag1 = sqrt(n1x*n1x + n1y*n1y);
				mag2 = sqrt(n2x*n2x + n2y*n2y);
				n1x /= mag1;
				n1y /= mag1;
				n2x /= mag2;
				n2y /= mag2;

				dp = n1x*n2x + n1y*n2y;

				if(dp < cos(cangle))
				{
					std::vector<int> c(4);
					c[0] = p2;						// global point number
					c[1] = ipoin;					// boundary point number
					c[2] = m->gbpointsb(ipoin,1);
					c[3] = m->gbpointsb(ipoin,2);
					corners[ipart].push_back(c);
					ncorners[ipart]++;
				}
			}
		}
		std::cout << "BoundaryReconstruction2d: detect_corners(): Part " << ipart << " contains " << ncorners[ipart] << " corner(s)" << std::endl;
		std::cout << "  at points ";
		for(int i = 0; i < ncorners[ipart]; i++)
			std::cout << corners[ipart][i][0] << ", ";
		std::cout << std::endl;
	}
}

void BoundaryReconstruction2d::read_corners(std::string cname)
{
	// read corner data from file
	std::ifstream fin(cname);
	
	std::cout << "BoundaryReconstruction2d: read_corners(): Reading file" << std::endl;
	int dum;
	std::vector<int> c(4,0);								// std::vector of 4 ints with value 0

	for(int ipart = 0; ipart < nparts; ipart++)
	{
		fin >> dum;									// read part number
		fin >> ncorners[ipart];						// read number of corners
		for(int i = 0; i < ncorners[ipart]; i++)
		{
			fin >> c[0];							// get global point number of this corner point
			corners[ipart].push_back(c);
		}
	}

	fin.close();

	std::cout << "BoundaryReconstruction2d: read_corners(): iterating over boundary points" << std::endl;
	
	// also store bpointsb number by searching in bpointsb for each corner point
	for(int ipart = 0; ipart < nparts; ipart++)
	{
		for(int i = 0; i < ncorners[ipart]; i++)
		{
			for(int ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
				if(m->gbpointsb(ipoin,0) == corners[ipart][i][0])
					corners[ipart][i][1] = ipoin;
		}
	}
}

void BoundaryReconstruction2d::split_parts()
{
	nnparts = 0;		// update to get total number of parts after all splits
	std::cout << "BoundaryReconstruction2d: split_parts(): Splitting parts" << std::endl;
	
	// Split parts containing corners
	for(int ipart = 0; ipart < nparts; ipart++)
	{
		std::vector<int> facelist;

		if(ncorners[ipart] == 0)
		{
			int nextface = startface[ipart], curface, nextpoint;

			// traverse through faces in this part
			while(true)
			{
				curface = nextface;
				facelist.push_back(curface);
				facepart(curface,0) = nnparts;
				facepart(curface,1) = facelist.size()-1;

				nextpoint = m->gbfacebp(curface,1);			// get boundary point number of next point (NOT global point number)
				nextface = m->gbpointsb(nextpoint,2);

				if(toRec(ipart,nextface) == 0 || nextface == startface[ipart]) 
				{	break;
				}
			}

			partfaces.push_back(facelist);
			isSplitClosed.push_back(isClosed[ipart]);
			nnparts++;
			continue;
		}

		if(isClosed[ipart])
		{
			int nsplits = ncorners[ipart];
			std::cout << "BoundaryReconstruction2d: split_parts(): Part " << ipart << " is closed and contains " << nsplits << " corners." << std::endl;

			// Change startface for this part to get a corner point face as startface
			// that is, set startface as the face to the right of the first corner point.
			startface[ipart] = m->gbpointsb(corners[ipart][0][1],2);

			std::vector<int> startfac(nsplits);
			startfac[0] = startface[ipart];
			for(int i = 1; i < ncorners[ipart]; i++)
				startfac[i] = m->gbpointsb(corners[ipart][i][1], 2);
			
			// now traverse the parent part in order, looking for startfacs
			int nextface = startface[ipart], curface, nextpoint; 
			//std::cout << "Startface is " << nextface << " " << m->gbface(nextface,0)  << std::endl;
			int nspl = 0;

			//std::cout << "BoundaryReconstruction2d: split_parts(): Starting loop" << std::endl;

			while(true)
			{
				curface = nextface;							// update current face
				facelist.push_back(curface);				// add face to list of faces in current split
				facepart(curface,0) = nnparts;				// store part number that curface belongs to
				facepart(curface,1) = facelist.size()-1;	// position of curface in the facelist of this part

				nextpoint = m->gbfacebp(curface,1);			// get boundary point number of next point (NOT global point number)
				nextface = m->gbpointsb(nextpoint,2);

				int i = 0;
				for(i = 0; i < nsplits; i++)
					if(nextface == startfac[i])				// if we've reached the beginning of the next part 
					{	
						// add the new part to partfces std::vector
						partfaces.push_back(facelist);

						// new split parts are always open
						isSplitClosed.push_back(false);

						// reset std::vector facelist
						facelist.clear();

						nnparts++;									// new number of parts
						nspl++;										// increment number of splits
						break;
					}

				if(i == 0)		// if we've reached the first corner again, we're done
					break;
				// check if we've reached the end of this parent part, either by going beyond the set of markers for this part, or reaching the startface again
				/*if(nextface == startface[ipart]) 
				{
					// add last new part to partfaces
					partfaces.push_back(facelist);
					isSplitClosed.push_back(false);
					nnparts++;
					nspl++;

					// exit from loop
					std::cout << "BoundaryReconstruction2d: split_parts(): Part " << ipart << " split into " << nspl << " part(s) (out of " << nsplits << " part(s))."  << std::endl;
					break;
				}*/
			}
		}

		else
		{
			// just split from startface to corner face 1, and then from corner face 2 to end face, etc
			// Note, however, that the corners are not necessarily in order.
			int nsplits = ncorners[ipart]+1;
			std::vector<int> startfac(nsplits);
			startfac[0] = startface[ipart];
			for(int i = 1; i < ncorners[ipart]+1; i++)
				startfac[i] = m->gbpointsb(corners[ipart][i-1][1], 2);
			
			// now traverse the parent part in order, looking for startfacs
			int nextface = startface[ipart], curface, nextpoint; 
			int nspl = 0;

			while(true)
			{
				curface = nextface;							// update current face
				facelist.push_back(curface);				// add face to list of faces in current split
				facepart(curface,0) = nnparts;
				facepart(curface,1) = facelist.size()-1;	// position of curface in the facelist of this part

				nextpoint = m->gbfacebp(curface,1);			// get boundary point number of next point (NOT global point number)
				nextface = m->gbpointsb(nextpoint,2);
				for(int i = nspl+1; i < ncorners[ipart]+1; i++)
					if(nextface == startfac[i]) 
					{	// if we've reached the beginning of the next part

						// add the new part to partfces std::vector
						partfaces.push_back(facelist);
						isSplitClosed.push_back(false);

						// reset std::vector facelist
						facelist.clear();

						nnparts++;					// new number of parts
						nspl++;										// increment number of splits
					}

				// check if we've reached the end of this parent part, either by going beyond the set of markers for this part, or reaching the startface again
				if(toRec(ipart,nextface) == 0 || nextface == startface[ipart]) 
				{
					// add last new part to partfaces
					partfaces.push_back(facelist);
					isSplitClosed.push_back(false);

					nnparts++;
					nspl++;

					// exit from loop
					std::cout << "BoundaryReconstruction2d: split_parts(): Part " << ipart << "split into " << nspl << "parts (out of " << nsplits << " parts)."  << std::endl;
					break;
				}
			}
		}
	}
}

void BoundaryReconstruction2d::compute_splines(double tol, int maxiter)
{
	sparts = new CSpline[nnparts];

	for(int ipart = 0; ipart < nnparts; ipart++)
	{
		sparts[ipart].setup(m,partfaces[ipart],isSplitClosed[ipart],true,tol,maxiter);
		sparts[ipart].sequence();
		sparts[ipart].compute();
	}
	std::cout << "BoundaryReconstruction2d: compute_splines(): Computed all spline pieces." << std::endl;
}

double BoundaryReconstruction2d::getcoords(int iface, int idim, double u)
{
	return sparts[facepart(iface,0)].getspline(iface,idim,u);
}

// ---------------------------- End of class BoundaryReconstruction2d ---------------------------------------------------------------------//

} // end namespace amc
