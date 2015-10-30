/** Classes to reconstruct a C^1 (or C^2) piecewise polynomial (cubic) boundary from a piecewise linear boundary given by a linear mesh.
   Aditya Kashi
   September 4, 2015
*/

#ifndef __AMESH2DGENERAL_H
#include "amesh2d.hpp"
#endif

#ifndef __ALINALG_H
#include "alinalg.hpp"
#endif

#define __AGEOMETRY_H 1

using namespace std;
using namespace amat;
using namespace acfd;

namespace acfd {

/** Class CSpline constructs a piecewise cubic C^2 spline curve interpolating the boundary points contained by the boundary faces either (a) having marker in rfl,
	or (b) listed in facelist. An overloaded setup() function is provided to distinguish the two situations.
	Whether the curve is open or closed needs to be specified in isClosed.
	NOTE: The boundary markers specified must form a continuous boundary.
*/

class CSpline
{
	UMesh2d* m;
	Matrix<int> rfl;					//< Contains boundary markers to be reconstructed to a spline
	vector<int> facelist;				//< Alternative to rfl; stores an ordered list of faces to be reconstructed
	int ndf;							//< number of degrees of freedom per spline - for cubic this is 4
	int dim;
	int nseg;							//< number of spline pieces
	int nspoin;							//< number of control points
	Matrix<int> seq_spoin;				//< stores sequence of global point numbers for use in spline construction
	Matrix<int> seq_bface;				//< for each bface, stores an order number indicating its occurrence order according to contiguity
	Matrix<int> segface;				//< inverse of seq_bface; stores bface number for each segment
	Matrix<double>* scf;
	Matrix<int> toRec;					//< stores for each bface face whether that face is to be reconstricted
	bool isClosed;						//< is the spline curve open or closed?
	bool issequenced;					//< is the list of faces already in sequence?
	bool face_list_available;			//< true if face list is available, false if rfl is available
	double tol;
	int maxiter;

	Matrix<double>* D;					//< D[idim](i) will contain the slope at point 0 of the ith spline piece
	SpMatrix slhs;						//< LHS of the system which is solved for D
	Matrix<double>* srhs;				//< RHS for each dimention

public:
	
	void setup(UMesh2d* mesh, Matrix<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter);
	/**< Use if you want to supply boundary markers for faces to reconstruct. */

	void setup(UMesh2d* mesh, vector<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter);
	/**< Use if you provide a pre-ordered list of faces to reconstruct. */

	~CSpline();
	
	void sequence();
	///< This function arranges the faces to be reconstructed in their sequence of contiguity.
	///< It calculates seq_poin and seq_bface, such that seq_bface(iface) contains the iface-th bface in geometrical order, and seq_poin(ibpoin)
	///<  is the first bpointsb point number of seq_bface(iface).

	void compute();
	///< This function computes the spline coeffs and stores them in scf. Depends on sequenced bfaces and points.

	double getspline(int iface, int idim, double t);
	///< returns idim-coordinate of iface-th spline segement with parameter t
};

// --- CSpline implementation ---

void CSpline::setup(UMesh2d* mesh, Matrix<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter)
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
		cout << "CSpline: setup(): Closed curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << endl;
	else
		cout << "CSpline: setup(): Open curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << endl;

	// allocate stuff
	
	scf = new Matrix<double>[dim];
	D = new Matrix<double>[dim];
	srhs = new Matrix<double>[dim];
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

void CSpline::setup(UMesh2d* mesh, vector<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter)
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
		cout << "CSpline: setup(): Closed curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << endl;
	else
		cout << "CSpline: setup(): Open curve. Number of pieces = " << nseg << ", number of control points = " << nspoin << endl;

	// allocate stuff
	
	scf = new Matrix<double>[dim];
	D = new Matrix<double>[dim];
	srhs = new Matrix<double>[dim];
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
				if(nextface == -1) cout << "UMesh2d: sequence(): ! Point not found!" << endl;
				if(nextface == startface) break;
			}
			seq_spoin(np) = seq_spoin(0);		// close the point list
		}
		else		// curve is open
		{
			nextface = startface;

			//first find the last element
			cout << "UMesh2d: sequence(): Finding first face of open boundary." << endl;
			while(true)
			{
				curface = nextface;
				nextpoin = m->gbface(curface,0);
				
				// search for the other face which contains nextpoin
				nextface = -1;
				for(int iface = 0; iface < m->gnface(); iface++)
					if(toRec(iface) == 1 && m->gbface(iface,1) == nextpoin) nextface = iface;
				
				if(nextface == -1) {
					cout << "UMesh2d: sequence(): Reached start of boundary to be reconstructed." << endl;
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
					cout << "UMesh2d: sequence(): " << nf << " boundary faces and " << np << " boundary points sequenced."  << endl;
					if(nf != nseg) cout << "UMesh2d: sequence(): ! All faces are not accounted for!" << endl;
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
			cout << "CSpline: sequence(): Boundary part is sequenced, but facelist is not available!!" << endl;
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
	//seq_bface.mprint();
	//seq_spoin.mprint();
}

void CSpline::compute()
{
	if(isClosed)
	{
		// Assemble RHS vector for each dimension
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
		// Assemble RHS vector for each dimension
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
	cout << "CSpline: compute(): Solving linear systems" << endl;
	Matrix<double> d0(nseg+1,1);
	d0.zeros();
	for(int idim = 0; idim < dim; idim++)
		D[idim] = sparseCG_d(&slhs, srhs[idim], d0, tol, maxiter);
	
	// get coeffs
	//cout << "CSpline: compute(): Getting spline coefficients" << endl;
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
	cout << "CSpline: compute(): Done." << endl;
}

double CSpline::getspline(int iface, int idim, double t)
{
	return scf[idim].get(seq_bface.get(iface),0) + scf[idim].get(seq_bface.get(iface),1)*t + scf[idim].get(seq_bface.get(iface),2)*t*t + scf[idim].get(seq_bface.get(iface),3)*t*t*t;
}

// --------------------- End of class CSpline --------------------------------------------------------------------------------------------//


/** Class BoundaryReconstruction2d handles spline reconstruction of multiple parts of the boundary into seperate c-splines.
	It accepts an arbitrary number of boundary parts (BPs) to be reconstructed independently of each other, each consisting of an arbitrary number
	of boundary markers.
	It scans each boundary part for corners, and splits them at the corners to get several boundary parts with no corners.
	NOTE: This process of splitting parts at corners modifies the original mesh file!
	The class can also store and retrieve spline coefficients for such multi-boundary-part meshes.
*/

class BoundaryReconstruction2d
{
	UMesh2d* m;								///< NOTE: make sure bpointsb has been computed!
	vector<vector<int>> marks;				///< to hold boundary markers of all parts
	double cangle;							///< minimum corner angle, above which an intersection is considered a corner
	int nparts;
	int nnparts;
	CSpline* sparts;
	vector<int> ncorners;
	vector<bool> isClosed;					///< contins true if a (parent) part is closed.
	vector<bool> isSplitClosed;				///< contains true if a split part is closed.
	vector<int> startface;
	Matrix<int> toRec;						///< nparts x nface array that stores 1 if a face belongs to a part.
	vector<vector<vector<int>>> corners;	///< contains a list of point number and two containing bfaces for each corner point in each part.
	vector<vector<int>> partfaces;			///< stores a list of ordered faces for each part
	Matrix<int> facepart;					///< stores part no. and local face number in that part, for each boundary face

public:
	void setup(UMesh2d* mesh, int num_parts, vector<vector<int>> boundary_markers, double angle_threshold);

	~BoundaryReconstruction2d();
	
	void preprocess();
	/**< Determines whether each part is open or closed, and stores a starting bface number for each part 
		NOTE: Make sure amesh2d::compute_boundary_points() has been executed.*/
	
	void detect_corners();
	/**< Detect corners in each part based on dot-product of normals of adjacent faces becoming too small. */

	void read_corners(string cname);
	/**< Can accept corners from a file rather than trying to detect them */

	void split_parts();
	/**< Splits parts based on corner points. */
	
	void compute_splines(double tol, int maxiter);
	/**< Calls the compute() function of class CSpline to compute spline coefficients of all parts. */
	
	/**	Function to return coordinates of the curve.
		NOTE: the argument iface must correspond to a face which was reconstructed!!
	*/
	double getcoords(int iface, int idim, double u);
	
	//void writeCoeffs(string fname);
	
	//void readCoeffs(string fname);
};

void BoundaryReconstruction2d::setup(UMesh2d* mesh, int num_parts, vector<vector<int>> boundary_markers, double angle_threshold)
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
	cout << "BoundaryReconstruction2d: preprocess(): toRec populated." << endl;
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
		//cout << "BoundaryReconstruction2d: preprocess(): Determining whether part " << ipart << " is open or closed" << endl;
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
	cout << "BoundaryReconstruction2d: preprocess(): Done." << endl;
}

/**	Detects corner points by computing dot-products of normals of the two faces associated with a point and comparing it with cos(theta) for each point,
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

	cout << "BoundaryReconstruction2d: detect_corners(): Searching each part for corners" << endl;
	
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
					vector<int> c(4);
					c[0] = p2;						// global point number
					c[1] = ipoin;					// boundary point number
					c[2] = m->gbpointsb(ipoin,1);
					c[3] = m->gbpointsb(ipoin,2);
					corners[ipart].push_back(c);
					ncorners[ipart]++;
				}
			}
		}
		cout << "BoundaryReconstruction2d: detect_corners(): Part " << ipart << " contains " << ncorners[ipart] << " corner(s)" << endl;
		cout << "  at points ";
		for(int i = 0; i < ncorners[ipart]; i++)
			cout << corners[ipart][i][0] << ", ";
		cout << endl;
	}
}

void BoundaryReconstruction2d::read_corners(string cname)
{
	// read corner data from file
	ifstream fin(cname);
	
	cout << "BoundaryReconstruction2d: read_corners(): Reading file" << endl;
	int dum;
	vector<int> c(4,0);								// vector of 4 ints with value 0

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

	cout << "BoundaryReconstruction2d: read_corners(): iterating over boundary points" << endl;
	
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
	cout << "BoundaryReconstruction2d: split_parts(): Splitting parts" << endl;
	
	// Split parts containing corners
	for(int ipart = 0; ipart < nparts; ipart++)
	{
		vector<int> facelist;

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
			cout << "BoundaryReconstruction2d: split_parts(): Part " << ipart << " is closed and contains " << nsplits << " corners." << endl;

			// Change startface for this part to get a corner point face as startface
			// that is, set startface as the face to the right of the first corner point.
			startface[ipart] = m->gbpointsb(corners[ipart][0][1],2);

			vector<int> startfac(nsplits);
			startfac[0] = startface[ipart];
			for(int i = 1; i < ncorners[ipart]; i++)
				startfac[i] = m->gbpointsb(corners[ipart][i][1], 2);
			
			// now traverse the parent part in order, looking for startfacs
			int nextface = startface[ipart], curface, nextpoint; 
			//cout << "Startface is " << nextface << " " << m->gbface(nextface,0)  << endl;
			int nspl = 0;

			//cout << "BoundaryReconstruction2d: split_parts(): Starting loop" << endl;

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
						// add the new part to partfces vector
						partfaces.push_back(facelist);

						// new split parts are always open
						isSplitClosed.push_back(false);

						// reset vector facelist
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
					cout << "BoundaryReconstruction2d: split_parts(): Part " << ipart << " split into " << nspl << " part(s) (out of " << nsplits << " part(s))."  << endl;
					break;
				}*/
			}
		}

		else
		{
			// just split from startface to corner face 1, and then from corner face 2 to end face, etc
			// Note, however, that the corners are not necessarily in order.
			int nsplits = ncorners[ipart]+1;
			vector<int> startfac(nsplits);
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

						// add the new part to partfces vector
						partfaces.push_back(facelist);
						isSplitClosed.push_back(false);

						// reset vector facelist
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
					cout << "BoundaryReconstruction2d: split_parts(): Part " << ipart << "split into " << nspl << "parts (out of " << nsplits << " parts)."  << endl;
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
	cout << "BoundaryReconstruction2d: compute_splines(): Computed all spline pieces." << endl;
}

double BoundaryReconstruction2d::getcoords(int iface, int idim, double u)
{
	return sparts[facepart(iface,0)].getspline(iface,idim,u);
}

// ---------------------------- End of class BoundaryReconstruction2d ---------------------------------------------------------------------//


/** Class HermiteSpline2d constructs one spline curve from all the boundary faces of a mesh which have markers contained in rfl.
   Currently reconstructs a C^1 boundary.
   But we need to ensure that out of the two possible tangents for each face, the correct one is chosen for consistency. For this, we calculate normal based on intfac,
    as face nodes in intfac are ordered to point outwards.
   NOTE: Make sure the boundaries with markers in rfl are contiguous.
   CAUTION: There are issues with the functionality of this implmenetation. Spurious loops are created at control points, strangely.
*/

class HermiteSpline2d
{
	typedef double(HermiteSpline2d::*Basispointer)(double);
	UMesh2d* m;												// NOTE: make sure compute_topological() has been done for this mesh before passing!
	Matrix<int> rfl;										// the markers of the boundaries to be reconstructed
	int nseg;												// number of curve segments - equal to number of boundary faces to reconstruct
	int ndeg;												// degree of spline interpolation - usually 3
	int ndf;												// Number of 'degrees of freedom' per face - this 1+ndeg.
	Matrix<double>* cf;										// to store geometric coefficients
	Matrix<double> gallfa;									// Normals etc of boundary faces in bfaces
	Matrix<int> toRec;										// Contains 1 if the corresponding intfac is to be reconstructed, otherwise contains 0
	Matrix<int> isCorner;									// contains 1 of this point is a corner point
	Matrix<double> ptangents;								// Contains average tangent at each boundary point. Uses UMesh2d::bpoints.
	Basispointer* F;										// Array of function pointers for Hermite basis functions
	double angle_threshold;									// Minimum angle (in radians) between two tangents for point to qualify as corner point

	bool store_intfac;										// Whether cf is stored according to intfac (true) or bface (false)
	int nndim;
	int nnbface;

	// definitions of the 4 Hermite basis functions
	double f0(double u) { return 2*pow(u,3) - 3*u*u + 1; }
	double f1(double u) { return -2*pow(u,3) + 3*u*u; }
	double f2(double u) { return pow(u,3) - 2*u*u + u; }
	double f3(double u) { return u*u*u - u*u; }

public:
	void setup(UMesh2d* mesh, Matrix<int> rflags, double angle)
	{
		m = mesh;
		rfl = rflags;
		ndeg = 3;
		ndf = 1+ndeg;

		angle_threshold = angle;

		F = new Basispointer[ndeg+1];

		F[0] = &HermiteSpline2d::f0;
		F[1] = &HermiteSpline2d::f1;
		F[2] = &HermiteSpline2d::f2;
		F[3] = &HermiteSpline2d::f3;
		//cout << (this->*F[0])(2.1);

		toRec.setup(m->gnbface(),1);
		toRec.zeros();
		for(int iface = 0; iface < m->gnbface(); iface++)
		{
			for(int im = 0; im < rfl.msize(); im++)
			{
				if(m->gbface(iface,m->gnnofa()) == rfl(im)) toRec(m->gifbmap(iface)) = 1;
			}
		}

		isCorner.setup(m->gnpoin(),1);
		isCorner.zeros();

		ptangents.setup(m->gnpoin(), m->gndim());			// stores components of unit tangent at each boundary point by averaging tangents of faces around that point.
		ptangents.zeros();

		cf = new Matrix<double>[m->gndim()];
		for(int idim = 0; idim < m->gndim(); idim++) {
			cf[idim].setup(m->gnbface(),ndeg+1);
			cf[idim].zeros();
		}

		gallfa.setup(m->gnbface(),4);		// gallfa(iseg,0) contains normalized tangent x-component at point 0 of corresponding intfac segment
		gallfa.zeros();

		store_intfac = true;
	}

	~HermiteSpline2d()
	{
		delete [] F;
		delete [] cf;
	}

	void compute_splines()
	{
		// get point tangents first
		cout << "HermiteSpline2D: compute_splines(): Getting tangents of boundary points." << endl;
		for(int ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
		{
			// in 2D, bpoints is not an array of stl vectors - we know there are only two faces surrounding each boundary point
			// if edge is bounded by (x1,y1) and (x2,y2), tangent to edge is (x2-x1)i + (y2-y1)j

			// first check if this point belongs to a face that needs to be reconstructed
			if(toRec(m->gbpoints(ipoin,1))==1  || toRec(m->gbpoints(ipoin,2))==1)
			{

				vector<int> en(2);

				// get intfac faces containing this point
				en[0] = m->gbpoints(ipoin,1);
				en[1] = m->gbpoints(ipoin,2);
				//cout << en[0] << " " << en[1] << endl;

				double xt, yt, dotx = 1, doty = 1, mag;

				// iterate over the two faces containing ipoin
				for(int i = 0; i < 2; i++)
				{
					xt = m->gcoords(m->gintfac(en[i],3),0) - m->gcoords(m->gintfac(en[i],2),0);
					yt = m->gcoords(m->gintfac(en[i],3),1) - m->gcoords(m->gintfac(en[i],2),1);
					mag = sqrt(xt*xt + yt*yt);
					ptangents(m->gbpoints(ipoin,0),0) += xt;
					ptangents(m->gbpoints(ipoin,0),1) += yt;
					dotx *= xt/mag;
					doty *= yt/mag;
				}
				ptangents(m->gbpoints(ipoin,0),0) /= 2.0;
				ptangents(m->gbpoints(ipoin,0),1) /= 2.0;

				// now normalize the averaged tangent vector
				mag = sqrt(ptangents(m->gbpoints(ipoin,0),0)*ptangents(m->gbpoints(ipoin,0),0) + ptangents(m->gbpoints(ipoin,0),1)*ptangents(m->gbpoints(ipoin,0),1));
				ptangents(m->gbpoints(ipoin,0),0) /= mag;
				ptangents(m->gbpoints(ipoin,0),1) /= mag;

				// check if this is a corner point
				double dotp = dotx+doty;			// dot product of tangents of the two faces
				if(dotp < cos(angle_threshold)) { 
					isCorner(m->gbpoints(ipoin,0)) = 1;				// set isCorner for global point number corresponding to boundary point ipoin
					cout << "Boundary point " << m->gbpoints(ipoin,0) << " is a corner point!" << endl;
				}
			}
		}

		// iterate over boundary faces of the mesh
		cout << "HermiteSpline2D: compute_splines(): Iterating over boundary faces to get tangents." << endl;
		vector<int> pnts(m->gnnofa());		// to store point numbers of points making up iface
		for(int iface = 0; iface < m->gnbface(); iface++)
		{
			if(toRec(iface) == 1)
			{
				// get tangent for this face
				for(int inofa = 0; inofa < m->gnnofa(); inofa++)
					pnts[inofa] = m->gintfac(iface,inofa+2);

				// magnitude of tangent (x2-x1)i + (y2-y1)j
				double mag = sqrt((m->gcoords(pnts[1],0)-m->gcoords(pnts[0],0))*(m->gcoords(pnts[1],0)-m->gcoords(pnts[0],0)) + (m->gcoords(pnts[1],1)-m->gcoords(pnts[0],1))*(m->gcoords(pnts[1],1)-m->gcoords(pnts[0],1)));
				
				for(int inofa = 0; inofa < m->gnnofa(); inofa++)
				{
					if(isCorner(pnts[inofa]) == 0)		// if not a corner
					{
						gallfa(iface,2*inofa) = ptangents(pnts[inofa],0);
						gallfa(iface,2*inofa+1) = ptangents(pnts[inofa],1);
					}
					else								// if pnts[inofa] is a corner, use this face's tangent at that point
					{
						gallfa(iface,2*inofa) = (m->gcoords(pnts[1],0) - m->gcoords(pnts[0],0))/mag;
						gallfa(iface,2*inofa+1) = (m->gcoords(pnts[1],1) - m->gcoords(pnts[0],1))/mag;
					}
				}

				/* Note that gallfa(*,0) contains x-component of tangent at point 1, gallfa(*,1) contains y-component of tangent at point 1,
				    gallfa(*,2) contains x-component of tangent at point 2, and gallfa(*,3) contains y-component of tangent at points 2. */
			}
		}

		// iterate over boundary faces of mesh again to calculate geometric spline coeffs from gallfa and intfac
		cout << "HermiteSpline2d: compute_splines(): Iterating over boundary faces again to compute geometric coefficients of the splines" << endl;
		for(int iface = 0; iface < m->gnbface(); iface++)
			if(toRec(iface) == 1)
			{
				for(int idim = 0; idim < m->gndim(); idim++)
				{
					cf[idim](iface,0) = m->gcoords(m->gintfac(iface,2),idim);
					cf[idim](iface,1) = m->gcoords(m->gintfac(iface,3),idim);
					cf[idim](iface,2) = gallfa(iface,idim);
					cf[idim](iface,3) = gallfa(iface,idim+2);
				}
			}
	}

	double spline(int iface, int idim, double u)
	// returns idim-coordinate of the spline-reconstruction of iface (of intfac or bface) corresponding to parameter u.
	{
		double val = 0;
		for(int i = 0; i < ndf; i++)
			val += (this->*F[i])(u)*cf[idim](iface,i);
		return val;
	}

	void writeGeomCoeffsToFile(string fname)
	{
		ofstream ofs(fname);
		if(store_intfac == true)
		{
			ofs << "Dims " << m->gndim() <<  " DFs " << ndf << " nbface " << m->gnbface() <<  '\n';
			for(int idim = 0; idim < m->gndim(); idim++)
			{
				ofs << idim << '\n';
				for(int inface = 0; inface < m->gnbface(); inface++)
				{
					ofs << m->gbifmap(inface);
					for(int idf = 0; idf < ndf; idf++)
						ofs << " " << cf[idim](inface,idf);
					ofs << '\n';
				}
			}
		}
		else
		{
			ofs << "Dims " << nndim  << " DFs " << ndf << " nbface " << nnbface  <<  '\n';
			for(int idim = 0; idim < nndim; idim++)
			{
				ofs << idim << '\n';
				for(int inface = 0; inface < nnbface; inface++)
				{
					ofs << inface;
					for(int idf = 0; idf < ndf; idf++)
						ofs << " " << cf[idim](inface,idf);
					ofs << '\n';
				}
			}
		}
		ofs.close();
	}

	void readGeomCoeffsFromFile(string fname)
	// Reads a file written by writeGeomCoeffsToFile() and stores the data in cf.
	// NOTE: Faces are arranged in cf according to BFACE numbers now.
	{
		ifstream fin(fname);
		string dum; int ddim; int nbface;
		fin >> dum; fin >> nndim; fin >> dum; fin >> ndf; fin >> dum; fin >> nnbface;
		F = new Basispointer[ndf];

		F[0] = &HermiteSpline2d::f0;
		F[1] = &HermiteSpline2d::f1;
		F[2] = &HermiteSpline2d::f2;
		F[3] = &HermiteSpline2d::f3;

		cf = new Matrix<double>[nndim];
		for(int idim = 0; idim < nndim; idim++) {
			cf[idim].setup(nnbface,ndf);
			cf[idim].zeros();
		}
		
		int iface;		// to hold bface number of each face

		for(int idim = 0; idim < nndim; idim++)
		{
			fin >> dum;
			for(int inface = 0; inface < nnbface; inface++)
			{
				fin >> iface;
				for(int idf = 0; idf < ndf; idf++)
					fin >> cf[idim](iface,idf);
			}
		}

		fin.close();
		store_intfac = false;
	}
};

}
