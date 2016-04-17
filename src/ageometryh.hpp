/** @file ageometryh.hpp
 * @brief Classes to reconstruct a C^1 (or C^2) piecewise polynomial (cubic) boundary from a piecewise linear boundary given by a hybrid linear mesh.
 * @author Aditya Kashi
 * @date November 2, 2015
 */

#ifndef __AGEOMETRYH_H

#ifndef __AMESH2DHYBRID_H
#include "amesh2dh.hpp"
#endif

#ifndef __ALINALG_H
#include "alinalg.hpp"
#endif

#define __AGEOMETRYH_H 1

namespace amc {

/** Class CSpline constructs a piecewise cubic C^2 spline curve interpolating the boundary points contained by the boundary faces either 
 * (a) having marker in rfl, or
 * (b) listed in facelist. 
 * An overloaded setup() function is provided to distinguish the two situations.
 * Whether the curve is open or closed needs to be specified in isClosed.
 * \note NOTE: The boundary markers specified must form a continuous boundary.
 */

class CSpline
{
	const UMesh2dh* m;
	amat::Matrix<int> rfl;					///< Contains boundary markers to be reconstructed to a spline
	std::vector<int> facelist;				///< Alternative to rfl; stores an ordered list of faces to be reconstructed
	int ndf;							///< number of degrees of freedom per spline - for cubic this is 4
	int dim;
	int nseg;							///< number of spline pieces
	int nspoin;							///< number of control points
	amat::Matrix<int> seq_spoin;				///< stores sequence of global point numbers for use in spline construction
	amat::Matrix<int> seq_bface;				///< for each bface, stores an order number indicating its occurrence order according to contiguity
	amat::Matrix<int> segface;				///< inverse of seq_bface; stores bface number for each segment
	amat::Matrix<double>* scf;
	amat::Matrix<int> toRec;					///< stores for each bface face whether that face is to be reconstricted
	bool isClosed;						///< is the spline curve open or closed?
	bool issequenced;					///< is the list of faces already in sequence?
	bool face_list_available;			///< true if face list is available, false if rfl is available
	double tol;
	int maxiter;

	amat::Matrix<double>* D;					///< D[idim](i) will contain the slope at point 0 of the ith spline piece
	amat::SpMatrix slhs;						///< LHS of the system which is solved for D
	amat::Matrix<double>* srhs;				///< RHS for each dimention

public:
	
	void setup(const UMesh2dh* const mesh, amat::Matrix<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter);
	/**< Use if you want to supply boundary markers for faces to reconstruct. */

	void setup(const UMesh2dh* const mesh, std::vector<int> recmarkers, bool closed, bool sequenced, double _tol, int _maxiter);
	/**< Use if you provide a pre-ordered list of faces to reconstruct. */

	~CSpline();
	
	/// This function arranges the faces to be reconstructed in their sequence of contiguity.
	/** It calculates seq_poin and seq_bface, such that seq_bface(iface) contains the iface-th bface in geometrical order, 
	 * and seq_poin(ibpoin) is the first bpointsb point number of seq_bface(iface).
	 */
	void sequence();

	void compute();
	///< This function computes the spline coeffs and stores them in scf. Depends on sequenced bfaces and points.

	double getspline(int iface, int idim, double t);
	///< returns idim-coordinate of iface-th spline segement with parameter t
};


/** @brief Class BoundaryReconstruction2d handles spline reconstruction of multiple parts of the boundary into seperate c-splines.
 * 
 * It accepts an arbitrary number of boundary parts (BPs) to be reconstructed independently of each other, each consisting of an arbitrary number of boundary markers.
 * It scans each boundary part for corners, and splits them at the corners to get several boundary parts with no corners.
 * The class can also store and retrieve spline coefficients for such multi-boundary-part meshes.
 * \note The mesh supplied must have bpointsb computed.
 */

class BoundaryReconstruction2d
{
	const UMesh2dh* m;								///< NOTE: make sure bpointsb has been computed!
	std::vector<std::vector<int>> marks;				///< to hold boundary markers of all parts
	double cangle;							///< minimum corner angle, above which an intersection is considered a corner
	int nparts;
	int nnparts;
	CSpline* sparts;
	std::vector<int> ncorners;
	std::vector<bool> isClosed;					///< contins true if a (parent) part is closed.
	std::vector<bool> isSplitClosed;				///< contains true if a split part is closed.
	std::vector<int> startface;
	amat::Matrix<int> toRec;						///< nparts x nface array that stores 1 if a face belongs to a part.
	std::vector<std::vector<std::vector<int>>> corners;	///< contains a list of point number and two containing bfaces for each corner point in each part.
	std::vector<std::vector<int>> partfaces;			///< stores a list of ordered faces for each part
	amat::Matrix<int> facepart;					///< stores part no. and local face number in that part, for each boundary face

public:
	void setup(const UMesh2dh* const mesh, int num_parts, std::vector<std::vector<int>> boundary_markers, double angle_threshold);

	~BoundaryReconstruction2d();
	
	void preprocess();
	/**< Determines whether each part is open or closed, and stores a starting bface number for each part 
		NOTE: Make sure amesh2d::compute_boundary_points() has been executed.*/
	
	void detect_corners();
	/**< Detect corners in each part based on dot-product of normals of adjacent faces becoming too small. */

	void read_corners(std::string cname);
	/**< Can accept corners from a file rather than trying to detect them */

	void split_parts();
	/**< Splits parts based on corner points. */
	
	void compute_splines(double tol, int maxiter);
	/**< Calls the compute() function of class CSpline to compute spline coefficients of all parts. */
	
	/**	Function to return coordinates of the curve.
		NOTE: the argument iface must correspond to a face which was reconstructed!!
	*/
	double getcoords(int iface, int idim, double u);
	
	//void writeCoeffs(std::string fname);
	
	//void readCoeffs(std::string fname);
};


// ---------------------------- End of class BoundaryReconstruction2d ---------------------------------------------------------------------//


}
#endif
