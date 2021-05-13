/** \file ageometry3d.hpp
 * \brief Header for implementation of the C^0 surface reconstruction of Jiao and Wang,
 * "Reconstructing high-order surfaces for meshing".
 * Other references:
 * Jiao and Zha, "Consistent computation of first- and second-order differential
 * quantities for surface meshes".
 *
 * \date March 14, 2016
 * \author Aditya Kashi
 */

#ifndef AMC_GEOMETRY3D_H
#define AMC_GEOMETRY3D_H

#include "amesh3d.hpp"
#include "alinalg.hpp"

namespace amc {

/// Recursively computes the factorial of an integer
int factorial(int x);

/// C1-discontinuity detection in boundaries of 3D linear meshes
class DiscontinuityDetection
{
protected:
	const UMesh* const m;
	const amat::Matrix<amc_real>* const fnormals;		///< Normals to faces
	amat::Matrix<amc_real> etangents;					///< tangents of edges
	/// Maximum angle between two faces to consider them as C1 continuous
	const double maxangle;
	/// Max angle between two edges to consider them as part of the same feature curve
	const double maxedgeangle;
	/// number of feature curves in the boundary
	int ncurves;
	/// Stores an ordered list of edges in each feature curve
	std::vector<std::vector<amc_int>> fecurve;
	/// Stores the feature curve that a boundary-edge (b-edge) belongs to, for each b-edge
	std::vector<int> febedge;
	/// Stores the feature-curve number that each boundary point belongs to
	std::vector<int> febpoint;

	/// For each boundary point, stores 0 if it's not a corner,
	///  and an integer indicating the type of corner if it is one
	std::vector<int> cornerpoint;

public:
	DiscontinuityDetection(const UMesh* const mesh,
						   const amat::Matrix<amc_real>* const fnormal,
						   const double max_angle, const double max_edge_angle);

	virtual void detect_C1_discontinuities();

	// getter functions
	int gncurves() const { return ncurves; }
	amc_int gfecurve(int icurve, int ied) const { return fecurve[icurve][ied]; }
	int gfebedge(amc_int ied) const { return febedge[ied]; }
	int gfebpoint(amc_int ipoin) const { return febpoint[ipoin]; }
	int gcornerpoint(amc_int ipoin) const { return cornerpoint[ipoin]; }
};

/// Abstract class for implementing a surface reconstruction using local Taylor expansion fittings
/** The local coordinate frames are decided as described in
 * "X. Jiao and H. Zha, Consistent computation of first- and second-order differential
 * quantities for surface meshes, 2008. In: ACM solid and physical modeling symposium,
 * pp 159-170. ACM."
 * We first compute the normal at a point somehow. This will be the local 'w' axis.
 * The u and v axes are decided as follows. u2 and u3 are selected arbitrarily.
 * u1 is set such that \f$ \mathbf{u} \cdot \mathbf{n} = 0 \f$.
 * \f$ \mathbf{u} \f$ is then normalized.
 * Finally, we set
 * \f$ \mathbf{v} = \mathbf{n} \times \mathbf{u} / ||\mathbf{n} \times \mathbf{u}|| \f$.
 */
class BoundaryReconstruction
{
protected:
	const UMesh* m;
	int degree;							///< Polynomial degree of reconstructed surface
	amat::Matrix<amc_real> fnormals;	///< Face normals
	amat::Matrix<amc_real> pnormals;	///< normals at each point
	amat::Matrix<amc_real> face_center;	///< contains coordinates of center of each face
	std::vector<amc_real> farea;		///< Areas of boundary faces

	/// Vandermonde matrices for surface points
	amat::Matrix<amc_real>* V;
	/// Unknowns (various derivatives)
	amat::Matrix<amc_real>* D;
	/// RHS of least-squares problem, consisting of point heights or coordinates
	amat::Matrix<amc_real>* F;
	/// A string describing the type of stencil to use:
	///  "half" or "full". "full" results in a more extended stencil
	std::string stencilType;

	/// Any number (to use for deciding the local coordinate frames)
	const amc_real s1;
	/// Any number (to use for deciding the local coordinate frames)
	const amc_real s2;

	/// Starting index for Taylor polynomials
	/** 0 for allowing a constant term in the Taylor series and
	 * 1 for starting the series with first-order terms.
	 */
	const int istart;
	
	/// computes vertex normals by using inverse distance to face-centers as weights
	void computePointNormalsInverseDistance();

	/// computes vertex normals by using area of faces as weights
	void computePointNormalsArea();

public:
	/// constructor; also computes face-normals for each b-face
	BoundaryReconstruction(const UMesh* mesh, int deg, std::string stencil_type, int i_start);
	virtual ~BoundaryReconstruction() { }

	virtual void preprocess();
	virtual void solve();
	virtual void getEdgePoint(const amc_real ratio, const amc_int edgenum,
							  std::vector<amc_real>& point) const = 0;
	virtual void getFacePoint(const std::vector<amc_real>& areacoords,
							  const amc_int facenum, std::vector<amc_real>& point) const = 0;
};

/// Implements WALF reconstruction ie, local fittings are calculated at each surface vertex
/** According to Jiao and Wang's paper.
 * The reconstructed surface at each point passes through that point.
 * One source of error is that we need point normals for computing the local u,v,w vectors.
 * Currently point normals are computed as average of surrounding face normals.
 */
class VertexCenteredBoundaryReconstruction : public BoundaryReconstruction
{
	/// Coordinate transformation (rotation) matrix for each point
	amat::Matrix<amc_real>* Q;

	bool isalloc;

	/// True if Jiao and Zha's safeguarded solution of least-squares is to be used
	bool safeguard;
	/// 1-norm upper limit for order-downgrade to not be done
	double normlimit;
	/// Flag containing the reconstruction order at each surface point
	std::vector<int> rec_order;

	/// Number of unknowns for the least-squares problems
	int nders;
	/// Number of points in stencil for each surface point
	std::vector<int> mpo;
	/// List of bpoint indices of points lying in the stencil of each surface point
	std::vector<int>* stencil;

	/// convert a point from local coord system of point ibpoin to the global xyz coord system
	void xyz_from_uvw(const amc_int ibpoin, const std::vector<amc_real>& uvwpoint,
					  std::vector<amc_real>& xyzpoint) const;

	/// convert a point from global coord system to the local uvw coord system of point ibpoin
	void uvw_from_xyz(const amc_int ibpoin, const std::vector<amc_real>& xyzpoint,
					  std::vector<amc_real>& uvwpoint) const;

public:
	VertexCenteredBoundaryReconstruction(const UMesh* mesh, int deg, std::string stencilsize,
										 bool _safeguard, double norm_limit);
	~VertexCenteredBoundaryReconstruction();

	/// compute normal, rotation matrix and stencil for each point
	void preprocess();
	/// solve linear least-squares problem at each point
	void solve();
	
	/** Returns coords of a point lying on the 'edgenum' edge,
	 * having length coordinate 'ratio' along the edge from point 0 to point 1.
	 *
	 * Uses WALF surfaces from the two edge points only.
	 */
	void getEdgePoint(const amc_real ratio, const amc_int edgenum,
					  std::vector<amc_real>& point) const;

	/** Returns coords of a point lying on the face 'facenum' and
	 * having area coordinates given by 'areacoords'
	 */
	void getFacePoint(const std::vector<amc_real>& areacoords, const amc_int facenum,
					  std::vector<amc_real>& point) const;
};

/** \brief Computes reconstructed surface using WALF with local Taylor polynomials
 * fitted to each face center, based on Jiao and Wang.
 *
 * This has the advantage that the normals at face centers are known better than
 * normals at vertices.
 * \note NOTE one important difference between this and VertexCenteredBoundaryReconstruction :
 * At each face, the value of the height function is NOT required to be zero.
 * So we have one more unknown compared to VertexCenteredBoundaryReconstruction.
 */
class FaceCenteredBoundaryReconstruction : public BoundaryReconstruction
{
	/// Coordinate transformation (rotation) matrix for each point
	amat::Matrix<amc_real>* Q;

	bool isalloc;

	/// true if Jiao and Zha's safeguarded solution of least-squares is to be used
	bool safeguard;
	/// 1-norm upper limit for order-downgrade to not be done
	double normlimit;
	/// Flag containing the reconstruction order at each surface point
	std::vector<int> rec_order;

	/// number of unknowns for the least-squares problems
	int nders;
	/// number of points in stencil for each surface point
	std::vector<int> mpo;
	/// List of bpoint indices of points lying in the stencil of each surface point
	std::vector<int>* stencil;

	/// Number of reconstructions to do, using normals from previous iteration
	/// (not used currently)
	int niter;

	/// convert a point from local coord system of point ibpoin to the global xyz coord system
	void xyz_from_uvw(const amc_int ibpoin, const std::vector<amc_real>& uvwpoint,
					  std::vector<amc_real>& xyzpoint) const;

	/// convert a point from global coord system to the local uvw coord system of point ibpoin
	void uvw_from_xyz(const amc_int ibpoin, const std::vector<amc_real>& xyzpoint,
					  std::vector<amc_real>& uvwpoint) const;
	
	/** Computes vertex normals by using a weighted average of face normals of
	 * faces surrounding the vertex
	 *
	 * Needed for computing weights of the weighted least-squares procedure.
	 */
	void computePointNormals();

public:
	FaceCenteredBoundaryReconstruction(const UMesh* mesh, int deg, std::string stencilsize,
									   bool _safeguard, double norm_limit);

	~FaceCenteredBoundaryReconstruction();

	/// Rotation matrix and stencil for each face
	/** For each boundary face, the stencil is computed as follows.
	 * First, the vertices of the face are added to the stencil.
	 * Then, all vertices surrounding the initially-added vertices are also added.
	 * NOTE: Just adding all vertices of faces surrounding the face should also suffice;
	 * perhaps we should be doing that.
	 */
	void preprocess();

	/// Solve linear least-squares problem at each point
	void solve();
	
	/** Returns coords of a point lying on the 'edgenum' edge, having
	 * length coordinate 'ratio' along the edge from point 0 to point 1.
	 *
	 * Uses WALF surfaces from the two edge points only.
	 */
	void getEdgePoint(const amc_real ratio, const amc_int edgenum,
					  std::vector<amc_real>& point) const;

	/** Returns coords of a point lying on the face 'facenum' and having area coordinates
	 * given by 'areacoords'
	 *
	 * There is no averaging involved in this case, unlike the Vertex-centered reconstruction.
	 * The height at the point specified by areacoords is determined by the reconstruction
	 * at the face that contains the point.
	 */
	void getFacePoint(const std::vector<amc_real>& areacoords, const amc_int facenum,
					  std::vector<amc_real>& point) const;
};

}
#endif
