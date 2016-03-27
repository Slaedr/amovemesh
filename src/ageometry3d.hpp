/** \brief Header for implementation of the C^0 surface reconstruction of Jiao and Wang, "Reconstructing high-order surfaces for meshing".
 * \date March 14, 2016
 * \author Aditya Kashi
 */

#ifndef __AGEOMETRY3D_H

#ifndef __AMESH3D_H
#include <amesh3d.hpp>
#endif

#ifndef __ALINALG_H
#include <alinalg.hpp>
#endif

#define __AGEOMETRY3D_H

namespace amc {

int factorial(int x);

/// Abstract class for implementing a surface reconstruction using local Taylor expansion fittings
/** The local coordinate frames are decided as described in
 * "X. Jiao and H. Zha, Consistent computation of first- and second-order differential quantities for surface meshes, 2008. In: ACM solid and physical modeling symposium, pp 159-170. ACM."
 * We first compute the normal at a point somehow. This will be the local 'w' axis.
 * The u and v axes are decided as follows. u2 and u3 are selected arbitrarily. u1 is set such that \f$ \mathbf{u} \cdot \mathbf{n} = 0 \f$. \f$ \mathbf{u} \f$ is then normalized.
 * Finally, we set \f$ \mathbf{v} = \mathbf{n} \times \mathbf{u} / ||\mathbf{n} \times \mathbf{u}|| \f$.
 */
class BoundaryReconstruction
{
protected:
	const UMesh* m;
	int degree;							///< Polynomial degree of reconstructed surface
	amat::Matrix<amc_real> fnormals;	///< Face normals
	amat::Matrix<amc_real>* V;			///< Vandermonde matrices for surface points
	amat::Matrix<amc_real>* D;			///< unknowns (various derivatives)
	amat::Matrix<amc_real>* F;			///< RHS of least-squares problem, consisting of point heights or coordinates
	const amc_real s1;					///< Any number (to use for deciding the local coordinate frames)
	const amc_real s2;					///< Any number (to use for deciding the local coordinate frames)

public:
	BoundaryReconstruction(const UMesh* mesh, int deg);
	virtual ~BoundaryReconstruction() { }

	virtual void preprocess();
	virtual void solve();
	virtual void getEdgePoint(const amc_real ratio, const amc_int edgenum, std::vector<amc_real>& point) const = 0;
	virtual void getFacePoint(const std::vector<amc_real>& areacoords, const amc_int facenum, std::vector<amc_real>& point) const = 0;
};

/// Implements WALF reconstruction according to Jiao and Wang's paper, ie, local fittings are calculated at each surface vertex
/** The reconstructed surface at each point passes through that point.
 *
 * Once source of error could be that we need point normals for computing the local u,v,w vectors. Currently point normals are computed as average of surrounding face normals.
 */
class VertexCenteredBoundaryReconstruction : public BoundaryReconstruction
{
	amat::Matrix<amc_real> pnormals;			///< normals at each point, calculated as average of face normals surrounding the point

	amat::Matrix<amc_real>* Q;					///< coordinate transformation (rotation) matrix for each point

	bool isalloc;

	int nders;									///< number of unknowns for the least-squares problems
	std::vector<int> mpo;						///< number of points in stencil for each surface point
	std::vector<int>* stencil;					///< List of bpoint indices of points lying in the stencil of each surface point

	/// convert a point from local coord system of point ibpoin to the global xyz coord system
	void xyz_from_uvw(const amc_int ibpoin, const std::vector<amc_real>& uvwpoint, std::vector<amc_real>& xyzpoint) const;

	/// convert a point from global coord system to the local uvw coord system of point ibpoin
	void uvw_from_xyz(const amc_int ibpoin, const std::vector<amc_real>& xyzpoint, std::vector<amc_real>& uvwpoint) const;

public:
	VertexCenteredBoundaryReconstruction(const UMesh* mesh, int deg);
	~VertexCenteredBoundaryReconstruction();

	/// compute normal, rotation matrix and stencil for each point
	void preprocess();
	/// solve linear least-squares problem at each point
	void solve();
	
	/// Returns coords of a point lying on the 'edgenum' edge, having length coordinate 'ratio' along the edge from point 0 to point 1.
	/** Uses WALF surfaces from the two edge points only.
	 */
	void getEdgePoint(const amc_real ratio, const amc_int edgenum, std::vector<amc_real>& point) const;
	/// Returns coords of a point lying on the face 'facenum' and having area coordinates given by 'areacoords'
	void getFacePoint(const std::vector<amc_real>& areacoords, const amc_int facenum, std::vector<amc_real>& point) const;
};

}
#endif
