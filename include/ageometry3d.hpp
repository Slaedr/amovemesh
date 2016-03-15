/** \brief Implementation of the C^0 surface reconstruction of Jiao and Wang, "Reconstructing high-order surfaces for meshing".
 * \date March 14, 2016
 * \author Aditya Kashi
 */

#ifndef __AMESH3D_H
#include <amesh3d.hpp>
#endif

#define __AGEOMETRY3D_H

namespace amc {

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
	const UMesh3d* m;
	int degree;							///< Polynomial degree of reconstructed surface
	amat::Matrix<amc_real> fnormals;	///< Face normals
	amat::Matrix<amc_real>* V;			///< Vandermonde matrices for surface points
	amat::Matrix<amc_real>* X;			///< unknowns (various derivatives)
	amat::Matrix<amc_real>* F;			///< RHS of least-squares problem, consisting of point heights or coordinates
	const amc_real u2;					///< Any number (to use for deciding the local coordinate frames)
	const amc_real u3;					///< Any number (to use for deciding the local coordinate frames)

public:
	BoundaryReconstruction(const UMesh3d* mesh, int deg);
	~BoundaryReconstruction();

	virtual std::vector<amc_real> getEdgePoint(const amc_real ratio, const amc_int edgenum) = 0;
	virtual std::vector<amc_real> getFacePoint(const std::vector<amc_real> areacoords, const amc_int facenum) = 0;
};

BoundaryReconstruction(const UMesh3d* mesh, int deg) 
	: m(mesh), degree(deg), s2(1.0), s3(2.0)
{
	V = new amat::Matrix<amc_real>[m->gnbpoin()];
	X = new amat::Matrix<amc_real>[m->gnbpoin()];
	F = new amat::Matrix<amc_real>[m->gnbpoin()];
	fnormals.setup(m->gnface(), m->gndim());

	// compute unit face normals
	int iface;
	amc_real x1,y1,z1,x2,y2,z2,mag;
	for(iface = 0; iface < m->gnface(); iface++)
	{
		x1 = m->gcoords(m->gbface(iface,1),0) - m->gcoords(m->gbface(iface,0),0);
		y1 = m->gcoords(m->gbface(iface,1),1) - m->gcoords(m->gbface(iface,0),1);
		z1 = m->gcoords(m->gbface(iface,1),2) - m->gcoords(m->gbface(iface,0),2);
		x2 = m->gcoords(m->gbface(iface,2),0) - m->gcoords(m->gbface(iface,0),0);
		y2 = m->gcoords(m->gbface(iface,2),1) - m->gcoords(m->gbface(iface,0),1);
		z2 = m->gcoords(m->gbface(iface,2),2) - m->gcoords(m->gbface(iface,0),2);
		fnormals(iface,0) = y1*z2 - y2*z1;
		fnormals(iface,1) = -(x1*z2 - x2*z1);
		fnormals(iface,2) = x1*y2 - x2*y1;
		mag = sqrt(fnormals(iface,0)*fnormals(iface,0) + fnormals(iface,1)*fnormals(iface,1) + fnormals(iface,2)*fnormals(iface,2));
		fnormals(iface,0) /= mag;
		fnormals(iface,1) /= mag;
		fnormals(iface,2) /= mag;
	}
}

BoundaryReconstruction::~BoundaryReconstruction()
{
	delete [] V;
	delete [] X;
	delete [] F;
}


/// Implements WALF reconstruction according to Jiao and Wang's paper, ie, local fittings are calculated at each surface vertex
class VertexCenteredBoundaryReconstruction : public BoundaryReconstruction
{
	/// convert a point from local coord system of point ibpoin to the global xyz coord system
	std::vector<amc_real> xyz_from_uvw(const amc_int ibpoin, const std::vector<amc_real>& uvwpoint);

	/// convert a point from global coord system to the local uvw coord system of point ibpoin
	std::vector<amc_real> uvw_from_xyz(amc_int ibpoin, const std::vector<amc_real>& xyzpoint);

public:
	/// compute Vandermonde matrix and stencils for each point
	void preprocess();
	/// solve linear least-squares problem at each point
	void solve();
	
	std::vector<amc_real> getEdgePoint(const amc_real ratio, const amc_int edgenum);
	std::vector<amc_real> getFacePoint(const std::vector<amc_real> areacoords, const amc_int facenum);
};

std::vector<amc_real> VertexCenteredBoundaryReconstruction::xyz_from_uvw(const amc_int ibpoin, std::vector<amc_real>& uvwpoint)

}
