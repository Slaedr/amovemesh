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
	amat::Matrix<amc_real>* D;			///< unknowns (various derivatives)
	amat::Matrix<amc_real>* F;			///< RHS of least-squares problem, consisting of point heights or coordinates
	const amc_real s1;					///< Any number (to use for deciding the local coordinate frames)
	const amc_real s2;					///< Any number (to use for deciding the local coordinate frames)

public:
	BoundaryReconstruction(const UMesh* mesh, int deg);
	~BoundaryReconstruction();

	virtual std::vector<amc_real> getEdgePoint(const amc_real ratio, const amc_int edgenum) = 0;
	virtual std::vector<amc_real> getFacePoint(const std::vector<amc_real> areacoords, const amc_int facenum) = 0;
};

BoundaryReconstruction(const UMesh* mesh, int deg) 
	: m(mesh), degree(deg), s1(1.0), s2(2.0)
{
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


/// Implements WALF reconstruction according to Jiao and Wang's paper, ie, local fittings are calculated at each surface vertex
class VertexCenteredBoundaryReconstruction : public BoundaryReconstruction
{
	amat::Matrix<amc_real> pnormals;			///< normals at each point, calculated as average of face normals surrounding the point

	amat::Matrix<amc_real>* Q;					///< coordinate transformation (rotation) matrix for each point

	bool isalloc;

	int nders;									///< number of unknowns for the least-squares problems
	std::vector<int> m;							///< number of points in stencil for each surface point
	std::vector<std::vector<int>> stencil;		///< List of bpoint indices of points lying in the stencil of each surface point

	/// convert a point from local coord system of point ibpoin to the global xyz coord system
	void xyz_from_uvw(const amc_int ibpoin, const std::vector<amc_real>& uvwpoint, std::vector<amc_real>& xyzpoint);

	/// convert a point from global coord system to the local uvw coord system of point ibpoin
	void uvw_from_xyz(const amc_int ibpoin, const std::vector<amc_real>& xyzpoint, std::vector<amc_real>& uvwpoint);

public:
	VertexCenteredBoundaryReconstruction(const UMesh* mesh, int deg);
	~VertexCenteredBoundaryReconstruction();

	/// compute normal, rotation matrix and stencil for each point
	void preprocess();
	/// solve linear least-squares problem at each point
	void solve();
	
	std::vector<amc_real> getEdgePoint(const amc_real ratio, const amc_int edgenum);
	std::vector<amc_real> getFacePoint(const std::vector<amc_real> areacoords, const amc_int facenum);
};

VertexCenteredBoundaryReconstruction::VertexCenteredBoundaryReconstruction(const UMesh* mesh, int deg) : BoundaryReconstruction(mesh, deg)
{
	V = new amat::Matrix<amc_real>[m->gnbpoin()];
	D = new amat::Matrix<amc_real>[m->gnbpoin()];
	F = new amat::Matrix<amc_real>[m->gnbpoin()];
	Q = new amat::Matrix<amc_real>[m->gnbpoin()];
	m.resize(m->gnbpoin());
	for(i = 0; i < m->gnbpoin(); i++)
	{
		Q[i].setup(m->gndim(), m->gndim());
		
		if(degree == 2)
			nders = 5;
		else
			nders = 9;

		D[i].setup(nders,1);
	}
	stencil.resize(m->gnbpoin());
}

VertexCenteredBoundaryReconstruction::~VertexCenteredBoundaryReconstruction()
{
	delete [] V;
	delete [] D;
	delete [] F;
	delete [] Q;
}

void VertexCenteredBoundaryReconstruction::preprocess()
{
	pnormals.setup(m->gnbpoin(), m->gndim());
	int ipoin, iface, idim, face, numfaces, inode;

	amat::Matrix<amc_real>* pnormals = &(this->pnormals);
	
	// get point normals and rotation matrices
	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		numfaces = 0;
		for(idim = 0; idim < m->gndim(); idim++)
			(*pnormals)(ipoin,idim) = 0.0;
		for(iface = m->gbfsubp_p(ipoin); iface < m->gbfsubp_p(ipoin+1); iface++)
		{
			face = m->gbfsubp(iface);
			for(idim = 0; idim < m->gndim(); idim++)
				(*pnormals)(ipoin,idim) += fnormals.get(face,idim);
			numfaces++;
		}
		for(idim = 0; idim < m->gndim(); idim++)
		{
			(*pnormals)(ipoin,idim) /= numfaces;
			
			Q[ipoin](idim,2) = pnormals->get(ipoin,idim);
		}

		if(fabs(Q[ipoin](0,2)) > ZERO_TOL)
		{
			Q[ipoin](1,0) = s1;
			Q[ipoin](2,0) = s2;
			Q[ipoin](0,0) = (-s1*Q[ipoin](1,2)-s2*Q[ipoin](2,2))/Q[ipoin](0,2);
		}
		else if(fabs(Q[ipoin](1,2)) > ZERO_TOL)
		{
			Q[ipoin](0,0) = s1;
			Q[ipoin](2,0) = s2;
			Q[ipoin](1,0) = (-s1*Q[ipoin](0,2) - s2*Q[ipoin](2,2))/Q[ipoin](1,2);
		}
		else
		{
			Q[ipoin](0,0) = s1;
			Q[ipoin](1,0) = s2;
			Q[ipoin](2,0) = (-s1*Q[ipoin](0,2) - s2*Q[ipoin](1,2))/Q(2,2);
		}

		Q[ipoin](0,1) = Q[ipoin](1,2)*Q[ipoin](2,0) - Q[ipoin](2,2)*Q[ipoin](1,0);
		Q[ipoin](1,1) = -( Q[ipoin](0,2)*Q[ipoin](2,0) - Q[ipoin](2,2)*Q[ipoin](0,0) );
		Q[ipoin](2,1) = Q[ipoin](0,2)*Q[ipoin](1,0) - Q[ipoin](1,2)*Q[ipoin](0,0);
	}

	// compute reconstruction stencils of each point and store
	std::vector<int> fflags(m->gnface(), 0);
	std::vector<amc_int> sfaces;
	std::vector<amc_int> facepo;		// for storing local node number of ipoin in each neighboring face
	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		if(m->gnnofa() == 3)
		{
			if(degree == 2)
			{
				for(iface = m->gbfsubp_p(ipoin); iface < m->gbfsubp_p(ipoin+1); iface++)
				{
					face = m->gbfsubp(iface);
					for(inode = 0; inode != m->gnnofa(); inode++)
					{
						if(m->gbface(face,inode) != m->gbpoints(ipoin))
							stencil[ipoin].push_back(m->gbpointsinv(m->gbface(face,inode)));
						else
							facepo.push_back(inode);
					}
					sfaces.push_back(face);
				}
				// TODO: 1-ring points added, now add points for the 1.5-ring
			}
			else if(degree == 3)
			{
			}
		}
		else if(m->gnnofa() == 4)
		{
		}
	}
}

void VertexCenteredBoundaryReconstruction::xyz_from_uvw(const amc_int ibpoin, const std::vector<amc_real>& uvwpoint, std::vector<amc_real>& xyzpoint)
{
	// local coordinate directions are the columns of Q
	int i,j;
	for(i = 0; i < m->gndim(); i++)
	{
		xyzpoint[i] = m->gcoords(m->gbpoints(ibpoin), i);
		for(j = 0; j < m->gndim(); j++)
			xyzpoint[i] += Q[ibpoin].get(i,j)*uvwpoint[j];
	}
}

void VertexCenteredBoundaryReconstruction::uvw_from_xyz(const amc_int ibpoin, const std::vector<amc_real>& xyzpoint, std::vector<amc_real>& uvwpoint)
{
	int i,j;
	for(i = 0; i < m->gndim(); i++)
	{
		uvwpoint[i] = 0;
		for(j = 0; j < m->gndim(); j++)
			uvwpoint[i] += Q[ibpoin].get(j,i) * (xyzpoint[j] - m->gcoords(m->gbpoints(ibpoin),j));
	}
}

}
