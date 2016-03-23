/** \brief Implementation of the C^0 surface reconstruction of Jiao and Wang, "Reconstructing high-order surfaces for meshing".
 * \date March 14, 2016
 * \author Aditya Kashi
 */

#ifndef __AMESH3D_H
#include <amesh3d.hpp>
#endif

#define __AGEOMETRY3D_H

namespace amc {

int factorial(int x)
{
	if(x == 0) 
		return 1;
	else
		return x * factorial(x-1);
}


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

BoundaryReconstruction::BoundaryReconstruction(const UMesh* mesh, int deg) 
	: m(mesh), degree(deg), s1(1.0), s2(2.0)
{
	fnormals.setup(m->gnface(), m->gndim());

	// compute unit face normals of triangular faces
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

void BoundaryReconstruction::preprocess() { }
void BoundaryReconstruction::solve() { }


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
	void getEdgePoint(const amc_real ratio, const amc_int edgenum, std::vector<amc_real>& point) const;
	/// Returns coords of a point lying on the face 'facenum' and having area coordinates given by 'areacoords'
	void getFacePoint(const std::vector<amc_real>& areacoords, const amc_int facenum, std::vector<amc_real>& point) const;
};

VertexCenteredBoundaryReconstruction::VertexCenteredBoundaryReconstruction(const UMesh* mesh, int deg) : BoundaryReconstruction(mesh, deg)
{
	D = new amat::Matrix<amc_real>[m->gnbpoin()];
	Q = new amat::Matrix<amc_real>[m->gnbpoin()];
	mpo.resize(m->gnbpoin());
	for(int i = 0; i < m->gnbpoin(); i++)
	{
		Q[i].setup(m->gndim(), m->gndim());
		
		if(degree == 2)
			nders = 5;
		else
			nders = 9;

		D[i].setup(nders,1);
	}
	stencil = new std::vector<int>[m->gnbpoin()];
}

VertexCenteredBoundaryReconstruction::~VertexCenteredBoundaryReconstruction()
{
	delete [] D;
	delete [] Q;
	delete [] stencil;
}

void VertexCenteredBoundaryReconstruction::preprocess()
{
	pnormals.setup(m->gnbpoin(), m->gndim());
	int ipoin, iface, idim, face, numfaces, inode, i,j, jed;
	amc_real normmag;

	amat::Matrix<amc_real>* pnormals = &(this->pnormals);
	
	// get point normals and rotation matrices
	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		numfaces = 0;
		normmag = 0;
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
			normmag += pnormals->get(ipoin,idim)*pnormals->get(ipoin,idim);
		}
		normmag = sqrt(normmag);

		for(idim = 0; idim < m->gndim(); idim++)
		{
			(*pnormals)(ipoin,idim) /= normmag;				// normalize the normal vector
			Q[ipoin](idim,2) = pnormals->get(ipoin,idim);
		}

		normmag = 0;

		if(fabs(Q[ipoin](0,2)) > ZERO_TOL)
		{
			Q[ipoin](1,0) = s1;
			Q[ipoin](2,0) = s2;
			Q[ipoin](0,0) = (-s1*Q[ipoin](1,2)-s2*Q[ipoin](2,2))/Q[ipoin](0,2);
			normmag = sqrt(Q[ipoin].get(1,0)*Q[ipoin].get(1,0) + Q[ipoin].get(2,0)*Q[ipoin].get(2,0) + Q[ipoin].get(0,0)*Q[ipoin].get(0,0));
			Q[ipoin](1,0) /= normmag;
			Q[ipoin](2,0) /= normmag;
			Q[ipoin](0,0) /= normmag;
		}
		else if(fabs(Q[ipoin](1,2)) > ZERO_TOL)
		{
			Q[ipoin](0,0) = s1;
			Q[ipoin](2,0) = s2;
			Q[ipoin](1,0) = (-s1*Q[ipoin](0,2) - s2*Q[ipoin](2,2))/Q[ipoin](1,2);
			normmag = sqrt(Q[ipoin].get(1,0)*Q[ipoin].get(1,0) + Q[ipoin].get(2,0)*Q[ipoin].get(2,0) + Q[ipoin].get(0,0)*Q[ipoin].get(0,0));
			Q[ipoin](1,0) /= normmag;
			Q[ipoin](2,0) /= normmag;
			Q[ipoin](0,0) /= normmag;
		}
		else
		{
			Q[ipoin](0,0) = s1;
			Q[ipoin](1,0) = s2;
			Q[ipoin](2,0) = (-s1*Q[ipoin](0,2) - s2*Q[ipoin](1,2))/Q[ipoin](2,2);
			normmag = sqrt(Q[ipoin].get(1,0)*Q[ipoin].get(1,0) + Q[ipoin].get(2,0)*Q[ipoin].get(2,0) + Q[ipoin].get(0,0)*Q[ipoin].get(0,0));
			Q[ipoin](1,0) /= normmag;
			Q[ipoin](2,0) /= normmag;
			Q[ipoin](0,0) /= normmag;
		}

		Q[ipoin](0,1) = Q[ipoin](1,2)*Q[ipoin](2,0) - Q[ipoin](2,2)*Q[ipoin](1,0);
		Q[ipoin](1,1) = -( Q[ipoin](0,2)*Q[ipoin](2,0) - Q[ipoin](2,2)*Q[ipoin](0,0) );
		Q[ipoin](2,1) = Q[ipoin](0,2)*Q[ipoin](1,0) - Q[ipoin](1,2)*Q[ipoin](0,0);
	}

	// compute reconstruction stencils of each point and store
	std::vector<int> pflags(m->gnbpoin());
	std::vector<amc_int> sfaces;
	std::vector<amc_int> facepo;		// for storing local node number of ipoin in each surrounding face
	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		pflags.assign(m->gnbpoin(),0);
		pflags[ipoin] = 1;
		sfaces.clear();
		facepo.clear();

		if(m->gnnofa() == 3)
		{
			if(degree == 2)
			{
				for(iface = m->gbfsubp_p(ipoin); iface < m->gbfsubp_p(ipoin+1); iface++)
				{
					face = m->gbfsubp(iface);
					for(inode = 0; inode != m->gnnofa(); inode++)
					{
						if(pflags[m->gbpointsinv(m->gbface(face,inode))] != 1)	
							stencil[ipoin].push_back(m->gbpointsinv(m->gbface(face,inode)));

						if(m->gbpointsinv(m->gbface(face,inode)) == ipoin)
							facepo.push_back(inode);

						pflags[m->gbpointsinv(m->gbface(face,inode))] = 1;
					}
					sfaces.push_back(face);
				}
				// 1-ring points added, now add points for the 1.5-ring
				for(i = 0; i < sfaces.size(); i++)
				{
					jed = (facepo[i]+1) % m->gnnofa();										// get the edge opposite to ipoin
					face = m->gbfsubf(sfaces[i],jed);										// get the face adjoining that edge
					for(j = 0; j < m->gnnofa(); j++)										// add nodes of that face to stencil provided they have not already been added
						if(pflags[m->gbpointsinv(m->gbface(face,j))] != 1)
							stencil[ipoin].push_back(m->gbpointsinv(m->gbface(face,j)));
				}
			}
			else if(degree == 3)
			{
			}
		}
		else if(m->gnnofa() == 4)
		{
		}

		mpo[ipoin] = stencil[ipoin].size();
	}
}

void VertexCenteredBoundaryReconstruction::xyz_from_uvw(const amc_int ibpoin, const std::vector<amc_real>& uvwpoint, std::vector<amc_real>& xyzpoint) const
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

void VertexCenteredBoundaryReconstruction::uvw_from_xyz(const amc_int ibpoin, const std::vector<amc_real>& xyzpoint, std::vector<amc_real>& uvwpoint) const
{
	int i,j;
	for(i = 0; i < m->gndim(); i++)
	{
		uvwpoint[i] = 0;
		for(j = 0; j < m->gndim(); j++)
			uvwpoint[i] += Q[ibpoin].get(j,i) * (xyzpoint[j] - m->gcoords(m->gbpoints(ibpoin),j));
	}
}

void VertexCenteredBoundaryReconstruction::solve()
{
	const UMesh* m = this->m;
	std::vector<int>* stencil = this->stencil;
	amat::Matrix<amc_real>* Q = this->Q;
	amat::Matrix<amc_real>* pnormals = &(this->pnormals);
	amat::Matrix<amc_real>* D = this->D;
	std::vector<int>* mpo = &(this->mpo);
	int nders = this->nders;
	int degree = this->degree;

	int ipoin;

	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		//std::cout << "VertexCenteredBoundaryReconstruction: solve(): Point " << m->gbpoints(ipoin) << " : ";
		int isp, i, j, idim, k, l;
		amc_int pno;
		amc_real weight, wd;

		// assemble V and F
		amat::Matrix<amc_real> V(mpo->at(ipoin), nders);			// least-squares LHS, Vandermonde matrix
		amat::Matrix<amc_real> F(mpo->at(ipoin),1);					// least-squares RHS, height values of stencil points
		std::vector<amc_real> xyzp(m->gndim()), uvwp(m->gndim());

		for(isp = 0; isp < mpo->at(ipoin); isp++)
		{
			pno = stencil[ipoin][isp];
			for(idim = 0; idim < m->gndim(); idim++)
				xyzp[idim] = m->gcoords(m->gbpoints(pno),idim);
			uvw_from_xyz(ipoin, xyzp, uvwp);
			l = 0;
			for(i = 1; i <= degree; i++)
			{
				for(j = i, k = 0; j >= 0, k <= i; j--, k++)
				{
					V(isp,l) = pow(uvwp[0],j)*pow(uvwp[1],k)/factorial(j)*factorial(k);
					l++;
				}
			}

			// for debug
			if(l != nders) std::cout << "VertexCenteredBoundaryReconstruction: solve(): ! LHS computation is wrong!!" << std::endl;

			F(isp) = uvwp[2];
			
			// compute weights
			weight = 0; wd = 0;
			for(i = 0; i < m->gndim(); i++)
			{
				weight += pnormals->get(ipoin,i)*pnormals->get(pno,i);
				wd += (m->gcoords(m->gbpoints(ipoin),i) - m->gcoords(m->gbpoints(pno),i))*(m->gcoords(m->gbpoints(ipoin),i) - m->gcoords(m->gbpoints(pno),i));
			}
			if(weight < A_SMALL_NUMBER) weight = A_SMALL_NUMBER;
			wd = pow(sqrt(wd + A_SMALL_NUMBER),degree/2.0);
			weight = weight/wd;

			for(i = 0; i < nders; i++)
				V(isp,i) *= weight;
			F(isp) *= weight;
		}

		leastSquares_NE(V, F, D[ipoin]);
		//std::cout << std::endl;
	}
}

void VertexCenteredBoundaryReconstruction::getEdgePoint(const amc_real ratio, const amc_int edgenum, std::vector<amc_real>& point) const
{
	int ipoin, jpoin, ibp,jbp, idim;
	ipoin = m->gintedge(edgenum,0);
	jpoin = m->gintedge(edgenum,1);
	ibp = m->gbpointsinv(ipoin);
	jbp = m->gbpointsinv(jpoin);

	std::vector<amc_real> xyzp(m->gndim()), xyzq(m->gndim()), uvw0(m->gndim()), uvw1(m->gndim());

	for(idim = 0; idim < m->gndim(); idim++)
		xyzp[idim] = m->gcoords(ipoin,idim) + ratio*(m->gcoords(jpoin,idim) - m->gcoords(ipoin,idim));

	uvw_from_xyz(ibp,xyzp,uvw0);
	uvw_from_xyz(jbp,xyzp,uvw1);

	// evaluate 2D Taylor polynomial for each point
	amc_real h1 = 0, h2 = 0;
	int l = 0, i,j,k, fj, fk;
	for(i = 1; i <= degree; i++)
	{
		for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
		{
			fj = factorial(j);
			fk = factorial(k);
			h1 += pow(uvw0[0],j)*pow(uvw0[1],k)/fj*fk * D[ibp].get(l);
			h2 += pow(uvw1[0],j)*pow(uvw1[1],k)/fj*fk * D[jbp].get(l);
			l++;
		}
	}

	uvw0[2] = h1;
	uvw1[2] = h2;
	xyz_from_uvw(ibp,uvw0,xyzp);
	xyz_from_uvw(jbp,uvw1,xyzq);

	for(idim = 0; idim < m->gndim(); idim++)
		point[idim] = (1.0-ratio)*xyzp[idim] + ratio*xyzq[idim];
}
	
void VertexCenteredBoundaryReconstruction::getFacePoint(const std::vector<amc_real>& areacoords, const amc_int facenum, std::vector<amc_real>& point) const
{
	std::vector<int> spo(m->gnnofa()), sbpo(m->gnnofa());
	int i,idim;
	
	for(i = 0; i < m->gnnofa(); i++)
	{
		spo[i] = m->gbface(facenum,i);
		sbpo[i] = m->gbpointsinv(m->gbface(facenum,i));
	}

	std::vector<std::vector<amc_real>> xyzp(m->gnnofa()), uvwp(m->gnnofa());
	for(i = 0; i < m->gnnofa(); i++)
	{
		xyzp[i].resize(m->gndim());
		uvwp[i].resize(m->gndim());
	}
	xyzp[0].assign(m->gndim(),0.0);

	for(i = 0; i < m->gnnofa(); i++)
		for(idim = 0; idim < m->gndim(); idim++)
			xyzp[0][idim] += areacoords[i]*m->gcoords(spo[i],idim);

	// get local coordinates of the point in the local frames of the three vertices
	for(i = 0; i < m->gnnofa(); i++)
		uvw_from_xyz(sbpo[i],xyzp[0],uvwp[i]);
	
	std::vector<amc_real> height(m->gnnofa(),0);
	int l = 0,j,k,inofa;
	for(i = 1; i <= degree; i++)
	{
		for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
		{
			for(inofa = 0; inofa < m->gnnofa(); inofa++)
				height[inofa] += pow(uvwp[inofa][0],j)*pow(uvwp[inofa][1],k)/factorial(j)*factorial(k) * D[sbpo[inofa]].get(l);
			l++;
		}
	}

	for(inofa = 0; inofa < m->gnnofa(); inofa++)
		uvwp[inofa][2] = height[inofa];

	for(inofa = 0; inofa < m->gnnofa(); inofa++)
		xyz_from_uvw(sbpo[inofa], uvwp[inofa], xyzp[inofa]);

	point.assign(m->gnnofa(),0.0);
	for(inofa = 0; inofa < m->gnnofa(); inofa++)
		for(idim = 0; idim < m->gndim(); idim++)
			point[idim] += areacoords[inofa]*xyzp[inofa][idim];
}

}
