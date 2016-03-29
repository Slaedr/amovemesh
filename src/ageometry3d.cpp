/** \brief Implementation of the C^0 surface reconstruction of Jiao and Wang, "Reconstructing high-order surfaces for meshing".
 * \date March 14, 2016
 * \author Aditya Kashi
 */

#include "ageometry3d.hpp"

namespace amc {

inline int factorial(int x)
{
	if(x == 0) 
		return 1;
	else
		return x * factorial(x-1);
}

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


VertexCenteredBoundaryReconstruction::VertexCenteredBoundaryReconstruction(const UMesh* mesh, int deg, bool _safeguard, double norm_limit) 
	: safeguard(_safeguard), normlimit(norm_limit), BoundaryReconstruction(mesh, deg)
{
	std::cout << "VertexCenteredBoundaryReconstruction: Computing with safeguard - " << safeguard << std::endl;
	D = new amat::Matrix<amc_real>[m->gnbpoin()];
	Q = new amat::Matrix<amc_real>[m->gnbpoin()];
	mpo.resize(m->gnbpoin());
	rec_order.resize(m->gnbpoin(), degree);
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
	std::cout << "VertexCenteredBoundaryReconstruction: solve(): Computing slopes, curvatures etc at each point" << std::endl;
	const UMesh* m = this->m;
	std::vector<int>* stencil = this->stencil;
	amat::Matrix<amc_real>* Q = this->Q;
	amat::Matrix<amc_real>* pnormals = &(this->pnormals);
	amat::Matrix<amc_real>* D = this->D;
	std::vector<int>* mpo = &(this->mpo);
	int nders = this->nders;
	int degree = this->degree;

	int ipoin;
	amc_real norm1;
	std::vector<amc_real>* v = new std::vector<amc_real>[nders];

	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		//std::cout << "VertexCenteredBoundaryReconstruction: solve(): Point " << m->gbpoints(ipoin) << " : ";
		int isp, i, j, idim, k, l, mp;
		mp = mpo->at(ipoin);
		amc_int pno;
		amc_real weight, wd = 0, csum;

		// assemble V and F
		amat::Matrix<amc_real> V(mp, nders);			// least-squares LHS, Vandermonde matrix
		amat::Matrix<amc_real> F(mp,1);					// least-squares RHS, height values of stencil points
		std::vector<amc_real> xyzp(m->gndim()), uvwp(m->gndim());
		std::vector<amc_real> weightsn(mp);				// numerators of row-weights for weighted least-squares
		std::vector<amc_real> weightsd(mp);				// denominators of row-weights

		for(isp = 0; isp < mp; isp++)
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
			weightsn[isp] = 0; weightsd[isp] = 0;
			for(i = 0; i < m->gndim(); i++)
			{
				weightsn[isp] += pnormals->get(ipoin,i)*pnormals->get(pno,i);
				//weightsd[isp] += (m->gcoords(m->gbpoints(ipoin),i) - m->gcoords(m->gbpoints(pno),i))*(m->gcoords(m->gbpoints(ipoin),i) - m->gcoords(m->gbpoints(pno),i));
				weightsd[isp] += uvwp[i]*uvwp[i];
			}
			if(weightsn[isp] < ZERO_TOL) weightsn[isp] = ZERO_TOL;
			//wd = pow(sqrt(wd + A_SMALL_NUMBER),degree/2.0);
			//weight = weight/wd;
			wd += weightsd[isp];
		}
	
		wd = wd / (100.0*mp);
		for(isp = 0; isp < mp; isp++)
		{
			weightsd[isp] += wd;
			weightsd[isp] = pow( sqrt(weightsd[isp]), degree/2.0 );
		}
		for(isp = 0; isp < mp; isp++)
		{
			for(i = 0; i < nders; i++)
				V(isp,i) *= weightsn[isp]/weightsd[isp];
			F(isp) *= weightsn[isp]/weightsd[isp];
		}

		//leastSquares_NE(V, F, D[ipoin]);
		//leastSquares_QR(V, F, D[ipoin]);
		leastSquares_SVD(V, F, D[ipoin]);

		/*std::vector<amc_real> scale(nders);		// for scaling the column of A
		
		// get norms of column-vectors of A and scale them
		for(j = 0; j < nders; j++)
		{
			csum = 0;
			for(i = 0; i < mp; i++)
				csum += V.get(i,j)*V.get(i,j);
			scale[j] = 1.0/sqrt(csum);
			for(i = 0; i < mp; i++)
				V(i,j) *= scale[j];
		}

		for(i = 0; i < nders; i++)
			v[i].resize(mp-i);

		// get QR decomposition (R is stored in V, Q is determined by v)
		qr(V,v);

		norm1 = V.matrixNorm_1();

		// if safeguard is on and the norm exceeds some limit, reduce the degree of reconstruction from 2 to 1 for this vertex
		if(safeguard && norm1 > normlimit)
		{
			std::cout << "VertexCenteredBoundaryReconstruction: solve(): Point " << ipoin << " demoted!" << std::endl;
			int n = nders-3;
			amat::Matrix<amc_real> R(mp-3,n); R.zeros();
			amat::Matrix<amc_real> b(mp-3,1);
			for(i = 0; i < mp-3; i++)
			{
				for(j = 0; j < n; j++)
					R(i,j) = V.get(i,j);
				b(i) = F.get(i);
			}
			solve_QR(v,R,b,D[ipoin]);
			for(i = 0; i < n; i++)
				D[ipoin](i) *= scale[i];

			rec_order[ipoin] = 1;
		}
		else
		{
			solve_QR(v,V,F,D[ipoin]);
			for(i = 0; i < nders; i++)
				D[ipoin](i) *= scale[i];
		}*/
	}

	delete [] v;
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
			if(rec_order[ibp] >= i)
				h1 += pow(uvw0[0],j)*pow(uvw0[1],k)/fj*fk * D[ibp].get(l);
			if(rec_order[jbp] >= i)
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
	int l = 0,j,k,inofa, fj,fk;
	for(i = 1; i <= degree; i++)
	{
		for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
		{
			fj = factorial(j);
			fk = factorial(k);
			for(inofa = 0; inofa < m->gnnofa(); inofa++)
				if(rec_order[sbpo[inofa]] >= i)
					height[inofa] += pow(uvwp[inofa][0],j)*pow(uvwp[inofa][1],k)/fj*fk * D[sbpo[inofa]].get(l);
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
