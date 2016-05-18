/** \file ageometry3d.cpp
 * \brief Implementation of the C^0 surface reconstruction of Jiao and Wang, "Reconstructing high-order surfaces for meshing".
 * 
 * See the corresponding header file for documentation of functions.
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


DiscontinuityDetection::DiscontinuityDetection(const UMesh* const mesh, const amat::Matrix<amc_real>* const fnormal, const double max_angle) : m(mesh), fnormals(fnormal), maxangle(max_angle)
{
	febedge.resize(m->gnbedge(),-1);
	febpoint.resize(m->gnbpoin(),-1);
}
	
void DiscontinuityDetection::detect_C1_discontinuities()
{
	int idim;
	amc_int ied, iface, jface, ibpoin, jbpoin;
	amc_real dotproduct;
	for(ied = 0; ied < m->gnbedge(); ied++)
	{
		iface = m->gintbedge(ied,0);
		jface = m->gintbedge(ied,1);
		ibpoin = m->gbpointsinv(m->gintbedge(ied,2));
		jbpoin = m->gbpointsinv(m->gintbedge(ied,3));
		dotproduct = 0;
		for(idim = 0; idim < 3; idim++)
			dotproduct += fnormals->get(iface,idim)*fnormals->get(jface,idim);
		if(dotproduct < cos(maxangle))
		{
			if(
			febedge[ied] = -2;
		}
	}
}

BoundaryReconstruction::BoundaryReconstruction(const UMesh* mesh, int deg, std::string stencil_type, int i_start) 
	: m(mesh), degree(deg), stencilType(stencil_type), s1(1.0), s2(2.0), istart(i_start)
{
	fnormals.setup(m->gnface(), m->gndim());
	std::cout << "BoundaryReconstruction: Stencil type is " << stencilType << std::endl;
	farea.resize(m->gnface());
	face_center.setup(m->gnface(),m->gndim());
	face_center.zeros();

	// compute unit face normals (by cross product) and areas (by Heron's formula) of triangular faces
	int iface, inode, idim;
	amc_real x1,y1,z1,x2,y2,z2,mag, a,b,c,s;
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
		
		a = b = c = 0;
		for(idim = 0; idim < m->gndim(); idim++)
		{
			a += (m->gcoords(m->gbface(iface,1),idim) - m->gcoords(m->gbface(iface,0),idim))*(m->gcoords(m->gbface(iface,1),idim) - m->gcoords(m->gbface(iface,0),idim));
			b += (m->gcoords(m->gbface(iface,2),idim) - m->gcoords(m->gbface(iface,1),idim))*(m->gcoords(m->gbface(iface,2),idim) - m->gcoords(m->gbface(iface,1),idim));
			c += (m->gcoords(m->gbface(iface,0),idim) - m->gcoords(m->gbface(iface,2),idim))*(m->gcoords(m->gbface(iface,0),idim) - m->gcoords(m->gbface(iface,2),idim));
		}
		a = sqrt(a); b = sqrt(b); c = sqrt(c);
		s = (a+b+c)*0.5;
		farea[iface] = sqrt(s*(s-a)*(s-b)*(s-c));
		
		// get face centers
		for(inode = 0; inode < m->gnnofa(); inode++)
			for(idim = 0; idim < m->gndim(); idim++)
				face_center(iface,idim) += m->gcoords(m->gbface(iface,inode),idim);

		for(idim = 0; idim < m->gndim(); idim++)
			face_center(iface,idim) /= m->gnnofa();

		// impose true normals
		/*mag = 0;
		for(idim=0; idim < m->gndim(); idim++)
		{
			fnormals(iface,idim) = face_center.get(iface,idim);
			mag += fnormals(iface,idim)*fnormals(iface,idim);
		}
		mag = sqrt(mag);
		for(idim = 0; idim < m->gndim(); idim++)
			fnormals(iface,idim) /= mag;*/
	}
}

void BoundaryReconstruction::preprocess() { }
void BoundaryReconstruction::solve() { }

// currently only for triangular surface mesh!
void BoundaryReconstruction::computePointNormalsInverseDistance()
{
	amc_int ipoin, ifa, iface;
	int idim;
	amc_real weight, weightsum, normmag, exponent = 1.0, a,b,c,s;
	pnormals.setup(m->gnbpoin(),m->gndim());
	pnormals.zeros();

	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		weightsum = 0;
		normmag = 0;
		for(ifa = m->gbfsubp_p(ipoin); ifa < m->gbfsubp_p(ipoin+1); ifa++)
		{
			iface = m->gbfsubp(ifa);
			weight = 0;
			for(idim = 0; idim < m->gndim(); idim++)
				weight += (face_center.get(iface,idim) - m->gcoords(m->gbpoints(ipoin),idim)) * (face_center.get(iface,idim) - m->gcoords(m->gbpoints(ipoin),idim));
			weight = 1.0/pow(sqrt(weight),exponent);
			weightsum += weight;
			for(idim = 0; idim < m->gndim(); idim++)
				pnormals(ipoin,idim) += fnormals(iface,idim)*weight;
		}

		for(idim = 0; idim < m->gndim(); idim++)
		{
			pnormals(ipoin,idim) /= weightsum;
			// also, for normalizing the normals
			normmag += pnormals.get(ipoin,idim)*pnormals.get(ipoin,idim);
		}
		normmag = sqrt(normmag);
		
		pnormals(ipoin,0) /= normmag;
		pnormals(ipoin,1) /= normmag;
		pnormals(ipoin,2) /= normmag;
	}
}

// currently only for triangular surface mesh!
void BoundaryReconstruction::computePointNormalsArea()
{
	amc_int ipoin, ifa, iface;
	int idim;
	amc_real weight, weightsum, normmag, exponent = 1.0;
	pnormals.setup(m->gnbpoin(),m->gndim());
	pnormals.zeros();

	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		weightsum = 0;
		normmag = 0;
		for(ifa = m->gbfsubp_p(ipoin); ifa < m->gbfsubp_p(ipoin+1); ifa++)
		{
			iface = m->gbfsubp(ifa);
			weight = 0;
			weight = farea[iface];
			weightsum += weight;
			for(idim = 0; idim < m->gndim(); idim++)
				pnormals(ipoin,idim) += fnormals(iface,idim)*weight;
		}

		for(idim = 0; idim < m->gndim(); idim++)
		{
			pnormals(ipoin,idim) /= weightsum;
			// also, for normalizing the normals
			normmag += pnormals.get(ipoin,idim)*pnormals.get(ipoin,idim);
		}
		normmag = sqrt(normmag);
		
		pnormals(ipoin,0) /= normmag;
		pnormals(ipoin,1) /= normmag;
		pnormals(ipoin,2) /= normmag;

		// set true normals for unit ball
		/*normmag = 0;
		for(idim = 0; idim < m->gndim(); idim++)
		{
			pnormals(ipoin,idim) = m->gcoords(m->gbpoints(ipoin),idim);
			normmag += pnormals(ipoin,idim)*pnormals(ipoin,idim);
		}
		normmag = sqrt(normmag);
		for(idim = 0; idim < m->gndim(); idim++)
			pnormals(ipoin,idim) /= normmag;*/
	}
}

VertexCenteredBoundaryReconstruction::VertexCenteredBoundaryReconstruction(const UMesh* mesh, int deg, std::string stencil_type, bool _safeguard, double norm_limit) 
	: safeguard(_safeguard), normlimit(norm_limit), BoundaryReconstruction(mesh, deg, stencil_type, 0)
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
			nders = 5+(1-istart);
		else
			nders = 9;

		D[i].setup(nders, m->gndim());
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
	int ipoin, jpoin, kpoin, iface, idim, face, inode, i,j,k, jed;
	amc_real normmag;

	amat::Matrix<amc_real>* pnormals = &(this->pnormals);
	
	// get point normals and rotation matrices
	
	computePointNormalsInverseDistance();

	for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
	{
		// true normals for a sphere centered at the origin
		/*for(idim = 0; idim < m->gndim(); idim++)
			normmag += m->gcoords(m->gbpoints(ipoin),idim)*m->gcoords(m->gbpoints(ipoin),idim);
		normmag = sqrt(normmag);
		for(idim = 0; idim < m->gndim(); idim++)
			(*pnormals)(ipoin,idim) = m->gcoords(m->gbpoints(ipoin),idim)/normmag;*/
		
		// next, get the rotation matrix to for the local coord system at this vertex
		
		for(idim = 0; idim < m->gndim(); idim++)
			Q[ipoin](idim,2) = pnormals->get(ipoin,idim);

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
	if(stencilType == "half")
	{
		for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
		{
			pflags.assign(m->gnbpoin(),0);
			pflags[ipoin] = 1;
			sfaces.clear();
			facepo.clear();

			// NOTE: adding the point itself in its stencil
			stencil[ipoin].push_back(ipoin);

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
	else if (stencilType == "full")
	{
		// currently only for a 2-ring (refer Jiao and Wang) in a triangular surface mesh
		// but can be generalized without much difficulty to n-ring stencils, by putting the loop over surpoints in a n-loop.
		std::vector<int> surpoints;
		amc_int fa;
		for(ipoin = 0; ipoin < m->gnbpoin(); ipoin++)
		{
			if(m->gnnofa() == 3)
			{
				surpoints.clear();
				pflags.assign(m->gnbpoin(),0);
				pflags[ipoin] = 1;

				// 1-ring
				for(j = m->gbpsubp_p(ipoin); j < m->gbpsubp_p(ipoin+1); j++)
				{
					jpoin = m->gbpsubp(j);
					if(pflags[jpoin] == 0)
					{
						pflags[jpoin] = 1;
						stencil[ipoin].push_back(jpoin);
						surpoints.push_back(jpoin);
					}
				}
				
				// 2-ring
				for(j = 0; j != surpoints.size(); j++)
				{
					jpoin = surpoints[j];
					for(k = m->gbpsubp_p(jpoin); k < m->gbpsubp_p(jpoin+1); k++)
					{
						kpoin = m->gbpsubp(k);
						if(pflags[kpoin] != 1) 
						{
							stencil[ipoin].push_back(kpoin);
							pflags[kpoin] = 1;
						}
					}
				}
				mpo[ipoin] = stencil[ipoin].size();
			}
			else
			{
				std::cout << "VertexCenteredBoundaryReconstruction: Not implemented for this type of face!" << std::endl;
			}
		}
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

/** \todo After implementation of templated Matrix class, use column-major arrays for storing V, F and D
 */
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
	const int istart = this->istart;

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
		amat::Matrix<amc_real> V(mp, nders);				// least-squares LHS, Vandermonde matrix
		amat::Matrix<amc_real> F(mp, m->gndim());			// least-squares RHS, coords of stencil points
		std::vector<amc_real> xyzp(m->gndim()), uvwp(m->gndim());
		std::vector<amc_real> weightsn(mp);					// numerators of row-weights for weighted least-squares
		std::vector<amc_real> weightsd(mp);					// denominators of row-weights

		for(isp = 0; isp < mp; isp++)
		{
			pno = stencil[ipoin][isp];
			for(idim = 0; idim < m->gndim(); idim++)
				xyzp[idim] = m->gcoords(m->gbpoints(pno),idim);
			uvw_from_xyz(ipoin, xyzp, uvwp);

			l = 0;
			for(i = istart; i <= degree; i++)
			{
				for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
				{
					V(isp,l) = pow(uvwp[0],j)*pow(uvwp[1],k)/factorial(j)*factorial(k);
					l++;
				}
			}

			// for debug
			if(l != nders) std::cout << "VertexCenteredBoundaryReconstruction: solve(): ! LHS computation is wrong!!" << std::endl;

			for(idim = 0; idim < m->gndim(); idim++)
				F(isp,idim) = xyzp[idim];
			
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
			for(idim = 0; idim < m->gndim(); idim++)
				F(isp,idim) *= weightsn[isp]/weightsd[isp];
		}

		//leastSquares_NE(V, F, D[ipoin]);
		leastSquares_QR(V, F, D[ipoin]);
		//leastSquares_SVD(V, F, D[ipoin]);

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

		solve_QR(v,V,F,D[ipoin]);
		for(i = 0; i < nders; i++)
			for(idim = 0; idim < m->gndim(); idim++)
				D[ipoin](idim,i) *= scale[i];	*/
	}

	delete [] v;
}

void VertexCenteredBoundaryReconstruction::getEdgePoint(const amc_real ratio, const amc_int edgenum, std::vector<amc_real>& point) const
{
	int ipoin, jpoin, ibp,jbp, idim;
	ipoin = m->gintbedge(edgenum,2);
	jpoin = m->gintbedge(edgenum,3);
	ibp = m->gbpointsinv(ipoin);
	jbp = m->gbpointsinv(jpoin);

	std::vector<amc_real> xyzp(m->gndim()), xyzq(m->gndim()), uvw0(m->gndim()), uvw1(m->gndim());

	for(idim = 0; idim < m->gndim(); idim++)
		xyzp[idim] = m->gcoords(ipoin,idim) + ratio*(m->gcoords(jpoin,idim) - m->gcoords(ipoin,idim));

	uvw_from_xyz(ibp,xyzp,uvw0);
	uvw_from_xyz(jbp,xyzp,uvw1);

	// evaluate 2D Taylor polynomial for each point
	xyzp.assign(m->gndim(),0); xyzq.assign(m->gndim(),0);
	int l, i,j,k, fj, fk;
	for(idim = 0; idim < m->gndim(); idim++)
	{
		l = 0;
		for(i = istart; i <= degree; i++)
		{
			for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
			{
				fj = factorial(j);
				fk = factorial(k);
				if(rec_order[ibp] >= i)
					xyzp[idim] += pow(uvw0[0],j)*pow(uvw0[1],k)/fj*fk * D[ibp].get(l,idim);
				if(rec_order[jbp] >= i)
					xyzq[idim] += pow(uvw1[0],j)*pow(uvw1[1],k)/fj*fk * D[jbp].get(l,idim);
				l++;
			}
		}
	}

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
		xyzp[i].assign(m->gndim(),0);
		uvwp[i].resize(m->gndim());
	}

	for(i = 0; i < m->gnnofa(); i++)
		for(idim = 0; idim < m->gndim(); idim++)
			xyzp[0][idim] += areacoords[i]*m->gcoords(spo[i],idim);

	// get local coordinates of the point in the local frames of the three vertices
	for(i = 0; i < m->gnnofa(); i++)
		uvw_from_xyz(sbpo[i],xyzp[0],uvwp[i]);
	
	xyzp[0].assign(m->gndim(),0);
	int l,j,k,inofa, fj,fk;
	for(idim = 0; idim < m->gndim(); idim++)
	{
		l = 0;
		for(i = istart; i <= degree; i++)
		{
			for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
			{
				fj = factorial(j);
				fk = factorial(k);
				for(inofa = 0; inofa < m->gnnofa(); inofa++)
					if(rec_order[sbpo[inofa]] >= i)
						xyzp[inofa][idim] += pow(uvwp[inofa][0],j)*pow(uvwp[inofa][1],k)/fj*fk * D[sbpo[inofa]].get(l,idim);
				l++;
			}
		}
	}

	point.assign(m->gnnofa(),0.0);
	for(inofa = 0; inofa < m->gnnofa(); inofa++)
		for(idim = 0; idim < m->gndim(); idim++)
			point[idim] += areacoords[inofa]*xyzp[inofa][idim];
}


// Implementation of face-centered reconstruction follows

FaceCenteredBoundaryReconstruction::FaceCenteredBoundaryReconstruction(const UMesh* mesh, int deg, std::string stencil_type, bool _safeguard, double norm_limit) 
	: safeguard(_safeguard), normlimit(norm_limit), BoundaryReconstruction(mesh, deg, stencil_type, 0)
{
	std::cout << "FaceCenteredBoundaryReconstruction: Computing with safeguard - " << safeguard << std::endl;
	D = new amat::Matrix<amc_real>[m->gnface()];
	Q = new amat::Matrix<amc_real>[m->gnface()];
	mpo.resize(m->gnface());
	rec_order.resize(m->gnface(), degree);
	for(int i = 0; i < m->gnface(); i++)
	{
		Q[i].setup(m->gndim(), m->gndim());
		
		if(degree == 2)
			nders = 6;
		else
			nders = 8;

		D[i].setup(nders,1);
	}
	std::cout << "FaceCenteredBoundaryReconstruction: Number of unknowns per face = " << nders << std::endl;
	stencil = new std::vector<int>[m->gnface()];
}

FaceCenteredBoundaryReconstruction::~FaceCenteredBoundaryReconstruction()
{
	delete [] D;
	delete [] Q;
	delete [] stencil;
}

/** For face-centered reconstruction, the stencil for a b-face consists of all b-points contained in all vertex-neighbors of the b-face.
 * This is for a P2 reconstruction; we may need a second layer for P3 (and we don't care beyond that).
 */
void FaceCenteredBoundaryReconstruction::preprocess()
{
	amc_int ipoin, iface, jface;
	int idim, face, numfaces, inode, i,j, jed;
	amc_real normmag;
	amat::Matrix<amc_real>* fnormals = &(this->fnormals);
	
	computePointNormalsInverseDistance();

	// get rotation matrices
	for(iface = 0; iface < m->gnface(); iface++)
	{
		for(idim = 0; idim < m->gndim(); idim++)
			Q[iface](idim,2) = fnormals->get(iface,idim);

		normmag = 0;

		if(fabs(Q[iface](0,2)) > ZERO_TOL)
		{
			Q[iface](1,0) = s1;
			Q[iface](2,0) = s2;
			Q[iface](0,0) = (-s1*Q[iface](1,2)-s2*Q[iface](2,2))/Q[iface](0,2);
			normmag = sqrt(Q[iface].get(1,0)*Q[iface].get(1,0) + Q[iface].get(2,0)*Q[iface].get(2,0) + Q[iface].get(0,0)*Q[iface].get(0,0));
			Q[iface](1,0) /= normmag;
			Q[iface](2,0) /= normmag;
			Q[iface](0,0) /= normmag;
		}
		else if(fabs(Q[iface](1,2)) > ZERO_TOL)
		{
			Q[iface](0,0) = s1;
			Q[iface](2,0) = s2;
			Q[iface](1,0) = (-s1*Q[iface](0,2) - s2*Q[iface](2,2))/Q[iface](1,2);
			normmag = sqrt(Q[iface].get(1,0)*Q[iface].get(1,0) + Q[iface].get(2,0)*Q[iface].get(2,0) + Q[iface].get(0,0)*Q[iface].get(0,0));
			Q[iface](1,0) /= normmag;
			Q[iface](2,0) /= normmag;
			Q[iface](0,0) /= normmag;
		}
		else
		{
			Q[iface](0,0) = s1;
			Q[iface](1,0) = s2;
			Q[iface](2,0) = (-s1*Q[iface](0,2) - s2*Q[iface](1,2))/Q[iface](2,2);
			normmag = sqrt(Q[iface].get(1,0)*Q[iface].get(1,0) + Q[iface].get(2,0)*Q[iface].get(2,0) + Q[iface].get(0,0)*Q[iface].get(0,0));
			Q[iface](1,0) /= normmag;
			Q[iface](2,0) /= normmag;
			Q[iface](0,0) /= normmag;
		}

		Q[iface](0,1) = Q[iface](1,2)*Q[iface](2,0) - Q[iface](2,2)*Q[iface](1,0);
		Q[iface](1,1) = -( Q[iface](0,2)*Q[iface](2,0) - Q[iface](2,2)*Q[iface](0,0) );
		Q[iface](2,1) = Q[iface](0,2)*Q[iface](1,0) - Q[iface](1,2)*Q[iface](0,0);
	}

	// compute reconstruction stencils of each point and store
	std::vector<int> pflags(m->gnbpoin());
	std::vector<amc_int> sfaces;
	int poin, jpoin;

	for(iface = 0; iface < m->gnface(); iface++)
	{
		pflags.assign(m->gnbpoin(),0);

		if(m->gnnofa() == 3)
		{
			// triangular face: for each point of the face, add the point and get surrounding boundary points
			for(inode = 0; inode < m->gnnofa(); inode++)
			{
				poin = m->gbpointsinv(m->gbface(iface,inode));
				if(pflags[poin] == 0)
				{
					stencil[iface].push_back(poin);
					pflags[poin] = 1;
				}
				for(j = m->gbpsubp_p(poin); j < m->gbpsubp_p(poin+1); j++)
				{
					jpoin = m->gbpsubp(j);
					if(pflags[jpoin] == 0)
					{
						pflags[jpoin] = 1;
						stencil[iface].push_back(jpoin);
					}
				}
			}
		}
		else
		{
			// \todo TODO: implement stencil for quad faces using only points of face-neighbors
			std::cout << "FaceCenteredBoundaryReconstruction: preprocess(): ! Not implemented for quad faces yet!" << std::endl;
		}
		mpo[iface] = stencil[iface].size();
	}
}

void FaceCenteredBoundaryReconstruction::xyz_from_uvw(const amc_int iface, const std::vector<amc_real>& uvwpoint, std::vector<amc_real>& xyzpoint) const
{
	// local coordinate directions are the columns of Q
	int i,j;
	for(i = 0; i < m->gndim(); i++)
	{
		xyzpoint[i] = face_center.get(iface, i);
		for(j = 0; j < m->gndim(); j++)
			xyzpoint[i] += Q[iface].get(i,j)*uvwpoint[j];
	}
}

void FaceCenteredBoundaryReconstruction::uvw_from_xyz(const amc_int iface, const std::vector<amc_real>& xyzpoint, std::vector<amc_real>& uvwpoint) const
{
	int i,j;
	for(i = 0; i < m->gndim(); i++)
	{
		uvwpoint[i] = 0;
		for(j = 0; j < m->gndim(); j++)
			uvwpoint[i] += Q[iface].get(j,i) * (xyzpoint[j] - face_center.get(iface,j));
	}
}

void FaceCenteredBoundaryReconstruction::solve()
{
	std::cout << "FaceCenteredBoundaryReconstruction: solve(): Computing slopes, curvatures etc at each face" << std::endl;
	const UMesh* m = this->m;
	std::vector<int>* stencil = this->stencil;
	amat::Matrix<amc_real>* Q = this->Q;
	amat::Matrix<amc_real>* fnormals = &(this->fnormals);
	amat::Matrix<amc_real>* pnormals = &(this->pnormals);
	amat::Matrix<amc_real>* face_center = &(this->face_center);
	amat::Matrix<amc_real>* D = this->D;
	std::vector<int>* mpo = &(this->mpo);
	int nders = this->nders;
	int degree = this->degree;

	amc_int iface;
	amc_real norm1;
	std::vector<amc_real>* v = new std::vector<amc_real>[nders];

	for(iface = 0; iface < m->gnface(); iface++)
	{
		//std::cout << "FaceCenteredBoundaryReconstruction: solve(): Face " << iface << " : ";
		int isp, i, j, idim, k, l, mp;
		mp = mpo->at(iface);
		amc_int pno;
		amc_real weight, wd = 0, csum;
		//std::cout << "FaceCenteredBoundaryReconstruction: solve(): Face " << iface << ": m = " << mp << ", n = " << nders << std::endl;

		// assemble V and F
		amat::Matrix<amc_real> V(mp, nders);			// least-squares LHS, Vandermonde matrix
		amat::Matrix<amc_real> F(mp,1);					// least-squares RHS, height values of stencil points
		std::vector<amc_real> xyzp(m->gndim()), uvwp(m->gndim());
		std::vector<amc_real> weightsn(mp);				// numerators of row-weights for weighted least-squares
		std::vector<amc_real> weightsd(mp);				// denominators of row-weights

		for(isp = 0; isp < mp; isp++)
		{
			pno = stencil[iface][isp];
			
			//std::cout << " " << m->gbpoints(pno);

			for(idim = 0; idim < m->gndim(); idim++)
				xyzp[idim] = m->gcoords(m->gbpointsinv(pno),idim);
			
			uvw_from_xyz(iface, xyzp, uvwp);

			l = 0;
			for(i = 0; i <= degree; i++)
			{
				// disregard slope terms
				//if(i == 1) continue;

				for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
				{
					V(isp,l) = pow(uvwp[0],j)*pow(uvwp[1],k)/factorial(j)*factorial(k);
					l++;
				}
			}

			// for debug
			if(l != nders) std::cout << "FaceCenteredBoundaryReconstruction: solve(): ! LHS computation is wrong!!" << std::endl;

			F(isp) = uvwp[2];
			
			// compute weights
			weightsn[isp] = 0; weightsd[isp] = 0;
			for(i = 0; i < m->gndim(); i++)
			{
				weightsn[isp] += fnormals->get(iface,i)*pnormals->get(pno,i);
				//weightsd[isp] += (m->gcoords(m->gbpoints(ipoin),i) - m->gcoords(m->gbpoints(pno),i))*(m->gcoords(m->gbpoints(ipoin),i) - m->gcoords(m->gbpoints(pno),i));
				weightsd[isp] += uvwp[i]*uvwp[i];
			}
			if(weightsn[isp] < ZERO_TOL) weightsn[isp] = ZERO_TOL;
			//wd = pow(sqrt(wd + A_SMALL_NUMBER),degree/2.0);
			//weight = weight/wd;
			wd += weightsd[isp];
		}
		//std::cout << std::endl;
	
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

		//leastSquares_NE(V, F, D[iface]);
		leastSquares_QR(V, F, D[iface]);
		//leastSquares_SVD(V, F, D[iface]);

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
			std::cout << "FaceCenteredBoundaryReconstruction: solve(): Point " << ipoin << " demoted!" << std::endl;
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

void FaceCenteredBoundaryReconstruction::getEdgePoint(const amc_real ratio, const amc_int edgenum, std::vector<amc_real>& point) const
{
	int ipoin, jpoin, ifa, jfa, idim;
	ipoin = m->gintbedge(edgenum,2);
	jpoin = m->gintbedge(edgenum,3);
	ifa = m->gintbedge(edgenum,0);
	jfa = m->gintbedge(edgenum,1);

	amc_real disti, distj;

	//std::cout << "Getting edge point for edge " << edgenum << " with points " << ipoin << ", " << jpoin << " and faces " << ifa << ", " << jfa << std::endl;

	std::vector<amc_real> xyzp(m->gndim()), xyzq(m->gndim()), uvw0(m->gndim()), uvw1(m->gndim());

	for(idim = 0; idim < m->gndim(); idim++)
		xyzp[idim] = m->gcoords(ipoin,idim) + ratio*(m->gcoords(jpoin,idim) - m->gcoords(ipoin,idim));

	disti = sqrt((xyzp[0]-face_center.get(ifa,0))*(xyzp[0]-face_center.get(ifa,0))
			+ (xyzp[1]-face_center.get(ifa,1))*(xyzp[1]-face_center.get(ifa,1)) + (xyzp[2]-face_center.get(ifa,2))*(xyzp[2]-face_center.get(ifa,2)));
	distj = sqrt((xyzp[0]-face_center.get(jfa,0))*(xyzp[0]-face_center.get(jfa,0))
			+ (xyzp[1]-face_center.get(jfa,1))*(xyzp[1]-face_center.get(jfa,1)) + (xyzp[2]-face_center.get(jfa,2))*(xyzp[2]-face_center.get(jfa,2)));

	uvw_from_xyz(ifa,xyzp,uvw0);
	uvw_from_xyz(jfa,xyzp,uvw1);

	// evaluate 2D Taylor polynomial for each point
	amc_real h1 = 0, h2 = 0;
	int l = 0, i,j,k, fj, fk;
	for(i = 0; i <= degree; i++)
	{
		//if(i == 1) continue;
		for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
		{
			fj = factorial(j);
			fk = factorial(k);
			//if(rec_order[ifa] >= i)
				h1 += pow(uvw0[0],j)*pow(uvw0[1],k)/fj*fk * D[ifa].get(l);
			//if(rec_order[jfa] >= i)
				h2 += pow(uvw1[0],j)*pow(uvw1[1],k)/fj*fk * D[jfa].get(l);
			l++;
		}
	}

	uvw0[2] = h1;
	uvw1[2] = h2;
	xyz_from_uvw(ifa,uvw0,xyzp);
	xyz_from_uvw(jfa,uvw1,xyzq);

	/*for(idim = 0; idim < m->gndim(); idim++)
		point[idim] = (xyzp[idim] + xyzq[idim])/2.0;*/
	/*for(idim = 0; idim < m->gndim(); idim++)
		point[idim] = (xyzp[idim]*farea[ifa] + xyzq[idim]*farea[jfa]) / (farea[ifa]+farea[jfa]);*/
	for(idim = 0; idim < m->gndim(); idim++)
		point[idim] = (xyzp[idim]*disti + xyzq[idim]*distj) / (disti+distj);
}
	
void FaceCenteredBoundaryReconstruction::getFacePoint(const std::vector<amc_real>& areacoords, const amc_int facenum, std::vector<amc_real>& point) const
{
	std::vector<int> spo(m->gnnofa()), sbpo(m->gnnofa());
	int i,idim;
	
	for(i = 0; i < m->gnnofa(); i++)
	{
		spo[i] = m->gbface(facenum,i);
	}

	std::vector<amc_real> xyzp, uvwp(m->gndim());
	xyzp.assign(m->gndim(),0.0);

	for(i = 0; i < m->gnnofa(); i++)
		for(idim = 0; idim < m->gndim(); idim++)
			xyzp[idim] += areacoords[i]*m->gcoords(spo[i],idim);

	// get local coordinates of the point in the local frames of the three vertices
	uvw_from_xyz(facenum, xyzp, uvwp);
	
	amc_real height = 0;
	int l = 0,j,k, fj,fk;
	for(i = 0; i <= degree; i++)
	{
		if(i == 1) continue;
		for(j = i, k = 0; j >= 0 && k <= i; j--, k++)
		{
			fj = factorial(j);
			fk = factorial(k);
			if(rec_order[facenum] >= i)
				height += pow(uvwp[0],j)*pow(uvwp[1],k)/fj*fk * D[facenum].get(l);
			l++;
		}
	}

	uvwp[2] = height;
	xyz_from_uvw(facenum, uvwp, point);
}

}
