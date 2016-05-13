/** \brief Code to construct curved mesh using a quadratic mesh that has its boundary curved, and therefore probably has invalid elements.
 * \author Aditya Kashi
 * \date March 3, 2016
 */

#ifndef __AMESH3D_H
#include <amesh3d.hpp>
#endif

#ifndef __ARBF_H
#include <arbf.hpp>
#endif

#define __ACURVEDMESHGEN_H

namespace amc {

/// Class to generate a fully curved mesh, using as input a high-order mesh with boundaries already curved
/** Usually, we assume that displacements of all boundaries have already been imposed, and that only displacements of interior nodes need to be computed by RBF.
 * However, sometimes we may need to apply RBF on a boundary, some of whose end edges have been curved according to the geometry. 
 * For this situation, set rbf_boundaries in CurvedMeshGen::CurvedMeshGen().
 */
class CurvedMeshGen
{
	UMesh* m;							///< the mesh to curve
	RBFmove* move;						///< mesh movement context
	amat::Matrix<amc_real> inpoints;	///< interior points of the mesh
	amat::Matrix<amc_real> bounpoints;	///< boundary points
	amat::Matrix<amc_real> boundisps;	///< boundary displacements
	std::vector<int> curvemarkers;		///< marker numbers of boundaries that need to be curved
	std::vector<int> rbfboundaries;		///< Marker numbers of boundaries in which RBF needs to be applied to curve the edges
	std::vector<int> torec;				///< contains 1 if the corresponding boundary node 
	amc_int nbpoin;						///< number of boundary points
	amc_int ninpoin;					///< number of interior points
	amat::Matrix<int> rbfbpoin;			///< Stores 1 if the corresponding node is part of a boundary which needs to be processed by RBF
	
public:

	/// constructor
	/** \param[in|out] mesh is the mesh to curve
	 * \param[in] rbf_boundaries contains physical-entity (marker) numbers of boundaries which need processing by RBF
	 * \param[in] choice corresponding to the mesh method method - for RBF, this is the type of RBF to use
	 * \param[in] param1 parameter required by mesh movement method - for RBF, this is the support radius
	 * \param[in] tol the tolerance to use in linear solvers, for instance
	 * \param[in] maxiter maximum iterations for linear solver
	 * \param[in] solver is a string describing the linear solver to use ("CG", "LU", "BICGSTAB")
	 */
	CurvedMeshGen(UMesh* mesh, const std::vector<int> rbf_boundaries, const int choice, const double param1, const double tol, const int maxiter, const std::string solver);

	~CurvedMeshGen();

	void generateCurvedMesh();
};

CurvedMeshGen::CurvedMeshGen(UMesh* mesh, const std::vector<int> rbf_boundaries, const int choice, const double param1, const double tol, const int maxiter, const std::string solver)
{
	m = mesh;
	
	amc_int ipoin, iface, inode, jnode, j, k, l, nexbpoin = 0;
	int idim;
	bool inrflag;

	rbfbpoin.setup(m->gnpoin(),1);
	rbfbpoin.ones();
	for(int i = 0; i < m->gnpoin(); i++)
		if(m->gflag_bpoin(i) == 0)
			rbfbpoin(i) = 0;

	amat::Matrix<amc_real>* midpoints = new amat::Matrix<amc_real>[m->gnface()];
	for(int i = 0; i < m->gnface(); i++)
	{
		if(m->gnnofa()==6)
			midpoints[i].setup(3,m->gndim());
		else if(m->gnnofa()==9)
			midpoints[i].setup(5,m->gndim());
		else if(m->gnnofa()==8)
			midpoints[i].setup(4,m->gndim());
		// TODO: add allocation for for other types of faces, if needed

		// now check for faces which are not in rbf_boundaries
		inrflag = false;
		for(j = 0; j < rbf_boundaries.size(); j++)
			if(m->gbface(i,m->gnnofa()) == rbf_boundaries[j])
			{
				inrflag = true;
				break;
			}
		if( ! inrflag )
			for(inode = 0; inode < m->gnnofa(); inode++)
				rbfbpoin(m->gbface(i,inode)) = 0;
	}

	// get number of  boundary points required to be moved by RBF
	int nrbfcoun = 0;
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		nrbfcoun += rbfbpoin(ipoin);

	nbpoin = 0;
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		nbpoin += m->gflag_bpoin(ipoin);

	nbpoin -= nrbfcoun;

	ninpoin = m->gnpoin() - nbpoin;
	inpoints.setup(ninpoin,m->gndim());
	bounpoints.setup(nbpoin,m->gndim());
	boundisps.setup(nbpoin,m->gndim());

	amat::Matrix<amc_real> disps;			// displacement of each point in the mesh
	disps.setup(m->gnpoin(),m->gndim());
	disps.zeros();

	// get displacements for each boundary point by iterating over faces
	
	for(iface = 0; iface < m->gnface(); iface++)
	{
		if(m->gnnofa() == 6)
		{
			for(inode = 0; inode < 3; inode++)
			{
				jnode = (inode+1) % 3;		// next local node
				for(idim = 0; idim < m->gndim(); idim++)
				{
					midpoints[iface](inode,idim) = (m->gcoords(m->gbface(iface,inode),idim) + m->gcoords(m->gbface(iface,jnode),idim)) / 2.0;
					
					disps(m->gbface(iface,inode+3),idim) = m->gcoords(m->gbface(iface,inode+3),idim) - midpoints[iface].get(inode,idim);
				}
			}
		}
		else if(m->gnnofa() == 9)
		{
			for(idim = 0; idim < m->gndim(); idim++)
				midpoints[iface](4,idim) = 0;
			
			for(inode = 0; inode < 4; inode++)
			{
				jnode = (inode+1)%4;
				for(idim = 0; idim < m->gndim(); idim++)
				{
					// edge-center nodes
					midpoints[iface](inode,idim) = (m->gcoords(m->gbface(iface,inode),idim) + m->gcoords(m->gbface(iface,jnode),idim)) / 2.0;

					// for face-center node
					midpoints[iface](4,idim) += m->gcoords(m->gbface(iface,inode),idim);

					// displacements of edge-center nodes
					disps(m->gbface(iface,inode+4),idim) = m->gcoords(m->gbface(iface,inode+4),idim) - midpoints[iface].get(inode,idim);
				}
			}

			// displacements of face-center node
			for(idim = 0; idim < m->gndim(); idim++)
			{
				midpoints[iface](4,idim) /= 4;
				disps(m->gbface(iface,8),idim) = m->gcoords(m->gbface(iface,8),idim) - midpoints[iface].get(4,idim);
			}
		}
		else if(m->gnnofa() == 8)
		{
			for(inode = 0; inode < 4; inode++)
			{
				jnode = (inode+1)%4;
				for(idim = 0; idim < m->gndim(); idim++)
				{
					// edge-center nodes
					midpoints[iface](inode,idim) = (m->gcoords(m->gbface(iface,inode),idim) + m->gcoords(m->gbface(iface,jnode),idim)) / 2.0;

					// displacements of edge-center nodes
					disps(m->gbface(iface,inode+4),idim) = m->gcoords(m->gbface(iface,inode+4),idim) - midpoints[iface].get(inode,idim);
				}
			}
		}
		
		// TODO: if necessary, add handlers for additional types of boundary faces here

		else
		{
			std::cout << "CurvedMeshGen: ! Unknown face type!" << std::endl;
		}
	}

	// next, make boundary edges and faces straight by setting the high-order nodes to corresponding midpoints computed above
	
	for(iface = 0; iface < m->gnface(); iface++)
	{
		if(m->gnnofa() == 6)
		{
			for(inode = 0; inode < 3; inode++)
			{
				for(idim = 0; idim < m->gndim(); idim++)
				{
					m->scoords(m->gbface(iface,inode+3),idim, midpoints[iface](inode,idim));
				}
			}
		}
		else if(m->gnnofa() == 9)
		{
			// handler for hexahedral mesh
			for(inode = 0; inode < 4; inode++)
			{
				for(idim = 0; idim < m->gndim(); idim++)
				{
					// set coordinate of the point as the midpoint, making the edge straight
					m->scoords(m->gbface(iface,inode+4),idim, midpoints[iface].get(inode,idim));
				}
			}
			for(idim = 0; idim < m->gndim(); idim++)
				m->scoords(m->gbface(iface,8),idim, midpoints[iface].get(4,idim));
		}
		else if(m->gnnofa() == 8)
		{
			for(inode = 0; inode < 4; inode++)
			{
				for(idim = 0; idim < m->gndim(); idim++)
				{
					// set coordinate of the point as the midpoint, making the edge straight
					m->scoords(m->gbface(iface,inode+4),idim, midpoints[iface].get(inode,idim));
				}
			}
		}
		
		// TODO: if necessary, add handlers for additional types of boundary faces here
		else
		{
			std::cout << "CurvedMeshGen: !> Unknown face type!" << std::endl;
		}
	}
	delete [] midpoints;
	
	// estimate support radius
	amat::Matrix<amc_real> supportradius(m->gnpoin(),1), srad(nbpoin,1);
	amc_real sradius, dist;
	supportradius.zeros();
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		if(m->gflag_bpoin(ipoin)==1 && rbfbpoin.get(ipoin)==0)
		{
			dist = sqrt(disps.get(ipoin,0)*disps.get(ipoin,0) + disps.get(ipoin,1)*disps.get(ipoin,1) + disps.get(ipoin,2)*disps.get(ipoin,2));
			supportradius(ipoin) = 100.0*dist;
		}

	k = 0; l = 0;

	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
	{
		if(m->gflag_bpoin(ipoin) == 1 && rbfbpoin.get(ipoin)==0)
		{
			for(idim = 0; idim < m->gndim(); idim++)
			{
				bounpoints(k,idim) = m->gcoords(ipoin,idim);
				boundisps(k,idim) = disps.get(ipoin,idim);
			}
			srad(k) = supportradius.get(ipoin);
			k++;
		}
		else
		{
			for(idim = 0; idim < m->gndim(); idim++)
				inpoints(l,idim) = m->gcoords(ipoin,idim);
			l++;
		}
	}

	std::cout << "CurvedMeshGen: Assembled lists of boundary points, boundary displacements and interior points. k = " << k << ", nbpoin = " << nbpoin << "; l = " << l << ", ninpoin = " << ninpoin << std::endl;

	std::cout << "CurvedMeshGen: Norm of boundary displacements ";
	std::vector<amc_real> nrm(m->gndim(), 0.0);
	for(int i = 0; i < nbpoin; i++)
		for(idim = 0; idim < m->gndim(); idim++)
			nrm[idim] += boundisps(i,idim)*boundisps(i,idim);
	for(idim = 0; idim < m->gndim(); idim++)
	{
		nrm[idim] = sqrt(nrm[idim]);
		std::cout << nrm[idim] << ", ";
	}
	std::cout << std::endl;
	
	// get avg support radius
	amc_real ssum = 0, smin = 1.0, smax = 0;
	for(ipoin = 0; ipoin < nbpoin; ipoin++)
		ssum += srad.get(ipoin);
	ssum /= nbpoin;
	std::cout << "CurvedMeshGen: Average estimated support radius = " << ssum << std::endl;
	
	for(ipoin = 0; ipoin < nbpoin; ipoin++)
		if(smax < srad.get(ipoin))
			smax = srad.get(ipoin);
	std::cout << "CurvedMeshGen: Max estimated support radius = " << smax << std::endl;
	
	for(ipoin = 0; ipoin < nbpoin; ipoin++)
		sradius = (0.6*ssum+0.4*smax);
	std::cout << "CurvedMeshGen: Final estimated support radius = " << sradius << std::endl;
	
	if(param1 < ZERO_TOL) 
	{
		std::cout << "CurvedMeshGen: Using estimated support radius." << std::endl;
	}
	else
	{
		std::cout << "CurvedMeshGen: Using user-provided support radius.\n";
		sradius = param1;
	}

	move = new RBFmove(&inpoints, &bounpoints, &boundisps, choice, sradius, 1, tol, maxiter, solver );
}

CurvedMeshGen::~CurvedMeshGen()
{
	delete move;
}

void CurvedMeshGen::generateCurvedMesh()
{
	move->move();

	inpoints = move->getInteriorPoints();
	bounpoints = move->getBoundaryPoints();
	//amat::Matrix<amc_real> ncoords(m->gnpoin(), m->gndim());
	
	amc_int ipoin, idim, k = 0, l = 0;

	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		if(m->gflag_bpoin(ipoin) == 1 && rbfbpoin.get(ipoin)==0)
		{
			for(idim = 0; idim < m->gndim(); idim++)
				m->scoords(ipoin,idim, bounpoints.get(k,idim));
			k++;
		}
		else
		{
			for(idim = 0; idim < m->gndim(); idim++)
				m->scoords(ipoin,idim, inpoints.get(l,idim));
			l++;
		}
}

} // end namespace
