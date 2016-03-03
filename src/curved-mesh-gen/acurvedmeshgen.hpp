/** \brief Code to construct curved mesh using a quadratic mesh that has its boundary curved, and therefore probably has invalid elements.
 * \author Aditya Kashi
 * \date March 3, 2016
 */

#ifndef __AMESH3D_H
#include <amesh3d.hpp>
#endif

#define __ACURVEDMESHGEN_H

namespace amc {

/// Class to generate a fully curved mesh, using as input a high-order mesh with boundaries already curved
class CurvedMeshGen
{
	UMesh* m;							///< the mesh to curve
	RBFmove* move;						///< mesh movement context
	amat::Matrix<amc_real> inpoints;	///< interior points of the mesh
	amat::Matrix<amc_real> bounpoints;	///< boundary points
	std::vector<int> curvemarkers;		///< marker numbers of boundaries that need to be curved
	std::vector<int> torec;				///< contains 1 if the corresponding boundary node 
	amc_int nbpoin;						///< number of boundary points
	amc_int ninpoin;					///< number of interior points

public:

	/// constructor
	/** \param[in|out] mesh is the mesh to curve
	 * \param[in] choice corresponding to the mesh method method - for RBF, this is the type of RBF to use
	 * \param[in] param1 parameter required by mesh movement method - for RBF, this is the support radius
	 * \param[in] tol the tolerance to use in linear solvers, for instance
	 * \param[in] maxiter maximum iterations for linear solver
	 * \param[in] solver a string describing the linear solver to use ("CG", "LU")
	 */
	CurvedMeshGen(UMesh* mesh, const int choice, const double param1, const double tol, const int maxiter, const string solver);

	~CurvedMeshGen();

	void generateCurvedMesh();
};

CurvedMeshGen::CurvedMeshGen(UMesh* mesh, const int choice, const double param1, const double tol, const int maxiter, const string solver)
{
	m = mesh;
	
	amc_int ipoin, ielem, iface, inode, jnode, k;
	std::vector<amc_real> midpoint(m->gndim());
	
	nbpoin = 0;
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
		nbpoin += m->gflag_bpoin(ipoin);
	ninpoin = m->gnpoin() - nbpoin;
	inpoints.setup(ninpoin,m->gndim());
	bounpoints.setup(nbpoin,m->gndim());

	amat::Matrix<amc_real> disps;			// displacement of each point in the mesh
	disps.setup(m->gnpoin(),m->gndim());

	// get displacements for each boundary point by iterating over faces
	for(iface = 0; iface < m->gnface(); iface++)
	{
		if(m->gnnofa() == 6)
		{
			for(inode = 0; inode < 3; inode++)
			{
				jnode = (inode+1) % 3;		// next local node
				for(idim = 0; idim < m->gndim(); idim++)
					midpoint[idim] = (m->gcoords(m->gbface(iface,inode),idim) + m->gcoords(m->gbface(iface,jnode),idim)) / 2.0;
			}
		}
	}

	move = new RBFmove();
}

} // end namespace
