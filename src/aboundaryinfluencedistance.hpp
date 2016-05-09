#ifndef __ABOUNDARYINFLUENCEDISTANCE_H

#define __ABOUNDARYINFLUENCEDISTANCE_H

#ifndef __AMESH2DH_H
#include <amesh2dh.hpp>
#endif

namespace amc {

/// Computes a support radius for each boundary point of a mesh
/** The support radius is computed as
 * \f[
 * r_s = l/2 + g(h/l) h
 * \f]
 * where \f$ l \f$ is the length of the linear edge, \f$ h \f$ is the displacement of the edge and \f$ g \f$ is given by
 * \f[
 * g(y) = exp(y^2) \text{, if } y \leq 1/2 \\
 * g(y) = exp(1/4) y + 1/2 exp(1/4) \text{, if } y > 1/2
 * \f]
 *
 * \param[in] pointdisps contains displacements for all points in the mesh (only those corresponding to boundary points are used, however)
 * \param[in|out] radii will contain support radii on output, but its size is npoin by 1
 */
void boundaryInfluenceDist2D(const UMesh2dh* const m, const amat::Matrix<amc_real>* const pointdisps, amat::Matrix<amc_real>* const radii);

}

#endif
