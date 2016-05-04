#ifndef __ABOUNDARYINFLUENCEDISTANCE_H

#define __ABOUNDARYINFLUENCEDISTANCE_H

#ifndef __AMESH2DH_H
#include <amesh2dh.hpp>
#endif

namespace amc {

void boundaryInfluenceDist2D(const UMesh2dh* const m, const amat::Matrix<amc_real>* const pointdisps, amat::Matrix<amc_real>* const radii);

}

#endif
